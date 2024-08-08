import logging
import os
import gzip
import pysam
import itertools
import tempfile
import numpy as np
import pickle
import subprocess
import shutil
import scipy
from . import misc
from Bio import SeqIO
from collections import defaultdict, Counter
from contextlib import nullcontext
from glob import glob
from multiprocessing import Pool
from threading import BoundedSemaphore
from scipy.sparse import lil_matrix
from .bc_aligner import CustomBCAligner
from .bc_decoders import BCDecoder, SBCDecoder
from .umi import get_umi_maps_from_bam_file


log = logging.getLogger(__name__)
pysam.set_verbosity(0)


n_first_seqs = 10000  # n seqs for finding score threshold

def preprocess_fastqs(arguments):
    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)

    paired_fpaths = misc.find_paired_fastqs_in_dir(arguments.fastq_dir)
    log.info('Files to process:')
    for i, (fpath1, fpath2) in enumerate(paired_fpaths):
        log.info(f'  {fpath1}')
        log.info(f'  {fpath2}')
        if i < len(paired_fpaths)-1:
            log.info('  -')

    paired_align_fqs_and_tags_fpaths = []
    namepairidxs = []

    process_fastqs_func = serial_process_fastqs if arguments.threads == 1 else parallel_process_fastqs
    bcs_on_both_reads = misc.config_has_barcodes_on_both_reads(arguments.config)
    keep_r1 = arguments.config["barcode_struct_r1"]["keep_nonbarcode"]
    keep_r2 = arguments.config["barcode_struct_r2"]["keep_nonbarcode"] if bcs_on_both_reads else True
    if not bcs_on_both_reads:
        bc_fq_idx, paired_fq_idx = misc.determine_bc_and_paired_fastq_idxs(paired_fpaths, arguments.config["barcode_struct_r1"]["blocks"], arguments.config["unknown_read_orientation"])
        log.info(f'Detected barcodes in read{bc_fq_idx+1} files')
    else:
        bc_fq_idx, paired_fq_idx = 0, 1

    for fpath_tup in paired_fpaths:
        bc_fq1_fpath = fpath_tup[bc_fq_idx]
        bc_fq2_fpath = fpath_tup[paired_fq_idx]

        sans_bc_fq1_fpath = os.path.join(arguments.output_dir, f'sans_bc_{misc.file_prefix_from_fpath(bc_fq1_fpath)}.fq') if keep_r1 else None
        sans_bc_fq2_fpath = os.path.join(arguments.output_dir, f'sans_bc_{misc.file_prefix_from_fpath(bc_fq2_fpath)}.fq') if keep_r2 else None

        tags1_fpath = os.path.join(arguments.output_dir, f'rec_names_and_tags1_{misc.file_prefix_from_fpath(bc_fq1_fpath)}.txt')
        tags2_fpath = os.path.join(arguments.output_dir, f'rec_names_and_tags2_{misc.file_prefix_from_fpath(bc_fq2_fpath)}.txt')

        namepairidx_fpath = os.path.join(arguments.output_dir, f'namepairidx_{misc.file_prefix_from_fpath(bc_fq1_fpath)}.pkl')

        namepairidx = misc.get_namepair_index(bc_fq1_fpath, bc_fq2_fpath)
        namepairidxs.append(namepairidx)
        with open(namepairidx_fpath, "wb") as pkl:
            pickle.dump(namepairidx, pkl)


        log.info('Processing files:')
        if bcs_on_both_reads:
            log.info(f'  barcode fastq 1: {bc_fq1_fpath}')
            log.info(f'  barcode fastq 2: {bc_fq2_fpath}')
            if sans_bc_fq1_fpath:
                log.info(f'  sans-bc fastq 1: {sans_bc_fq1_fpath}')
            if sans_bc_fq2_fpath:
                log.info(f'  sans-bc fastq 2: {sans_bc_fq2_fpath}')
            log.info(f'  tags file 1:     {tags1_fpath}')
            log.info(f'  tags file 2:     {tags2_fpath}')
        else:
            log.info(f'  barcode fastq: {bc_fq1_fpath}')
            log.info(f'  paired fastq:  {bc_fq2_fpath}')
            if sans_bc_fq1_fpath:
                log.info(f'  sans-bc fastq: {sans_bc_fq1_fpath}')
            log.info(f'  sans-bc fastq: {sans_bc_fq2_fpath}')
            log.info(f'  tags file:     {tags1_fpath}')

        if all(os.path.exists(fpath) for fpath in [sans_bc_fq1_fpath, sans_bc_fq2_fpath, tags1_fpath, tags2_fpath] if fpath):
            log.info('Barcode output found. Skipping barcode detection')
        else:
            process_fastqs_func(
                    arguments,
                    bc_fq1_fpath,
                    bc_fq2_fpath,
                    sans_bc_fq1_fpath,
                    sans_bc_fq2_fpath,
                    tags1_fpath,
                    tags2_fpath,
                    bcs_on_both_reads)

        if bc_fq_idx == 0:
            paired_align_fqs_and_tags_fpaths.append((sans_bc_fq1_fpath, sans_bc_fq2_fpath, tags1_fpath, tags2_fpath))
        else:
            paired_align_fqs_and_tags_fpaths.append((sans_bc_fq2_fpath, sans_bc_fq1_fpath, tags2_fpath, tags1_fpath))
    return paired_align_fqs_and_tags_fpaths, namepairidxs

def process_fastqs(arguments):
    """
    Output single file with parsed bcs from bc_fastq in read names and seqs from paired_fastq.
    """
    star_w_bc_fpath = os.path.join(arguments.output_dir, 'with_bc.bam')
    star_w_bc_sorted_fpath = os.path.join(arguments.output_dir, 'with_bc.sorted.bam')
    star_w_bc_umi_fpath = os.path.join(arguments.output_dir, 'with_bc_umi.bam')
    star_w_bc_umi_sorted_fpath = os.path.join(arguments.output_dir, 'with_bc_umi.sorted.bam')
    if os.path.exists(star_w_bc_umi_sorted_fpath):
        log.info('Sorted bam with UMIs found. Skipping ahead...')
        completed = 4
    elif os.path.exists(star_w_bc_umi_fpath):
        log.info('Corrected UMI file found. Skipping ahead...')
        completed = 3
    elif os.path.exists(star_w_bc_sorted_fpath):
        log.info('Sorted STAR results found. Skipping ahead...')
        completed = 2
    elif os.path.exists(star_w_bc_fpath):
        log.info('STAR results found. Skipping ahead...')
        completed = 1
    else:
        completed = 0

    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)
    # find barcodes, write to file, remove from reads before mapping
    # map reads
    # Add barcodes etc to bam file
    # Sort, index
    if completed < 1:
        log.info('Writing output to:')
        log.info(f'  {star_w_bc_fpath}')
    
        if not arguments.star_output_path:
            paired_align_fqs_and_tags_fpaths, namepairidxs = preprocess_fastqs(arguments)

            log.info(f'Running STAR alignment...')
            star_out_dirs = set()
            star_bam_and_tags_fpaths = []
            for R1_fpath, R2_fpath, tags1_fpath, tags2_fpath in paired_align_fqs_and_tags_fpaths:
                if R1_fpath:
                    log.info(f'  {R1_fpath}')
                log.info(f'  {R2_fpath}')
                if R1_fpath != paired_align_fqs_and_tags_fpaths[-1][0]:
                    log.info('  -')
                star_out_dir, star_out_fpath = run_STAR(arguments, R1_fpath, R2_fpath)
                star_out_dirs.add(star_out_dir)
                star_bam_and_tags_fpaths.append((star_out_fpath, tags1_fpath, tags2_fpath))
        else:
            log.info(f'Using STAR results from {arguments.star_output_path}')
            bam_fpath = glob(os.path.join(arguments.star_output_path, "*Aligned.out.bam")).sort()
            tag_fpaths = glob(os.path.join(arguments.fastq_dir, "rec_names_and_tags_*.txt")).sort()
            namepairidx_fpaths = glob(os.path.join(arguments.fastq_dir, "rec_names_and_tags_*.pkl")).sort()

            star_bam_and_tags_fpaths = list(zip(bam_fpath, tag_fpaths[::2], tag_fpaths[1::2]))
            namepairidxs = [pickle.load(open(fpath, "rb")) for fpath in namepairidx_fpaths]

        log.info('Adding barcode tags to mapped reads...')
        template_bam_fpath = star_bam_and_tags_fpaths[0][0]
        with pysam.AlignmentFile(star_w_bc_fpath, 'wb', template=pysam.AlignmentFile(template_bam_fpath)) as bam_out:
            for (star_out_fpath, tags1_fpath, tags2_fpath), namepairidx in zip(star_bam_and_tags_fpaths, namepairidxs):
                if arguments.threads == 1 and not arguments.star_output_path:
                    readsiter = iter(serial_add_tags_to_reads(tags1_fpath, tags2_fpath, star_out_fpath))
                else:
                    readsiter = iter(parallel_add_tags_to_reads(tags1_fpath, tags2_fpath, star_out_fpath, namepairidx, arguments.threads))
                for read in readsiter:
                    bam_out.write(read)

        for fpaths in paired_align_fqs_and_tags_fpaths:
            for fpath in fpaths:
                if fpath and os.path.isfile(fpath):
                    os.remove(fpath)
        for star_out_dir in star_out_dirs:
            shutil.rmtree(star_out_dir)  # clean up intermediate STAR files
    
    if completed < 2:
        log.info('Sorting bam...')
        pysam.sort('-@', str(arguments.threads), '-o', star_w_bc_sorted_fpath, star_w_bc_fpath)
        os.remove(star_w_bc_fpath)  #clean up unsorted bam
        log.info('Indexing bam...')
        pysam.index(star_w_bc_sorted_fpath)

    if completed < 3:
        log.info('Correcting UMIs...')
        correct_UMIs(arguments, star_w_bc_sorted_fpath, star_w_bc_umi_fpath)

    if completed < 4:
        log.info('Sorting bam...')
        pysam.sort('-@', str(arguments.threads), '-o', star_w_bc_umi_sorted_fpath, star_w_bc_umi_fpath)
        log.info('Indexing bam...')
        pysam.index(star_w_bc_umi_sorted_fpath)
        os.remove(star_w_bc_umi_fpath)
        os.remove(star_w_bc_sorted_fpath)
        os.remove(star_w_bc_sorted_fpath + '.bai')

    count_matrix(arguments, star_w_bc_umi_sorted_fpath)
    log.info('Done')


def run_STAR(arguments, R1_fpath, R2_fpath):
    """
    Run STAR aligner for gDNA files.

    Returns STAR output directory and bam path.
    """
    star_out_dir = os.path.join(arguments.output_dir, 'STAR_files')
    R2_bname = misc.file_prefix_from_fpath(R2_fpath)
    out_prefix = os.path.join(star_out_dir, f'{R2_bname}_')
    cmd_star = [
        'STAR',
        f'--runThreadN {arguments.threads}',
        f'--genomeDir {arguments.star_ref_dir}',
        f'--readFilesIn {R1_fpath if R1_fpath else ""} {R2_fpath}',
        f'--outFileNamePrefix {out_prefix}',
        '--outFilterMultimapNmax 1',
        '--outSAMtype BAM Unsorted',
        '--outSAMattributes NH HI AS nM GX GN',
    ]
    if R2_fpath.endswith('.gz'):
        if R1_fpath and not R1_fpath.endswith('.gz'):
            raise ValueError('Paired read files must be both zipped or both unzipped')
        cmd_star.append('--readFilesCommand zcat')
    star_out_fpath = f'{out_prefix}Aligned.out.bam'
    if os.path.exists(star_out_fpath):
        log.info("STAR results found. Skipping alignment")
    else:
        subprocess.run(cmd_star, check=True)
    return star_out_dir, star_out_fpath


def process_bc_rec(blocks, bc_rec, aligners, decoders):
    """
    Find barcodes etc in bc_rec 
    """
    scores_pieces_end_pos = [al.find_norm_score_pieces_and_end_pos(bc_rec.seq, return_seq=True) for al in aligners]
    raw_score, raw_pieces, raw_end_pos, raw_seq = max(scores_pieces_end_pos)
    raw_bcs = [raw_piece.upper().replace('N', 'A') for raw_piece, block in zip(raw_pieces, blocks) if block["blocktype"] == "barcodeList"]
    bcs = [decoder.decode(raw_bc) for raw_bc, decoder in zip(raw_bcs, decoders)]

    best_aligner = next(al for al, (s, p, e, sq) in zip(aligners, scores_pieces_end_pos) if s == raw_score)
    commonseqs = [best_aligner.prefixes[i] for i, block in enumerate(blocks) if block["blocktype"] == "constantRegion"]

    bc_idx = commonseq_idx = 0
    corrected_pieces = []
    for block in blocks:
        if block["blocktype"] == "barcodeList":
            if bcs[bc_idx] is None:
                return raw_score, None, None
            corrected_pieces.append(bcs[bc_idx])
            bc_idx += 1
        elif block["blocktype"] == "constantRegion":
            if commonseqs[commonseq_idx] is None:
                return raw_score, None, None
            corrected_pieces.append(commonseqs[commonseq_idx])
            commonseq_idx += 1
        else:
            corrected_pieces.append('N' * block["length"])

    new_aligner = CustomBCAligner(*corrected_pieces)
    new_score, new_pieces = new_aligner.find_norm_score_and_pieces(raw_seq)

    bam_bcs = []
    bam_raw_bcs = []
    bam_umis = []
    bam_commonseqs = []

    bc_idx = commonseq_idx = 0
    for i, block in enumerate(blocks):
        if block["blockfunction"] != "discard":
            if block["blocktype"] == "barcodeList":
                bam_bcs.append(bcs[bc_idx])
                bam_raw_bcs.append(raw_pieces[i])
                bc_idx += 1
            elif block["blocktype"] == "constantRegion":
                bam_commonseqs.append(commonseqs[commonseq_idx])
                commonseq_idx += 1
            else:
                if block["blockfunction"] == "UMI":
                    bam_umis.append(new_pieces[i])
                else:
                    bam_bcs.append(new_pieces[i])
                    bam_raw_bcs.append(new_pieces[i])


    tags = []
    # Cell barcode
    tags.append(('CB', '.'.join(bam_bcs)))
    tags.append(('CR', '.'.join(bam_raw_bcs)))
    # Filler sequences
    tags.append(('FL', sum(len(seq) for seq in bam_commonseqs)))
    # And raw UMI
    tags.append(("UR", ".".join(bam_umis)))

    new_rec = bc_rec[raw_end_pos:]
    new_rec.id = bc_rec.id
    new_rec.description = bc_rec.description

    return raw_score, new_rec, tags

def umi_parallel_wrapper(ref_and_input_bam_fpath):
    ref, input_bam_fpath = ref_and_input_bam_fpath
    return ref, get_umi_maps_from_bam_file(input_bam_fpath, chrm=ref)

def correct_UMIs(arguments, input_bam_fpath, out_bam_fpath):
    with pysam.AlignmentFile(input_bam_fpath) as bamfile:
        reference_names = bamfile.references
    reference_names_with_input_bam = [(ref, input_bam_fpath) for ref in reference_names]

    with pysam.AlignmentFile(out_bam_fpath, 'wb', template=pysam.AlignmentFile(input_bam_fpath)) as bam_out, \
            Pool(arguments.threads) as pool:
        for i, (ref, umi_map_given_bc_then_feature) in enumerate(pool.imap_unordered(
                umi_parallel_wrapper,
                reference_names_with_input_bam)):
            log.info(f'  {ref}')
            for read in pysam.AlignmentFile(input_bam_fpath).fetch(ref):
                for gx_gn_tup in misc.gx_gn_tups_from_read(read):
                    corrected_umi = umi_map_given_bc_then_feature[read.get_tag('CB')][gx_gn_tup][read.get_tag('UR')]
                    read.set_tag('UB', corrected_umi)
                    bam_out.write(read)
                    break # use the first one. Ideally same across all

def build_tags_iter(tags_fpath):
    for line in open(tags_fpath):
        words = line.strip().split('\t')
        read_name = words[0]
        tags = [parse_tag_str(word) for word in words[1:]]
        yield read_name, tags

def merge_tags(tags1, tags2):
    if not len(tags2):
        return tags1
    tags1 = defaultdict(tags1, str)
    tags2 = dict(tags2)
    for name, val in tags2.items():
        tags1[name] = ".".join(tags1[name], val)
    return tags1.items()

def serial_add_tags_to_reads(tags1_fpath, tags2_fpath, bam_fpath):
    """
    Iterates bam reads and matching tags records and adds tags to reads.
    """
    tags1_iter = build_tags_iter(tags1_fpath)
    tags2_iter = build_tags_iter(tags2_fpath) if os.path.isfile(tags2_fpath) else itertools.repeat(("", ()))
    read_name = 'Lorem ipsum'
    for read in pysam.AlignmentFile(bam_fpath).fetch(until_eof=True):
        while not misc.names_pair(read_name, str(read.query_name)):
            read_name, tags1 = next(tags1_iter)
            _, tags2 = next(tags2_iter)
        for name, val in merge_tags(tags1, tags2):
            read.set_tag(name, val)
        yield read

def parallel_add_tags_to_reads(tags1_fpath, tags2_fpath, bam_fpath, namepairidx, threads):
    """
    Iterates bam reads and matching tags records and adds tags to reads.
    """
    sortedbampath = bam_fpath + "_readname_sorted.bam"
    bamidx = misc.sort_and_index_readname_bam(bam_fpath, sortedbampath, namepairidx, threads)
    tags1_iter = build_tags_iter(tags1_fpath)
    tags2_iter = build_tags_iter(tags2_fpath) if os.path.isfile(tags2_fpath) else itertools.repeat(("", ()))
    with pysam.AlignmentFile(sortedbampath) as bam:
        for (read_name, tags1), (_, tags2) in zip(tags1_iter, tags2_iter):
            for read in misc.get_bam_read_by_name(read_name, bam, bamidx):
                for name, val in merge_tags(tags1, tags2):
                    read.set_tag(name, val)
                yield read
    os.remove(sortedbampath)

def tag_type_from_val(val):
    if isinstance(val, str):
        return 'Z'
    elif isinstance(val, int):
        return 'i'
    else:
        # no other types implemented
        raise ValueError('Unexpected tag type')


def parse_tag_str(tag_str):
    tag, tag_type, val = tag_str.split(':')
    if tag_type == 'i':
        # only worry about strings and ints
        val = int(val)
    return tag, val


def output_rec_name_and_tags(rec, tags, out_fh):
    out_fh.write('\t'.join([f'{str(rec.id)}'] + [f'{tag}:{tag_type_from_val(val)}:{val}' for tag, val in tags]) + '\n')


def serial_process_fastqs(arguments, fq1_fpath, fq2_fpath, sans_bc_fq1_fpath, sans_bc_fq2_fpath, tags1_fpath, tags2_fpath, bcs_on_both_reads):
    blocks = [arguments.config["barcode_struct_r1"]["blocks"]]
    if bcs_on_both_reads:
        blocks.append(arguments.config["barcode_struct_r2"]["blocks"])

    log.info('Building aligners and barcode decoders')
    aligners = tuple(misc.build_bc_aligners(bblocks, arguments.config["unknown_read_orientation"]) for bblocks in blocks)
    decoders = tuple(misc.build_bc_decoders(bblocks) for bblocks in blocks)

    log.info(f'Processing first {n_first_seqs:,d} for score threshold...')

    thresh = [np.inf, -np.inf]
    for i, (fq_fpath, bblocks, baligners, bdecoders) in enumerate(zip((fq1_fpath, fq2_fpath), blocks, aligners, decoders, strict=False)):
        first_scores_recs_tags = []
        for j, bc_rec in enumerate(SeqIO.parse(misc.gzip_friendly_open(fq_fpath), 'fastq')):
            first_scores_recs_tags.append(process_bc_rec(bblocks, bc_rec, baligners, bdecoders))
            if j >= n_first_seqs:
                break

        scores = [score for score, rec, tags in first_scores_recs_tags]
        thresh[i] = np.average(scores) - 2 * np.std(scores)

    log.info(f'Score threshold: {thresh[0]:.2f}{f"{thresh[1]:.2f}" if thresh[1] != -np.inf else ""}')

    total_out = 0
    with (open(sans_bc_fq1_fpath, 'w') if sans_bc_fq1_fpath else nullcontext(None)) as bc_fq1_fh, \
            (open(sans_bc_fq2_fpath, 'w') if sans_bc_fq2_fpath else nullcontext(None)) as bc_fq2_fh, \
            open(tags1_fpath, 'w') as tag1_fh, \
            (open(tags2_fpath, 'w') if bcs_on_both_reads else nullcontext(None)) as tag2_fh:
        log.info('Continuing...')
        for i, (bc_rec1, bc_rec2) in enumerate(zip(
            SeqIO.parse(misc.gzip_friendly_open(fq1_fpath), 'fastq'),
            SeqIO.parse(misc.gzip_friendly_open(fq2_fpath), 'fastq'))):
            if i % 100000 == 0 and i > 0:
                log.info(f'  {i:,d}')

            scores = [0, 0]
            sans_bc_rec = [None, True]
            tags = [None, None]
            for j, (bblocks, bc_rec, baligners, bdecoders) in enumerate(zip(blocks, (bc_rec1, bc_rec2), aligners, decoders, strict=False)):
                scores[j], sans_bc_rec[j], tags[j] = process_bc_rec(bblocks, bc_rec, baligners, bdecoders)
            if all(score >= cthresh for score, cthresh in zip(scores, thresh)) and all(sans_bc_rec):
                total_out += 1
                for bc_rec, bc_fq_fh, csans_bc_rec, ctags, tag_fh in zip((bc_rec1, bc_rec2), (bc_fq1_fh, bc_fq2_fh), sans_bc_rec, tags, (tag1_fh, tag2_fh)):
                    if bc_fq_fh:
                        if csans_bc_rec and csans_bc_rec is not True:
                            SeqIO.write(csans_bc_rec, bc_fq_fh, 'fastq')
                        else:
                            SeqIO.write(bc_rec, bc_fq_fh, 'fastq')
                    if ctags:
                        output_rec_name_and_tags(bc_rec, ctags, tag_fh)
    log.info(f'{i+1:,d} barcode records processed')
    log.info(f'{total_out:,d} pairs of records output')

def worker_build_aligners(ablocks, unknown_read_orientation, athresh):
    global blocks, aligners, decoders, thresh
    blocks = ablocks
    thresh = athresh
    aligners = tuple(misc.build_bc_aligners(bblocks, unknown_read_orientation) for bblocks in blocks)
    decoders = tuple(misc.build_bc_decoders(bblocks) for bblocks in blocks)

def worker_process_read(bc_recs):
    bc_rec1, bc_rec2 = bc_recs

    scores = [0, 0]
    sans_bc_rec = [None, True]
    tags = [None, None]
    for i, (bblocks, bc_rec, baligners, bdecoders) in enumerate(zip(blocks, (bc_rec1, bc_rec2), aligners, decoders, strict=False)):
        scores[i], sans_bc_rec[i], tags[i] = process_bc_rec(bblocks, bc_rec, baligners, bdecoders)

    output_recs = [None, None]
    output_tags = [None, None]
    if all(score >= cthresh for score, cthresh in zip(scores, thresh)) and all(sans_bc_rec):
        # TODO: don't write second fastq if barcodes only in first. Need to do this currently in case input
        # fastqs are gzipped: STAR input must be either all gzipped or all uncompressed
        for i, (bc_rec, csans_bc_rec, ctags) in enumerate(zip((bc_rec1, bc_rec2), sans_bc_rec, tags)):
            output_recs[i] = csans_bc_rec if csans_bc_rec and csans_bc_rec is not True else bc_rec
            output_tags[i] = ctags
    return output_recs, output_tags

def read_iterator(fq1_fpath, fq2_fpath, sem):
    for bc_rec1, bc_rec2 in zip(SeqIO.parse(misc.gzip_friendly_open(fq1_fpath), 'fastq'), SeqIO.parse(misc.gzip_friendly_open(fq2_fpath), 'fastq')):
        sem.acquire()
        yield bc_rec1, bc_rec2

def parallel_process_fastqs(arguments, fq1_fpath, fq2_fpath, sans_bc_fq1_fpath, sans_bc_fq2_fpath, tags1_fpath, tags2_fpath, bcs_on_both_reads):
    chunksize = 100

    blocks = [arguments.config["barcode_struct_r1"]["blocks"]]
    if bcs_on_both_reads:
        blocks.append(arguments.config["barcode_struct_r2"]["blocks"])

    log.info('Building aligners and barcode decoders')
    aligners = tuple(misc.build_bc_aligners(bblocks, arguments.config["unknown_read_orientation"]) for bblocks in blocks)
    decoders = tuple(misc.build_bc_decoders(bblocks) for bblocks in blocks)

    log.info(f'Processing first {n_first_seqs:,d} for score threshold...')

    thresh = [np.inf, -np.inf]
    for i, (fq_fpath, bblocks, baligners, bdecoders) in enumerate(zip((fq1_fpath, fq2_fpath), blocks, aligners, decoders, strict=False)):
        first_scores_recs_tags = []
        for j, bc_rec in enumerate(SeqIO.parse(misc.gzip_friendly_open(fq_fpath), 'fastq')):
            first_scores_recs_tags.append(process_bc_rec(bblocks, bc_rec, baligners, bdecoders))
            if j >= n_first_seqs:
                break

        scores = [score for score, rec, tags in first_scores_recs_tags]
        thresh[i] = np.average(scores) - 2 * np.std(scores)

    log.info(f'Score threshold: {thresh[0]:.2f}{f"{thresh[1]:.2f}" if thresh[1] != -np.inf else ""}')

    # multiprocessing.pool.imap tends to prefetch too many values into memory. The semaphore prevents that.
    sem = BoundedSemaphore(3 * chunksize * arguments.threads)
    with (open(sans_bc_fq1_fpath, 'w') if sans_bc_fq1_fpath else nullcontext(None)) as bc_fq1_fh, \
            (open(sans_bc_fq2_fpath, 'w') if sans_bc_fq2_fpath else nullcontext(None)) as bc_fq2_fh, \
            open(tags1_fpath, 'w') as tag1_fh, \
            (open(tags2_fpath, 'w') if bcs_on_both_reads else nullcontext(None)) as tag2_fh, \
            Pool(arguments.threads, initializer = worker_build_aligners, initargs=(blocks, arguments.config["unknown_read_orientation"], thresh)) as pool:
        log.info('Continuing...')
        total_out = 0
        for i, (output_recs, output_tags) in enumerate(pool.imap(worker_process_read, read_iterator(fq1_fpath, fq2_fpath, sem), chunksize=chunksize)):
            if i % 100000 == 0 and i > 0:
                log.info(f'  {i:,d}')
            sem.release()

            processed = False
            for rec, bc_fh, tags, tag_fh in zip(output_recs, (bc_fq1_fh, bc_fq2_fh), output_tags, (tag1_fh, tag2_fh)):
                if rec and bc_fh:
                    processed = True
                    SeqIO.write(rec, bc_fh, 'fastq')
                if tags:
                    output_rec_name_and_tags(rec, tags, tag_fh)
            total_out += processed
    log.info(f'{i+1:,d} barcode records processed')
    log.info(f'{total_out:,d} pairs of records output')


def count_parallel_wrapper(ref_and_input_bam_fpath):
    ref, input_bam_fpath = ref_and_input_bam_fpath
    read_count_given_bc_then_feature_then_umi = defaultdict(lambda : defaultdict(Counter))
    for read in pysam.AlignmentFile(input_bam_fpath).fetch(ref):
        if not read.is_paired or read.is_read1 or (read.is_read2 and read.mate_is_unmapped):
            for gx_gn_tup in misc.gx_gn_tups_from_read(read): # count read toward all compatible genes
                read_count_given_bc_then_feature_then_umi[read.get_tag('CB')][gx_gn_tup][read.get_tag('UB')] += 1
    for k in read_count_given_bc_then_feature_then_umi.keys():
        read_count_given_bc_then_feature_then_umi[k] = dict(read_count_given_bc_then_feature_then_umi[k])
    read_count_given_bc_then_feature_then_umi = dict(read_count_given_bc_then_feature_then_umi)
    return ref, read_count_given_bc_then_feature_then_umi


def count_matrix(arguments, input_bam_fpath):
    """
    Counts the reads from the input bam file and outputs sparse matrices of read and UMI counts.
    """
    raw_reads_output_dir = os.path.join(arguments.output_dir, 'raw_reads_bc_matrix')
    raw_umis_output_dir = os.path.join(arguments.output_dir, 'raw_umis_bc_matrix')
    for out_dir in [raw_reads_output_dir, raw_umis_output_dir]:
        if os.path.exists(out_dir):
            log.info('Matrix output folder exists. Skipping count matrix build')
            return
        else:
            os.makedirs(out_dir)

    log.info('Finding all barcodes present...')
    sorted_complete_bcs, sorted_features = misc.get_bcs_and_features_from_bam(
            input_bam_fpath,
            arguments.threads-1
            )
    i_given_feature = {feat: i for i, feat in enumerate(sorted_features)}
    j_given_complete_bc = {comp_bc: j for j, comp_bc in enumerate(sorted_complete_bcs)}


    with pysam.AlignmentFile(input_bam_fpath) as bamfile:
        reference_names = bamfile.references
    reference_names_with_input_bam = [(ref, input_bam_fpath) for ref in reference_names]


    log.info('Counting reads...')
    M_reads = lil_matrix((len(sorted_features), len(sorted_complete_bcs)), dtype=int)
    M_umis = lil_matrix((len(sorted_features), len(sorted_complete_bcs)), dtype=int)
    with Pool(arguments.threads) as pool:
        for ref, read_count_given_bc_then_feature_then_umi in pool.imap_unordered(
                count_parallel_wrapper,
                reference_names_with_input_bam):
            for comp_bc, read_count_given_feature_then_umi in read_count_given_bc_then_feature_then_umi.items():
                j = j_given_complete_bc[comp_bc]
                for gx_gn_tup, umi_cntr in read_count_given_feature_then_umi.items():
                    i = i_given_feature[gx_gn_tup]
                    for umi, count in umi_cntr.items():
                        M_reads[i, j] += count
                        M_umis[i, j] += 1

    log.info('Writing raw read count matrix...')
    for out_dir, M in [(raw_reads_output_dir, M_reads), (raw_umis_output_dir, M_umis)]:
        misc.write_matrix(M, sorted_complete_bcs, sorted_features, out_dir)
