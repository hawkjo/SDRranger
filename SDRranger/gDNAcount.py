import logging
import os
import gzip
import pysam
import itertools
import tempfile
import numpy as np
import subprocess
import shutil
import scipy
from . import misc
from Bio import SeqIO
from collections import defaultdict, Counter
from multiprocessing import Pool
from scipy.sparse import lil_matrix
from .bc_aligner import CustomBCAligner
from .bc_decoders import BCDecoder, SBCDecoder
from .constants import commonseq1_options, commonseq2_gDNA 


log = logging.getLogger(__name__)
pysam.set_verbosity(0)


n_first_seqs = 10000  # n seqs for finding score threshold


def process_gDNA_fastqs(arguments):
    """
    Output single file with parsed bcs from bc_fastq in read names and seqs from paired_fastq.
    """
    star_w_bc_fpath = os.path.join(arguments.output_dir, 'gDNA_with_bc.bam')
    star_w_bc_sorted_fpath = os.path.join(arguments.output_dir, 'gDNA_with_bc.sorted.bam')
    if os.path.exists(star_w_bc_sorted_fpath):
        log.info('Sorted STAR results found. Skipping ahead...')
        completed = 2
    elif os.path.exists(star_w_bc_fpath):
        log.info('STAR results found. Skipping ahead...')
        completed = 1
    else:
        completed = 0

    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)

    paired_fpaths = misc.find_paired_fastqs_in_dir(arguments.fastq_dir)
    log.info('Files to process:')
    for i, (fpath1, fpath2) in enumerate(paired_fpaths):
        log.info(f'  {fpath1}')
        log.info(f'  {fpath2}')
        if i < len(paired_fpaths)-1:
            log.info('  -')

    # find barcodes, write to file, remove from reads before mapping
    # map reads
    # Add barcodes etc to bam file
    # Sort, index
    if completed < 1:
        bc_fq_idx, paired_fq_idx = misc.determine_bc_and_paired_fastq_idxs(paired_fpaths, 'gDNA') 
        log.info(f'Detected barcodes in read{bc_fq_idx+1} files')
        log.info('Writing output to:')
        log.info(f'  {star_w_bc_fpath}')
    
        paired_align_fqs_and_tags_fpaths = []
        process_fastqs_func = serial_process_gDNA_fastqs if arguments.threads == 1 else parallel_process_gDNA_fastqs
        for fpath_tup in paired_fpaths:
            bc_fq_fpath = fpath_tup[bc_fq_idx]
            paired_fq_fpath = fpath_tup[paired_fq_idx]
            sans_bc_fq_fname = f'sans_bc_{misc.file_prefix_from_fpath(bc_fq_fpath)}.fq'
            sans_bc_fq_fpath = os.path.join(arguments.output_dir, sans_bc_fq_fname)
            sans_bc_paired_fq_fname = f'sans_paired_{misc.file_prefix_from_fpath(paired_fq_fpath)}.fq'
            sans_bc_paired_fq_fpath = os.path.join(arguments.output_dir, sans_bc_paired_fq_fname)
            tags_fname = f'rec_names_and_tags_{misc.file_prefix_from_fpath(bc_fq_fpath)}.txt'
            tags_fpath = os.path.join(arguments.output_dir, tags_fname)
            if bc_fq_idx == 0:
                paired_align_fqs_and_tags_fpaths.append((sans_bc_fq_fpath, sans_bc_paired_fq_fpath, tags_fpath))
            else:
                paired_align_fqs_and_tags_fpaths.append((sans_bc_paired_fq_fpath, sans_bc_fq_fpath, tags_fpath))

            log.info('Processing files:')
            log.info(f'  barcode fastq: {bc_fq_fpath}')
            log.info(f'  paired fastq:  {paired_fq_fpath}')
            log.info(f'  sans-bc fastq: {sans_bc_fq_fpath}')
            log.info(f'  sans-bc fastq: {sans_bc_paired_fq_fpath}')
            log.info(f'  tags file:     {tags_fpath}')
            if all(os.path.exists(fpath) for fpath in [sans_bc_fq_fpath, sans_bc_paired_fq_fpath, tags_fpath]):
                log.info('Barcode output found. Skipping barcode detection')
            else:
                process_fastqs_func(
                        arguments, 
                        bc_fq_fpath, 
                        paired_fq_fpath, 
                        sans_bc_fq_fpath,
                        sans_bc_paired_fq_fpath, 
                        tags_fpath)

        log.info(f'Running STAR alignment...')
        star_out_dirs = set()
        star_bam_and_tags_fpaths = []
        for R1_fpath, R2_fpath, tags_fpath in paired_align_fqs_and_tags_fpaths:
            log.info(f'  {R1_fpath}')
            log.info(f'  {R2_fpath}')
            if R1_fpath != paired_align_fqs_and_tags_fpaths[-1][0]:
                log.info('  -')
            star_out_dir, star_out_fpath = run_STAR_gDNA(arguments, R1_fpath, R2_fpath)
            star_out_dirs.add(star_out_dir)
            star_bam_and_tags_fpaths.append((star_out_fpath, tags_fpath))

        log.info('Adding barcode tags to mapped reads...')
        template_bam_fpath = star_bam_and_tags_fpaths[0][0]
        with pysam.AlignmentFile(star_w_bc_fpath, 'wb', template=pysam.AlignmentFile(template_bam_fpath)) as bam_out:
            for star_out_fpath, tags_fpath in star_bam_and_tags_fpaths:
                for read in gDNA_add_tags_to_reads(tags_fpath, star_out_fpath):
                    bam_out.write(read)
        for fpath in (sans_bc_fq_fpath, sans_bc_paired_fq_fpath, tags_fpath):
            os.remove(fpath)
        for star_out_dir in star_out_dirs:
            shutil.rmtree(star_out_dir)  # clean up intermediate STAR files
    
    if completed < 2:
        log.info('Sorting bam...')
        pysam.sort('-@', str(arguments.threads), '-o', star_w_bc_sorted_fpath, star_w_bc_fpath)
        os.remove(star_w_bc_fpath)  #clean up unsorted bam
        log.info('Indexing bam...')
        pysam.index(star_w_bc_sorted_fpath)

    gDNA_count_matrix(arguments, star_w_bc_sorted_fpath)
    log.info('Done')


def run_STAR_gDNA(arguments, R1_fpath, R2_fpath):
    """
    Run STAR aligner for gDNA files.

    Returns STAR output directory and bam path.
    """
    star_out_dir = os.path.join(arguments.output_dir, 'STAR_files')
    R1_bname = misc.file_prefix_from_fpath(R1_fpath)
    out_prefix = os.path.join(star_out_dir, f'{R1_bname}_')
    cmd_star = [
        'STAR',
        f'--runThreadN 1', # required to keep order matching with fastq file
        f'--genomeDir {arguments.star_ref_dir}',
        f'--readFilesIn {R1_fpath} {R2_fpath}',
        f'--outFileNamePrefix {out_prefix}',
        '--outFilterMultimapNmax 1', 
        '--outSAMtype BAM Unsorted', 
        '--outSAMattributes NH HI AS nM GX GN',
    ]
    if R1_fpath.endswith('.gz'):
        if not R2_fpath.endswith('.gz'):
            raise ValueError('Paired read files must be both zipped or both unzipped')
        cmd_star.append('--readFilesCommand zcat')
    star_out_fpath = f'{out_prefix}Aligned.out.bam'
    if os.path.exists(star_out_fpath):
        log.info("STAR results found. Skipping alignment")
    else:
        subprocess.run(cmd_star, check=True)
    return star_out_dir, star_out_fpath


def process_bc_rec(bc_rec, aligners, bcd):
    """
    Find barcodes etc in bc_rec 
    """
    bc_seq = str(bc_rec.seq)
    scores_pieces_end_pos = [al.find_norm_score_pieces_and_end_pos(bc_seq) for al in aligners]
    raw_score, raw_pieces, raw_end_pos = max(scores_pieces_end_pos)
    raw_bc1, raw_bc2 = [raw_pieces[i].upper().replace('N', 'A') for i in [0, 2]]

    bc1 = bcd.decode(raw_bc1)
    bc2 = bcd.decode(raw_bc2)
    best_aligner = next(al for al, (s, p, e) in zip(aligners, scores_pieces_end_pos) if s == raw_score)
    commonseq1, commonseq2 = [best_aligner.prefixes[i] for i in [1, 3]]
    corrected_pieces = [bc1, commonseq1, bc2, commonseq2]
    if None in corrected_pieces:
        return raw_score, None, None

    tags = []
    # Cell barcode
    tags.append(('CB', f'{bc1}.{bc2}'))
    tags.append(('CR', f'{raw_pieces[0]}.{raw_pieces[2]}'))
    # Filler sequences
    tags.append(('FL', len(commonseq1) + len(commonseq2)))

    new_rec = bc_rec[raw_end_pos:]
    new_rec.id = bc_rec.id
    new_rec.description = bc_rec.description

    return raw_score, new_rec, tags


def build_tags_iter(tags_fpath):
    for line in open(tags_fpath):
        words = line.strip().split('\t')
        read_name = words[0]
        tags = [parse_tag_str(word) for word in words[1:]]
        yield read_name, tags


def gDNA_add_tags_to_reads(tags_fpath, bam_fpath):
    """
    Iterates bam reads and matching tags records and adds tags to reads.
    """
    tags_iter = build_tags_iter(tags_fpath)
    read_name = 'Lorem ipsum'
    for read in pysam.AlignmentFile(bam_fpath).fetch(until_eof=True):
        while not misc.names_pair(read_name, str(read.qname)):
            read_name, tags = next(tags_iter)
        for name, val in tags:
            read.set_tag(name, val)
        yield read


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


def serial_process_gDNA_fastqs(arguments, bc_fq_fpath, paired_fq_fpath, sans_bc_fq_fpath, sans_bc_paired_fq_fpath, tags_fpath):
    log.info('Building aligners and barcode decoders')
    aligners = misc.build_gDNA_bc_aligners()
    bcd = BCDecoder(arguments.barcode_whitelist, arguments.max_bc_err_decode)

    log.info(f'Processing first {n_first_seqs:,d} for score threshold...')
    first_scores_recs_tags = []
    for i, bc_rec in enumerate(SeqIO.parse(misc.gzip_friendly_open(bc_fq_fpath), 'fastq')):
        first_scores_recs_tags.append(process_bc_rec(bc_rec, aligners, bcd))
        if i >= n_first_seqs:
            break

    scores = [score for score, rec, tags in first_scores_recs_tags]
    thresh = np.average(scores) - 2 * np.std(scores)
    log.info(f'Score threshold: {thresh:.2f}')

    total_out = 0
    with open(sans_bc_fq_fpath, 'w') as bc_fq_fh, \
            open(sans_bc_paired_fq_fpath, 'w') as paired_fq_fh, \
            open(tags_fpath, 'w') as tag_fh:
        log.info('Continuing...')
        for i, (bc_rec, paired_rec) in enumerate(zip(
            SeqIO.parse(misc.gzip_friendly_open(bc_fq_fpath), 'fastq'),
            SeqIO.parse(misc.gzip_friendly_open(paired_fq_fpath), 'fastq'))):
            if i % 100000 == 0 and i > 0:
                log.info(f'  {i:,d}')
            score, sans_bc_rec, tags = process_bc_rec(bc_rec, aligners, bcd)
            if score >= thresh and sans_bc_rec:
                total_out += 1
                SeqIO.write(sans_bc_rec, bc_fq_fh, 'fastq')
                SeqIO.write(paired_rec, paired_fq_fh, 'fastq')
                output_rec_name_and_tags(sans_bc_rec, tags, tag_fh)
    log.info(f'{i+1:,d} barcode records processed')
    log.info(f'{total_out:,d} pairs of records output')


def write_chunk(arguments, tmpdirname, template_bam_fpath, i, bc_chunk, p_chunk):
    """
    Writes chunks to files
    """
    tmp_bc_fq_fpath = os.path.join(tmpdirname, f'{i}_bc.fq')
    tmp_paired_fq_fpath = os.path.join(tmpdirname, f'{i}_paired.fq')
    tmp_out_bc_fq_fpath = os.path.join(tmpdirname, f'{i}_bc.processed.fq')
    tmp_out_paired_fq_fpath = os.path.join(tmpdirname, f'{i}_paired.processed.fq')
    tmp_out_tags_fpath = os.path.join(tmpdirname, f'{i}_tags.fq')
    with open(tmp_bc_fq_fpath, 'w') as fq_out:
        SeqIO.write(bc_chunk, fq_out, 'fastq')
    with open(tmp_paired_fq_fpath, 'w') as fq_out:
        SeqIO.write(p_chunk, fq_out, 'fastq')
    return (tmp_bc_fq_fpath,
            tmp_paired_fq_fpath,
            tmp_out_bc_fq_fpath,
            tmp_out_paired_fq_fpath,
            tmp_out_tags_fpath,
            template_bam_fpath)


def chunked_gDNA_paired_recs_tmp_files_iterator(arguments, thresh, bc_fq_fpath, paired_fq_fpath, tmpdirname, chunksize):
    """
    Breaks pairs into chunks and writes to files.
    """
    bc_chunk, p_chunk = [], []
    for i, (bc_rec, paired_rec) in enumerate(zip(
        SeqIO.parse(misc.gzip_friendly_open(bc_fq_fpath), 'fastq'),
        SeqIO.parse(misc.gzip_friendly_open(paired_fq_fpath), 'fastq'))):
        bc_chunk.append(bc_rec)
        p_chunk.append(paired_rec)
        if i % chunksize == 0 and i > 0:
            yield arguments, thresh, write_chunk(arguments, tmpdirname, paired_fq_fpath, i, bc_chunk, p_chunk)
            bc_chunk, p_chunk = [], []
    if i % chunksize:
        yield arguments, thresh, write_chunk(arguments, tmpdirname, paired_fq_fpath, i, bc_chunk, p_chunk)


def process_chunk_of_reads(args_and_fpaths):
    """
    Processing chunks of reads. Required to build aligners in each parallel process.
    """
    arguments, thresh, (tmp_bc_fq_fpath,
            tmp_paired_fq_fpath,
            tmp_out_bc_fq_fpath,
            tmp_out_paired_fq_fpath,
            tmp_out_tags_fpath,
            template_bam_fpath) = args_and_fpaths
    aligners = misc.build_gDNA_bc_aligners()
    bcd = BCDecoder(arguments.barcode_whitelist, arguments.max_bc_err_decode)
    with open(tmp_out_bc_fq_fpath, 'w') as bc_fq_fh, \
            open(tmp_out_paired_fq_fpath, 'w') as paired_fq_fh, \
            open(tmp_out_tags_fpath, 'w') as tag_fh:
        for i, (bc_rec, paired_rec) in enumerate(zip(
            SeqIO.parse(misc.gzip_friendly_open(tmp_bc_fq_fpath), 'fastq'),
            SeqIO.parse(misc.gzip_friendly_open(tmp_paired_fq_fpath), 'fastq'))):
            score, sans_bc_rec, tags = process_bc_rec(bc_rec, aligners, bcd)
            if score >= thresh and sans_bc_rec:
                SeqIO.write(sans_bc_rec, bc_fq_fh, 'fastq')
                SeqIO.write(paired_rec, paired_fq_fh, 'fastq')
                output_rec_name_and_tags(sans_bc_rec, tags, tag_fh)
    os.remove(tmp_bc_fq_fpath)
    os.remove(tmp_paired_fq_fpath)
    return tmp_out_bc_fq_fpath, tmp_out_paired_fq_fpath, tmp_out_tags_fpath

def parallel_process_gDNA_fastqs(arguments, bc_fq_fpath, paired_fq_fpath, sans_bc_fq_fpath, sans_bc_paired_fq_fpath, tags_fpath):
    """
    Parallel version of serial process.
    """
    chunksize=100000
    aligners = misc.build_gDNA_bc_aligners()
    bcd = BCDecoder(arguments.barcode_whitelist, arguments.max_bc_err_decode)
    with Pool(arguments.threads) as pool, \
            tempfile.TemporaryDirectory(prefix='/dev/shm/') as tmpdirname:
        log.info(f'Processing first {n_first_seqs:,d} for score threshold...')
        first_scores_recs_tags = []
        for i, bc_rec in enumerate(SeqIO.parse(misc.gzip_friendly_open(bc_fq_fpath), 'fastq')):
            first_scores_recs_tags.append(process_bc_rec(bc_rec, aligners, bcd))
            if i >= n_first_seqs:
                break

        scores = [score for score, rec, tags in first_scores_recs_tags]
        thresh = np.average(scores) - 2 * np.std(scores)
        log.info(f'Score threshold: {thresh:.2f}')

        log.info(f'Using temporary directory {tmpdirname}')
        total_out = 0

        # See RNAcount for explanation of architecture choice. Could possibly be restructured here
        # in a more normal fashion due to not passing around pysam objects. That would speed things
        # up moderately, but since the expected gain is low we keep the same architecture here.
        chunk_iter = chunked_gDNA_paired_recs_tmp_files_iterator(
                arguments,
                thresh,
                bc_fq_fpath,
                paired_fq_fpath,
                tmpdirname,
                chunksize=chunksize)

        with open(sans_bc_fq_fpath, 'w') as sans_bc_fh, \
                open(sans_bc_paired_fq_fpath, 'w') as sans_bc_paired_fh, \
                open(tags_fpath, 'w') as tags_fh:
            i = 0
            while True:
                chunk_of_args_and_fpaths = list(itertools.islice(chunk_iter, arguments.threads))
                if not chunk_of_args_and_fpaths:
                    break
                for j, (tmp_out_bc_fq_fpath, tmp_out_paired_fq_fpath, tmp_out_tags_fpath) in enumerate(pool.imap(
                    process_chunk_of_reads,
                    chunk_of_args_and_fpaths)):
                    it_idx = i*arguments.threads+j
                    log.info(f'  {it_idx*chunksize:,d}-{(it_idx+1)*chunksize:,d}')
                    SeqIO.write(SeqIO.parse(tmp_out_bc_fq_fpath, 'fastq'), sans_bc_fh, 'fastq')
                    SeqIO.write(SeqIO.parse(tmp_out_paired_fq_fpath, 'fastq'), sans_bc_paired_fh, 'fastq')
                    for line in open(tmp_out_tags_fpath):
                        tags_fh.write(line)
                        total_out += 1
                    for fpath in (tmp_out_bc_fq_fpath, tmp_out_paired_fq_fpath, tmp_out_tags_fpath):
                        os.remove(fpath)
                i += 1
    
    nrecs = int(misc.file_prefix_from_fpath(tmp_out_bc_fq_fpath).split('_')[0]) 
    log.info(f'{nrecs:,d} barcode records processed')
    log.info(f'{total_out:,d} pairs of records output')


def build_complete_bc(read):
    return read.get_tag('CB')

def count_parallel_wrapper(ref_and_input_bam_fpath):
    ref, input_bam_fpath = ref_and_input_bam_fpath
    read_count_given_bc_then_feature = defaultdict(Counter)
    for read in pysam.AlignmentFile(input_bam_fpath).fetch(ref):
        if read.is_read1 or (read.is_read2 and read.mate_is_unmapped):
            for gx_gn_tup in misc.gx_gn_tups_from_read(read): # count read toward all compatible genes
                read_count_given_bc_then_feature[build_complete_bc(read)][gx_gn_tup] += 1
    return ref, read_count_given_bc_then_feature

def gDNA_count_matrix(arguments, input_bam_fpath):
    """
    Counts the reads from the input bam file and outputs a sparse matrix of read counts.
    """
    raw_reads_output_dir = os.path.join(arguments.output_dir, 'raw_reads_bc_matrix')
    if os.path.exists(raw_reads_output_dir):
        log.info('Matrix output folder exists. Skipping count matrix build')
        return 
    else:
        os.makedirs(raw_reads_output_dir)

    log.info('Finding all barcodes present...')
    sorted_complete_bcs, sorted_features = misc.get_bcs_and_features_from_bam(
            input_bam_fpath,
            build_complete_bc,
            arguments.threads-1
            )
    i_given_feature = {feat: i for i, feat in enumerate(sorted_features)}
    j_given_complete_bc = {comp_bc: j for j, comp_bc in enumerate(sorted_complete_bcs)}


    with pysam.AlignmentFile(input_bam_fpath) as bamfile:
        reference_names = bamfile.references
    reference_names_with_input_bam = [(ref, input_bam_fpath) for ref in reference_names]

    log.info('Counting reads...')
    M_reads = lil_matrix((len(sorted_features), len(sorted_complete_bcs)), dtype=int)
    with Pool(arguments.threads) as pool:
        for ref, read_count_given_bc_then_feature in pool.imap_unordered(
                count_parallel_wrapper,
                reference_names_with_input_bam):
            for comp_bc, read_count_given_feature in read_count_given_bc_then_feature.items():
                j = j_given_complete_bc[comp_bc]
                for gx_gn_tup, read_count in read_count_given_feature.items():
                    i = i_given_feature[gx_gn_tup]
                    M_reads[i, j] += read_count

    log.info('Writing raw read count matrix...')
    misc.write_matrix(M_reads, sorted_complete_bcs, sorted_features, raw_reads_output_dir)
