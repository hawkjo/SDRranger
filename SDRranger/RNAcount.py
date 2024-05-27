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
from functools import partial
from multiprocessing import Pool
from scipy.sparse import lil_matrix
from .bc_aligner import CustomBCAligner
from .umi import get_umi_maps_from_bam_file


log = logging.getLogger(__name__)
pysam.set_verbosity(0)


n_first_seqs = 10000  # n seqs for finding score threshold


def process_RNA_fastqs(arguments):
    """
    Output single file with parsed bcs from bc_fastq in read names and seqs from paired_fastq.
    """
    star_w_bc_fpath = os.path.join(arguments.output_dir, 'RNA_with_bc.bam')
    star_w_bc_sorted_fpath = os.path.join(arguments.output_dir, 'RNA_with_bc.sorted.bam')
    star_w_bc_umi_fpath = os.path.join(arguments.output_dir, 'RNA_with_bc_umi.bam')
    star_w_bc_umi_sorted_fpath = os.path.join(arguments.output_dir, 'RNA_with_bc_umi.sorted.bam')
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

    paired_fpaths = misc.find_paired_fastqs_in_dir(arguments.fastq_dir)
    log.info('Files to process:')
    for i, (fpath1, fpath2) in enumerate(paired_fpaths):
        log.info(f'  {fpath1}')
        log.info(f'  {fpath2}')
        if i < len(paired_fpaths)-1:
            log.info('  -')

    if completed < 1:
        bc_fq_idx, paired_fq_idx = misc.determine_bc_and_paired_fastq_idxs(paired_fpaths, arguments.config)
        log.info(f'Detected barcodes in read{bc_fq_idx+1} files')

        paired_fq_bam_fpaths = []
        namepairidxs = []
        if not arguments.star_output_path:
            log.info(f'Running STAR alignment...')
            for tup_fastq_fpaths in paired_fpaths:
                bc_fq_fpath = tup_fastq_fpaths[bc_fq_idx]
                paired_fq_fpath = tup_fastq_fpaths[paired_fq_idx]

                namepairidxs.append(misc.get_namepair_index(bc_fq_fpath, paired_fq_fpath))

                log.info(f'  {paired_fq_fpath}')
                star_out_dir, star_out_fpath = run_STAR_RNA(arguments, paired_fq_fpath)
                paired_fq_bam_fpaths.append((bc_fq_fpath, star_out_fpath))
        else:
            log.info(f'Using STAR results from {arguments.star_output_path}')
            for fpaths, bampath in zip(paired_fpaths, arguments.star_output_path, strict=True):
                paired_fq_bam_fpaths.append((fpaths[bc_fq_idx], bampath))
    
        log.info('Writing output to:')
        log.info(f'  {star_w_bc_fpath}')
    
        template_bam = paired_fq_bam_fpaths[0][1]
        process_fastqs_func = partial(serial_process_RNA_fastqs, sameorder=not arguments.star_output_path) if arguments.threads == 1 else parallel_process_RNA_fastqs
        with pysam.AlignmentFile(star_w_bc_fpath, 'wb', template=pysam.AlignmentFile(template_bam)) as star_w_bc_fh:
            for (bc_fq_fpath, star_raw_fpath), namepairidx in zip(paired_fq_bam_fpaths, namepairidxs):
                log.info('Processing files:')
                log.info(f'  barcode fastq: {bc_fq_fpath}')
                log.info(f'  paired bam:    {star_raw_fpath}')
                process_fastqs_func(arguments, bc_fq_fpath, star_raw_fpath, star_w_bc_fh, namepairidx)
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

    RNA_count_matrix(arguments, star_w_bc_umi_sorted_fpath)
    log.info('Done')


def run_STAR_RNA(arguments, fastq_fpath):
    """
    Run STAR aligner for RNA files.

    Returns STAR output directory and bam path.
    """
    star_out_dir = os.path.join(arguments.output_dir, 'STAR_files')
    fastq_bname = misc.file_prefix_from_fpath(fastq_fpath)
    out_prefix = os.path.join(star_out_dir, f'{fastq_bname}_')
    cmd_star = [
        'STAR',
        f'--runThreadN {arguments.threads}',
        f'--genomeDir {arguments.star_ref_dir}',
        f'--readFilesIn {fastq_fpath}',
        f'--outFileNamePrefix {out_prefix}',
        '--outFilterMultimapNmax 1',
        '--outSAMtype BAM Unsorted',
        '--outSAMattributes NH HI AS nM GX GN',
    ]
    if fastq_fpath.endswith('gz'):
        cmd_star.append('--readFilesCommand zcat')
    star_out_fpath = f'{out_prefix}Aligned.out.bam'
    if os.path.exists(star_out_fpath):
        log.info("STAR results found. Skipping alignment")
    else:
        subprocess.run(cmd_star, check=True)
    return star_out_dir, star_out_fpath


def process_bc_rec_and_p_read(config, bc_rec, p_read, aligners, decoders):
    """
    Find barcodes etc in bc_rec and add them as tags to p_read
    """
    blocks = config["barcode_struct_r1"]["blocks"]
    scores_and_pieces = [al.find_norm_score_and_pieces(bc_rec.seq, return_seq=True) for al in aligners]
    raw_score, raw_pieces, raw_seq = max(scores_and_pieces)
    raw_bcs = [raw_piece.upper().replace('N', 'A') for raw_piece, block in zip(raw_pieces, blocks) if block["blocktype"] == "barcodeList"]
    bcs = [decoder.decode(raw_bc) for raw_bc, decoder in zip(raw_bcs, decoders)]

    best_aligner = next(al for al, (s, p, sq) in zip(aligners, scores_and_pieces) if s == raw_score)
    commonseqs = [best_aligner.prefixes[i] for i, block in enumerate(blocks) if block["blocktype"] == "constantRegion"]

    bc_idx = commonseq_idx = 0
    corrected_pieces = []
    for block in blocks:
        if block["blocktype"] == "barcodeList":
            if bcs[bc_idx] is None:
                return raw_score, None
            corrected_pieces.append(bcs[bc_idx])
            bc_idx += 1
        elif block["blocktype"] == "constantRegion":
            if commonseqs[commonseq_idx] is None:
                return raw_score, None
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
                if block["blocktype"] == "UMI":
                    bam_umis.append(new_pieces[i])
                else:
                    bam_bcs.append(new_pieces[i])
                    bam_raw_bcs.append(new_pieces[i])

    # Add tags for corrected and raw:
    # Cell barcode
    p_read.set_tag("CB", ".".join(bam_bcs))
    p_read.set_tag("CR", ".".join(bam_raw_bcs))
    # Filler sequences
    p_read.set_tag("FL", sum(len(seq) for seq in bam_commonseqs))
    # And raw UMI
    p_read.set_tag("UR", ".".join(bam_umis))
    return raw_score, p_read


def RNA_paired_recs_iterator(bc_fq_fpath, p_bam_fpath):
    """
    Iterates bc fastq reads with matching paired bam records.
    """
    bc_fq_iter = iter(SeqIO.parse(misc.gzip_friendly_open(bc_fq_fpath), 'fastq'))
    for p_read in pysam.AlignmentFile(p_bam_fpath).fetch(until_eof=True):
        bc_rec = next(rec for rec in bc_fq_iter if misc.names_pair(str(rec.id), str(p_read.qname)))
        yield bc_rec, p_read

def RNA_unpaired_recs_iterator(bc_fq_fpath, p_bam_fpath, bamidx):
    """
    Iterates bc fastq reads with matching bam records when order of reads is not the same.
    """
    bam = pysam.AlignmentFile(p_bam_fpath)
    for bc_rec in SeqIO.parse(misc.gzip_friendly_open(bc_fq_fpath), 'fastq'):
        for p_read in misc.get_bam_read_by_name(bc_rec.id, bam, bamidx):
            yield bc_rec, p_read


def serial_process_RNA_fastqs(arguments, bc_fq_fpath, star_raw_fpath, star_w_bc_fh, namepairidx, sameorder=True):
    log.info('Building aligners and barcode decoders')
    aligners = misc.build_bc_aligners(arguments.config)
    decoders = misc.build_bc_decoders(arguments.config)

    if not sameorder:
        star_readname_sorted_fpath = star_raw_fpath + "_readname_sorted.bam"
        bamidx = misc.sort_and_index_readname_bam(star_raw_fpath, star_readname_sorted_fpath, namepairidx, arguments.threads)
        iterator = RNA_unpaired_recs_iterator(bc_fq_fpath, star_readname_sorted_fpath, bamidx)
    else:
        iterator = RNA_paired_recs_iterator(bc_fq_fpath, star_raw_fpath)

    log.info(f'Processing first {n_first_seqs:,d} for score threshold...')
    first_scores_and_reads = []
    for i, (bc_rec, p_read) in enumerate(iterator):
        first_scores_and_reads.append(process_bc_rec_and_p_read(arguments.config, bc_rec, p_read, aligners, decoders))
        if i >= n_first_seqs:
            break

    scores = [score for score, read in first_scores_and_reads]
    thresh = np.average(scores) - 2 * np.std(scores)
    log.info(f'Score threshold: {thresh:.2f}')
    out_reads = [read for score, read in first_scores_and_reads if score >= thresh and read]
    total_out = len(out_reads)
    for read in out_reads:
        star_w_bc_fh.write(read)

    log.info('Continuing...')
    for i, (bc_rec, p_read) in enumerate(iterator):
        if i <= n_first_seqs:
            continue
        if i % 100000 == 0 and i > 0:
            log.info(f'  {i:,d}')
        score, read = process_bc_rec_and_p_read(arguments.config, bc_rec, p_read, aligners, decoders)
        if score >= thresh and read:
            total_out += 1
            star_w_bc_fh.write(read)

    if not sameorder:
        os.remove(star_readname_sorted_fpath)

    log.info(f'{i+1:,d} records processed')
    log.info(f'{total_out:,d} records output')


def write_chunk(tmpdirname, i, bc_chunk):
    """
    Writes chunks to files
    """
    tmp_fq_fpath = os.path.join(tmpdirname, f'{i}.fq')
    tmp_out_bam_fpath = os.path.join(tmpdirname, f'{i}.parsed.bam')
    with open(tmp_fq_fpath, 'w') as fq_out:
        SeqIO.write(bc_chunk, fq_out, 'fastq')
    return tmp_fq_fpath, tmp_out_bam_fpath


def chunked_RNA_recs_tmp_files_iterator(bc_fq_fpath, tmpdirname, chunksize):
    """
    Breaks pairs into chunks and writes to files.
    """
    bc_chunk = []
    for i, bc_req in enumerate(SeqIO.parse(misc.gzip_friendly_open(bc_fq_fpath), 'fastq')):
        bc_chunk.append(bc_req)
        if i % chunksize == 0 and i > 0:
            yield write_chunk(tmpdirname, i, bc_chunk)
            bc_chunk, p_chunk = [], []
    if i % chunksize:
        yield write_chunk(tmpdirname, i, bc_chunk)


def process_chunk_of_reads(args_and_fpaths):
    """
    Processing chunks of reads. Required to build aligners in each parallel process.
    """
    (tmp_fq_fpath, tmp_out_bam_fpath), (config, thresh, sorted_bam_fpath, sorted_bam_idx) = args_and_fpaths
    aligners = misc.build_bc_aligners(config)
    decoders = misc.build_bc_decoders(config)
    with pysam.AlignmentFile(tmp_out_bam_fpath, 'wb', template=pysam.AlignmentFile(sorted_bam_fpath)) as out:
        for bc_rec, p_read in RNA_unpaired_recs_iterator(tmp_fq_fpath, sorted_bam_fpath, sorted_bam_idx):
                score, read = process_bc_rec_and_p_read(config, bc_rec, p_read, aligners, decoders)
                if score >= thresh and read:
                    out.write(read)
    os.remove(tmp_fq_fpath)
    return tmp_out_bam_fpath


def parallel_process_RNA_fastqs(arguments, bc_fq_fpath, star_raw_fpath, star_w_bc_fh, namepairidx):
    """
    Parallel version of serial process.

    Rather more involved. pysam doesn't parallelize well. AlignedSegment's don't pickle. So one
    must create a large number of temporary files and process things that way, with multiple levels
    of helper functions
    """
    chunksize=100000
    aligners = misc.build_bc_aligners(arguments.config)
    decoders = misc.build_bc_decoders(arguments.config)
    with Pool(arguments.threads) as pool, \
            tempfile.TemporaryDirectory(prefix='/dev/shm/') as tmpdirname:
        log.info(f'Processing first {n_first_seqs:,d} for score threshold...')
        first_scores_and_reads = []
        for i, (bc_rec, p_read) in enumerate(
                RNA_paired_recs_iterator(bc_fq_fpath, star_raw_fpath)
                ):
            first_scores_and_reads.append(process_bc_rec_and_p_read(arguments.config, bc_rec, p_read, aligners, decoders))
            if i >= n_first_seqs:
                break

        scores = [score for score, read in first_scores_and_reads]
        thresh = np.average(scores) - 2 * np.std(scores)
        log.info(f'Score threshold: {thresh:.2f}')

        star_readname_sorted_fpath = star_raw_fpath + "_readname_sorted.bam"
        bamidx = misc.sort_and_index_readname_bam(star_raw_fpath, star_readname_sorted_fpath, namepairidx, arguments.threads)

        log.info(f'Using temporary directory {tmpdirname}')
        total_out = 0

        # The following iteration architecture is designed to overcome a few limitaitons.
        #
        # First: pysam AlignedSegment objects are not picklable, so they cannot be sent to parallel
        # processes directly. Hence, each chunk of reads must be written to an intermediate file
        # and processed as a file. The only objects passed around are arguments and filenames
        #
        # Second: pysam does not easily append to an existing bam file, so we pass a filehandle
        # opened in the parent function that does not close until the full, combined bam file is
        # produced.
        #
        # Third: We use the imap method of the Pool object from the multiprocessing library.
        # The original intention was to hand it an iterator that generated the intermediate files
        # as they were lazily loaded by the imap function. However, it turns out the imap function
        # actually does construct all the objects of the iterator first, generating all files
        # instead of just the currently needed ones. So we use a construction suggested on
        # stackoverflow to created it first as an iterator, use itertools.islice to break off a
        # chunk of arguments at a time (generating only those intermediate files), and then call
        # imap separately for each chunk thus created. This works as desired, though obviously
        # sacrifices a bit in performance as the chunks are hard-separated before parallelism.
        # However, for this problem all pieces in each chunk take very similar amount of time to
        # process, so the performance sacrifice is not too significant.
        chunk_iter = chunked_RNA_recs_tmp_files_iterator(
                bc_fq_fpath,
                tmpdirname,
                chunksize=chunksize)
        i = 0
        while True:
            chunk = list(itertools.islice(chunk_iter, arguments.threads))
            if not chunk:
                break
            chunk_of_args_and_fpaths = zip(chunk,
                                           itertools.repeat((arguments.config, thresh, star_readname_sorted_fpath, bamidx)))
            for j, tmp_out_bam_fpath in enumerate(pool.imap(
                process_chunk_of_reads,
                chunk_of_args_and_fpaths)):
                it_idx = i*arguments.threads+j
                log.info(f'  {it_idx*chunksize:,d}-{(it_idx+1)*chunksize:,d}')
                for read in pysam.AlignmentFile(tmp_out_bam_fpath).fetch(until_eof=True):
                    total_out += 1
                    star_w_bc_fh.write(read)
                os.remove(tmp_out_bam_fpath)
            i += 1

    os.remove(star_readname_sorted_fpath)
    nrecs = int(misc.file_prefix_from_fpath(tmp_out_bam_fpath).split('.')[0])
    log.info(f'{nrecs:,d} records processed')
    log.info(f'{total_out:,d} records output')


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


def count_parallel_wrapper(ref_and_input_bam_fpath):
    ref, input_bam_fpath = ref_and_input_bam_fpath
    read_count_given_bc_then_feature_then_umi = defaultdict(lambda : defaultdict(Counter))
    for read in pysam.AlignmentFile(input_bam_fpath).fetch(ref):
        for gx_gn_tup in misc.gx_gn_tups_from_read(read): # count read toward all compatible genes
            read_count_given_bc_then_feature_then_umi[read.get_tag('CB')][gx_gn_tup][read.get_tag('UB')] += 1
    for k in read_count_given_bc_then_feature_then_umi.keys():
        read_count_given_bc_then_feature_then_umi[k] = dict(read_count_given_bc_then_feature_then_umi[k])
    read_count_given_bc_then_feature_then_umi = dict(read_count_given_bc_then_feature_then_umi)
    return ref, read_count_given_bc_then_feature_then_umi

def RNA_count_matrix(arguments, input_bam_fpath):
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
