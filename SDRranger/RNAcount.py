import logging
import os
import pysam
import itertools
import tempfile
import numpy as np
import subprocess
import shutil
from Bio import SeqIO
from multiprocessing import Pool
from .bc_aligner import CustomBCAligner
from .bc_decoders import BCDecoder, SBCDecoder
from .misc import gzip_friendly_open, names_pair, find_paired_fastqs_in_dir, file_prefix_from_fpath
from .constants import commonseq1_options, commonseq2_RNA 
from .umi import get_umi_maps_from_bam_file


log = logging.getLogger(__name__)
pysam.set_verbosity(0)


def process_RNA_fastqs(arguments):
    """
    Output single file with parsed bcs from bc_fastq in read names and seqs from paired_fastq.
    """
    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)
    paired_fpaths = find_paired_fastqs_in_dir(arguments.fastq_dir)
    log.info('Files to process:')
    for fpath1, fpath2 in paired_fpaths:
        log.info(f'  {fpath1}')
        log.info(f'  - {fpath2}')
    bc_fq_idx, paired_fq_idx = determine_bc_and_paired_fastq_idxs(paired_fpaths) # determine which paired end has bcs
    log.info(f'Detected barcodes in read{bc_fq_idx+1} files')
    log.info(f'Running STAR alignment...')
    paired_fq_bam_fpaths = []
    for tup_fastq_fpaths in paired_fpaths:
        bc_fq_fpath = tup_fastq_fpaths[bc_fq_idx]
        paired_fq_fpath = tup_fastq_fpaths[paired_fq_idx]
        log.info(f'  {paired_fq_fpath}')
        star_out_dir, star_out_fpath = run_STAR_RNA(arguments, paired_fq_fpath)
        paired_fq_bam_fpaths.append((bc_fq_fpath, star_out_fpath))

    star_w_bc_fname = 'RNA_with_bc.bam'
    star_w_bc_fpath = os.path.join(arguments.output_dir, star_w_bc_fname)
    log.info('Writing output to:')
    log.info(f'  {star_w_bc_fpath}')

    template_bam = paired_fq_bam_fpaths[0][1]
    process_fastqs_func = serial_process_RNA_fastqs if arguments.threads == 1 else parallel_process_RNA_fastqs
    with pysam.AlignmentFile(star_w_bc_fpath, 'wb', template=pysam.AlignmentFile(template_bam)) as star_w_bc_fh:
        for bc_fq_fpath, star_raw_fpath in paired_fq_bam_fpaths:
            log.info('Processing files:')
            log.info(f'  barcode fastq: {bc_fq_fpath}')
            log.info(f'  paired bam:    {star_raw_fpath}')
            process_fastqs_func(arguments, bc_fq_fpath, star_raw_fpath, star_w_bc_fh)
    shutil.rmtree(star_out_dir)  # clean up intermediate STAR files

    star_w_bc_sorted_fname = 'RNA_with_bc.sorted.bam'
    star_w_bc_sorted_fpath = os.path.join(arguments.output_dir, star_w_bc_sorted_fname)
    log.info('Sorting bam...')
    pysam.sort('-o', star_w_bc_sorted_fpath, star_w_bc_fpath)
    os.remove(star_w_bc_fpath)  #clean up unsorted bam
    log.info('Indexing bam...')
    pysam.index(star_w_bc_sorted_fpath)

    log.info('Correcting UMIs...')
    star_w_bc_umi_sorted_fname = 'RNA_with_bc_umi.sorted.bam'
    star_w_bc_umi_sorted_fpath = os.path.join(arguments.output_dir, star_w_bc_umi_sorted_fname)
    correct_UMIs(star_w_bc_sorted_fpath, star_w_bc_umi_sorted_fpath)
    log.info('Indexing bam...')
    pysam.index(star_w_bc_umi_sorted_fpath)
    os.remove(star_w_bc_sorted_fpath)
    os.remove(star_w_bc_sorted_fpath + '.bai')
    log.info('Done')


def run_STAR_RNA(arguments, fastq_fpath):
    """
    Run STAR aligner for RNA files.

    Returns STAR output directory and bam path.
    """
    star_out_dir = os.path.join(arguments.output_dir, 'STAR_files')
    fastq_bname = file_prefix_from_fpath(fastq_fpath)
    out_prefix = os.path.join(star_out_dir, f'{fastq_bname}_')
    cmd_star = [
        'STAR',
        f'--runThreadN 1', # required to keep order matching with fastq file
        f'--genomeDir {arguments.star_ref_dir}',
        f'--readFilesIn {fastq_fpath}',
        f'--outFileNamePrefix {out_prefix}',
        '--outFilterMultimapNmax 1', 
        '--outSAMtype BAM Unsorted', 
    ]
    if fastq_fpath.endswith('gz'):
        cmd_star.append('--readFilesCommand zcat')
    subprocess.run(cmd_star, check=True)
    star_out_fpath = f'{out_prefix}Aligned.out.bam'
    return star_out_dir, star_out_fpath


def process_bc_rec_and_p_read(bc_rec, p_read, aligners, bcd, sbcd):
    """
    Find barcodes etc in bc_rec and add them as tags to p_read
    """
    bc_seq = str(bc_rec.seq)
    scores_and_pieces = [al.find_norm_score_and_pieces(bc_seq) for al in aligners]
    raw_score, raw_pieces = max(scores_and_pieces)
    raw_bc1, raw_bc2, raw_sbc = [raw_pieces[i] for i in [0, 2, 4]]

    bc1 = bcd.decode(raw_bc1)
    bc2 = bcd.decode(raw_bc2)
    sbc = sbcd.decode(raw_sbc)
    best_aligner = next(al for al, (s, p) in zip(aligners, scores_and_pieces) if s == raw_score)
    commonseq1, commonseq2 = [best_aligner.prefixes[i] for i in [1, 3]]
    corrected_pieces = [bc1, commonseq1, bc2, commonseq2, sbc]
    if None in corrected_pieces:
        return raw_score, None

    new_aligner = CustomBCAligner(*corrected_pieces, 'N'*8)
    new_score, new_pieces = new_aligner.find_norm_score_and_pieces(bc_seq)

    # Add tags for corrected and raw:
    # Cell barcode
    p_read.set_tag('CB', f'{bc1}.{bc2}')
    p_read.set_tag('CR', f'{raw_pieces[0]}.{raw_pieces[2]}')
    # Sample barcode
    p_read.set_tag('SB', sbc)
    p_read.set_tag('SR', raw_pieces[4])
    # Filler sequences
    p_read.set_tag('FB', f'{commonseq1}.{commonseq2}')
    p_read.set_tag('FR', f'{raw_pieces[1]}.{raw_pieces[3]}')
    # And raw UMI
    p_read.set_tag('UR', new_pieces[-1])
    return raw_score, p_read


def RNA_paired_recs_iterator(bc_fq_fpath, p_bam_fpath):
    """
    Iterates bc fastq reads with matching paired bam records.
    """
    bc_fq_iter = iter(SeqIO.parse(gzip_friendly_open(bc_fq_fpath), 'fastq'))
    for p_read in pysam.AlignmentFile(p_bam_fpath).fetch(until_eof=True):
        bc_rec = next(rec for rec in bc_fq_iter if names_pair(str(rec.id), str(p_read.qname)))
        yield bc_rec, p_read


def determine_bc_and_paired_fastq_idxs(paired_fpaths):
    """Returns (bc_fq_idx, paired_fq_idx), based on first pair of fastqs"""
    fpath0, fpath1 = paired_fpaths[0]
    avg_score0 = average_align_score_of_first_recs(fpath0)
    avg_score1 = average_align_score_of_first_recs(fpath1)
    return (0, 1) if avg_score0 > avg_score1 else (1, 0)  

def build_RNA_bc_aligners():
    """Helper function, kept as function for parallelism"""
    return [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_RNA, 'N'*8, 'N'*8) for cso in commonseq1_options]

def average_align_score_of_first_recs(fastq_fpath, n_first_seqs=500):
    """Return average alignment score of first n_first_seqs records"""
    aligners = build_RNA_bc_aligners()
    first_scores = []
    for i, rec in enumerate(SeqIO.parse(gzip_friendly_open(fastq_fpath), 'fastq')):
        scores_and_pieces = [al.find_norm_score_and_pieces(str(rec.seq)) for al in aligners]
        raw_score, raw_pieces = max(scores_and_pieces)
        first_scores.append(raw_score)
        if i >= n_first_seqs:
            break
    return np.average(first_scores)


def serial_process_RNA_fastqs(arguments, bc_fq_fpath, star_raw_fpath, star_w_bc_fh):
    log.info('Building aligners and barcode decoders')
    aligners = build_RNA_bc_aligners()
    bcd = BCDecoder(arguments.barcode_whitelist, arguments.max_bc_err_decode)
    sbcd = SBCDecoder(arguments.sample_barcode_whitelist, arguments.max_sbc_err_decode, arguments.sbc_reject_delta)

    n_first_seqs = 1000
    log.info(f'Processing first {n_first_seqs:,d} for score threshold...')
    first_scores_and_reads = []
    for i, (bc_rec, p_read) in enumerate(
            RNA_paired_recs_iterator(bc_fq_fpath, star_raw_fpath)
            ):
        first_scores_and_reads.append(process_bc_rec_and_p_read(bc_rec, p_read, aligners, bcd, sbcd))
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
    for i, (bc_rec, p_read) in enumerate(
            RNA_paired_recs_iterator(bc_fq_fpath, star_raw_fpath)
            ):
        if i <= n_first_seqs:
            continue
        if i % 10000 == 0 and i > 0:
            log.info(f'  {i:,d}')
        score, read = process_bc_rec_and_p_read(bc_rec, p_read, aligners, bcd, sbcd)
        if score >= thresh and read:
            total_out += 1
            star_w_bc_fh.write(read)
    log.info(f'{i+1:,d} records processed')
    log.info(f'{total_out:,d} records output')


def write_chunk(arguments, tmpdirname, template_bam_fpath, i, bc_chunk, p_chunk):
    """
    Writes chunks to files
    """
    tmp_fq_fpath = os.path.join(tmpdirname, f'{i}.fq')
    tmp_bam_fpath = os.path.join(tmpdirname, f'{i}.bam')
    tmp_out_bam_fpath = os.path.join(tmpdirname, f'{i}.parsed.bam')
    with open(tmp_fq_fpath, 'w') as fq_out:
        SeqIO.write(bc_chunk, fq_out, 'fastq')
    with pysam.AlignmentFile(tmp_bam_fpath, 'wb', template=pysam.AlignmentFile(template_bam_fpath)) as bam_out:
        for p_read in p_chunk:
            bam_out.write(p_read)
    return tmp_fq_fpath, tmp_bam_fpath, tmp_out_bam_fpath, template_bam_fpath


def chunked_RNA_paired_recs_tmp_files_iterator(arguments, thresh, bc_fq_fpath, p_bam_fpath, tmpdirname, chunksize):
    """
    Breaks pairs into chunks and writes to files.
    """
    bc_chunk, p_chunk = [], []
    for i, (bc_rec, p_read) in enumerate(RNA_paired_recs_iterator(bc_fq_fpath, p_bam_fpath)):
        bc_chunk.append(bc_rec)
        p_chunk.append(p_read)
        if i % chunksize == 0 and i > 0:
            yield arguments, thresh, write_chunk(arguments, tmpdirname, p_bam_fpath, i, bc_chunk, p_chunk)
            bc_chunk, p_chunk = [], []
    if i % chunksize:
        yield arguments, thresh, write_chunk(arguments, tmpdirname, p_bam_fpath, i, bc_chunk, p_chunk)


def process_chunk_of_reads(args_and_fpaths):
    """
    Processing chunks of reads. Required to build aligners in each parallel process.
    """
    arguments, thresh, (tmp_fq_fpath, tmp_bam_fpath, tmp_out_bam_fpath, template_bam_fpath) = args_and_fpaths
    aligners = build_RNA_bc_aligners()
    bcd = BCDecoder(arguments.barcode_whitelist, arguments.max_bc_err_decode)
    sbcd = SBCDecoder(arguments.sample_barcode_whitelist, arguments.max_sbc_err_decode, arguments.sbc_reject_delta)
    with pysam.AlignmentFile(tmp_out_bam_fpath, 'wb', template=pysam.AlignmentFile(template_bam_fpath)) as out:
        for bc_rec, p_read in RNA_paired_recs_iterator(tmp_fq_fpath, tmp_bam_fpath):
            score, read = process_bc_rec_and_p_read(bc_rec, p_read, aligners, bcd, sbcd)
            if score >= thresh and read:
                out.write(read)
    os.remove(tmp_fq_fpath)
    os.remove(tmp_bam_fpath)
    return tmp_out_bam_fpath


def parallel_process_RNA_fastqs(arguments, bc_fq_fpath, star_raw_fpath, star_w_bc_fh):
    """
    Parallel version of serial process.
    
    Rather more involved. pysam doesn't parallelize well. AlignedSegment's don't pickle. So one
    must create a large number of temporary files and process things that way, with multiple levels
    of helper functions
    """
    n_first_seqs = 1000
    chunksize=2000
    aligners = build_RNA_bc_aligners()
    bcd = BCDecoder(arguments.barcode_whitelist, arguments.max_bc_err_decode)
    sbcd = SBCDecoder(arguments.sample_barcode_whitelist, arguments.max_sbc_err_decode, arguments.sbc_reject_delta)
    with Pool(arguments.threads) as pool, \
            tempfile.TemporaryDirectory(prefix='/dev/shm/') as tmpdirname:
        log.info(f'Processing first {n_first_seqs:,d} for score threshold...')
        first_scores_and_reads = []
        for i, (bc_rec, p_read) in enumerate(
                RNA_paired_recs_iterator(bc_fq_fpath, star_raw_fpath)
                ):
            first_scores_and_reads.append(process_bc_rec_and_p_read(bc_rec, p_read, aligners, bcd, sbcd))
            if i >= n_first_seqs:
                break

        scores = [score for score, read in first_scores_and_reads]
        thresh = np.average(scores) - 2 * np.std(scores)
        log.info(f'Score threshold: {thresh:.2f}')

        log.info(f'Using temporary directory {tmpdirname}')
        total_out = 0
        chunk_iter = chunked_RNA_paired_recs_tmp_files_iterator(
                arguments,
                thresh,
                bc_fq_fpath,
                star_raw_fpath,
                tmpdirname,
                chunksize=chunksize)
        i = 0
        while True:
            chunk_of_args_and_fpaths = list(itertools.islice(chunk_iter, arguments.threads))
            if not chunk_of_args_and_fpaths:
                break
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
    
    nrecs = int(file_prefix_from_fpath(tmp_out_bam_fpath).split('.')[0])
    log.info(f'{nrecs:,d} records processed')
    log.info(f'{total_out:,d} records output')

def correct_UMIs(input_bam_fpath, out_bam_fpath):
    with pysam.AlignmentFile(input_bam_fpath) as bamfile:
        reference_names = bamfile.references

    with pysam.AlignmentFile(out_bam_fpath, 'wb', template=pysam.AlignmentFile(input_bam_fpath)) as bam_out:
        for ref in reference_names:
            log.info(f'  {ref}')
            umi_map_given_bc = get_umi_maps_from_bam_file(input_bam_fpath, chrm=ref)
            for read in pysam.AlignmentFile(input_bam_fpath).fetch(ref):
                corrected_umi = umi_map_given_bc[read.get_tag('CB')][read.get_tag('UR')]
                read.set_tag('UB', corrected_umi)
                bam_out.write(read)
