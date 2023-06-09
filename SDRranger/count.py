import logging
import warnings
import os
import pysam
import numpy as np
from Bio import SeqIO
from collections import Counter, defaultdict
from multiprocessing import Pool
from .bc_aligner import CustomBCAligner
from .bc_decoders import BCDecoder, SBCDecoder
from .misc import gzip_friendly_open

log = logging.getLogger(__name__)
warnings.filterwarnings('ignore', '.*SQ')


commonseq1_options = ['GTCAGTACGTACGAGTC'[i:] for i in range(4)]
commonseq2_RNA = 'GTACTCGCAGTAGTCGACACGTC'
commonseq2_gDNA = 'GTACTCGCAGTAGTC'

def process_RNA_fastqs(arguments):
    """
    Output single file with parsed bcs from bc_fastq in read names and seqs from paired_fastq.
    """
    paired_fpath = arguments.paired_fastq_file
    paired_bname = os.path.splitext(paired_fpath[:-3] if paired_fpath.endswith('.gz') else paired_fpath)[0]
    out_fname = f'{os.path.basename(paired_bname)}_parsed.bam'
    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)
    out_fpath = os.path.join(arguments.output_dir, out_fname)

    log.info('Processing fastqs:')
    log.info(f'  barcode fastq: {arguments.bc_fastq_file}')
    log.info(f'  paired fastq:  {arguments.paired_fastq_file}')
    log.info('Writing output to:')
    log.info(f'  {out_fpath}')

    if arguments.threads == 1:
        serial_process_RNA_fastqs(arguments, out_fpath)
    else:
        parallel_process_RNA_fastqs(arguments, out_fpath)


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
    p_read.setTag('CB', bc1 + bc2)
    p_read.setTag('CR', raw_pieces[0] + raw_pieces[2])
    # Sample barcode
    p_read.setTag('SB', sbc)
    p_read.setTag('SR', raw_pieces[4])
    # Filler sequences
    p_read.setTag('FB', commonseq1 + commonseq2)
    p_read.setTag('FR', raw_pieces[1] + raw_pieces[3])
    # And raw UMI
    p_read.setTag('UR', new_pieces[-1])
    return raw_score, p_read


def RNA_paired_recs_iterator(bc_fastq_fpath, p_bam_fpath):
    """
    Iterates bc fastq reads with matching paired bam records.
    """
    def reads_match(r1, r2):
        return all(c1 == c2 or set([c1, c2]) == set('12') for c1, c2 in zip(r1, r2))
    bc_fq_iter = iter(SeqIO.parse(gzip_friendly_open(bc_fastq_fpath), 'fastq'))
    for p_read in pysam.AlignmentFile(p_bam_fpath).fetch(until_eof=True):
        bc_rec = next(rec for rec in bc_fq_iter if reads_match(str(rec.id), str(p_read.qname)))
        yield bc_rec, p_read


def chunked_RNA_paired_recs_iterator(arguments, bc_fastq_fpath, p_bam_fpath, chunksize=10000):
    """
    Breaks into pairs into chunks for parallelism
    """
    chunk = []
    for i, (bc_rec, p_read) in enumerate(RNA_paired_recs_iterator(bc_fastq_fpath, p_bam_fpath)):
        chunk.append((bc_rec, p_read.to_string()))
        if i % chunksize == 0 and i > 0:
            yield arguments, chunk
            chunk = []
    if i % chunksize:
        yield arguments, chunk


def serial_process_RNA_fastqs(arguments, out_fpath):
    log.info('Building aligners and barcode decoders')
    aligners = [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_RNA, 'N'*8, 'N'*8) for cso in commonseq1_options]
    bcd = BCDecoder(arguments.barcode_whitelist, arguments.max_bc_err_decode)
    sbcd = SBCDecoder(arguments.sample_barcode_whitelist, arguments.max_sbc_err_decode, arguments.sbc_reject_delta)

    n_first_seqs = 1000
    with pysam.AlignmentFile(out_fpath, 'wb', template=pysam.AlignmentFile(arguments.paired_fastq_file)) as out:
        log.info(f'Processing first {n_first_seqs:,d} for score threshold...')
        first_scores_and_reads = []
        for i, (bc_rec, p_read) in enumerate(
                RNA_paired_recs_iterator(arguments.bc_fastq_file, arguments.paired_fastq_file)
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
            out.write(read)

        log.info('Continuing...')
        for i, (bc_rec, p_read) in enumerate(
                RNA_paired_recs_iterator(arguments.bc_fastq_file, arguments.paired_fastq_file)
                ):
            if i <= n_first_seqs:
                continue
            if i % 10000 == 0 and i > 0:
                log.info(f'  {i:,d}')
            score, read = process_bc_rec_and_p_read(bc_rec, p_read, aligners, bcd, sbcd)
            if score >= thresh and read:
                total_out += 1
                out.write(read)
    log.info(f'{i+1:,d} records processed')
    log.info(f'{total_out:,d} records output')


def process_chunk_of_reads(args_and_chunk):
    """
    Processing chunks of reads. Required to build aligners in each parallel process.
    """
    arguments, chunk = args_and_chunk
    aligners = [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_RNA, 'N'*8, 'N'*8) for cso in commonseq1_options]
    bcd = BCDecoder(arguments.barcode_whitelist, arguments.max_bc_err_decode)
    sbcd = SBCDecoder(arguments.sample_barcode_whitelist, arguments.max_sbc_err_decode, arguments.sbc_reject_delta)
    header = pysam.AlignmentHeader()
    chunk = [(bc_rec, pysam.AlignedSegment.fromstring(p_read, header=header)) for bc_rec, p_read in chunk]
    return [process_bc_rec_and_p_read(bc_rec, p_read, aligners, bcd, sbcd) for bc_rec, p_read in chunk]


def parallel_process_RNA_fastqs(arguments, out_fpath):
    n_first_seqs = 4000
    chunksize=1000
    with pysam.AlignmentFile(out_fpath, 'wb', template=pysam.AlignmentFile(arguments.paired_fastq_file)) as out:
        with Pool(arguments.threads) as pool:
            log.info(f'Processing first {n_first_seqs:,d} for score threshold...')
            first_scores_and_reads = []
            for i, score_and_read_chunk in enumerate(pool.imap(
                    process_chunk_of_reads,
                    chunked_RNA_paired_recs_iterator(
                        arguments,
                        arguments.bc_fastq_file,
                        arguments.paired_fastq_file,
                        chunksize=chunksize),
                    )):
                first_scores_and_reads.extend(score_and_read_chunk)
                if (i+1)*chunksize >= n_first_seqs:
                    break
            scores = [score for score, read in first_scores_and_reads]
            thresh = np.average(scores) - 2 * np.std(scores)
            log.info(f'Score threshold: {thresh:.2f}')

            total_out = 0
            for i, score_and_read_chunk in enumerate(pool.imap(
                    process_chunk_of_reads,
                    chunked_RNA_paired_recs_iterator(
                        arguments,
                        arguments.bc_fastq_file,
                        arguments.paired_fastq_file,
                        chunksize=chunksize)
                    )):
                log.info(f'  {(i+1)*chunksize:,d}')
                for score, read in score_and_read_chunk:
                    if score >= thresh and read:
                        total_out += 1
                        out.write(read)

    log.info(f'{i*chunksize+len(score_and_read_chunk):,d} records processed')
    log.info(f'{total_out:,d} records output')


def process_gDNA_fastq(arguments):
    aligners = [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_gDNA) for cso in commonseq1_options]
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(gzip_friendly_open(bc_fastq_fpath), 'fastq'):
            seq = str(rec.seq)
            norm_score, pieces, end_pos = max([al.find_norm_score_pieces_and_end_pos(seq) for al in aligners])
            pieces_str = ','.join(pieces)
            name = f'{pieces_str}/{norm_score}/{rec.id}'
            rec = rec[end_pos:]
            rec.id = name
            rec.name = name
            rec.description = ''
            SeqIO.write(rec, out, 'fastq')


def norm_score_from_rec(rec):
    return float(str(rec.id).split('/')[1])


def bc_and_sbc_counter_from_fastq(bc_fastq_fpath, bc_parser):
    # Determine minimum norm_score threshold = mean-2sigma. See notebook for figures
    norm_scores_sample = []
    for i, rec in enumerate(SeqIO.parse(open(bc_fastq_fpath), 'fastq')):
        if i >= 100000:
            break
        norm_scores_sample.append(norm_score_from_rec(rec))
    thresh = np.average(norm_scores_sample) - 2 * np.std(norm_scores_sample)





def process_gDNA_fastqs(arguments):
    aligners = [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_gDNA) for cso in commonseq1_options]
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(gzip_friendly_open(bc_fastq_fpath), 'fastq'):
            seq = str(rec.seq)
            norm_score, pieces, end_pos = max([al.find_norm_score_pieces_and_end_pos(seq) for al in aligners])
            pieces_str = ','.join(pieces)
            name = f'{pieces_str}/{norm_score}/{rec.id}'
            rec = rec[end_pos:]
            rec.id = name
            rec.name = name
            rec.description = ''
            SeqIO.write(rec, out, 'fastq')


def norm_score_from_rec(rec):
    return float(str(rec.id).split('/')[1])


def bc_and_sbc_counter_from_fastq(bc_fastq_fpath, bc_parser):
    # Determine minimum norm_score threshold = mean-2sigma. See notebook for figures
    norm_scores_sample = []
    for i, rec in enumerate(SeqIO.parse(open(bc_fastq_fpath), 'fastq')):
        if i >= 100000:
            break
        norm_scores_sample.append(norm_score_from_rec(rec))
    thresh = np.average(norm_scores_sample) - 2 * np.std(norm_scores_sample)
    log.info(f'Norm score threshold: {thresh:.2f} for {bc_fastq_fpath}')

    bc_sbc_cntr = defaultdict(Counter)
    for rec in SeqIO.parse(open(bc_fastq_fpath), 'fastq'):
        norm_score = norm_score_from_rec(rec)
        if norm_score < thresh:
            continue
        bcs, sbc = bc_parser.bc_and_sbc_from_rec(rec)
        if None not in bcs and sbc is not None:
            bc_sbc_cntr[bcs][sbc] += 1
    return bc_sbc_cntr


