import os
import numpy as np
import logging
from Bio import SeqIO
from collections import Counter, defaultdict
from .bc_aligner import CustomBCAligner
from .bc_decoders import BCDecoder, SBCDecoder
from .misc import gzip_friendly_open

log = logging.getLogger(__name__)


commonseq1_options = ['GTCAGTACGTACGAGTC'[i:] for i in range(4)]
commonseq2_RNA = 'GTACTCGCAGTAGTCGACACGTC'
commonseq2_gDNA = 'GTACTCGCAGTAGTC'

def process_RNA_fastqs(arguments):
    """
    Output single file with parsed bcs from bc_fastq in read names and seqs from paired_fastq.
    """
    paired_fpath = arguments.paired_fastq_file
    paired_bname = os.path.splitext(paired_fpath[:-3] if paired_fpath.endswith('.gz') else paired_fpath)[0]
    out_fname = f'{os.path.basename(paired_bname)}_parsed.fq'
    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)
    out_fpath = os.path.join(arguments.output_dir, out_fname)

    log.info('Building aligners and barcode decoders')
    aligners = [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_RNA, 'N'*8, 'N'*8) for cso in commonseq1_options]
    bcd = BCDecoder(arguments.barcode_whitelist, arguments.max_bc_err_decode)
    sbcd = SBCDecoder(arguments.sample_barcode_whitelist, arguments.max_sbc_err_decode, arguments.sbc_reject_delta)
    def parse_recs(bc_rec, p_rec):
        if not all(c1 == c2 or set(c1, c2) == set('12') for c1, c2 in zip(str(bc_rec.id), str(p_rec.id))):
            raise ValueError('Non-matching records:\n{bc_rec.id}\n{p_rec.id}')
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
        corrected_pieces.append(new_pieces[-1])

        raw_pieces_str = ':'.join(raw_pieces)
        corrected_pieces_str = ':'.join(corrected_pieces)
        name = f'{corrected_pieces_str}/{new_score:.2f}/{raw_pieces_str}/{raw_score:.2f}/{bc_rec.id}'
        p_rec.id = name
        p_rec.name = name
        p_rec.description = ''
        return raw_score, p_rec

    log.info('Processing fastqs:')
    log.info(f'  barcode fastq: {arguments.bc_fastq_file}')
    log.info(f'  paired fastq:  {arguments.paired_fastq_file}')
    log.info('Writing output to:')
    log.info(f'  {out_fpath}')

    n_first_seqs = 1000
    with open(out_fpath, 'w') as out:
        log.info(f'Processing first {n_first_seqs:,d} for score threshold...')
        first_scores_and_recs = []
        for i, (bc_rec, p_rec) in enumerate(zip(
            SeqIO.parse(gzip_friendly_open(arguments.bc_fastq_file), 'fastq'),
            SeqIO.parse(gzip_friendly_open(arguments.paired_fastq_file), 'fastq')
            )):
            first_scores_and_recs.append(parse_recs(bc_rec, p_rec))
            if i >= n_first_seqs:
                break

        scores = [score for score, rec in first_scores_and_recs]
        thresh = np.average(scores) - 2 * np.std(scores)
        log.info(f'Score threshold: {thresh:.2f}')
        out_recs = [rec for score, rec in first_scores_and_recs if score >= thresh and rec]
        total_out = len(out_recs)
        SeqIO.write(out_recs, out, 'fastq')

        log.info('Continuing...')
        for i, (bc_rec, p_rec) in enumerate(zip(
            SeqIO.parse(gzip_friendly_open(arguments.bc_fastq_file), 'fastq'),
            SeqIO.parse(gzip_friendly_open(arguments.paired_fastq_file), 'fastq')
            )):
            if i <= n_first_seqs:
                continue
            if i % 10000 == 0 and i > 0:
                log.info(f'  {i:,d}')
            score, rec = parse_recs(bc_rec, p_rec)
            if score >= thresh and rec:
                total_out += 1
                SeqIO.write(rec, out, 'fastq')
        log.info(f'{i+1:,d} records processed')
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
    print(f'Norm score threshold: {thresh:.2f} for {bc_fastq_fpath}')

    bc_sbc_cntr = defaultdict(Counter)
    for rec in SeqIO.parse(open(bc_fastq_fpath), 'fastq'):
        norm_score = norm_score_from_rec(rec)
        if norm_score < thresh:
            continue
        bcs, sbc = bc_parser.bc_and_sbc_from_rec(rec)
        if None not in bcs and sbc is not None:
            bc_sbc_cntr[bcs][sbc] += 1
    return bc_sbc_cntr


