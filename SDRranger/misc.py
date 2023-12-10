import os
import gzip
import glob
import logging
import pysam
import scipy
import numpy as np
from Bio import SeqIO
from .bc_aligner import CustomBCAligner
from .constants import commonseq1_options, commonseq2_RNA, commonseq2_gDNA


log = logging.getLogger(__name__)


def gzip_friendly_open(fpath):
    return gzip.open(fpath, 'rt') if fpath.endswith('gz') else open(fpath)


def load_bc_list(fpath):
    bc_list = [line.strip() for line in gzip_friendly_open(fpath)]
    if not bc_list:
        raise RuntimeError(f'No barcodes found in {arguments.barcode_file}')
    bc_list = [bc[:bc.index('-')] if '-' in bc else bc for bc in bc_list] # clean cellranger bcs
    return bc_list


def write_stats_file_from_cntr(cntr, fpath):
    with open(fpath, 'w') as out:
        for k, v in sorted(cntr.items()):
            s = f'{k}: {v:,d}'
            log.info(s)
            out.write(s + '\n')


def find_paired_fastqs_in_dir(fq_dir):
    raw_prompts = ['*.fastq', '*.fq', '*.txt']
    glob_prompts = [os.path.join(fq_dir, prompt) for rp in raw_prompts for prompt in [rp, f'{rp}.gz']]
    fq_fpaths = [fpath for glob_prompt in glob_prompts for fpath in glob.glob(glob_prompt)]
    fq_fpaths.sort()

    paired_names = []
    while fq_fpaths:
        fpath1 = fq_fpaths.pop(0)
        fpath2 = fq_fpaths.pop(0)
        if not names_pair(fpath1, fpath2):
            raise ValueError(f'Unexpected input in fastq directory. Following files do not pair:\n{fpath1}\n{fpath2}')
        paired_names.append((fpath1, fpath2))
    return paired_names


def names_pair(s1, s2):
    """
    A very general test for whether two names are pair-compatible.
    
    Intentionally permissive. Flags true any pair of strings identical except for possibly 1/2.
    """
    return all(c1 == c2 or set([c1, c2]) == set('12') for c1, c2 in zip(s1, s2))


def file_prefix_from_fpath(fpath):
    """Strip away directories and extentions from file path"""
    bname = os.path.splitext(fpath[:-3] if fpath.endswith('.gz') else fpath)[0]
    return os.path.basename(bname)


def determine_bc_and_paired_fastq_idxs(paired_fpaths, RNA_or_gDNA):
    """Returns (bc_fq_idx, paired_fq_idx), based on first pair of fastqs"""
    fpath0, fpath1 = paired_fpaths[0]
    avg_score0 = average_align_score_of_first_recs(fpath0, RNA_or_gDNA)
    avg_score1 = average_align_score_of_first_recs(fpath1, RNA_or_gDNA)
    return (0, 1) if avg_score0 > avg_score1 else (1, 0)


def build_gDNA_bc_aligners():
    """Helper function, kept as function for parallelism"""
    return [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_gDNA) for cso in commonseq1_options]


def build_RNA_bc_aligners():
    """Helper function, kept as function for parallelism"""
    return [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_RNA, 'N'*8, 'N'*8) for cso in commonseq1_options]


def average_align_score_of_first_recs(fastq_fpath, RNA_or_gDNA, n_seqs=500):
    """Return average alignment score of first n_seqs records"""
    if RNA_or_gDNA == 'RNA':
        aligners = build_RNA_bc_aligners()
    elif RNA_or_gDNA == 'gDNA':
        aligners = build_gDNA_bc_aligners()
    else:
        raise ValueError('Invalid input for RNA_or_gDNA')
    first_scores = []
    for i, rec in enumerate(SeqIO.parse(gzip_friendly_open(fastq_fpath), 'fastq')):
        scores_and_pieces = [al.find_norm_score_and_pieces(str(rec.seq)) for al in aligners]
        raw_score, raw_pieces = max(scores_and_pieces)
        first_scores.append(raw_score)
        if i >= n_seqs:
            break
    return np.average(first_scores)


def gx_gn_tups_from_read(read):
    try:
        gxs = read.get_tag('GX').split(';')
        gns = read.get_tag('GN').split(';')
        return list(zip(gxs, gns))
    except:
        return []


def get_bcs_and_features_from_bam(input_bam_fpath, build_complete_bc, threads=1):
    complete_bcs, features = set(), set()
    for read in pysam.AlignmentFile(input_bam_fpath, threads=threads).fetch():
        complete_bcs.add(build_complete_bc(read))
        features.update(gx_gn_tups_from_read(read))
    sorted_complete_bcs = sorted(complete_bcs)
    sorted_features = sorted(features)
    return sorted_complete_bcs, sorted_features


def write_matrix(M, bcs, features, out_dir):
    matrix_fpath = os.path.join(out_dir, 'matrix.mtx.gz')
    with gzip.open(matrix_fpath, 'wb') as out:
        scipy.io.mmwrite(out, M)

    rows_fpath = os.path.join(out_dir, 'barcodes.tsv.gz')
    with gzip.open(rows_fpath, 'wt') as out:
        out.write('\n'.join(bcs))

    cols_fpath = os.path.join(out_dir, 'features.tsv.gz')
    with gzip.open(cols_fpath, 'wt') as out:
        out.write('\n'.join([f'{gx}\t{gn}\tGene Expression' for gx, gn in features]))
