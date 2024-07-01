import os
import gzip
import glob
import logging
import pysam
import json
from itertools import product
from bisect import bisect

import scipy
import numpy as np
from Bio import SeqIO
from pywfa import WavefrontAligner


log = logging.getLogger(__name__)

class DistanceThresh:
    def __init__(self, dist_type, max_dist):
        self._freediv = False
        if dist_type == "hamming":
            self._aligner = WavefrontAligner(distance="linear", mismatch=1, gap_extension=1000, span="end-to-end", scope="score", max_steps=max_dist + 1)
        elif dist_type == "levenshtein":
            self._aligner = WavefrontAligner(distance="levenshtein", scope="score", max_steps=max_dist + 1)
        elif dist_type == "freediv":
            self._aligner = WavefrontAligner(distance="linear", mismatch=1, gap_extension=1, span="ends-free", scope="score", pattern_begin_free=0, text_begin_free=0, pattern_end_free=1000, text_end_free=1000, max_steps=max_dist + 1)
            self._freediv = True
        else:
            raise ValueError('dist_type must be either hamming, levenshtein, or freediv')

    def __call__(self, s1, s2):
        lendiff = abs(len(s1) - len(s2))
        if lendiff  >= self._aligner.max_steps:
            return False
        else:
            if self._freediv:
                self._aligner.pattern_end_free = len(s2)
                self._aligner.text_end_free = len(s1)
            score = -self._aligner.wavefront_align(s1, s2)
            if self._aligner.status != 0:
                return False
            else:
                return score + self._freediv * lendiff

def gzip_friendly_open(fpath):
    with open(fpath, "rb") as f:
        header = f.read(2)
    return gzip.open(fpath, 'rt') if header == b"\x1f\x8b" else open(fpath)


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

def config_has_barcodes_on_both_reads(config):
    return "barcode_struct_r2" in config and "blocks" in config["barcode_struct_r2"] and len(config["barcode_struct_r2"]["blocks"])

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

def namepair_differences(r1, r2):
    return [i for i, (c1, c2) in enumerate(zip(r1, r2)) if c1 != c2]

def make_paired_name(readname, idx):
    name = ""
    lastidx = 0
    for i in idx:
        name += readname[lastidx:i]
        lastidx = i + 1
    name += readname[lastidx:]
    return name

def get_namepair_index(R1_fpath, R2_fpath):
    r1 = next(SeqIO.parse(gzip_friendly_open(R1_fpath), 'fastq'))
    r2 = next(SeqIO.parse(gzip_friendly_open(R2_fpath), 'fastq'))
    return namepair_differences(r1.id, r2.id)


def file_prefix_from_fpath(fpath):
    """Strip away directories and extentions from file path"""
    bname = os.path.splitext(fpath[:-3] if fpath.endswith('.gz') else fpath)[0]
    return os.path.basename(bname)


def determine_bc_and_paired_fastq_idxs(paired_fpaths, blocks, unknown_read_orientation):
    """Returns (bc_fq_idx, paired_fq_idx), based on first pair of fastqs"""
    fpath0, fpath1 = paired_fpaths[0]
    avg_score0 = average_align_score_of_first_recs(fpath0, blocks, unknown_read_orientation)
    avg_score1 = average_align_score_of_first_recs(fpath1, blocks, unknown_read_orientation)
    return (0, 1) if avg_score0 > avg_score1 else (1, 0)

def build_bc_aligners(blocks, unknown_read_orientation):
    from .bc_aligner import CustomBCAligner

    naligners = []
    barcodelengths = []
    for i, block in enumerate(blocks):
        if block["blocktype"] == "constantRegion":
            naligners.append(block["sequence"])
            barcodelengths.append(None)
        elif block["blocktype"] == "barcodeList":
            barcodelengths.append(len(block["sequence"][0]))
        else:
            barcodelengths.append(block["length"])

    aligners = []
    for aln in product(*naligners):
        alnidx = 0
        prefixes = []
        for ln in barcodelengths:
            if ln is not None:
                prefixes.append("N" * ln)
            else:
                prefixes.append(aln[alnidx])
                alnidx += 1
        aligners.append(CustomBCAligner(*prefixes, unknown_read_orientation=unknown_read_orientation))
    return aligners

def build_bc_decoders(blocks):
    from .bc_decoders import BCDecoder, SBCDecoder

    decoders = []
    lastblock = None
    lastblockidx = -1
    for i, block in enumerate(blocks):
        if block["blocktype"] == "barcodeList":
            if lastblock is not None:
                decoders.append(BCDecoder(lastblock["sequence"], lastblock["maxerrors"]))
            lastblock = block
            lastblockidx = i
    if lastblock is not None:
        if lastblockidx == len(blocks) - 1 or blocks[lastblockidx + 1]["blocktype"] == "randomBarcode":
            decoders.append(SBCDecoder(lastblock["sequence"], lastblock["maxerrors"], 0)) # TODO: expose max_reject_delta?
        else:
            decoders.append(BCDecoder(lastblock["sequence"], lastblock["maxerrors"]))
    return decoders

def average_align_score_of_first_recs(fastq_fpath, blocks, unknown_read_orientation, n_seqs=500):
    """Return average alignment score of first n_seqs records"""
    aligners = build_bc_aligners(blocks, unknown_read_orientation)
    first_scores = []
    for i, rec in enumerate(SeqIO.parse(gzip_friendly_open(fastq_fpath), 'fastq')):
        scores_and_pieces = [al.find_norm_score_and_pieces(rec.seq) for al in aligners]
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


def get_bcs_and_features_from_bam(input_bam_fpath, threads=1):
    complete_bcs, features = set(), set()
    for read in pysam.AlignmentFile(input_bam_fpath, threads=threads).fetch():
        complete_bcs.add(read.get_tag('CB'))
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

def sort_and_index_readname_bam(input_bam_fpath, output_bam_fpath, namepairidx, threads=1):
    pysam.sort("-N", "-@", str(threads), "-o", output_bam_fpath, input_bam_fpath)
    with pysam.AlignmentFile(output_bam_fpath, "r", threads=threads) as bam:
        i_readnames = []
        i_offsets = []
        lastblock = -1
        offset = bam.tell()
        lastreadname = ""
        lastreadoffset = offset
        for read in bam.fetch(until_eof=True):
            block = offset >> 16
            readname = make_paired_name(read.query_name, namepairidx)
            if block > lastblock:
                i_readnames.append(readname)
                i_offsets.append(offset if readname != lastreadname else lastreadoffset)
                lastblock = block
            if readname != lastreadname:
                lastreadname = readname
                lastreadoffset = offset
            offset = bam.tell()

    return i_readnames, i_offsets, namepairidx

def get_bam_read_by_name(name, bam, index, threads=1):
    i_readnames, i_offsets, namepairidx = index
    name = make_paired_name(name, namepairidx)

    idx = bisect(i_readnames, make_paired_name(name, namepairidx)) - 1
    needclose = False
    if not isinstance(bam, pysam.AlignmentFile):
        bam = pysam.AlignmentFile(bam, "r", threads=threads)
        needclose = True
    bam.seek(i_offsets[idx])
    startblock = None
    haveread = False
    for bread in bam.fetch(until_eof=True):
        if startblock is None:
            startblock = bam.tell() >> 16
        if make_paired_name(bread.query_name, namepairidx) == name:
            haveread = True
            yield bread
        elif haveread or bam.tell() >> 16 > startblock:
            break
    if needclose:
        bam.close()
