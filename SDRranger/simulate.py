#!/usr/bin/env python3
import sys
import os
import gzip
from collections import defaultdict

import numpy as np

alphabet = ("A", "T", "G", "C")
complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

alphabet_subsets = {"A": ("T", "G", "C"),
                    "T": ("A", "G", "C"),
                    "G": ("A", "T", "C"),
                    "C": ("A", "T", "G")}

def one_mutation(seq: str, position: int, substitutionprob: float, insertionprob: float, rng: np.random.Generator):
    if rng.random() < substitutionprob:
        seq = seq[:position] + rng.choice(alphabet_subsets[seq[position]]) + seq[position + 1:]
    else:
        if rng.random() < insertionprob:
            seq = seq[:position] + rng.choice(alphabet) + seq[position:] # 1bp insertion
        else:
            seq = seq[:position] + seq[position + 1:] # 1bp deletion
    return seq

def mutate_read_probabilistic(seq: str, errorprob: float, substitutionprob: float, insertionprob: float, rng: np.random.Generator):
    pos = 0
    while pos < len(seq):
        if rng.random() < errorprob:
            seq = one_mutation(seq, pos, substitutionprob, insertionprob, rng)
        pos += 1
    return seq

def mutate_read_maxerrors(seq: str, maxerrors: int, substitutionprob: float, insertionprob: float, rng: np.random.Generator):
    positions = rng.choice(len(seq), size=maxerrors, replace=False)
    for pos in positions:
        seq = one_mutation(seq, pos, substitutionprob, insertionprob, rng)
    return seq

def make_read(blocks: dict, name_prefix: str, mutate_umi: bool, mutatereadfun, error_probability: float, substitution_probability, insertion_probability: float, umis: dict[list[str]], rng: np.random.Generator):
    read = ""
    readname = name_prefix

    for blockidx, block in enumerate(blocks):
        if block["blocktype"] == "randomBarcode":
            if not mutate_umi:
                blocksequence = "".join(rng.choice(alphabet, block["length"]))
                blockname = blocksequence
                umis[blockidx].append(blocksequence)
            else:
                blockname = rng.choice(umis[blockidx])
                blocksequence = mutatereadfun(blockname, block["maxerrors"], error_probability, substitution_probability, insertion_probability, rng)
        elif block["blocktype"] == "constantRegion" or block["blocktype"] == "barcodeList":
            seqidx = rng.choice(len(block["sequence"]))
            blockname = block["name"][seqidx] if "name" in block else block["sequence"][seqidx]
            blocksequence = mutatereadfun(block["sequence"][seqidx], block["maxerrors"], error_probability, substitution_probability, insertion_probability, rng)
        readname += f"_{blockname}"
        read += blocksequence

    nonbc_sequence = "".join(rng.choice(alphabet, rng.poisson(20)))
    read += nonbc_sequence

    return readname, read

def simulate_reads(arguments):
    rng = np.random.default_rng(arguments.seed)

    if "barcode_struct_r2" in arguments.config:
        fq1name = f"{arguments.fastq_prefix}_1.txt.gz"
        fq2name = f"{arguments.fastq_prefix}_2.txt.gz"
        fq2 = gzip.open(fq2name, "wt")
    else:
        fq1name = f"{arguments.fastq_prefix}.txt.gz"

    os.makedirs(os.path.dirname(arguments.fastq_prefix), exist_ok=True)

    fq1 = gzip.open(fq1name, "wt")

    read1_umis = defaultdict(list)
    read2_umis = defaultdict(list)

    if arguments.error_probability <= 0:
        mutatereadfun = lambda seq, maxerrors, errorprob, substitutionprob, insertionprob, rng: mutate_read_maxerrors(seq, maxerrors, substitutionprob, insertionprob, rng)
    else:
        mutatereadfun = lambda seq, maxerrors, errorprob, substitutionprob, insertionprob, rng: mutate_read_probabilistic(seq, errorprob, substitutionprob, insertionprob, rng)

    for readidx in range(arguments.nreads):
        read1name, read1 = make_read(arguments.config["barcode_struct_r1"]["blocks"], str(readidx), readidx > arguments.unique_umis * arguments.nreads,  mutatereadfun, arguments.error_probability, arguments.substitution_probability, arguments.insertion_probability, read1_umis, rng)
        revcomp = False
        if arguments.config["unknown_read_orientation"] and rng.random() > 0.5:
            read1 = read1.translate(complement)[::-1]
            revcomp = True
        print(read1name, read1, "+", "E" * len(read1), sep="\n", file=fq1)

        if "barcode_struct_r2" in arguments.config:
            read2name, read2 = make_read(arguments.config["barcode_struct_r2"]["blocks"], str(readidx), readidx > arguments.unique_umis * arguments.nreads, mutatereadfun, arguments.error_probability, arguments.substitution_probability, arguments.insertion_probability, read1_umis, rng)
            if revcomp:
                read2 = read2.translate(complement)[::-1]
            print(read2name, read2, "+", "E" * len(read2), sep="\n", file=fq2)

    fq1.close()

    if "barcode_struct_r2" in arguments.config:
        fq2.close()
