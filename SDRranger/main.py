"""
SDRranger: Process SDR-seq data 

Usage:
  SDRranger count            <fastq_dir> (--STAR-ref-dir=<> | --STAR-output=<>...) --config=<> [--output-dir=<>] [--threads=<>] [-v | -vv | -vvv]
  SDRranger preprocess       <fastq_dir> --config=<> [--output-dir=<>] [--threads=<>] [-v | -vv | -vvv]
  SDRranger count_matrix     <SDR_bam_file> --output-dir=<> [--threads=<>] [-v | -vv | -vvv]
  SDRranger simulate_reads   --config=<> --fastq-prefix=<> --nreads=<> [--unique-umis=<>] [--seed=<>] [--error-probability=<>] [--substitution-probability=<>] [--insertion-probability=<>] [-v | -vv | -vvv]

Options:
  --STAR-ref-dir=<>:              Path to directory with STAR index.
  --STAR-output=<>:               Path to STAR output file (BAM/SAM). Can be repeated multiple times,
                                    in which case the order must correspond to the lexicographic ordering
                                    of paired FASTQ files in <fastq_dir>.
  --config=<>:                    Path to JSON configuration.
  --output-dir=<>:                Path to output directory [default: .].
  --threads=<>:                   Number of threads [default: 1].
  -v:                             Verbose output.
  --fastq-prefix=<>:              Prefix for output FASTQ files.
  --nreads=<>:                    Number of reads to simulate.
  --unique-umis=<>:               Fraction of all reads that have unique UMIs [default: 0.5].
  --seed=<>:                      Random seed [default: 42].
  --error-probability=<>:         Probability of an error occurring [default: 0.1]. Set to a negative number to
                                    always introduce as many errors as allowed by the configuration.
  --substitution-probability=<>:  Probability of generating a substitution as opposed to an indel [default: 0.7].
  --insertion-probability=<>:     Probability of generating an insertion as opposed to a deletion when generating an indel [default: 0.5].
  -h --help                       Show this screen.
  --version                       Show version.

Commands:
  preprocess       Preprocess files such that STAR can be run on the output.
  count            Process and count input files.
  count_matrix     Build a count matrix (or matrices) from an existing bam file.
  simulate_reads   Generate synthetic sequencing reads given a barcode configuration.
"""
import logging
import os
from docopt import docopt
from .__init__ import __version__
from .config import CommandLineArguments
from .count import process_fastqs, preprocess_fastqs
from .count_matrix import build_count_matrices_from_bam
from .simulate import simulate_reads

def main(**kwargs):
    docopt_args = docopt(__doc__, version=__version__)
    arguments = CommandLineArguments(docopt_args)

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    commands = {
        'count': process_fastqs,
        'preprocess': preprocess_fastqs,
        'count_matrix': build_count_matrices_from_bam,
        'simulate_reads': simulate_reads
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
