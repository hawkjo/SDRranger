"""
SDRranger: Processing Genomic-Transcriptomic TAP-Seq data.

Usage:
  SDRranger count_RNA        <bc_fastq_file> <paired_fastq_file> [--barcode-whitelist=<barcode_whitelist>] [--max-bc-err-decode=<max_bc_err_decode>] [--sample-barcode-whitelist=<sample_barcode_whitelist>] [--max-sbc-err-decode=<max_sbc_err_decode>] [--sbc-reject-delta=<sbc_reject_delta>] [--output-dir=<output_dir>] [-v | -vv | -vvv]
  SDRranger count_gDNA       <bc_fastq_file> <paired_fastq_file> [--barcode-whitelist=<barcode_whitelist>] [--max-bc-err-decode=<max_bc_err_decode>] [--output-dir=<output_dir>] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  gDNA          Process and count Genomic gDNA files
  RNA           Process and count Transcriptomic RNA files
"""
import logging
import os
from docopt import docopt
from .__init__ import __version__
from .config import CommandLineArguments
from .count import process_RNA_fastqs, process_gDNA_fastqs 


def main(**kwargs):
    docopt_args = docopt(__doc__, version=__version__)
    arguments = CommandLineArguments(docopt_args, os.getcwd())

    log = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(arguments.log_level)
    log.debug(docopt_args)

    commands = {
        'count_RNA': process_RNA_fastqs,
        'count_gDNA': process_gDNA_fastqs,
    }

    commands[arguments.command](arguments)


if __name__ == '__main__':
    main()
