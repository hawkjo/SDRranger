"""
SDRranger: Processing Genomic-Transcriptomic TAP-Seq data.

Usage:
  SDRranger count_RNA        <fastq_dir> [--STAR-ref-dir=<>] [--barcode-whitelist=<>] [--max-bc-err-decode=<>] [--sample-barcode-whitelist=<>] [--max-sbc-err-decode=<>] [--sbc-reject-delta=<>] [--output-dir=<>] [--threads=<>] [-v | -vv | -vvv]
  SDRranger count_gDNA       <fastq_dir> [--STAR-ref-dir=<>] [--barcode-whitelist=<>] [--max-bc-err-decode=<>] [--output-dir=<>] [--threads=<>] [-v | -vv | -vvv]

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
from .RNAcount import process_RNA_fastqs
from .gDNAcount import process_gDNA_fastqs 


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
