import logging
import os
import json

from .misc import gzip_friendly_open


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.
    """
    def __init__(self, arguments):
        self._arguments = arguments

        try:
            with gzip_friendly_open(arguments["--config"]) as f:
                self._config = json.load(f)
        except KeyError:
            self._config = None

    def _comma_delimited_arg(self, key):
        if self._arguments[key]:
            return self._arguments[key].split(',')
        return []

    @property
    def command(self):
        # We have to do this weird loop to deal with the way docopt stores the command name
        for possible_command in ('count_gDNA', 'count_RNA', 'count_matrix'):
            if self._arguments.get(possible_command):
                return possible_command

    @property
    def fastq_dir(self):
        return self._arguments['<fastq_dir>']

    @property
    def SDR_bam_file(self):
        return self._arguments['<SDR_bam_file>']

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)

    @property
    def star_ref_dir(self):
        return self._arguments['--STAR-ref-dir']

    @property
    def star_output_path(self):
        return self._arguments['--STAR-output']

    @property
    def config(self):
        return self._config

    @property
    def threads(self):
        return int(self._arguments['--threads'])

    @property
    def output_dir(self):
        return self._arguments['--output-dir']
