import logging
import os
import json

from .misc import gzip_friendly_open

class CommandLineArguments:
    def __new__(cls, arguments):
        if arguments["simulate_reads"]:
            return SimulationCommandLineArguments(arguments)
        else:
            return AnalysisCommandLineArguments(arguments)

class CommandLineArgumentsBase:
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
        for possible_command in ('count', 'preprocess', 'count_matrix', 'simulate_reads'):
            if self._arguments.get(possible_command):
                return possible_command
    @property
    def config(self):
        return self._config

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)

class AnalysisCommandLineArguments(CommandLineArgumentsBase):
    """
    Wraps the raw arguments provided by docopt.
    """
    def __init__(self, arguments):
        super().__init__(arguments)

    @property
    def fastq_dir(self):
        return self._arguments['<fastq_dir>']

    @property
    def SDR_bam_file(self):
        return self._arguments['<SDR_bam_file>']

    @property
    def star_ref_dir(self):
        return self._arguments['--STAR-ref-dir']

    @property
    def star_output_path(self):
        return self._arguments['--STAR-output']

    @property
    def threads(self):
        return int(self._arguments['--threads'])

    @property
    def output_dir(self):
        return self._arguments['--output-dir']

class SimulationCommandLineArguments(CommandLineArgumentsBase):
    def __init__(self, arguments):
        super().__init__(arguments)

    @property
    def fastq_prefix(self):
        return self._arguments['--fastq-prefix']

    @property
    def nreads(self):
        return int(self._arguments['--nreads'])

    @property
    def unique_umis(self):
        return float(self._arguments['--unique-umis'])

    @property
    def seed(self):
        return int(self._arguments['--seed'])

    @property
    def error_probability(self):
        return float(self._arguments['--error-probability'])

    @property
    def substitution_probability(self):
        return float(self._arguments['--substitution-probability'])

    @property
    def insertion_probability(self):
        return float(self._arguments['--insertion-probability'])
