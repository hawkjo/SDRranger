import logging
import os


class CommandLineArguments(object):
    """
    Wraps the raw arguments provided by docopt.
    """
    def __init__(self, arguments, current_directory):
        self._arguments = arguments
        self._current_directory = current_directory

    def _comma_delimited_arg(self, key):
        if self._arguments[key]:
            return self._arguments[key].split(',')
        return []

    @property
    def command(self):
        # We have to do this weird loop to deal with the way docopt stores the command name
        for possible_command in ('count_gDNA', 'count_RNA'):
            if self._arguments.get(possible_command):
                return possible_command

    @property
    def fastq_dir(self):
        return self._arguments('<fastq_dir>')

    @property
    def log_level(self):
        log_level = {0: logging.ERROR,
                     1: logging.WARN,
                     2: logging.INFO,
                     3: logging.DEBUG}
        # default to silent if the user supplies no verbosity setting
        return log_level.get(self._arguments['-v'], logging.ERROR)

    @property
    def max_bc_err_decode(self):
        return int(self._arguments['--max-bc-err-decode'] or 1)

    @property
    def max_sbc_err_decode(self):
        return int(self._arguments['--max-sbc-err-decode'] or 1)

    @property
    def sbc_reject_delta(self):
        return int(self._arguments['--sbc-reject-delta'] or 0)

    @property
    def barcode_whitelist(self):
        return os.path.expanduser(self._arguments['--barcode-whitelist']) if self._arguments['--barcode-whitelist'] else None

    @property
    def sample_barcode_whitelist(self):
        return os.path.expanduser(self._arguments['--sample-barcode-whitelist']) if self._arguments['--sample-barcode-whitelist'] else None

    @property
    def threads(self):
        return int(self._arguments['--threads'] or 1)

    @property
    def output_dir(self):
        return self._arguments['--output-dir'] or '.'
