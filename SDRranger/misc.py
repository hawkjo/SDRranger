import gzip
import logging


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


def names_pair(s1, s2):
    return all(c1 == c2 or set([c1, c2]) == set('12') for c1, c2 in zip(s1, s2))
