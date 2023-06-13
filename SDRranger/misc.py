import os
import gzip
import glob
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


def find_paired_fastqs_in_dir(fq_dir):
    raw_prompts = ['*.fastq', '*.fq', '*.txt']
    glob_prompts = [prompt for rs in raw_prompts for prompt in [rs, f'{rs}.gz']]
    fq_fpaths = [
        fpath
        for glob_prompt in [os.path.join(fq_dir, s) for s in glob_prompts]
        for fpath in glob.glob(glob_prompt)
    ]
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
