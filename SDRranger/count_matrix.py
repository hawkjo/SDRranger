import logging
import pysam
from .count import count_matrix

log = logging.getLogger(__name__)
pysam.set_verbosity(0)

def build_count_matrices_from_bam(arguments):
    for read in pysam.AlignmentFile(arguments.SDR_bam_file).fetch():
        break
    if read.has_tag('UB'):
        log.info('Detected bam file with UMIs...')
        RNA_count_matrix(arguments, arguments.SDR_bam_file)
    else:
        log.info('Detected bam file without UMIs...')
        gDNA_count_matrix(arguments, arguments.SDR_bam_file)
