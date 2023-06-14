import editdistance
import pysam
import numpy as np
from Bio import SeqIO
from typing import Tuple, Dict
from collections import Counter, defaultdict
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import connected_components


def simple_hamming_distance(s1, s2):
    return sum(1 for c1, c2 in zip(s1, s2) if c1 != c2)

def get_connected_components(
        umis: list, 
        max_dist: int = 1,
        dist_type: str = 'hamming'
        ) -> Tuple[int, np.array]:
    """
    Builds adjacency matrix of umis given max_dist and calls connected_componenets

    returns
        :int:       n_vals
        :int_array: component_array
    """
    if dist_type == 'hamming':
        dist_func = simple_hamming_distance
    elif dist_type == 'edit':
        dist_func = editdistance.distance
    else:
        raise ValueError('dist type must be either hamming or edit')

    adj_mat = lil_matrix((len(umis), len(umis)), dtype=np.uint8)
    for i, umi_i in enumerate(umis):
        for j in range(i+1, len(umis)):
            umi_j = umis[j]
            if dist_func(umi_i, umi_j) <= max_dist:
                adj_mat[i, j] = 1
                adj_mat[j, i] = 1
    return connected_components(adj_mat)

def get_umi_map_from_cntr(umi_cntr: Counter, max_dist: int = 1, dist_type: str = 'hamming') -> dict:
    """
    Builds a dict from observed umis to connected component umi with max count.
    """
    umi_list = list(umi_cntr.keys())
    n_vals, component_array = get_connected_components(umi_list, max_dist=max_dist, dist_type=dist_type)
    component_max_umi = [None for _ in range(n_vals)]
    for umi, component in zip(umi_list, component_array):
        if component_max_umi[component] is None or umi_cntr[umi] > umi_cntr[component_max_umi[component]]:
            component_max_umi[component] = umi
    umi_map = {umi: component_max_umi[component] for umi, component in zip(umi_list, component_array)}
    return umi_map

def get_umi_maps_from_bam_file(
        bam_fpath: str,
        chrm: str = None,
        start: int = None,
        end: int = None,
        max_dist: int = 1,
        dist_type: str = 'hamming'
        ) -> Dict[str, Counter]:
    """
    Builds umi_map_given_bc for reads from (specified region of) a given file.

    returns
        :dict: umi_map_given_bc
    """
    umi_cntr_given_bc = defaultdict(Counter)
    for read in pysam.AlignmentFile(bam_fpath).fetch(chrm, start, end):
        bc = read.get_tag('CB')
        umi = read.get_tag('UR')
        umi_cntr_given_bc[bc][umi] += 1
    umi_cntr_given_bc = dict(umi_cntr_given_bc)

    umi_map_given_bc = {}
    for bc, umi_cntr in umi_cntr_given_bc.items():
        umi_map_given_bc[bc] = get_umi_map_from_cntr(umi_cntr, max_dist=max_dist, dist_type=dist_type)

    return umi_map_given_bc
