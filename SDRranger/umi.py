import pysam
import numpy as np
from .misc import DistanceThresh, gx_gn_tups_from_read
from Bio import SeqIO
from typing import Tuple, Dict
from collections import Counter, defaultdict
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import connected_components


def get_umi_map_from_cntr(
        umi_cntr: Counter, 
        max_dist: int = 1, 
        dist_type: str = 'freediv',
        connection_type: str = 'directional'
        ) -> dict:
    """
    Builds a dict from observed umis to connected component umi with max count.
    """
    umi_list = list(umi_cntr.keys())
    if connection_type == 'undirected':
        n_vals, component_array = get_connected_components(umi_list, max_dist=max_dist, dist_type=dist_type)
    else:
        assert connection_type == 'directional', connection_type
        n_vals, component_array = get_directional_connected_components(umi_list, umi_cntr, max_dist=max_dist, dist_type=dist_type)
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
        dist_type: str = 'freediv',
        connection_type: str = 'directional'
        ) -> Dict[str, Counter]:
    """
    Builds umi_map_given_bc_then_feature for reads from (specified region of) a given file.

    returns
        :dict: umi_map_given_bc_then_feature
    """
    umi_cntr_given_bc_then_feature = defaultdict(lambda: defaultdict(Counter))
    for read in pysam.AlignmentFile(bam_fpath).fetch(chrm, start, end):
        bc = read.get_tag('CB')
        umi = read.get_tag('UR')
        for gx_gn_tup in gx_gn_tups_from_read(read):
            umi_cntr_given_bc_then_feature[bc][gx_gn_tup][umi] += 1

    umi_map_given_bc_then_feature = defaultdict(dict)
    for bc, umi_cntr_given_feature in umi_cntr_given_bc_then_feature.items():
        for feature, umi_cntr in umi_cntr_given_feature.items():
            umi_map_given_bc_then_feature[bc][feature] = get_umi_map_from_cntr(umi_cntr, max_dist=max_dist, dist_type=dist_type, connection_type=connection_type)

    umi_map_given_bc_then_feature = dict(umi_map_given_bc_then_feature)
    return umi_map_given_bc_then_feature

def get_connected_components(
        umis: list, 
        max_dist: int = 1,
        dist_type: str = 'freediv'
        ) -> Tuple[int, np.array]:
    """
    Builds adjacency matrix of umis given max_dist and calls connected_componenets

    returns
        :int:       n_vals
        :int_array: component_array
    """
    dist_func = DistanceThresh(dist_type, max_dist)

    adj_mat = lil_matrix((len(umis), len(umis)), dtype=np.uint8)
    for i, umi_i in enumerate(umis):
        for j in range(i+1, len(umis)):
            umi_j = umis[j]
            if dist_func(umi_i, umi_j) is not False:
                adj_mat[i, j] = 1
                adj_mat[j, i] = 1
    return connected_components(adj_mat)

def get_directional_connected_components(
        umis: list,
        umi_cntr: Counter,
        max_dist: int = 1,
        dist_type: str = 'freediv'
        ) -> Tuple[int, np.array]:
    """
    Builds adjacency matrix of umis given max_dist and calls connected_componenets

    returns
        :int:       n_vals
        :int_array: component_array
    """
    dist_func = DistanceThresh(dist_type, max_dist)

    # Build direction-based adjacency matrix: i->j if cnt_i >= 2*cnt_j - 1 as in umitools
    adj_mat = lil_matrix((len(umis), len(umis)), dtype=np.uint8)
    nonzero_rows_given_col = defaultdict(list)
    for i, umi_i in enumerate(umis):
        umi_i_cnt = umi_cntr[umi_i]
        for j in range(i+1, len(umis)):
            umi_j = umis[j]
            umi_j_cnt = umi_cntr[umi_j]
            if dist_func(umi_i, umi_j) is not False:
                if umi_i_cnt >= 2 * umi_j_cnt - 1:
                    adj_mat[i, j] = 1
                    nonzero_rows_given_col[j].append(i)
                elif umi_j_cnt >= 2 * umi_i_cnt - 1:
                    adj_mat[j, i] = 1
                    nonzero_rows_given_col[i].append(j)

    # Only allow one incoming edge per node. Keep max count or random if equal.
    for j, i_list in nonzero_rows_given_col.items():
        if len(i_list) > 1:
            max_i = max(i_list, key=lambda i: umi_cntr[umis[i]])
            for i in i_list:
                if i != max_i:
                    adj_mat[i, j] = 0

    return connected_components(adj_mat)
