import logging
from Bio import Align
from Bio.Seq import Seq

log = logging.getLogger(__name__)

class CustomBCAligner:
    aligner = Align.PairwiseAligner()
    aligner.wildcard = 'N'
    aligner.mismatch = -1
    aligner.gap_score = -1.1
    aligner.target_left_gap_score = -1.9
    aligner.query_left_gap_score = -1.9
    aligner.target_right_gap_score = 0
        
    def __init__(self, *args, unknown_read_orientation=False):
        """
        Input a list of prefixes that are strings of interest. Ns wild.
        """
        self._unknown_read_orientation = unknown_read_orientation
        self.prefixes = args
        self.full_prefix = ''.join(self.prefixes)
        self.prefix_ends = [sum(len(p) for p in self.prefixes[:i+1]) for i in range(len(self.prefixes))]
        self.prefix_all_Ns = [set(prefix) == set('N') for prefix in self.prefixes]
        self.max_query_len = int(1.5*len(self.full_prefix))
        
    def find_norm_score_and_key_boundaries(self, seq: Seq):
        """
        Find best alignment and return norm_score=score/alignment_length and boundary positions.
        """
        alignment = self.aligner.align(self.full_prefix, str(seq[:self.max_query_len]))[0]
        if self._unknown_read_orientation:
            revcompseq = seq.reverse_complement()
            alignment2 = self.aligner.align(self.full_prefix, str(revcompseq[:self.max_query_len]))[0]
            if alignment2.score > alignment.score:
                alignment = alignment2
                seq = revcompseq

        obs_ends = [None for _ in range(len(self.prefixes))]
        obs_idx = 0
        for i in range(len(alignment.aligned[0])):
            tstart, tend = alignment.aligned[0][i]
            qstart, qend = alignment.aligned[1][i]

            for obs_idx, prefix_end in enumerate(self.prefix_ends[:-1]):
                if obs_ends[obs_idx] is None:
                    if tstart <= prefix_end <= tend:
                        obs_ends[obs_idx] = qstart + prefix_end - tstart
                    elif tstart > prefix_end:
                        if i == 0 and obs_idx == 0:  # bizarre alignment. discard
                            return None
                        obs_ends[obs_idx] = alignment.aligned[1][i-1][1]  # prev_qend

            if obs_ends[-2] is not None:
                break
                    
                
        tstart, tend = alignment.aligned[0][-1]
        qstart, qend = alignment.aligned[1][-1]
        obs_ends[-1] = qend + len(self.full_prefix) - tend
        
        return alignment.score/len(self.full_prefix), obs_ends, seq

        
    def find_norm_score_and_pieces(self, seq: Seq, return_seq=False):
        norm_score, obs_ends, seq = self.find_norm_score_and_key_boundaries(seq)
        pieces = [str(seq[:obs_ends[0]])] + [str(seq[obs_ends[i]:obs_ends[i+1]]) for i in range(len(obs_ends)-1)]
        return (norm_score, pieces) if not return_seq else (norm_score, pieces, seq)

    def find_norm_score_pieces_and_end_pos(self, seq: Seq, return_seq=False):
        norm_score, obs_ends, seq = self.find_norm_score_and_key_boundaries(seq)
        pieces = [str(seq[:obs_ends[0]])] + [str(seq[obs_ends[i]:obs_ends[i+1]]) for i in range(len(obs_ends)-1)]
        return (norm_score, pieces, obs_ends[-1]) if not return_seq else (norm_score, pieces, obs_ends[-1], seq)
