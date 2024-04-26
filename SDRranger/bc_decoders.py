import Levenshtein
import freebarcodes
import freebarcodes.decode
import logging

log = logging.getLogger(__name__)

class BCDecoder:
    def __init__(self, bc_whitelist, bc_maxdist):
        self.bcs = bc_whitelist
        self.bcs_set = set(self.bcs)
        self.bc_maxdist = bc_maxdist
        self.bc_len = len(self.bcs[0])
        assert all(len(bc) == self.bc_len for bc in self.bcs)

    def decode(self, raw_bc):
        if raw_bc in self.bcs_set:
            return raw_bc
        dists_and_scores = [(Levenshtein.distance(raw_bc, bc, score_cutoff=self.bc_maxdist), bc) for bc in self.bcs]
        min_dist, bc = min(dists_and_scores)
        if min_dist > self.bc_maxdist:
            return None
        if sum(1 for dist, bc in dists_and_scores if dist == min_dist) > 1:
            return None
        return bc


class SBCDecoder:
    def __init__(self, sbc_whitelist, sbc_maxdist, sbc_reject_delta):
        self.sbcs = sbc_whitelist
        self.sbc_len = len(self.sbcs[0])
        assert all(len(sbc) == self.sbc_len for sbc in self.sbcs)
        self.sbc_maxdist = sbc_maxdist
        self.sbc_reject_delta = sbc_reject_delta
        self.sbcd = freebarcodes.decode.FreeDivBarcodeDecoder()
        self.sbcd.build_codebook_from_random_codewords(self.sbcs, self.sbc_maxdist, self.sbc_reject_delta)

    def decode(self, raw_sbc):
        sbc = self.sbcd.decode(raw_sbc)
        return sbc if isinstance(sbc, str) else None

