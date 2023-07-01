

def process_gDNA_fastq(arguments):
    aligners = [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_gDNA) for cso in commonseq1_options]
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(gzip_friendly_open(bc_fq_fpath), 'fastq'):
            seq = str(rec.seq)
            norm_score, pieces, end_pos = max([al.find_norm_score_pieces_and_end_pos(seq) for al in aligners])
            pieces_str = ','.join(pieces)
            name = f'{pieces_str}/{norm_score}/{rec.id}'
            rec = rec[end_pos:]
            rec.id = name
            rec.name = name
            rec.description = ''
            SeqIO.write(rec, out, 'fastq')


def norm_score_from_rec(rec):
    return float(str(rec.id).split('/')[1])


def bc_and_sbc_counter_from_fastq(bc_fq_fpath, bc_parser):
    # Determine minimum norm_score threshold = mean-2sigma. See notebook for figures
    norm_scores_sample = []
    for i, rec in enumerate(SeqIO.parse(open(bc_fq_fpath), 'fastq')):
        if i >= 100000:
            break
        norm_scores_sample.append(norm_score_from_rec(rec))
    thresh = np.average(norm_scores_sample) - 2 * np.std(norm_scores_sample)





def process_gDNA_fastqs(arguments):
    aligners = [CustomBCAligner('N'*9, cso, 'N'*9, commonseq2_gDNA) for cso in commonseq1_options]
    with open(out_fpath, 'w') as out:
        for rec in SeqIO.parse(gzip_friendly_open(bc_fq_fpath), 'fastq'):
            seq = str(rec.seq)
            norm_score, pieces, end_pos = max([al.find_norm_score_pieces_and_end_pos(seq) for al in aligners])
            pieces_str = ','.join(pieces)
            name = f'{pieces_str}/{norm_score}/{rec.id}'
            rec = rec[end_pos:]
            rec.id = name
            rec.name = name
            rec.description = ''
            SeqIO.write(rec, out, 'fastq')


def norm_score_from_rec(rec):
    return float(str(rec.id).split('/')[1])


def bc_and_sbc_counter_from_fastq(bc_fq_fpath, bc_parser):
    # Determine minimum norm_score threshold = mean-2sigma. See notebook for figures
    norm_scores_sample = []
    for i, rec in enumerate(SeqIO.parse(open(bc_fq_fpath), 'fastq')):
        if i >= 100000:
            break
        norm_scores_sample.append(norm_score_from_rec(rec))
    thresh = np.average(norm_scores_sample) - 2 * np.std(norm_scores_sample)
    log.info(f'Norm score threshold: {thresh:.2f} for {bc_fq_fpath}')

    bc_sbc_cntr = defaultdict(Counter)
    for rec in SeqIO.parse(open(bc_fq_fpath), 'fastq'):
        norm_score = norm_score_from_rec(rec)
        if norm_score < thresh:
            continue
        bcs, sbc = bc_parser.bc_and_sbc_from_rec(rec)
        if None not in bcs and sbc is not None:
            bc_sbc_cntr[bcs][sbc] += 1
    return bc_sbc_cntr



