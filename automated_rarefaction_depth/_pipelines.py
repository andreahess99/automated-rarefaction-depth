from ._methods import _nw_align_defaults


#example from the tutorial
def align_and_summarize(
        ctx, seq1, seq2,
        gap_open_penalty=_nw_align_defaults['gap_open_penalty'],
        gap_extend_penalty=_nw_align_defaults['gap_extend_penalty'],
        match_score=_nw_align_defaults['match_score'],
        mismatch_score=_nw_align_defaults['mismatch_score']):
    nw_align_action = ctx.get_action('dwq2', 'nw_align')
    summarize_alignment_action = ctx.get_action('dwq2', 'summarize_alignment')

    msa, = nw_align_action(
                    seq1, seq2, gap_open_penalty=gap_open_penalty,
                    gap_extend_penalty=gap_extend_penalty,
                    match_score=match_score, mismatch_score=mismatch_score)
    msa_summary, = summarize_alignment_action(msa)

    return (msa, msa_summary)