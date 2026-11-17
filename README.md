# Whippet: a shorter tail for lattice-based polynomial commitments

Everything needed to regenerate the paper `out/whippet.pdf`.

    model/params.py      size model of LaBRADOR (baseline), ranks, coding and rounding formulas
    model/wmodel.py      size model of the linear tail (levels, terminal level, beam search)
    model/gh.py          first fold of Greyhound as the source of the tail instance
    model/run_final.py   all model numbers of the paper -> model/final.json   (about 1 h on 2 cores, resumable)
    proto/ring.py        R_q = Z_q[X]/(X^64+1), q = 2^32-99; the field F = Z_q[X^8+X^-8]; trace; challenges
    proto/coder.py       range coder with a discretised Gaussian model; fixed-radix packing
    proto/tail.py        reference prover and verifier of the linear tail (non-interactive)
    proto/run_tests.py   small instances, tamper tests            -> proto/tests.json
    proto/run_full.py    full-size instances                      -> proto/full_{C,A,Q}.json (+ .log)
    proto/chk_pairs.py   challenge pairs with difference 3(1-X) and 4(1-X) (Proposition 8.1)
    proto/run_retries.py numbers of attempts of the honest prover, random and structured witnesses -> proto/retries.json
    CLAUDE_CODE_PROMPT.md  instructions for the next implementation step (first fold, end-to-end experiments)
    tex/*.tex            sources of the paper
    build.py             numbers, tables, figures; assembles out/whippet.tex and runs pdflatex

Usage:

    cd model && python3 run_final.py            # optional: recompute the model numbers
    cd proto && python3 run_tests.py            # optional: small prover/verifier runs
    cd proto && python3 run_full.py whip_gh30_A_438_full full_A.json    # about 10 min, 3 GB
    python3 build.py

Requirements: Python 3 with NumPy; pdflatex with TikZ and pgfplots.

Settings.  'C' conventional and 'A' proof-based ranks (model/params.py, wmodel.py); WOpt(gc=True) selects the
guaranteed-completeness norm bounds (deterministic bounds for digits and low parts, joint Markov bound for
(z, e), projection margin); keys in model/final.json: whip_gh30_A_483_full, whip_gh30_A_483_gc, ...
The full-size runs of the paper are
    python3 run_full.py whip_gh30_C_438_full full_C.json
    python3 run_full.py whip_gh30_A_438_full full_A.json
    python3 run_full.py whip_gh30_A_483_full full_Q.json
    python3 run_full.py whip_gh30_A_483_gc   full_G.json
    python3 run_full.py whip_gh30_A_483_gc   full_Gc.json const

Status of the numbers.  Proof sizes are computed by the model.  The reference implementation covers the
tail (not the first fold) and is used to check the model and the completeness of the protocol; it is not
optimised and its running times are not meaningful.  Baselines are models of LaBRADOR and Greyhound.
