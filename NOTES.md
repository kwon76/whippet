# Whippet — working notes (keep updated; survives context compaction)

User request (Korean, 2026-11-17): aim at a clear improvement over latest SOTA, CRYPTO-level paper, not
restricted to LaBRADOR/Greyhound; deliver tex + pdf, no AI tells, figures + protocol diagrams.
Honesty constraints: no fabricated experiments; sizes from model + reference implementation; say what is
not achieved. Reply to user in Korean. Do not spawn agents. Memory writes only if user asks.

## Idea (decided)
New tail replacing LaBRADOR for *linear* relations (what Greyhound/Hachi/Akita hand to their tail):
 1. Linear fold: all constraints after the first fold are linear -> no quadratic garbage g_ij.
 2. Trace lifting: F = Z[eta]/q, eta = X^8+X^-8, is the field F_{q^4} inside R_q (q = 5 mod 8 inert in
    Q(zeta_16)^+). Theta = normalised trace R_q -> F, ct = tau_0 o Theta. Lifting a Z_q-constraint costs
    3 coefficients (not 63); R-constraints are aggregated via Theta(alpha * f), alpha in R_q.
    Theta(a) = a[0] + sum_{k=1..3} (a[8k]-a[64-8k])/2 * eta_k; emb(c) = c0 + sum c_k (X^{8k} - X^{64-8k}).
 3. One sum-check over F (2 elements of F per round, log2 n rounds over the position inside a block) turns
    the aggregated constraint Theta(<Phi,w>) = b into partial evaluations W_i = <eq(rho), s_i> in R_q
    (r of them instead of (r^2+r)/2 garbage), exactly the shape of Greyhound's first fold.
    Final claim: Theta(sum_i Phi'_i W_i) = claim, Phi'_i = sum_k eq(rho,k) Phi_{i,k}.
 4. Rounded commitments everywhere (inner at all levels, outer), entropy coding (Akita already Rice-codes z).
 5. JL projection only on s (not on low parts); binding needed: outer 2*SL*beta', inner (normal form)
    8T*zf*SL*beta' (recursion) / 8T*sqrt(gamma^2+gamma_e^2) (terminal). Only "same c-bar" weak binding used.
Soundness error per level ~ 2^-128 + q^-4 + (2 mu+1)/|F| + 2/|C|.
Name: Whippet ("a smaller greyhound").

## SOTA facts (verified by WebFetch of Akita PDF, eprint 2026/1983)
- Akita Table 8 (KB, excl. commitment; bits committed 2^27..2^35 = N 2^22..2^30 of 32-bit coeffs):
  Akita direct 61.3/61.8/63.1/64.5/64.6; offloaded -/66.2/66.3/66.9/67.3; Greyhound recalibrated
  60.1/60.6/63.3/64.9/67.2; RoKoko 109-115. Recalibration: 128-bit *quantum* ADPS16 (BKZ ~483), ranks for
  every matrix role, sparse-ternary JL slack sqrt(128/29)?, grinding. Verify 8.1-15.9 ms vs Greyhound 155-1423 ms.
- Akita terminal follows LaBRADOR 5.6 (reveals z, checks A z = sum c_i t_i), Rice-codes z; relaxed-opening
  radius^2 64*Gamma^2*S_max (i.e. 8*Gamma); future work (Sec. 14): "cheaper folds and a smaller tail",
  rigorous completeness (response model heuristic, eq. (185)), QROM, ZK. Commitments 128 B via 2-stage
  binary Ajtai chain. Hachi: ~55 KB incl. estimated Greyhound tail (secondary source).
- Greyhound relation handed to LaBRADOR is purely linear (eq. (12): P z = h, ||z|| <= gamma).

## Files
- model/params.py (old LaBRADOR model), model/run_eval.py (greyhound_instance), model/wmodel.py (new tail
  model: wlevel/wlast/search, set_bkz), model/run_w.py -> weval.json (results).
- proto/ring.py (R_q, F, Theta, transcript, challenges), proto/coder.py (range coder, radix packing),
  proto/tail.py (prover/verifier/statements), proto/run_small.py, proto/run_tests.py -> tests.json.
- old paper: /home/claude/paper (tex/*.tex, build.py, chal_results.json, ccdf.dat) — reuse preamble,
  completeness section (Thm moments), rounding lemmas.

## Results so far (model, KB; Greyhound N=2^26/28/30)
- acc C (reproduces published LaBRADOR sizes): ~29.6/30.6/31.0 (baseline model 50.9/52.6/53.7; published 46/53/53)
- acc A (proof-based), BKZ 438: ~38-39 (baseline 60.5/62.9/64.6).  BKZ 483: see weval.json.
- prototype: 1500-ring-element instance: verified; measured 28865 B vs model 28856 B.

## TODO
- finish tests (tamper), larger prototype run (target m=304383, beta=2^16.54 = Greyhound 2^30 tail input)
- converse-of-weak-binding proposition + numeric mu for challenge pairs
- baseline under BKZ 483; ablations; figures; paper tex (build script); verification pass; Korean summary


## State on 2026-11-17 (after context compaction) — PAPER WRITTEN
- tex/00..08, 90_refs, 91_appendix written; build.py generates out/numbers.tex, tab_*.tex, fig_*.tex and
  out/whippet.tex (single file) and compiles (25 pages).
- model/final.json complete (90 entries incl. linear-specialised LaBRADOR baseline 'lin', 'linall').
  Baselines for Greyhound are reduced by 0.25 KB in build.py (y' is not sent, [NS24, Remark 4.3]).
- Headline (KB, N=2^26/28/30): conv/438 Whippet 29.6/30.6/31.0 vs 50.6/52.4/53.5 (lin 47.2/47.7/48.7);
  proof-based/483 Whippet 39.0/40.0/40.7 vs 64.3/65.8/68.5 (lin 60.9/61.9/64.3). Savings 39-42% / 35-38%.
- Full-scale prototype: A/438 38,509 B (model 38,478), C/438 30,208 B, A/483 39,787 B; all verified.
- Proposition 8.1 (non-uniqueness of weak openings below ~6 nu) with delta=3(1-X), delta'=4(1-X);
  existence of such challenge pairs checked by proto/chk_pairs.py.
- Remaining: final visual check, deliver pdf+tex (+code zip), save tex to Project, Korean summary.


## Revision after the first review (2026-11-17, afternoon)
Review (Korean, uploaded by the user): borderline / weak reject; main requests: completeness for the
parameters, experiments on real first-fold witnesses, a final PCS security theorem with a clear ROM status,
strengthened baseline at block size 483, Akita revision/units, small fixes in Appendix A.
Done in this revision:
- model: WOpt(gc=True) = guaranteed-completeness setting (wmodel._wlevel_gc/_wlast_gc, gh.first_fold_gc,
  instance_gc); new jobs in run_final.py (gc for gh26/28/30 at A/483, linall baselines for all N and settings,
  ablation variants at A/483).  Results: guaranteed 45.8 / 46.7 / 47.4 KB; rounded-and-coded LaBRADOR-linear
  baseline at A/483: 50.7 / 50.9 / 51.7 KB (after removing y').
- proto/tail.py: bind_statement (statement + schedule absorbed into the transcript), projection bound 'bjl',
  joint bound 'qmax'/'joint' with retry at every level, STATS (norm ratios, attempts per level).
- proto/run_retries.py (random / all-equal / sparse witnesses, optimised vs guaranteed), run_full.py stats and
  'const' witness; full_G.json, full_Gc.json.
- paper: abstract and contributions reorganised (two folds on two axes first; 40.7 KB first), Theorem 5.2
  (PCS: binding, knowledge error 2/|C| + lambda*eps0), Section 5.3 (status of Fiat-Shamir: no theorem),
  Section 6.6 (Lemma 6.2 Cantelli for projections, Theorem 6.3 completeness with guaranteed parameters),
  Section 7 rewritten (scope, sizes for 4 settings, attempts and norms, repeated runs), Appendix A (extractor
  with restart, uniqueness of s*, padding), Appendix C (explicit challenge pairs), k_q = 4, units.
- NOT done (see CLAUDE_CODE_PROMPT.md): first fold implementation and end-to-end runs on polynomials,
  timings/memory, fast ring arithmetic.
- Could not re-read the Akita PDF (robots); abstract page of 2026-10-15 says 58-68 KB.  The numbers 64.6 (Akita
  direct, 2^30) agree between my earlier reading and the reviewer; recalibrated Greyhound is 67.2 in my earlier
  reading and 66.1 according to the reviewer -> paper says "66 to 67", units of 1000 bytes (reviewer).  CHECK.
