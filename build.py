#!/usr/bin/env python3
"""Generates the numbers, tables and figures of the paper from model/final.json, proto/tests.json,
proto/retries.json and proto/full_*.json, assembles out/whippet.tex (one file) and compiles it."""
import json, os, re, subprocess, sys
from math import log2, floor, ceil, sqrt

ROOT = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.join(ROOT, "tex")
OUT = os.path.join(ROOT, "out")
os.makedirs(OUT, exist_ok=True)
F = json.load(open(os.path.join(ROOT, "model", "final.json")))
TESTS = json.load(open(os.path.join(ROOT, "proto", "tests.json")))
RET = json.load(open(os.path.join(ROOT, "proto", "retries.json")))
PAIRS = json.load(open(os.path.join(ROOT, "proto", "pairs.json")))
FULL = {}
for tag in ("C", "A", "Q", "G", "Gc"):
    fn = os.path.join(ROOT, "proto", "full_%s.json" % tag)
    if os.path.exists(fn):
        FULL[tag] = json.load(open(fn))

if "--draft" in sys.argv:           # layout check while experiments are still running: placeholders, not for release
    pass
    for tag in ("C", "A", "Q", "G", "Gc"):
        if tag not in FULL:
            FULL[tag] = dict(FULL["G" if "G" in FULL else "Q"])
        FULL[tag].setdefault("stats", [dict(fill=0.99, jl_tries=1, c_tries=1)] * FULL[tag]["levels"])

YPRIME = 64 * 32 / 8192.0        # Greyhound does not send y' [NS24, Remark 4.3]: remove it from the baselines
SET = {"C": ("C", 438), "A": ("A", 438), "Q": ("A", 483)}
SETNAME = {"C": "conventional, $438$", "A": "proof-based, $438$", "Q": "proof-based, $483$"}
NS = (26, 28, 30)
ABC = dict(zip(NS, "abc"))
comma = lambda x: "{:,}".format(int(round(x))).replace(",", "{,}")


def ent(fam, kind, size, s, var):
    acc, bkz = SET[s]
    return F["%s_%s%d_%s_%d_%s" % (fam, kind, size, acc, bkz, var)]


def tot(fam, kind, size, s, var):
    e = ent(fam, kind, size, s, var)
    return e["total"] - (YPRIME if (fam == "base" and kind == "gh") else 0.0)


W = lambda n, s, var="full": tot("whip", "gh", n, s, var)
B = lambda n, s, var="base": tot("base", "gh", n, s, var)
G = lambda n: F["whip_gh%d_A_483_gc" % n]
f1 = lambda x: "%.1f" % x
f2 = lambda x: "%.2f" % x
pct = lambda new, old: 100.0 * (1 - new / old)
num = {}

# ---------------------------------------------------------------- headline numbers
for s in "CAQ":
    for n in NS:
        num["W%s%s" % (s, ABC[n])] = f1(W(n, s))
        num["B%s%s" % (s, ABC[n])] = f1(B(n, s))
        num["L%s%s" % (s, ABC[n])] = f1(B(n, s, "lin"))
        num["S%s%s" % (s, ABC[n])] = f1(B(n, s, "linall"))
for n in NS:
    num["GQ%s" % ABC[n]] = f1(G(n)["total"])
sav = [pct(W(n, s), B(n, s)) for s in "CAQ" for n in NS]
savl = [pct(W(n, s), B(n, s, "lin")) for s in "CAQ" for n in NS]
savs = [pct(W(n, s), B(n, s, "linall")) for s in "CAQ" for n in NS]
num["PctLo"], num["PctHi"] = str(floor(min(sav))), str(ceil(max(sav)))
num["PctLinLo"], num["PctLinHi"] = str(floor(min(savl))), str(ceil(max(savl)))
num["PctStrongLo"], num["PctStrongHi"] = str(floor(min(savs))), str(ceil(max(savs)))
num["PctStrongQ"] = str(round(pct(W(30, "Q"), B(30, "Q", "linall"))))
gp = [100.0 * (G(n)["total"] / W(n, "Q") - 1) for n in NS]
num["GPctLo"], num["GPctHi"] = str(floor(min(gp))), str(ceil(max(gp)))
num["BQcDec"] = f1(B(30, "Q") * 1.024)
e = ent("base", "gh", 30, "C", "base")
num["TailShare"] = str(round(100 * e["tail"] / (e["total"] - YPRIME)))
coms = [ent("whip", "gh", n, s, "full")["com"] for s in "CAQ" for n in NS]
num["ComLo"], num["ComHi"] = str(round(min(coms))), str(round(max(coms)))
num["GComHi"] = str(round(max(G(n)["com"] for n in NS)))
eq = ent("whip", "gh", 30, "Q", "full")
num["LiftW"], num["ScW"] = f2(eq["comp"]["lift"] + 12 / 1024), f2(eq["comp"]["sc"])
tr = eq["trace"]
num["SchedM"] = comma(eq["W0"])
num["SchedMtwo"] = comma(tr[0]["m2"])
num["SchedLastKB"] = f1(tr[-1]["bits"] / 8192)
num["SchedLastPct"] = str(round(100 * tr[-1]["bits"] / 8192 / eq["total"]))

# ---------------------------------------------------------------- ablation
abl = lambda fam, var, s: round(tot(fam, "gh", 30, s, var), 1)      # differences of the displayed values
num["AblG"] = f1(abl("base", "base", "Q") - abl("base", "lin", "Q"))
num["AblS"] = f1(abl("base", "lin", "Q") - abl("whip", "plain", "Q"))
num["AblE"] = f1(abl("whip", "plain", "Q") - abl("whip", "code", "Q"))
num["AblR"] = f1(abl("whip", "code", "Q") - abl("whip", "noinner", "Q"))
num["AblI"] = f1(abl("whip", "noinner", "Q") - abl("whip", "full", "Q"))
num["AblStrong"] = f1(abl("base", "linall", "Q") - abl("whip", "full", "Q"))
num["AblFeight"] = f1(abl("whip", "f8", "Q") - abl("whip", "full", "Q"))
num["AblP"] = f1(abl("whip", "P", "Q") - abl("whip", "full", "Q"))
rows = [("Greyhound with LaBRADOR", "base", "base"),
        ("\\quad without quadratic garbage $g_{ij}$", "base", "lin"),
        ("\\quad linear tail (sum-check, lifting to $\\FF$)", "whip", "plain"),
        ("\\quad + entropy coding", "whip", "code"),
        ("\\quad + rounding (outer, first fold, terminal)", "whip", "noinner"),
        ("\\quad + rounding (inner, recursive levels): this work", "whip", "full"),
        None,
        ("LaBRADOR without $g_{ij}$, rounded and coded", "base", "linall"),
        ("this work without entropy coding", "whip", "round"),
        ("this work with $\\lvert\\FF\\rvert=q^8$", "whip", "f8"),
        ("this work with the bound $2\\tau\\beta^2$ for $\\vz$", "whip", "P")]
with open(os.path.join(OUT, "tab_ablation.tex"), "w") as fh:
    for r in rows:
        if r is None:
            fh.write("\\midrule\n")
            continue
        fh.write("%s & %s\\\\\n" % (r[0], " & ".join(f1(abl(r[1], r[2], s)) for s in "CAQ")))
    fh.write("this work with guaranteed completeness & -- & %s & %s\\\\\n" % (
        f1(F["whip_gh30_A_438_gc"]["total"]), f1(G(30)["total"])))

# ---------------------------------------------------------------- main table and headline table
with open(os.path.join(OUT, "tab_main.tex"), "w") as fh:
    for s in "CAQ":
        fh.write("%s & Greyhound & %s\\\\\n" % (SETNAME[s], " & ".join(f1(B(n, s)) for n in NS)))
        fh.write(" & \\quad without $g_{ij}$ & %s\\\\\n" % " & ".join(f1(B(n, s, "lin")) for n in NS))
        fh.write(" & \\quad rounded and coded & %s\\\\\n" % " & ".join(f1(B(n, s, "linall")) for n in NS))
        fh.write(" & \\textbf{this work} & %s\\\\\n" % " & ".join("\\textbf{%s}" % f1(W(n, s)) for n in NS))
        fh.write(" & saving (\\%%) & %s\\\\\n" % " & ".join(
            "%d / %d / %d" % (round(pct(W(n, s), B(n, s))), round(pct(W(n, s), B(n, s, "lin"))),
                              round(pct(W(n, s), B(n, s, "linall")))) for n in NS))
        fh.write(" & commitment (bytes) & %s\\\\\n" % " & ".join(
            str(round(ent("whip", "gh", n, s, "full")["com"])) for n in NS))
        fh.write("\\midrule\n")
    fh.write("guaranteed, $483$ & \\textbf{this work} & %s\\\\\n" % " & ".join("\\textbf{%s}" % f1(G(n)["total"]) for n in NS))
    fh.write(" & commitment (bytes) & %s\\\\\n" % " & ".join(str(round(G(n)["com"])) for n in NS))
with open(os.path.join(OUT, "tab_head.tex"), "w") as fh:
    fh.write("\\multicolumn{3}{@{}l}{\\emph{Size model}}\\\\\n")
    fh.write("Greyhound with LaBRADOR & proof-based & $%s$\\\\\n" % f1(B(30, "Q")))
    fh.write("\\quad without $g_{ij}$, rounded and coded & proof-based & $%s$\\\\\n" % f1(B(30, "Q", "linall")))
    fh.write("\\textbf{Greyhound with the linear tail (this work)} & proof-based & $\\mathbf{%s}$\\\\\n" % f1(W(30, "Q")))
    fh.write("\\textbf{this work} & guaranteed & $\\mathbf{%s}$\\\\\n" % f1(G(30)["total"]))
    fh.write("\\addlinespace\n")
    fh.write("Greyhound with LaBRADOR & conventional & $%s$\\\\\n" % f1(B(30, "C")))
    fh.write("\\textbf{this work} & conventional & $\\mathbf{%s}$\\\\\n" % f1(W(30, "C")))

# ---------------------------------------------------------------- schedule
SHAPE = {26: (3156, 333, 6, 5, 7, 5), 28: (6312, 665, 5, 6, 6, 5), 30: (12625, 1329, 4, 8, 6, 5)}
with open(os.path.join(OUT, "tab_sched.tex"), "w") as fh:
    m0, r0, lb0, d0, lb, _ = SHAPE[30]
    fh.write("0 & $%s$ & -- & %d & $%s$ & -- & %d & %d & %d & %d & %d & %s\\\\\n" % (
        comma(d0 * m0 * r0), r0, comma(d0 * m0), eq["n"], eq["n1"], eq["Dt"], eq["Du"], 2 ** lb, f2(eq["own"])))
    m = int(eq["W0"])
    for i, lv in enumerate(tr):
        last = lv["kind"] == "last"
        mark = ""
        if not last:
            mk = ("*" if lv["wmode"] == "clear" else "") + ("\\dagger" if not lv["dec"] else "")
            mark = "^{%s}" % mk if mk else ""
        fh.write("$%d%s$ & $%s$ & %.1f & %d & %d & %d & %d & %s & %d & %s & %s & %s\\\\\n" % (
            i + 1, mark, comma(m), log2(lv["beta"]), lv["r"], lv["n"], lv["mu"],
            lv["kap"], "--" if last else str(lv["k1"]), lv["di"], "--" if last else str(lv["do"]),
            "--" if (last or not lv["dec"]) else str(lv["b"]), f2(lv["bits"] / 8192)))
        if not last:
            m = int(lv["m2"])

# ---------------------------------------------------------------- composition figure
def comp_whip(s):
    e = ent("whip", "gh", 30, s, "full")
    c = e["comp"]
    return [e["own"], c["outer"], c["t"], c["z"], c["garb"], c["jl"], c["lift"], c["sc"], c["chal"]]


def comp_base(s, var):
    e = ent("base", "gh", 30, s, var)
    c = e["comp"]
    return [e["own"] - YPRIME, c["outer"], c["t"], c["z"], c["gh"], c["jl"], c["lift"], 0.0, c["chal"]]


bars = [("W", comp_whip("Q")), ("L$'$", comp_base("Q", "lin")), ("L", comp_base("Q", "base")),
        ("W", comp_whip("C")), ("L$'$", comp_base("C", "lin")), ("L", comp_base("C", "base"))]
names = ["first fold", "outer commitments", "inner commitments", "opening $\\vz$", "cross terms", "projections",
         "lifting", "sum-check", "seeds"]
cols = ["ours!55", "cOuter", "cInner", "cZ", "cGarb", "cJL", "cLift", "cSC", "cChal"]
with open(os.path.join(OUT, "fig_bars.tex"), "w") as fh:
    fh.write("\\begin{tikzpicture}\n\\begin{axis}[xbar stacked, width=0.86\\linewidth, height=6.1cm, bar width=9.5pt,\n"
             "xmin=0, xmax=78, xlabel={KB}, ytick={0,1,2,3.6,4.6,5.6}, yticklabels={%s},\n"
             "ymin=-0.7, ymax=6.3, tick label style={font=\\small}, label style={font=\\small},\n"
             "legend style={at={(0.5,-0.24)}, anchor=north, legend columns=5, font=\\footnotesize, draw=none, legend cell align=left, /tikz/every even column/.append style={column sep=7pt}},\n"
             "xmajorgrids, grid style={gray!25}, axis line style={gray!60}, clip=false]\n" % ",".join(b[0] for b in bars))
    ys = [0, 1, 2, 3.6, 4.6, 5.6]
    for k in range(9):
        fh.write("\\addplot[fill=%s, draw=black!55, line width=0.2pt] coordinates {%s};\n" % (
            cols[k], " ".join("(%.3f,%s)" % (b[1][k], ys[i]) for i, b in enumerate(bars))))
    fh.write("\\legend{%s}\n" % ",".join(names))
    for i, b in enumerate(bars):
        fh.write("\\node[anchor=west,font=\\footnotesize] at (axis cs:%.2f,%s) {%s};\n" % (sum(b[1]) + 0.5, ys[i], f1(sum(b[1]))))
    fh.write("\\node[anchor=east,font=\\footnotesize,rotate=90] at (axis cs:-7.5,6.1) {conventional};\n")
    fh.write("\\node[anchor=east,font=\\footnotesize,rotate=90] at (axis cs:-7.5,2.55) {proof-based};\n")
    fh.write("\\end{axis}\n\\end{tikzpicture}\n")

# ---------------------------------------------------------------- linear relations
LW = (12, 14, 16, 18, 20, 22)
lin = lambda fam, lw, s, var: tot(fam, "lin", lw, s, var)
wl = [lin("whip", lw, "C", "full") for lw in LW]
num["LinLo"], num["LinHi"] = f1(min(wl)), f1(max(wl))
sl = [pct(lin("whip", lw, "C", "full"), lin("base", lw, "C", "lin")) for lw in LW]
num["LinPctLo"], num["LinPctHi"] = str(floor(min(sl))), str(ceil(max(sl)))
with open(os.path.join(OUT, "fig_lin.tex"), "w") as fh:
    fh.write("\\begin{tikzpicture}\n\\begin{axis}[width=0.74\\linewidth, height=6.4cm, xlabel={$\\log_2m$}, ylabel={KB},\n"
             "xmin=11.5, xmax=22.5, ymin=0, ymax=70, xtick={12,14,16,18,20,22}, grid=major, grid style={gray!25},\n"
             "tick label style={font=\\small}, label style={font=\\small},\n"
             "legend style={at={(1.03,0.5)}, anchor=west, font=\\footnotesize, draw=none}, legend cell align=left]\n")
    spec = [("base", "base", "A", "cOuter, dashed, mark=square*, mark options={solid}", "LaBRADOR (P)"),
            ("base", "lin", "A", "cZ, dashed, mark=triangle*, mark options={solid}", "without $g_{ij}$ (P)"),
            ("whip", "full", "A", "ours, dashed, mark=*, mark options={solid}", "this work (P)"),
            ("base", "base", "C", "cOuter, mark=square*", "LaBRADOR (C)"),
            ("base", "lin", "C", "cZ, mark=triangle*", "without $g_{ij}$ (C)"),
            ("whip", "full", "C", "ours, thick, mark=*", "this work (C)")]
    for fam, var, s, style, name in spec:
        fh.write("\\addplot[%s, mark size=1.6pt] coordinates {%s};\n\\addlegendentry{%s}\n" % (
            style, " ".join("(%d,%.2f)" % (lw, lin(fam, lw, s, var)) for lw in LW), name))
    fh.write("\\end{axis}\n\\end{tikzpicture}\n")

# ---------------------------------------------------------------- implementation: sizes
TAGS = ("C", "A", "Q", "G")
modb = lambda r: r["model_total"] - r["seeds"]
devs = [abs(100 * (FULL[t]["size"] / modb(FULL[t]) - 1)) for t in "CAQ"]
num["FullDev"] = "0.1" if max(devs) <= 0.1 else "%.2f" % max(devs)
fills = sorted(FULL[t]["fill"] for t in "CAQ")
num["FullFill"] = "%.2f" % fills[0] if "%.2f" % fills[0] == "%.2f" % fills[-1] else "%.2f--%.2f" % (fills[0], fills[-1])
num["FullGFill"] = "%.3f" % FULL["G"]["fill"]
num["FullGBytes"], num["FullGModel"] = comma(FULL["G"]["size"]), comma(modb(FULL["G"]))
assert all(FULL[t]["ok"] and FULL[t]["tampered_rejected"] for t in FULL)
comp_rows = [("outer commitments", "outer", "outer"), ("projections", "jl", "jl"), ("lifted values", "lift", "lift"),
             ("sum-check", "sumcheck", "sumcheck"), ("partial evaluations", "evals", "evals"),
             ("inner commitments", "inner", "inner"), ("opening $\\vz$", "z", "z"), ("nonces", "nonce", None)]
with open(os.path.join(OUT, "tab_impl.tex"), "w") as fh:
    for name, km, kmod in comp_rows:
        cells = []
        for tag in TAGS:
            r = FULL[tag]
            cells.append("$%s$" % comma(r["sizes"].get(km, 0)))
            cells.append("$%s$" % comma(r["model"].get(kmod, 0)) if kmod else "--")
        fh.write("%s & %s\\\\\n" % (name, " & ".join(cells)))
    fh.write("\\midrule\n")
    cells = []
    for tag in TAGS:
        r = FULL[tag]
        cells += ["$\\mathbf{%s}$" % comma(r["size"]), "$%s$" % comma(modb(r))]
    fh.write("total & %s\\\\\n" % " & ".join(cells))
    fh.write("witness (ring elements) & %s\\\\\n" % " & ".join("\\multicolumn{2}{c}{$%s$}" % comma(FULL[t]["m"]) for t in TAGS))
    fh.write("verified; flipped bit rejected & %s\\\\\n" % " & ".join("\\multicolumn{2}{c}{yes}" for t in TAGS))

# ---------------------------------------------------------------- implementation: attempts and norms
st = {t: FULL[t]["stats"] for t in ("Q", "G", "Gc")}
nl = max(len(v) for v in st.values())
with open(os.path.join(OUT, "tab_stats.tex"), "w") as fh:
    for i in range(nl):
        cells = []
        for t in ("Q", "G", "Gc"):
            if i < len(st[t]):
                s = st[t][i]
                cells += ["%.3f" % s["fill"], str(s["jl_tries"]), str(s["c_tries"])]
            else:
                cells += ["", "", ""]
        fh.write("%d & %s\\\\\n" % (i + 1, " & ".join(cells)))
opt_fills = [s["fill"] for t in "CAQ" for s in FULL[t]["stats"][1:]]
num["FillLo"], num["FillHi"] = "%.3f" % min(opt_fills), "%.3f" % max(opt_fills)
num["JlMaxOpt"] = str(max(s["jl_tries"] for t in "CAQ" for s in FULL[t]["stats"]))
num["CMaxOpt"] = str(max(s["c_tries"] for t in "CAQ" for s in FULL[t]["stats"]))
assert max(s["fill"] for t in ("G", "Gc") for s in FULL[t]["stats"]) <= 1.0

# ---------------------------------------------------------------- implementation: repeated runs
KIND = {"uniform": "random", "const": "all coefficients equal", "sparse": "one coefficient per ring element"}
runs0 = RET["opt_uniform"]["runs"]
num["RetRuns"] = str(len(runs0))
num["RetLevOpt"], num["RetLevGc"] = str(RET["opt_uniform"]["levels"]), str(RET["gc_uniform"]["levels"])
att = lambda r: sum(r["jl"]) + sum(r["c"])
with open(os.path.join(OUT, "tab_retries.tex"), "w") as fh:
    for mode, mname in (("opt", "optimised"), ("gc", "guaranteed")):
        for j, kind in enumerate(("uniform", "const", "sparse")):
            rr = RET["%s_%s" % (mode, kind)]["runs"]
            okr = [r for r in rr if r["ok"]]
            a = [att(r) for r in okr]
            fh.write("%s & %s & %d of %d & %.1f & %d\\\\\n" % (mname if j == 0 else "", KIND[kind], len(okr), len(rr),
                                                               sum(a) / len(a), max(a)))
        if mode == "opt":
            fh.write("\\midrule\n")
oc = RET["opt_const"]["runs"]
num["RetOptUni"] = "%.1f" % (sum(att(r) for r in runs0) / len(runs0))
num["RetOptConstMax"] = str(max(att(r) for r in oc if r["ok"]))
num["RetOptConstFail"] = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}.get(sum(1 for r in oc if not r["ok"]), str(sum(1 for r in oc if not r["ok"])))
gcr = [r for k in ("gc_uniform", "gc_const", "gc_sparse") for r in RET[k]["runs"]]
assert all(r["ok"] for r in gcr)
num["RetGcMax"] = str(max(att(r) for r in gcr))
num["RetM"], num["RetLogB"] = comma(3000), "13"

good = [t for t in TESTS if t["opts"].get("E", True) is not False]
assert all(t["ok"] for t in TESTS)
num["TestN"] = {5: "five", 6: "six", 7: "seven", 8: "eight"}[len(good)]
dev = max(abs(t["dev"]) for t in good)
num["TestDev"] = "0.1" if dev <= 0.1 else "%.1f" % dev

# ---------------------------------------------------------------- challenge pairs
with open(os.path.join(OUT, "tab_pairs.tex"), "w") as fh:
    for k in ("3", "4"):
        c = PAIRS[k]["c"]
        assert len(c) == 64 and c[0] - c[1] == int(k)
        fh.write("\\begin{center}\\footnotesize\\setlength{\\tabcolsep}{3.5pt}\n\\begin{tabular}{@{}*{16}{r}@{}}\n")
        for i0 in range(0, 64, 16):
            fh.write(" & ".join("$%d$" % v for v in c[i0:i0 + 16]) + "\\\\\n")
        fh.write("\\end{tabular}\\\\[4pt]\n$\\opn c=%.2f$,\\quad$\\opn{c'}=%.2f$,\\quad$c-c'=%s(1-X)$\n\\end{center}\n" % (
            PAIRS[k]["op"], PAIRS[k]["op2"], k))

# ---------------------------------------------------------------- assemble
with open(os.path.join(OUT, "numbers.tex"), "w") as fh:
    for k in sorted(num):
        fh.write("\\newcommand{\\n%s}{%s}\n" % (k, num[k]))
print("numbers: %d macros" % len(num))

order = ["01_intro", "02_prelim", "03_field", "04_tail", "05_pcs", "06_params", "07_impl", "08_discussion"]
body = open(os.path.join(TEX, "00_preamble.tex")).read()
body = body.replace("%%NUMBERS%%", open(os.path.join(OUT, "numbers.tex")).read())
for name in order:
    body += "\n" + open(os.path.join(TEX, name + ".tex")).read()
appendix = open(os.path.join(TEX, "91_appendix.tex")).read()
inline = lambda txt: re.sub(r"\\input\{gen/([a-z_]+)\.tex\}",
                            lambda mo: open(os.path.join(OUT, mo.group(1) + ".tex")).read().rstrip("\n"), txt)
text = inline(body + appendix)
used = set(re.findall(r"\\n([A-Z][A-Za-z]*)", text))
missing = sorted(u for u in used if u not in num)
assert not missing, missing
cited = set()
for grp in re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]*)\}", text):
    cited.update(k.strip() for k in grp.split(","))
refs = open(os.path.join(TEX, "90_refs.tex")).read()
head, rest = refs.split("\\bibitem", 1)
items = ("\\bibitem" + rest).split("\\end{thebibliography}")[0]
items = re.split(r"(?=\\bibitem)", items)
bkey = lambda it: re.match(r"\\bibitem(?:\[[^\]]*\])?\{([^}]*)\}", it.replace("\n", " ")).group(1)
keep = [it for it in items if it.strip() and bkey(it) in cited]
have = set(bkey(it) for it in keep)
assert cited <= have, sorted(cited - have)
refs = head + "".join(keep).rstrip() + "\n\n\\end{thebibliography}\n"
full = inline(body + "\n" + refs + "\n" + appendix + "\n\\end{document}\n")
open(os.path.join(OUT, "whippet.tex"), "w").write(full)
print("cited %d references, kept %d of %d" % (len(cited), len(keep), len([i for i in items if i.strip()])))
if "--nocompile" not in sys.argv:
    for _ in range(3):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "whippet.tex"], cwd=OUT,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if r.returncode != 0:
            print(r.stdout[-2500:])
            sys.exit(1)
    log = open(os.path.join(OUT, "whippet.log"), errors="replace").read()
    print("pages:", re.findall(r"Output written on whippet.pdf \((\d+) pages", log))
    print("overfull:", len(re.findall(r"Overfull \\hbox", log)), " undefined:", len(re.findall(r"undefined", log)))
    for mo in re.finditer(r"Overfull \\hbox \(([0-9.]+)pt too wide\)[^\n]*lines? (\d+)", log):
        if float(mo.group(1)) > 4:
            print("  overfull %spt at line %s" % (mo.group(1), mo.group(2)))
    for mo in re.finditer(r"LaTeX Warning: (Reference|Citation)[^\n]*", log):
        print(" ", mo.group(0))
