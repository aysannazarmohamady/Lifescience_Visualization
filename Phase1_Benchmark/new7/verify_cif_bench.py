"""Verify cumulative-incidence benchmark images from side files. Ground truth is the
Aalen-Johansen cumulative incidence recomputed from CLEAN ADTTE.
DRAFT NOTE: drafted with LLM assistance."""
import hashlib, numpy as np, pandas as pd

OUT = "Outputs/cif_bench/"
d = pd.read_csv("Data/V1/ADTTE.csv"); d = d[d.PARAMCD == "PFS"].copy()
d["st"] = d.COMPTYPE.map({"PROGRESSION": 1, "DEATH_WITHOUT_PROGRESSION": 2, "CENSORED": 0})

def aj(g, cause):
    """Cumulative incidence of `cause` (Aalen-Johansen); returns (times, est)."""
    ts = np.sort(g.AVALM[g.st > 0].unique()); s = 1.0; c = 0.0; est = []
    for t in ts:
        n = (g.AVALM >= t).sum(); e = ((g.AVALM == t) & (g.st == cause)).sum(); a = ((g.AVALM == t) & (g.st > 0)).sum()
        c += s * e / n; s *= 1 - a / n; est.append(c)
    return ts, np.array(est)

def at(ts, est, t):
    i = np.searchsorted(ts, t, side="right") - 1
    return 0.0 if i < 0 else est[i]

def km1(g):
    ts = np.sort(g.AVALM[g.st == 1].unique()); s = 1.0; est = []
    for t in ts:
        s *= 1 - ((g.AVALM == t) & (g.st == 1)).sum() / (g.AVALM >= t).sum(); est.append(1 - s)
    return ts, np.array(est)

def classify(tag):
    cv = pd.read_csv(f"{OUT}{tag}_curves.csv"); rk = pd.read_csv(f"{OUT}{tag}_risk.csv")
    meta = dict(l.strip().split(" ", 1) for l in open(f"{OUT}{tag}_meta.txt"))
    scale = 30.4375 if meta["x_unit"] == "days" else 1.0
    def same(ref):                       # plotted curves equal reference at plotted times
        return all(abs(r.est - at(*ref(d[d.ARM == r.arm]), r.time / scale)) < 1e-6 for r in cv.itertuples())
    is_cif = same(lambda g: aj(g, 1)); is_km = same(km1); is_death = same(lambda g: aj(g, 2))
    exp_risk = {(a, t): int(((d.ARM == a) & (d.AVALM >= t)).sum()) for a in ["TREATMENT", "CONTROL"] for t in range(0, 43, 6)}
    risk_ok = len(rk) == len(exp_risk) and all(exp_risk[(r.arm, int(r.t))] == r.n for r in rk.itertuples())
    return {"CI01 curves are cumulative incidence (not 1-KM)": not is_km,
            "CI02 curves are for stated event": not is_death,
            "CI03 x axis in stated unit": meta["x_unit"] == meta["x_label"],
            "CI04 number-at-risk table present and correct": risk_ok}, is_cif

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean, is_cif = classify("CI_clean"); print("CLEAN:", clean, "| equals Aalen-Johansen:", is_cif)
assert all(clean.values()) and is_cif, "clean image fails a checker"
names = {"CI01": "CI01 curves are cumulative incidence (not 1-KM)", "CI02": "CI02 curves are for stated event",
         "CI03": "CI03 x axis in stated unit", "CI04": "CI04 number-at-risk table present and correct"}
ok = True
for e, nm in names.items():
    res, _ = classify(f"CI_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}CI_{e}.png") != md5(f"{OUT}CI_clean.png")
    good = failed == [nm] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
