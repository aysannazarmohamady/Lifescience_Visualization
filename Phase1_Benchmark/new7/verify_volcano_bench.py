"""Verify volcano benchmark images from plotted data written beside each figure.
Ground truth recomputed from CLEAN ADBM (t-test of per-patient log2 ratio, BH).
DRAFT NOTE: drafted with LLM assistance."""
import hashlib, numpy as np, pandas as pd
from scipy import stats

OUT = "Outputs/volcano_bench/"
d = pd.read_csv("Data/V1/ADBM.csv")
d = d[(d.AVISIT != "BASELINE") & (d.AVAL > 0) & (d.BASE > 0)].copy()
d["l2"] = np.log2(d.AVAL / d.BASE)
d["feature"] = d.ARM.str[:3] + " " + d.PARAMCD
rows = []
for f, g in d.groupby("feature"):
    pat = g.groupby("USUBJID").l2.mean()
    if len(pat) >= 10:
        rows.append((f, pat.mean(), stats.ttest_1samp(pat, 0).pvalue))
t = pd.DataFrame(rows, columns=["feature", "log2fc", "p"])
t["sig"] = t.p < 0.05; t = t.set_index("feature")

def run(tag):
    pts = pd.read_csv(f"{OUT}{tag}_points.csv").set_index("feature").join(t)
    meta = {l.split(" ", 1)[0]: l.split(" ", 1)[1].strip() for l in open(f"{OUT}{tag}_meta.txt")}
    return {"VP01 threshold line at alpha": abs(float(meta["thr_y"]) + np.log10(0.05)) < 1e-6,
            "VP02 y is -log10(p)": np.allclose(pts.y, -np.log10(pts.p), atol=1e-6),
            "VP03 x is log2 fold change": np.allclose(pts.x, pts.log2fc, atol=1e-6),
            "VP04 highlight only significant": bool((pts.highlighted == pts.sig).all())}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
print("significant features:", int(t.sig.sum()), "of", len(t))
clean = run("VP_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
names = {"VP01": "VP01 threshold line at alpha", "VP02": "VP02 y is -log10(p)",
         "VP03": "VP03 x is log2 fold change", "VP04": "VP04 highlight only significant"}
ok = True
for e, nm in names.items():
    res = run(f"VP_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}VP_{e}.png") != md5(f"{OUT}VP_clean.png")
    good = failed == [nm] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
