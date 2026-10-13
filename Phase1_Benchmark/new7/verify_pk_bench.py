"""Verify PK benchmark images from side files. Ground truth recomputed from CLEAN ADPK.
DRAFT NOTE: drafted with LLM assistance."""
import hashlib, numpy as np, pandas as pd

OUT = "Outputs/pk_bench/"
d = pd.read_csv("Data/V1/ADPK.csv"); d = d[(d.PARAMCD == "CONC") & (d.NOMTPT > 0)]
true = {(f"Cycle {c}", t): np.exp(np.log(g.AVAL[g.AVAL > 0]).mean())
        for (c, t), g in d.groupby(["AVISITN", "NOMTPT"])}

def check(tag):
    p = pd.read_csv(f"{OUT}{tag}_points.csv")
    meta = dict(l.strip().split(" ", 1) for l in open(f"{OUT}{tag}_meta.txt"))
    # match each point to the true value of its stated label at its time; time may be
    # wrong for PK03, so resolve by sorted order within the label
    times = sorted(d.NOMTPT.unique())
    res_geo, res_lab = True, True
    for lab, g in p.groupby("label"):
        g = g.sort_values("x")
        for t, y in zip(times, g.y):
            if abs(y - true[(lab, t)]) > 1e-6:
                other = "Cycle 3" if lab == "Cycle 1" else "Cycle 1"
                # wrong label if values equal the other cycle's, wrong statistic otherwise
                if abs(y - true[(other, t)]) < 1e-6: res_lab = False
                else: res_geo = False
    return {"PK01 log10 concentration axis": meta["y_scale"] == "log10",
            "PK02 geometric mean plotted": res_geo,
            "PK03 x axis in hours": bool(p.x.isin(times).all()),
            "PK04 cycle labels correct": res_lab}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean = check("PK_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
names = {"PK01": "PK01 log10 concentration axis", "PK02": "PK02 geometric mean plotted",
         "PK03": "PK03 x axis in hours", "PK04": "PK04 cycle labels correct"}
ok = True
for e, nm in names.items():
    res = check(f"PK_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}PK_{e}.png") != md5(f"{OUT}PK_clean.png")
    good = failed == [nm] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
