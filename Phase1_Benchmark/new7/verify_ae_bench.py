"""Verify AE bar-chart benchmark images from the bar data written beside each figure.
Ground truth recomputed from CLEAN ADAE/ADSL. DRAFT NOTE: drafted with LLM assistance."""
import hashlib, numpy as np, pandas as pd

OUT = "Outputs/ae_bench/"
adae = pd.read_csv("Data/V1/ADAE.csv"); adsl = pd.read_csv("Data/V1/ADSL.csv")
n_arm = int((adsl.ARM == "TREATMENT").sum())
w = adae[adae.ARM == "TREATMENT"].groupby(["USUBJID", "AEPT"]).AETOXGR.max().reset_index()
true = w.groupby(["AEPT", "AETOXGR"]).size().rename("n").reset_index()
true["pct"] = 100 * true.n / n_arm
true_tot = true.groupby("AEPT").pct.sum()
true_top = set(true_tot.sort_values(ascending=False).index[:12])

def lum(h):
    h = h.lstrip("#"); r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4)); return 0.3 * r + 0.59 * g + 0.11 * b

def run(tag):
    b = pd.read_csv(f"{OUT}{tag}_bars.csv")
    meta = dict(l.strip().split(" ", 1) for l in open(f"{OUT}{tag}_meta.txt"))
    cols = dict(kv.split("=") for kv in meta["colours"].split(","))
    lums = [lum(cols[str(g)]) for g in range(1, 6)]
    shown = b.groupby("term").pct.sum()
    t = b.groupby("term").order.first()
    present = set(zip(b.term, b.grade)); terms = set(b.term)
    need = [(r.AEPT, r.AETOXGR) for r in true.itertuples() if r.AEPT in terms]
    tp = {(r.AEPT, r.AETOXGR): r.pct for r in true.itertuples()}
    by_freq = list(t.sort_values().index)               # bottom -> top
    exp_order = list(shown.sort_values().index)
    return {"TX01 grade colours darken with grade": all(x > y for x, y in zip(lums, lums[1:])),
            "TX02 all grades incl. 4-5 shown": all((t, g) in present for t, g in need),
            "TX03 percent uses arm N": all(abs(r.pct - tp.get((r.term, r.grade), -1)) < 1e-6
                                           for r in b.itertuples()),
            "TX04 terms ordered by frequency": bool(np.all(np.diff(shown[by_freq].values) >= -1e-9))}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean = run("TX_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
names = {"TX01": "TX01 grade colours darken with grade", "TX02": "TX02 all grades incl. 4-5 shown",
         "TX03": "TX03 percent uses arm N", "TX04": "TX04 terms ordered by frequency"}
ok = True
for e, nm in names.items():
    res = run(f"TX_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}TX_{e}.png") != md5(f"{OUT}TX_clean.png")
    good = failed == [nm] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
