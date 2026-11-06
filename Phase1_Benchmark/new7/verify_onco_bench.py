"""Verify OncoPrint benchmark images from side files. Ground truth recomputed from
CLEAN ADMUT/ADSL. DRAFT NOTE: drafted with LLM assistance."""
import hashlib, pandas as pd

OUT = "Outputs/onco_bench/"
mut = pd.read_csv("Data/V1/ADMUT.csv"); adsl = pd.read_csv("Data/V1/ADSL.csv")
N = adsl.USUBJID.nunique()
pairs = mut[["USUBJID", "HUGO_SYMBOL"]].drop_duplicates()
true_n = pairs.groupby("HUGO_SYMBOL").size()
COLS = {"Truncating": "#222222", "Missense": "#2E7D32", "Amplification": "#C62828",
        "Deletion": "#1565C0", "Other": "#9E9E9E"}

def run(tag):
    cells = pd.read_csv(f"{OUT}{tag}_cells.csv"); top = pd.read_csv(f"{OUT}{tag}_topbar.csv")
    meta = dict(l.strip().split(" ", 1) for l in open(f"{OUT}{tag}_meta.txt"))
    freq = {k: int(v) for k, v in (kv.split("=") for kv in meta["freq"].split(","))}
    leg = dict(kv.split("=") for kv in meta["legend"].split(","))
    genes = meta["gene_order"].split(",")
    true_alt = pairs[pairs.HUGO_SYMBOL.isin(genes)].groupby("USUBJID").size()
    # memo-sort expectation: lexicographic presence in gene order, then id
    pres = pd.DataFrame({g: adsl.USUBJID.isin(pairs[pairs.HUGO_SYMBOL == g].USUBJID).astype(int) for g in genes})
    pres["id"] = adsl.USUBJID.values
    exp = pres.sort_values([*genes, "id"], ascending=[False] * len(genes) + [True]).id.tolist()
    return {"OP01 frequency uses all patients": all(freq[g] == round(100 * true_n[g] / N) for g in genes),
            "OP02 legend colours match cells": leg == COLS and all(COLS[c] == col for c, col in zip(cells.cat, cells.colour)),
            "OP03 samples sorted": top.sort_values("x")["sample"].tolist() == exp,
            "OP04 top-bar counts correct": all(int(r.top) == int(true_alt.get(r.sample, 0)) for r in top.itertuples())}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean = run("OP_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
names = {"OP01": "OP01 frequency uses all patients", "OP02": "OP02 legend colours match cells",
         "OP03": "OP03 samples sorted", "OP04": "OP04 top-bar counts correct"}
ok = True
for e, nm in names.items():
    res = run(f"OP_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}OP_{e}.png") != md5(f"{OUT}OP_clean.png")
    good = failed == [nm] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
