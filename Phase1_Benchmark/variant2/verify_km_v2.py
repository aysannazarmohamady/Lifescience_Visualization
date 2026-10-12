"""Verify the variant-2 Kaplan-Meier images from the plotted data written next
to each figure. Ground truth: censored-patient counts from the CLEAN ADTTE, and
risk sets recomputed from the data each figure was drawn from. The clean image
must pass every check; each flawed image must fail exactly the check of its own
error. Run from the repository root after run_km_v2.R.
DRAFT NOTE: drafted with LLM assistance."""
import hashlib, sys
import pandas as pd
sys.path.insert(0, "Phase1_Benchmark/checkers")
import checkers as c

OUT = "Outputs/km_v2/"; MONTH = 30.4375
clean_os = pd.read_csv("Data/V1/ADTTE.csv"); clean_os = clean_os[clean_os.PARAMCD == "OS"]
true_cens = clean_os.groupby("ARM").CNSR.sum().to_dict()

def run(tag):
    curve = pd.read_csv(f"{OUT}{tag}_curve.csv"); marks = pd.read_csv(f"{OUT}{tag}_marks.csv")
    tab = pd.read_csv(f"{OUT}{tag}_risktable.csv"); used = pd.read_csv(f"{OUT}{tag}_data_used.csv")
    mono = all((g.surv.diff().dropna() <= 1e-12).all() for _, g in curve.groupby("arm"))
    start = all(c.check_km_starts_at_one(g.sort_values("t").surv.iloc[0]) for _, g in curve.groupby("arm"))
    # one tick per censored patient (ties at the same time share a tick position)
    n_marks = {a: int((marks.arm == a).sum()) for a in true_cens}
    cens_ok = all((n_marks[a] > 0) == (true_cens[a] > 0) and
                  (n_marks[a] == clean_os[(clean_os.ARM == a) & (clean_os.CNSR == 1)].AVAL.nunique())
                  for a in true_cens)
    true_rs = [int((used[used.ARM == r.arm].AVAL / MONTH >= r.t).sum()) for r in tab.itertuples()]
    risk_ok = c.check_km_at_risk_matches_data(tab.n_risk.tolist(), true_rs)
    return {"KM01 curve non-increasing": mono, "KM02 censor marks present": cens_ok,
            "KM03 curve starts at 1.0": start, "KM04 at-risk table matches data": risk_ok}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean = run("KM_v2_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
ok = True
for e, name in {"KM01": "KM01 curve non-increasing", "KM02": "KM02 censor marks present",
                "KM03": "KM03 curve starts at 1.0", "KM04": "KM04 at-risk table matches data"}.items():
    res = run(f"KM_v2_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}KM_v2_{e}.png") != md5(f"{OUT}KM_v2_clean.png")
    good = failed == [name] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
