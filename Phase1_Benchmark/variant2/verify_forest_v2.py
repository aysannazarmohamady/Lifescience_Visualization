"""Verify the variant-2 forest images: the clean image must pass every check
and each flawed image must fail exactly the check of its own error. Ground
truth for subgroup n is recomputed from ADTTE. Run from the repository root
after run_forest_v2.R. DRAFT NOTE: drafted with LLM assistance."""
import hashlib, sys
import pandas as pd
sys.path.insert(0, "Phase1_Benchmark/checkers")
import checkers as c

OUT = "Outputs/forest_v2/"
d = pd.read_csv("Data/V1/ADTTE.csv"); d = d[d.PARAMCD == "PFS"]
masks = {"All patients": d.index == d.index, "NSCLC": d.TUMORTYPE == "NSCLC", "BRCA": d.TUMORTYPE == "BRCA",
 "HCC": d.TUMORTYPE == "HCC", "CRC": d.TUMORTYPE == "CRC", "PDAC": d.TUMORTYPE == "PDAC",
 "Age <65": d.AGEGR1 == "<65", "Age >=65": d.AGEGR1 == ">=65", "Male": d.SEX == "M", "Female": d.SEX == "F",
 "ECOG 0": d.ECOG == 0, "ECOG >=1": d.ECOG != 0, "Stage I-III": d.STAGE != "STAGE IV", "Stage IV": d.STAGE == "STAGE IV",
 "Liver mets: yes": d.LIVERMETS == "Y", "Liver mets: no": d.LIVERMETS == "N", "TMB-High": d.TMBHIGH == "Y",
 "TMB-Low": d.TMBHIGH == "N", "PD-L1 High": d.PDL1GRP == "HIGH", "PD-L1 not High": d.PDL1GRP != "HIGH"}
true_n = {k: d[m].USUBJID.nunique() for k, m in masks.items()}

def run(tag):
    t = pd.read_csv(f"{OUT}{tag}_table.csv")
    meta = dict(l.split() for l in open(f"{OUT}{tag}_meta.txt"))
    f = t.rename(columns={"lo_c": "lo_c", "hi_c": "hi_c", "n_shown": "n"})
    tn = pd.Series([true_n[l] for l in t.label])
    return {"FP01 CI contains point": c.check_forest_ci_contains_point(f),
            "FP02 null line at 1": c.check_forest_null_line_correct(float(meta["null_x"])),
            "FP03 log axis": c.check_forest_axis_is_log(meta["scale"]),
            "FP04 n matches data": c.check_forest_subgroup_n_matches_data(f, tn)}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean = run("FP_v2_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
ok = True
for e, name in {"FP01": "FP01 CI contains point", "FP02": "FP02 null line at 1",
                "FP03": "FP03 log axis", "FP04": "FP04 n matches data"}.items():
    res = run(f"FP_v2_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}FP_v2_{e}.png") != md5(f"{OUT}FP_v2_clean.png")
    good = failed == [name] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
