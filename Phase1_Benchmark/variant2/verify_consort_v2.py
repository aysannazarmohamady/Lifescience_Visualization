"""Verify the variant-2 CONSORT images. Numbers are read back from the text
written in the figure's boxes (<name>_boxes.csv) and compared with ground
truth recomputed from the CLEAN ADRAND. The clean image must pass every check;
each flawed image must fail exactly the check of its own error. Run from the
repository root after run_consort_v2.R. DRAFT NOTE: drafted with LLM assistance."""
import hashlib, re, sys
import pandas as pd
sys.path.insert(0, "Phase1_Benchmark/checkers")
import checkers as c

OUT = "Outputs/consort_v2/"
adrand = pd.read_csv("Data/V1/ADRAND.csv")
true_reasons = adrand.SFREASN.dropna().unique().tolist()
true_arm = adrand[adrand.RANDFL == "Y"].ARM.value_counts().to_dict()

def first_n(text): return int(re.search(r"n = (\d+)", text).group(1))

def run(tag):
    b = dict(pd.read_csv(f"{OUT}{tag}_boxes.csv").values)
    n_screen, n_excl, n_rand = (first_n(b[k]) for k in ("assessed", "excluded", "randomized"))
    arms = {"TREATMENT": first_n(b["alloc_TREATMENT"]), "CONTROL": first_n(b["alloc_CONTROL"])}
    return {"CN01 screened = randomized + excluded": c.check_consort_stage_arithmetic(n_screen, n_excl, n_rand),
            "CN02 exclusion reasons shown": (c.check_consort_exclusion_reasons_present(b["excluded"])
                                              and all(r in b["excluded"] for r in true_reasons)),
            "CN03 arm counts match randomization": (c.check_consort_arm_totals_match_randomization(n_rand, arms)
                                                   and arms == true_arm),
            "CN04 discontinuation detail shown": all(c.check_consort_followup_reported(b[k])
                                                    for k in ("disc_TREATMENT", "disc_CONTROL"))}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean = run("CN_v2_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
ok = True
for e, name in {"CN01": "CN01 screened = randomized + excluded", "CN02": "CN02 exclusion reasons shown",
                "CN03": "CN03 arm counts match randomization", "CN04": "CN04 discontinuation detail shown"}.items():
    res = run(f"CN_v2_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}CN_v2_{e}.png") != md5(f"{OUT}CN_v2_clean.png")
    good = failed == [name] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
