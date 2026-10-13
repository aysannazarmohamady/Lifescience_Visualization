"""Verify the spider benchmark images from the plotted data written next to each
figure. Ground truth is recomputed from the CLEAN ADTR/ADRS/ADSL. The clean
image must pass every check; each flawed image must fail exactly its own check.
DRAFT NOTE: drafted with LLM assistance."""
import hashlib
import pandas as pd

OUT = "Outputs/spider_bench/"
adtr = pd.read_csv("Data/V1/ADTR.csv"); adrs = pd.read_csv("Data/V1/ADRS.csv"); adsl = pd.read_csv("Data/V1/ADSL.csv")
adtr["ADT"] = pd.to_datetime(adtr.ADT)
true_day = {(r.USUBJID, r.ADT): (r.ADT - adtr[adtr.USUBJID == r.USUBJID].ADT.min()).days
            for r in adtr.itertuples()}
pd_date = (adrs[(adrs.PARAMCD == "OVRLRESP") & (adrs.AVALC == "PD") & adrs.ADT.notna()]
           .assign(ADT=lambda d: pd.to_datetime(d.ADT)).groupby("USUBJID").ADT.min())
death = adsl[adsl.OSCR == 1].set_index("USUBJID").OSDTC.pipe(pd.to_datetime)
stop = pd.concat([pd_date, death], axis=1).min(axis=1)

def run(tag):
    s = pd.read_csv(f"{OUT}{tag}_series.csv", parse_dates=["date"])
    meta = {l.split(" ", 1)[0]: l.split(" ", 1)[1].strip() for l in open(f"{OUT}{tag}_meta.txt")}
    first = s.sort_values("date").groupby("USUBJID").pchg.first()
    last_date = s.groupby("USUBJID").date.max()
    over = [u for u, d in last_date.items() if u in stop.index and d > stop[u]]
    unit_ok = all(abs(true_day[(r.USUBJID, r.date)] - r.day) < 1e-6 for r in s.itertuples())
    return {"SP01 -30/+20 reference lines": meta["ref_lines"] == "-30,20",
            "SP02 trajectories start at 0": bool((first.abs() < 1e-9).all()),
            "SP03 lines stop at PD or death": len(over) == 0,
            "SP04 uniform time units": unit_ok}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean = run("SP_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
ok = True
for e, name in {"SP01": "SP01 -30/+20 reference lines", "SP02": "SP02 trajectories start at 0",
                "SP03": "SP03 lines stop at PD or death", "SP04": "SP04 uniform time units"}.items():
    res = run(f"SP_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}SP_{e}.png") != md5(f"{OUT}SP_clean.png")
    good = failed == [name] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
