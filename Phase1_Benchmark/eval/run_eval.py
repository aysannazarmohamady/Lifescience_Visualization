#!/usr/bin/env python3
"""One-file VLM evaluation of the ONCVIZ benchmark (60 images: 12 plot types x [1 clean + 4 flawed]). Standard library only.

  python3 run_eval.py            # uses ANTHROPIC_API_KEY and/or OPENAI_API_KEY from the env
  python3 run_eval.py --mock     # offline self-test, no keys, no network

Picks the newest Claude and GPT models from each provider's model list (override with
--claude-model / --openai-model). Runs two prompt modes (closed, open) per model, resumes
if interrupted, and writes results_*.jsonl and report.md. Keys are only read from the
environment and are never written anywhere.
DRAFT NOTE: drafted with LLM assistance."""
import argparse, base64, csv, collections, json, os, re, sys, time, urllib.request, urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
# Optional gateways/proxies (e.g. a company endpoint): set ANTHROPIC_BASE_URL / OPENAI_BASE_URL.
A_URL = os.environ.get("ANTHROPIC_BASE_URL", "https://api.anthropic.com").rstrip("/")
O_URL = os.environ.get("OPENAI_BASE_URL", "https://api.openai.com/v1").rstrip("/")
TAX = collections.defaultdict(list)
for r in csv.DictReader(open(os.path.join(HERE, "taxonomy.csv"))):
    TAX[r["plot_type"]].append((r["error_id"], r["description"]))

def http(url, headers, body=None, timeout=180):
    req = urllib.request.Request(url, data=None if body is None else json.dumps(body).encode(),
                                 headers={"content-type": "application/json", **headers})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r: return json.loads(r.read())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"HTTP {e.code}: {e.read().decode()[:300]}")

def newest_claude():
    j = http(A_URL + "/v1/models?limit=100",
             {"x-api-key": os.environ["ANTHROPIC_API_KEY"], "anthropic-version": "2023-06-01"})
    ids = [m["id"] for m in j["data"] if m["id"].startswith("claude")]   # API lists newest first
    big = [i for i in ids if not re.search(r"haiku", i)]                # prefer a flagship tier
    return (big or ids)[0]

def newest_gpt():
    j = http(O_URL + "/models", {"Authorization": "Bearer " + os.environ["OPENAI_API_KEY"]})
    cands = []
    for m in j["data"]:
        mt = re.fullmatch(r"gpt-(\d+(?:\.\d+)?)(-(?:pro|turbo))?", m["id"])
        if mt: cands.append((float(mt.group(1)), m["id"]))
    return max(cands)[1]

def prompt(plot, mode):
    base = (f"This figure is a '{plot}' from an oncology clinical-study dataset. Check whether the figure "
            "contains a visualization or statistical error (wrong axis, wrong denominator, mislabelled "
            "element, missing required element, etc.).\n")
    if mode == "closed":
        opts = "\n".join(f"- {i}: {d}" for i, d in TAX[plot])
        return base + f"At most one of these errors is present:\n{opts}\n- NONE: the figure is correct.\n" \
            'Answer with JSON only: {"has_error": true/false, "error_id": "<id or NONE>", "reason": "<one sentence>"}'
    return base + 'Answer with JSON only: {"has_error": true/false, "error_id": "NONE", "reason": "<what is wrong, one sentence>"}'

def parse(t):
    m = re.search(r"\{.*\}", t, re.S)
    try:
        j = json.loads(m.group(0)); return bool(j.get("has_error")), str(j.get("error_id", "NONE")), str(j.get("reason", ""))
    except Exception:
        return None, "PARSE_FAIL", t[:200]

def ask(provider, model, b64, text):
    if provider == "anthropic":
        j = http(A_URL + "/v1/messages",
                 {"x-api-key": os.environ["ANTHROPIC_API_KEY"], "anthropic-version": "2023-06-01"},
                 {"model": model, "max_tokens": 1500, "messages": [{"role": "user", "content": [
                  {"type": "image", "source": {"type": "base64", "media_type": "image/png", "data": b64}},
                  {"type": "text", "text": text}]}]})
        return "".join(b.get("text", "") for b in j["content"])
    if provider == "openai":
        j = http(O_URL + "/chat/completions", {"Authorization": "Bearer " + os.environ["OPENAI_API_KEY"]},
                 {"model": model, "max_completion_tokens": 16000, "messages": [{"role": "user", "content": [
                  {"type": "text", "text": text},
                  {"type": "image_url", "image_url": {"url": "data:image/png;base64," + b64}}]}]})
        return j["choices"][0]["message"].get("content") or ""
    return '{"has_error": false, "error_id": "NONE", "reason": "mock"}'

ap = argparse.ArgumentParser()
ap.add_argument("--mock", action="store_true"); ap.add_argument("--limit", type=int, default=0)
ap.add_argument("--claude-model"); ap.add_argument("--openai-model")
a = ap.parse_args()
runs = []
if a.mock: runs.append(("mock", "mock"))
else:
    if os.environ.get("ANTHROPIC_API_KEY"): runs.append(("anthropic", a.claude_model or newest_claude()))
    if os.environ.get("OPENAI_API_KEY"): runs.append(("openai", a.openai_model or newest_gpt()))
    if not runs: sys.exit("Set ANTHROPIC_API_KEY and/or OPENAI_API_KEY in the environment (or use --mock).")
print("Models:", runs)
rows = list(csv.DictReader(open(os.path.join(HERE, "manifest.csv"))))
if a.limit: rows = rows[:a.limit]
files = []
for prov, model in runs:
    for mode in ("closed", "open"):
        fn = os.path.join(HERE, f"results_{prov}_{model.replace('/', '_')}_{mode}.jsonl"); files.append((fn, prov, model, mode))
        done = {json.loads(l)["image"] for l in open(fn)} if os.path.exists(fn) else set()
        for r in rows:
            if r["image"] in done: continue
            b64 = base64.b64encode(open(os.path.join(HERE, r["image"]), "rb").read()).decode(); txt = ""
            for k in range(4):
                try: txt = ask(prov, model, b64, prompt(r["plot"], mode)); break
                except Exception as e: print("retry", r["image"], str(e)[:200], file=sys.stderr); time.sleep(5 * (k + 1))
            he, eid, why = parse(txt)
            with open(fn, "a") as f:
                f.write(json.dumps({**r, "provider": prov, "model": model, "mode": mode, "has_error": he,
                                    "error_id": eid, "reason": why, "raw": txt}) + "\n")
        print("done", fn)

lines = ["# ONCVIZ VLM evaluation report", "", "| model | mode | n | detection recall | false alarm on clean | exact error identified | parse failures |", "|---|---|---|---|---|---|---|"]
per_lines = []
for fn, prov, model, mode in files:
    rs = [json.loads(l) for l in open(fn)]
    fl = [r for r in rs if r["truth"] != "NONE"]; cl = [r for r in rs if r["truth"] == "NONE"]
    rec = sum(r["has_error"] is True for r in fl) / max(1, len(fl)); fpr = sum(r["has_error"] is True for r in cl) / max(1, len(cl))
    idn = sum(r["error_id"] == r["truth"] for r in fl) / max(1, len(fl)); pf = sum(r["error_id"] == "PARSE_FAIL" for r in rs)
    lines.append(f"| {model} | {mode} | {len(rs)} | {rec:.2f} | {fpr:.2f} | {(f"{idn:.2f}" if mode == "closed" else "n/a (judge reasons manually)")} | {pf} |")
    if mode == "closed":
        per = collections.defaultdict(list)
        for r in fl: per[r["plot"]].append(r["error_id"] == r["truth"])
        per_lines.append(f"\n**{model}, closed mode, identified per plot:** " + "; ".join(f"{p} {sum(v)}/{len(v)}" for p, v in per.items()))
open(os.path.join(HERE, "report.md"), "w").write("\n".join(lines + per_lines) + "\n")
print("\n".join(lines + per_lines))
