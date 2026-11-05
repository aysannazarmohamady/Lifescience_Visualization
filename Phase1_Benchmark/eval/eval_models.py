"""Ask a vision-language model to find the error in each benchmark image.

  export ANTHROPIC_API_KEY=...   # or OPENAI_API_KEY=...   (never commit keys)
  python3 Phase1_Benchmark/eval/eval_models.py --provider anthropic --model <model-id>
  python3 Phase1_Benchmark/eval/eval_models.py --provider openai    --model <model-id>
  python3 Phase1_Benchmark/eval/eval_models.py --provider mock      # offline self-test

Results are appended to Phase1_Benchmark/eval/results_<model>.jsonl (resumable).
Temperature 0. Two prompt modes: 'closed' lists the plot type's candidate errors,
'open' only asks whether the figure has a statistical/visual error and why.
DRAFT NOTE: drafted with LLM assistance."""
import argparse, base64, csv, json, os, re, sys, time
import requests

TAX = {}
for r in csv.DictReader(open("Phase1_Benchmark/new7/taxonomy_new7.csv")):
    TAX.setdefault(r["plot_type"], []).append((r["error_id"], r["description"]))

def prompt(plot, mode):
    base = (f"This figure is a '{plot}' from an oncology clinical-study dataset. "
            "Check whether the figure contains a visualization or statistical error "
            "(wrong axis, wrong denominator, mislabelled element, missing required element, etc.).\n")
    if mode == "closed":
        opts = "\n".join(f"- {i}: {d}" for i, d in TAX[plot])
        return base + f"At most one of these errors is present:\n{opts}\n- NONE: the figure is correct.\n" \
            'Answer with JSON only: {"has_error": true/false, "error_id": "<id or NONE>", "reason": "<one sentence>"}'
    return base + 'Answer with JSON only: {"has_error": true/false, "error_id": "NONE", "reason": "<what is wrong, one sentence>"}'

def parse(text):
    m = re.search(r"\{.*\}", text, re.S)
    try:
        j = json.loads(m.group(0)); return bool(j.get("has_error")), str(j.get("error_id", "NONE")), str(j.get("reason", ""))
    except Exception:
        return None, "PARSE_FAIL", text[:200]

def call_anthropic(model, img_b64, text):
    r = requests.post("https://api.anthropic.com/v1/messages", timeout=120,
        headers={"x-api-key": os.environ["ANTHROPIC_API_KEY"], "anthropic-version": "2023-06-01"},
        json={"model": model, "max_tokens": 400, "temperature": 0, "messages": [{"role": "user", "content": [
            {"type": "image", "source": {"type": "base64", "media_type": "image/png", "data": img_b64}},
            {"type": "text", "text": text}]}]})
    r.raise_for_status(); return r.json()["content"][0]["text"]

def call_openai(model, img_b64, text):
    body = {"model": model, "messages": [{"role": "user", "content": [
        {"type": "text", "text": text},
        {"type": "image_url", "image_url": {"url": f"data:image/png;base64,{img_b64}"}}]}]}
    r = requests.post("https://api.openai.com/v1/chat/completions", timeout=120,
        headers={"Authorization": "Bearer " + os.environ["OPENAI_API_KEY"]}, json=body)
    if r.status_code == 400 and "temperature" in r.text: pass
    r.raise_for_status(); return r.json()["choices"][0]["message"]["content"]

def call_mock(model, img_b64, text):
    return '{"has_error": false, "error_id": "NONE", "reason": "mock"}'

CALLS = {"anthropic": call_anthropic, "openai": call_openai, "mock": call_mock}

ap = argparse.ArgumentParser()
ap.add_argument("--provider", required=True, choices=CALLS); ap.add_argument("--model", default="mock")
ap.add_argument("--mode", default="closed", choices=["closed", "open"])
ap.add_argument("--manifest", default="Phase1_Benchmark/eval/manifest.csv")
ap.add_argument("--limit", type=int, default=0, help="only first N images (cheap trial)")
a = ap.parse_args()
out = f"Phase1_Benchmark/eval/results_{a.provider}_{a.model.replace('/', '_')}_{a.mode}.jsonl"
done = set()
if os.path.exists(out):
    done = {json.loads(l)["image"] for l in open(out)}
rows = list(csv.DictReader(open(a.manifest)))
if a.limit: rows = rows[:a.limit]
for r in rows:
    if r["image"] in done: continue
    b64 = base64.b64encode(open(r["image"], "rb").read()).decode()
    for attempt in range(4):
        try: txt = CALLS[a.provider](a.model, b64, prompt(r["plot"], a.mode)); break
        except Exception as e:
            print("retry", r["image"], e, file=sys.stderr); time.sleep(5 * (attempt + 1)); txt = ""
    he, eid, why = parse(txt)
    with open(out, "a") as f:
        f.write(json.dumps({**r, "mode": a.mode, "model": a.model, "has_error": he, "error_id": eid,
                            "reason": why, "raw": txt}) + "\n")
print("written", out)
