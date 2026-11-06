# VLM evaluation (seven new plot types, 35 images)

1. Create the images and manifest (repo root):
   `cd Data/V1 && unzip -o oncviz_data_v1.zip && cd ../..`
   `bash Phase1_Benchmark/new7/run_all.sh`
   `python3 Phase1_Benchmark/eval/build_manifest_new7.py`
2. Offline self-test (no key): `python3 Phase1_Benchmark/eval/eval_models.py --provider mock`
3. Real runs (key only in your shell, never in the repo or chat; `pip install requests`):
   `export ANTHROPIC_API_KEY=...` then `... --provider anthropic --model <latest Claude model id>`
   `export OPENAI_API_KEY=...` then `... --provider openai --model <latest GPT model id>`
   Add `--mode open` for the free-form variant; `--limit 5` for a cheap trial first.
4. Score: `python3 Phase1_Benchmark/eval/score_eval.py Phase1_Benchmark/eval/results_*.jsonl`

Runs are resumable and temperature 0 (if a model rejects `temperature`, remove it in the
call function). The five original plot types (variant 2 and originals) are not in the
manifest yet; add rows `image,plot,truth` and entries in `TAX` to include them.
Scripts are LLM-drafted; review before use.
