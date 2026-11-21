# ONCVIZ final evaluation run (reported in the manuscript)

Frozen 60-image test set: 12 plot types x (1 clean + 4 flawed) = 60 images, 48 injected errors (`manifest.csv`; error definitions in `taxonomy.csv`).

| item | value |
|---|---|
| Model identifier | `us.anthropic.claude-opus-5-5` |
| Run | October 2026, after renderers, verifiers, prompts and scorer were frozen |
| Input | rendered PNG only, plus the fixed prompt (no code, data or label) |
| Prompts | `closed` (four candidate error codes for the plot type) and `open` (no candidate list) |
| Responses | one per image per prompt (120 in total) |
| Sampling settings | `max_tokens` = 1500; no temperature set (provider default) |
| Unparseable responses | counted as failures, not repaired (closed: OP02, OP04; open: OP04) |

Files
- `results_*_closed.jsonl`, `results_*_open.jsonl`: raw model responses (`raw`) with parsed fields.
- `open_judgement.csv`: manual judgment of each open-ended response against the taxonomy definition (`correct` = 1/0, with a note).
- `run_eval.py`: the script that produced the responses (reads keys from the environment).

Metrics reported in the paper are recomputed directly from these files.
