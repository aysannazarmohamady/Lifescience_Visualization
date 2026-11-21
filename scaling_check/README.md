# Scaling check at 1,000 patients

This folder documents the trial-scale check reported in the ONCVIZ manuscript. The check tests the synthetic-data generator, the 12 renderers and the 12 verifiers. It does not test the vision-language model, which was not run on large-sample images.

| item | value |
|---|---|
| Patients | 1,000 (683 TREATMENT, 317 CONTROL; counted in `ADSL.csv`) |
| Dataset | 13 ADaM domains, in `oncviz001_n1000_data.zip` |
| Generator | `Data/V1/generate_adam_synthetic.R`, run with larger patient-count parameters (`n_phase1`, `n_phase2_trt`, `n_phase2_ctl`) |
| Renderers and verifiers | all 12 ran; all 12 verifiers returned OK after one fix to the PK verifier's handling of cycle labels |
| Swimmer plot | one row per patient; became unreadable at roughly 683 treatment rows |

The 80-patient catalogue is kept for visual review because patient-level figures (swimmer, spider, waterfall) must remain legible. The 1,000-patient version serves as a stress test for the generator and the verifiers.

## Reproduce

1. Make a scratch copy of the repository (do not overwrite the released 80-patient files).
2. Extract `oncviz001_n1000_data.zip` into `Data/V1/` of the scratch copy, replacing the CSVs there.
3. Run `bash run12.sh` from the scratch copy's root. It renders and verifies each of the 12 plot types and prints one line per plot.

Requires R with the packages listed in the root `README.md`, and Python 3.
