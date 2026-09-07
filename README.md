# CNV Pipeline — MMRF CoMMpass (GDC)

A GitHub-ready Jupyter notebook pipeline to **process**, **QC**, and **analyze** Copy Number Variation (CNV) segments from the **MMRF CoMMpass** cohort via the NCI **Genomic Data Commons (GDC)**.

This repository focuses on *downstream* CNV segment analytics (recurrence summaries, cytoband / fixed-bin aggregation, clinical integration, and survival/staging evaluation). The pipeline is designed to run in **Google Colab** or locally.

> Updated notebook version: **V7** (2026-02-05)

## Repository contents

- `CNV_MMRF_COMMPASS_Pipeline_GitHub_Documented.ipynb`  
  Main notebook (outputs stripped to keep diffs clean).

- `CNV_MMRF_COMMPASS_Pipeline_GitHub_Documented_NO_OUTPUTS.ipynb`  
  Same notebook, explicitly output-free (recommended for version control).

- `requirements.txt`  
  Best-effort dependency list for local execution.

- `docs/PIPELINE_OVERVIEW.md`  
  High-level description of stages, inputs, and outputs.

## Quickstart (local)

```bash
python -m venv .venv
# Linux/Mac:
source .venv/bin/activate
# Windows (PowerShell):
# .venv\Scripts\Activate.ps1

pip install -r requirements.txt
jupyter lab
```

Open the notebook and run it from top to bottom.

## Quickstart (Google Colab)

Upload the notebook to Colab and run top-to-bottom.  
If you store inputs in Google Drive, mount Drive and point the notebook to your `.tsv` inputs.

## Outputs

All artifacts are written under:

- `outputs/run_<RUN_ID>/raw/`
- `outputs/run_<RUN_ID>/processed/`
- `outputs/run_<RUN_ID>/results/`
- `outputs/run_<RUN_ID>/logs/`

The notebook ends with a concise list of deliverables and their paths.

## Notes / limitations

- **Exact-breakpoint recurrence** (`SegmentID = chr:start-end`) is included for comparability, but breakpoint variation can underestimate biological recurrence.
- **Cytoband overlap** and **fixed 1Mb bins** provide more stable recurrence summaries.
- When working with CNVs already called upstream, **purity/ploidy adjustment may be unavailable**; the pipeline is explicit about this limitation.

## License

MIT (see `LICENSE`).
