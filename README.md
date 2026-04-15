# mr-pipeline · v1.0.0

**Mendelian Randomization pipeline — Phase B of the CDC integrated genomics suite**  
Nadeem Khan · INRS-CAFSB · [github.com/nkhan119](https://github.com/nkhan119)

---

## Overview

Takes harmonised REGENIE summary statistics from `gwas-pipeline-v3` and runs full
Mendelian Randomization analysis across two cohorts in both directions.

| Step | Module | Description |
|------|--------|-------------|
| B1 | `BUILD_RSID_MAP` | rsID lookup from 1000G `.bim` files |
| B2 | `UNIVARIABLE_MR` | IVW · MR-Egger · Weighted Median |
| B3 | `HETEROGENEITY` | Cochran Q · I² · Egger intercept · Steiger · MR-PRESSO · LOO · Single-SNP |
| B4 | `MVMR` | MVMR-IVW · MVMR-Egger + conditional F-statistics |
| B5 | `MR_FIGURES` | 6 publication-quality figures (300 DPI PNG + PDF) |
| B6 | `MR_REPORT` | Interactive Plotly HTML report (8 tabs) |

---

## Requirements

### R packages (GitHub-only — not on Bioconda)

```r
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("mrcieu/ieugwasr")
remotes::install_github("WSpiller/MVMR")
remotes::install_github("rondolab/MR-PRESSO")
```

### Docker (recommended)

A `Dockerfile` is included that installs all R packages, `plink2`, and `plink 1.9`:

```bash
docker build -t nkhan119/mr-pipeline:1.0.0 .
docker push nkhan119/mr-pipeline:1.0.0  
```

---

## Directory structure

```
mr-pipeline-v1/
├── main.nf
├── nextflow.config
├── modules/
│   ├── build_rsid_map.nf
│   ├── univariable_mr.nf
│   ├── heterogeneity.nf
│   ├── mvmr.nf
│   ├── mr_figures.nf
│   └── mr_report.nf
├── bin/
│   ├── mr_figures.py
│   └── render_mr_report.py
├── assets/
│   └── mr_report_template.html
└── envs/
    └── mr_env.yaml
```

### Expected input layout (from `gwas-pipeline-v3`)

```
{base_dir}/GWAS_analysis/
  gwas/Cohort_A/sumstats/Cohort_A_<trait>_sumstats.tsv.gz
  gwas/Cohort_B/sumstats/Cohort_B_<trait>_sumstats.tsv.gz
  pca/1000G_dedup.{bim,bed,fam}
```

---

## Quick start

### 1. Set the base path in `nextflow.config`

```groovy
params.base_dir = "/path/to/1000G_phase3_common_norel"
```

Everything else derives from it automatically.

### 2. Run

| Environment | Command |
|-------------|---------|
| Local (Docker) | `nextflow run main.nf -profile local,docker -resume` |
| SLURM (Singularity) | `nextflow run main.nf -profile slurm,singularity -resume` |

Singularity pulls and caches the Docker Hub image automatically on first run.

---

## Skip flags

| Flag | Skips | Reads cached results from |
|------|-------|--------------------------|
| `--skip_uvmr` | B2 univariable MR | `results/uvmr/*_mr.tsv` |
| `--skip_het` | B3 heterogeneity | `results/heterogeneity/*_het.tsv` |
| `--skip_mvmr` | B4 MVMR | `results/mvmr/*_mvmr.tsv` |
| `--skip_report` | B5–B6 figures + HTML | — |

---

## Causal pairs (defaults)

**Univariable MR** — 5 pairs × 2 directions = 10 analyses

| Exposure | Outcome | p-threshold |
|----------|---------|-------------|
| LDL cholesterol | CAD | 5×10⁻⁸ |
| BMI | T2D | 5×10⁻⁸ |
| BMI | Hypertension | 5×10⁻⁸ |
| CRP (log) | CAD | 1×10⁻⁵ |
| LDL cholesterol | T2D | 5×10⁻⁸ |

**MVMR** — 3 pairs × 2 directions = 6 analyses

| Exposures | Outcome |
|-----------|---------|
| LDL + BMI | CAD |
| LDL + BMI | T2D |
| BMI + CRP | CAD |

Modify pairs in `nextflow.config` under `params.uvmr_pairs` and `params.mvmr_pairs`.

---

## Outputs

```
results/
  MR_Report.html
  figures/
    F1_forest_uvmr.{png,pdf}
    F2_heterogeneity.{png,pdf}
    F3_scatter_funnel.{png,pdf}
    F4_leave_one_out.{png,pdf}
    F5_single_snp.{png,pdf}
    F6_mvmr_forest.{png,pdf}
  uvmr/
    *_mr.tsv
    *_wald.tsv
  heterogeneity/
    *_het.tsv
    *_loo.tsv
    *_single.tsv
  mvmr/
    *_mvmr.tsv
  refs/
    rsid_map.tsv.gz
  logs/
    trace_*.txt
    nf_report_*.html
```
