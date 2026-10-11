---

editor_options: 
  markdown: 
    wrap: 72
---

# 16S rRNA analysis of diarrheic piglets 

Code for the 16S rRNA (V3–V4) part of the manuscript *"Integrated culture, microbiome profiling, and whole-genome sequencing reveal multilevel bacterial heterogeneity in diarrheic piglets in Kazakhstan"* (34 samples; age groups 14, 24 and 30 days).

## Layout

| Path | Purpose |
|------------------------------------|------------------------------------|
| `scripts/01_qiime2_dada2.sh` | QIIME 2 import → DADA2 → SILVA taxonomy → **chloroplast/mitochondria removal** → export of tables and level-6 counts |
| `scripts/02_alpha_beta.R` | rarefied alpha diversity, Bray–Curtis / Jaccard / Aitchison, PCoA, NMDS, PERMANOVA, PERMDISP, genus barplot, top-50 heatmap |
| `scripts/03_dmm.R` | Dirichlet-multinomial mixture models on genus-resolved taxa (top-15/20/30), ARI sensitivity, DMM vs age (Fisher) |
| `scripts/04_ancombc2.R` | ANCOM-BC2 (global + pairwise) at genus level |
| `config.R` | all paths, seeds, rarefaction depth, helper functions |
| `data/` | local inputs (feature tables, taxonomy, level-6, `metadata_pigs.csv`) — **git-ignored until publication**, see `data/README.md` |
| `results/` | script outputs (git-ignored) |

## Order of execution

``` bash
# 1. QIIME 2 (conda env qiime2-amplicon-2024.10), needs raw FASTQ + manifest
MANIFEST=/path/to/manifest.tsv CLASSIFIER=/path/silva-138.2-v3v4-classifier.qza bash scripts/01_qiime2_dada2.sh

# 2. R (run from the repository root)
Rscript scripts/install_R_packages.R
PIGLET_VARIANT=noorg bash run_all.sh      # corrected analysis (organelles removed)
PIGLET_VARIANT=all   bash run_all.sh      # first-submission table (organelles retained)
```

Step 01 writes `data/all/` (no organelle filter) and `data/noorg/` (chloroplast + mitochondria removed). Every R script reads `data/<variant>/` and writes to `results/<variant>/`.

## DADA2 / QIIME 2 parameters

| Step | Setting |
|------------------------------------|------------------------------------|
| Import | `PairedEndFastqManifestPhred33V2` |
| Primer removal | inside DADA2: `--p-trim-left-f 17 --p-trim-left-r 21` (no cutadapt) |
| Truncation | `--p-trunc-len-f 280 --p-trunc-len-r 220`, `--p-trunc-q 2` |
| Filtering | `--p-max-ee-f 2 --p-max-ee-r 2`, `--p-min-overlap 12` |
| Chimeras / pooling | `consensus` / `independent`, 8 threads |
| Taxonomy | `classify-sklearn`, SILVA 138.2 V3–V4 classifier, `--p-confidence 0.7` |
| Organelles | `taxa filter-table --p-exclude mitochondria,chloroplast` |
| Rarefaction | 4,871 reads/sample, 100 iterations (alpha diversity only) |

Read retention (unfiltered run): mean 70.3% filtered, 43.7% merged, 29.4% non-chimeric; final depth 4,871–27,763 reads/sample.

## License

MIT (see `LICENSE`)
