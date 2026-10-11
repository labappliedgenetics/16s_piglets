#!/usr/bin/env bash
# Full pipeline. Step 01 needs the QIIME 2 conda env and the raw reads; steps 02-04 need only R.
set -euo pipefail
# bash scripts/01_qiime2_dada2.sh            # (QIIME 2 env, raw FASTQ + manifest)
export PIGLET_VARIANT="${PIGLET_VARIANT:-noorg}"   # use "all" to reproduce the first-submission table
Rscript scripts/02_alpha_beta.R
Rscript scripts/03_dmm.R
Rscript scripts/04_ancombc2.R
Rscript scripts/session_info.R
