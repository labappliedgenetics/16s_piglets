#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════
# 01 — 16S rRNA (V3–V4) processing: QIIME 2 + DADA2 + SILVA, with
#      chloroplast / mitochondria removal.
#
# Usage:   conda activate qiime2-amplicon-2024.10
#          bash scripts/01_qiime2_dada2.sh
#
# Inputs : $MANIFEST      QIIME 2 paired-end manifest (sample-id, forward-absolute-filepath,
#                         reverse-absolute-filepath); raw FASTQ are NOT stored in this repo
#          $CLASSIFIER    SILVA 138.2 V3–V4 naive Bayes classifier (.qza)   [CONFIRM origin]
# Outputs: data/all/    tables WITHOUT organelle filtering (= first-submission analysis)
#          data/noorg/  tables AFTER removing chloroplast + mitochondria (recommended)
#          Each folder holds: feature-table.tsv, level-6.csv, taxonomy.tsv, dna-sequences.fasta
#          data/stats.tsv (DADA2 read tracking)
#
# Primers are removed INSIDE DADA2 (--p-trim-left-*), so no separate cutadapt step is run.
# Primer presence before trimming was checked separately (README: "Primer check").
# ══════════════════════════════════════════════════════════════════════════
set -euo pipefail

# ── settings (edit) ──────────────────────────────────────────────────────
MANIFEST="${MANIFEST:-manifest.tsv}"
CLASSIFIER="${CLASSIFIER:-silva-138.2-v3v4-classifier.qza}"   # [CONFIRM] trained in-house vs downloaded
WORK="${WORK:-qiime2_work}"
THREADS="${THREADS:-8}"

TRIM_LEFT_F=17;  TRIM_LEFT_R=21        # primer lengths removed by DADA2
TRUNC_LEN_F=280; TRUNC_LEN_R=220
MAX_EE_F=2;      MAX_EE_R=2
TRUNC_Q=2;       MIN_OVERLAP=12
CONFIDENCE=0.7                          # classify-sklearn
# ─────────────────────────────────────────────────────────────────────────
mkdir -p "$WORK" data/all data/noorg
exportdir="$WORK/export"; mkdir -p "$exportdir"

echo "[1/7] Import"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path "$WORK/demux.qza"

qiime demux summarize --i-data "$WORK/demux.qza" --o-visualization "$WORK/demux.qzv"

echo "[2/7] DADA2 denoise-paired (primers trimmed, chimeras removed)"
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$WORK/demux.qza" \
  --p-trim-left-f $TRIM_LEFT_F --p-trim-left-r $TRIM_LEFT_R \
  --p-trunc-len-f $TRUNC_LEN_F --p-trunc-len-r $TRUNC_LEN_R \
  --p-max-ee-f $MAX_EE_F --p-max-ee-r $MAX_EE_R \
  --p-trunc-q $TRUNC_Q --p-min-overlap $MIN_OVERLAP \
  --p-chimera-method consensus --p-pooling-method independent \
  --p-n-threads "$THREADS" \
  --p-hashed-feature-ids \
  --o-table "$WORK/table.qza" \
  --o-representative-sequences "$WORK/rep-seqs.qza" \
  --o-denoising-stats "$WORK/stats.qza"

qiime tools export --input-path "$WORK/stats.qza" --output-path "$exportdir/stats"
tr -d '\r' < "$exportdir/stats/stats.tsv" > data/stats.tsv

echo "[3/7] Taxonomy (SILVA 138.2, confidence $CONFIDENCE)"
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$WORK/rep-seqs.qza" \
  --p-confidence $CONFIDENCE \
  --p-n-jobs "$THREADS" \
  --o-classification "$WORK/taxonomy.qza"

echo "[4/7] Chloroplast / mitochondria removal"
# SILVA labels: o__Chloroplast, f__Mitochondria. Matching is case-insensitive substring.
qiime taxa filter-table \
  --i-table "$WORK/table.qza" --i-taxonomy "$WORK/taxonomy.qza" \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table "$WORK/table-noorg.qza"
qiime taxa filter-seqs \
  --i-sequences "$WORK/rep-seqs.qza" --i-taxonomy "$WORK/taxonomy.qza" \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences "$WORK/rep-seqs-noorg.qza"

# per-sample summary BEFORE/AFTER (report ASVs and reads removed in the manuscript)
qiime feature-table summarize --i-table "$WORK/table.qza"       --o-visualization "$WORK/table.qzv"
qiime feature-table summarize --i-table "$WORK/table-noorg.qza" --o-visualization "$WORK/table-noorg.qzv"

# Optional stricter filter (uncomment to also drop unassigned / archaeal / eukaryotic ASVs):
# qiime taxa filter-table --i-table "$WORK/table-noorg.qza" --i-taxonomy "$WORK/taxonomy.qza" \
#   --p-include d__Bacteria --o-filtered-table "$WORK/table-noorg.qza"

echo "[5/7] Export tables, sequences, taxonomy"
export_variant () {   # $1 = all|noorg   $2 = table.qza   $3 = rep-seqs.qza
  local v="$1" t="$2" r="$3" out="data/$1"
  qiime tools export --input-path "$t" --output-path "$exportdir/table-$v"
  biom convert -i "$exportdir/table-$v/feature-table.biom" \
               -o "$out/feature-table.tsv" --to-tsv
  qiime tools export --input-path "$r" --output-path "$exportdir/seqs-$v"
  cp "$exportdir/seqs-$v/dna-sequences.fasta" "$out/dna-sequences.fasta"
  qiime tools export --input-path "$WORK/taxonomy.qza" --output-path "$exportdir/tax-$v"
  cp "$exportdir/tax-$v/taxonomy.tsv" "$out/taxonomy.tsv"

  echo "[6/7] Genus-level (level 6) counts for $v"
  qiime taxa barplot --i-table "$t" --i-taxonomy "$WORK/taxonomy.qza" \
    --m-metadata-file "$MANIFEST" --o-visualization "$WORK/taxa-barplot-$v.qzv"
  # level-6.csv sits at <uuid>/data/level-6.csv inside the .qzv
  unzip -o -j "$WORK/taxa-barplot-$v.qzv" "*/data/level-6.csv" -d "$out"
}
export_variant all   "$WORK/table.qza"       "$WORK/rep-seqs.qza"
export_variant noorg "$WORK/table-noorg.qza" "$WORK/rep-seqs-noorg.qza"

echo "[7/7] Done."
qiime --version
