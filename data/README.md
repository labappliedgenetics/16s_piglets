# data/ (contents not versioned until publication)

Place locally:
- `metadata_pigs.csv` — columns `sample`, `group` (group_1 = 30 d, group_2 = 14 d, group_3 = 24 d)
- `all/` and `noorg/` — written by `scripts/01_qiime2_dada2.sh` (feature-table.tsv, level-6.csv, taxonomy.tsv, dna-sequences.fasta)
- `stats.tsv` — DADA2 read tracking

The QIIME 2 manifest (`sample-id`, `forward-absolute-filepath`, `reverse-absolute-filepath`; TSV) is kept outside the repository.
