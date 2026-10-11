# ══════════════════════════════════════════════════════════════════════════
# Shared configuration + helpers for all R scripts (source("config.R"))
# Run every script from the repository root.
#
#   PIGLET_VARIANT = "noorg"  -> table after removing chloroplast/mitochondria
#                                (recommended; produced by scripts/01_*.sh)
#   PIGLET_VARIANT = "all"    -> table WITHOUT organelle filtering
#                                (= the table used for the first submission)
# Override from the shell, e.g.:  PIGLET_VARIANT=all Rscript scripts/02_alpha_beta.R
# ══════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

DATA_DIR <- Sys.getenv("PIGLET_DATA", "data")       # inputs
OUT_DIR  <- Sys.getenv("PIGLET_OUT",  "results")    # outputs
VARIANT  <- Sys.getenv("PIGLET_VARIANT", "noorg")
stopifnot(VARIANT %in% c("all", "noorg"))

TABLE_DIR     <- file.path(DATA_DIR, VARIANT)       # feature-table.tsv, level-6.csv ...
METADATA_PATH <- file.path(DATA_DIR, "metadata_pigs.csv")
FT_PATH       <- file.path(TABLE_DIR, "feature-table.tsv")
L6_PATH       <- file.path(TABLE_DIR, "level-6.csv")

RARE_DEPTH <- 4871      # rarefaction depth used for alpha diversity (see README)
SEED_ALPHA <- 123; SEED_NMDS <- 42; SEED_DMM <- 42; SEED_ANCOM <- 42

out_path <- function(...) {
  d <- file.path(OUT_DIR, VARIANT)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, ...)
}

# metadata: columns `sample` and `group` (group_1/2/3 -> 30_d/14_d/24_d)
load_metadata <- function(path = METADATA_PATH) {
  md <- read.csv(path, stringsAsFactors = FALSE)
  md$group <- dplyr::recode(md$group,
                            "group_1" = "30_d", "group_2" = "14_d", "group_3" = "24_d")
  md$sample <- as.character(md$sample)
  md$group  <- factor(md$group, levels = c("14_d", "24_d", "30_d"))
  md$Age    <- md$group                       # alias used in the age-confounding analyses
  md
}

# QIIME 2 / biom TSV export starts with "# Constructed from biom file" and
# "#OTU ID": read it robustly (rows = ASVs, cols = samples).
read_feature_table <- function(path = FT_PATH) {
  first <- readLines(path, n = 1)
  skip  <- if (grepl("^# Constructed", first)) 1 else 0
  read.table(path, sep = "\t", header = TRUE, row.names = 1, skip = skip,
             comment.char = "", check.names = FALSE, quote = "")
}

# QIIME 2 taxa-barplot export (level-6.csv): rows = samples, cols = taxonomy
# strings + metadata columns (manifest file paths), which are dropped here.
read_level6 <- function(path = L6_PATH) {
  l6 <- read.csv(path, row.names = 1, check.names = FALSE)
  l6$`forward-absolute-filepath` <- NULL
  l6$`reverse-absolute-filepath` <- NULL
  l6
}

# last assigned taxonomic rank of a QIIME taxonomy string
shorten_taxa <- function(x) {
  parts <- trimws(unlist(strsplit(x, ";")))
  gsub("^[a-z]__", "", tail(parts[nzchar(parts)], 1))
}

# genus name only if g__ is resolved and not a generic placeholder
extract_genus <- function(tax_string) {
  parts  <- trimws(strsplit(tax_string, ";")[[1]])
  g_part <- parts[grepl("^g__", parts)]
  if (length(g_part) == 0) return(NA_character_)
  g_name <- gsub("^g__", "", g_part[1])
  if (g_name == "" || tolower(g_name) %in%
      c("uncultured", "unclassified", "unidentified", "metagenome")) return(NA_character_)
  g_name
}

cat(sprintf("[config] variant = %s | tables = %s | outputs = %s\n",
            VARIANT, TABLE_DIR, file.path(OUT_DIR, VARIANT)))
