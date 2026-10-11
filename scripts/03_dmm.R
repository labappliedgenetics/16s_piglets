# ══════════════════════════════════════════════════════════════════════════
# DMM RE-ANALYSIS — genus-resolved taxa only
#
# Основной DMM (модель отбора признаков): top-15 однозначно определённых
# до рода таксонов. Дополнительно (sensitivity): top-20 и top-30, с ARI
# между всеми парами (15 vs 20, 15 vs 30, 20 vs 30) для оценки устойчивости
# разбиения по числу признаков.
#
# Поскольку top-15 даёт best_k = 1 (нет разбиения), а top-20 и top-30
# устойчиво дают k = 2, дальнейшее описание групп DMM_1/DMM_2 и сравнение
# с возрастными группами (14/24/30 сут.) проводится по модели TOP-30.
#
# Возраст НЕ используется при построении DMM — только в постфактум-тесте.
# ══════════════════════════════════════════════════════════════════════════

source("config.R")   # paths, helpers (run from repository root)

# ── 0. Libraries ─────────────────────────────────────────────────────────
if (!requireNamespace("DirichletMultinomial", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DirichletMultinomial", update = FALSE, ask = FALSE)
}
if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")

library(DirichletMultinomial)
library(mclust)      # adjustedRandIndex()
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. Read genus-level count table ─────────────────────────────────────
# Edit path if your genus-level table lives elsewhere. Expected format:
# rows = samples, cols = full taxonomy strings (e.g. "d__...;g__Prevotella"),
# values = raw READ COUNTS (not relative abundance — DMM requires counts).
l6_counts <- read_level6(L6_PATH)
l6_counts <- round(as.matrix(l6_counts))
l6_counts[l6_counts < 0] <- 0

# ── 2. Keep only unambiguously genus-resolved taxa ──────────────────────
# Drops taxa with no g__ level (e.g. family-only assignments like
# f__Lachnospiraceae;g__ or f__Lachnospiraceae;__), plus generic/uninformative
# genus labels (uncultured, unclassified, unidentified, metagenome).
# extract_genus() is defined in config.R
genus_names_all <- vapply(colnames(l6_counts), extract_genus, character(1))
genus_resolved  <- !is.na(genus_names_all)

cat(sprintf("Excluded %d / %d taxa with no resolved genus-level assignment.\n",
            sum(!genus_resolved), length(genus_resolved)))

l6_genus_resolved <- l6_counts[, genus_resolved, drop = FALSE]
colnames(l6_genus_resolved) <- genus_names_all[genus_resolved]
# collapse duplicate genus columns (same genus name from different lineages)
l6_genus_resolved <- t(rowsum(t(l6_genus_resolved), group = colnames(l6_genus_resolved)))

# ── 3. Reusable DMM runner (identical logic for top-15 / 20 / 30) ──────
run_dmm_panel <- function(top_n, counts_full, k_max = 5, seed = 42) {
  rel <- sweep(counts_full, 1, rowSums(counts_full), "/") * 100
  top_taxa   <- names(sort(colMeans(rel), decreasing = TRUE))[1:min(top_n, ncol(counts_full))]
  counts_sub <- counts_full[, top_taxa, drop = FALSE]
  counts_sub <- counts_sub[rowSums(counts_sub) > 0, , drop = FALSE]
  
  set.seed(seed)   # SEED_DMM = 42
  k_range <- 1:min(k_max, nrow(counts_sub) - 1)
  fits <- lapply(k_range, function(k) dmn(counts_sub, k, verbose = FALSE))
  
  fit_stats <- data.frame(
    k       = k_range,
    Laplace = sapply(fits, laplace),
    AIC     = sapply(fits, AIC),
    BIC     = sapply(fits, BIC)
  )
  
  best_k   <- fit_stats$k[which.min(fit_stats$Laplace)]
  best_fit <- fits[[best_k]]
  prob     <- mixture(best_fit)
  group    <- apply(prob, 1, which.max)
  
  assign_df <- data.frame(
    Sample    = rownames(counts_sub),
    DMM_group = paste0("DMM_", group),
    DMM_prob  = apply(prob, 1, max)
  )
  
  list(top_taxa = top_taxa, counts = counts_sub, fit_stats = fit_stats,
       best_k = best_k, best_fit = best_fit, assign = assign_df)
}

# ── 4. MAIN analysis: top-15 genus-resolved panel, k = 1–5 ─────────────
res15 <- run_dmm_panel(15, l6_genus_resolved, k_max = 5)

cat("\n=== MAIN: top-15 genus-resolved panel ===\n")
cat("Genera in panel (ranked by mean relative abundance):\n")
print(res15$top_taxa)
cat("\nModel fit across k = 1..5:\n")
print(res15$fit_stats)
cat(sprintf("\nBest k by Laplace approximation = %d\n", res15$best_k))
print(table(res15$assign$DMM_group))

# per-group mean abundance of the 15 genera
rel15 <- sweep(res15$counts, 1, rowSums(res15$counts), "/") * 100
group_profile15 <- as.data.frame(rel15) %>%
  mutate(Sample = rownames(.)) %>%
  left_join(res15$assign[, c("Sample", "DMM_group")], by = "Sample") %>%
  group_by(DMM_group) %>%
  summarise(across(all_of(res15$top_taxa), mean), n_samples = n(), .groups = "drop")

write.csv(res15$assign,      out_path("dmm_assignments_top15_genus.csv"), row.names = FALSE)
write.csv(group_profile15,   out_path("dmm_group_profile_top15_genus.csv"), row.names = FALSE)
write.csv(res15$fit_stats,   out_path("dmm_fit_stats_top15_genus.csv"), row.names = FALSE)
cat("Saved: dmm_assignments_top15_genus.csv, dmm_group_profile_top15_genus.csv, ",
    "dmm_fit_stats_top15_genus.csv\n", sep = "")

# ── 5. Heatmaps — genus-resolved taxa only, ordered by DMM group ───────
# Reusable so top-15 (Figure 4), top-20, and top-30 all use identical logic.
make_dmm_heatmap <- function(res, label, top_n) {
  rel <- sweep(res$counts, 1, rowSums(res$counts), "/") * 100
  heat_long <- as.data.frame(rel) %>%
    mutate(Sample = rownames(.)) %>%
    left_join(res$assign[, c("Sample", "DMM_group")], by = "Sample") %>%
    arrange(DMM_group) %>%
    mutate(Sample = factor(Sample, levels = unique(Sample))) %>%
    pivot_longer(cols = all_of(res$top_taxa), names_to = "Taxon", values_to = "RelAbund")
  
  p <- ggplot(heat_long, aes(x = Sample, y = Taxon, fill = log1p(RelAbund))) +
    geom_tile() +
    facet_grid(~ DMM_group, scales = "free_x", space = "free_x") +
    scale_fill_viridis_c(name = "log1p(%)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          strip.background = element_rect(fill = "#4A90A4"),
          strip.text = element_text(color = "white", face = "bold")) +
    labs(title = sprintf("Top-%d genus-resolved taxa, ordered by DMM group", top_n),
         x = "Sample", y = "Genus")
  
  ggsave(out_path(sprintf("heatmap_%s.pdf", label)), p,
         width = max(10, nrow(res$assign) * 0.3),
         height = max(6, length(res$top_taxa) * 0.3))
  ggsave(out_path(sprintf("heatmap_%s.png", label)), p,
         width = max(10, nrow(res$assign) * 0.3),
         height = max(6, length(res$top_taxa) * 0.3), dpi = 300)
  p
}

# top-15 heatmap = "Figure 4" for the manuscript
p_fig4 <- make_dmm_heatmap(res15, "Figure4_top15_genus_resolved", 15)

# ── 6. Supplementary Figure — DMM model selection (Laplace & BIC) ──────
p_supp <- ggplot(res15$fit_stats, aes(x = k)) +
  geom_line(aes(y = Laplace, color = "Laplace")) +
  geom_point(aes(y = Laplace, color = "Laplace")) +
  geom_line(aes(y = BIC, color = "BIC")) +
  geom_point(aes(y = BIC, color = "BIC")) +
  theme_bw(base_size = 13) +
  labs(title = "Supplementary Figure. DMM model selection (top-15 genus-resolved taxa)",
       x = "Number of community types (k)", y = "Score", color = "Metric")

ggsave(out_path("SuppFig_DMM_model_selection.pdf"), p_supp, width = 7, height = 5)
ggsave(out_path("SuppFig_DMM_model_selection.png"), p_supp, width = 7, height = 5, dpi = 300)
cat("Saved: Figure4_heatmap_genus_resolved.*, SuppFig_DMM_model_selection.*\n")

# ── 7. SENSITIVITY: top-20 / top-30 vs top-15 (Adjusted Rand Index) ────
res20 <- run_dmm_panel(20, l6_genus_resolved, k_max = 5)
res30 <- run_dmm_panel(min(30, ncol(l6_genus_resolved)), l6_genus_resolved, k_max = 5)

for (nm in c("res20", "res30")) {
  r <- get(nm)
  cat(sprintf("\n=== SENSITIVITY: top-%d genus-resolved panel ===\n", length(r$top_taxa)))
  print(r$fit_stats)
  cat(sprintf("Best k by Laplace = %d\n", r$best_k))
  print(table(r$assign$DMM_group))
}

# top-30 is used going forward to describe DMM_1/DMM_2 and test vs age
# (top-15 collapses to a single component; top-20/30 consistently find k=2)
rel30 <- sweep(res30$counts, 1, rowSums(res30$counts), "/") * 100
group_profile30 <- as.data.frame(rel30) %>%
  mutate(Sample = rownames(.)) %>%
  left_join(res30$assign[, c("Sample", "DMM_group")], by = "Sample") %>%
  group_by(DMM_group) %>%
  summarise(across(all_of(res30$top_taxa), mean), n_samples = n(), .groups = "drop")

write.csv(res30$assign,    out_path("dmm_assignments_top30_genus.csv"), row.names = FALSE)
write.csv(group_profile30, out_path("dmm_group_profile_top30_genus.csv"), row.names = FALSE)
write.csv(res30$fit_stats, out_path("dmm_fit_stats_top30_genus.csv"), row.names = FALSE)
cat("Saved: dmm_assignments_top30_genus.csv, dmm_group_profile_top30_genus.csv, ",
    "dmm_fit_stats_top30_genus.csv\n", sep = "")

ari_15_20 <- {
  m <- inner_join(res15$assign, res20$assign, by = "Sample", suffix = c("_15", "_20"))
  adjustedRandIndex(m$DMM_group_15, m$DMM_group_20)
}
ari_15_30 <- {
  m <- inner_join(res15$assign, res30$assign, by = "Sample", suffix = c("_15", "_30"))
  adjustedRandIndex(m$DMM_group_15, m$DMM_group_30)
}
ari_20_30 <- {
  m <- inner_join(res20$assign, res30$assign, by = "Sample", suffix = c("_20", "_30"))
  adjustedRandIndex(m$DMM_group_20, m$DMM_group_30)
}

cat(sprintf("\nARI (top-15 vs top-20) = %.3f\n", ari_15_20))
cat(sprintf("ARI (top-15 vs top-30) = %.3f\n", ari_15_30))
cat(sprintf("ARI (top-20 vs top-30) = %.3f\n", ari_20_30))
cat("Note: ARI is only a meaningful 'agreement' number when both panels\n",
    "converge on best_k > 1. Check best_k for each panel before interpreting.\n", sep = "")

# heatmaps for the sensitivity panels (same style/ordering as Figure 4)
make_dmm_heatmap(res20, "top20_genus_resolved", 20)
make_dmm_heatmap(res30, "top30_genus_resolved", 30)
cat("Saved: heatmap_Figure4_top15_genus_resolved.*, heatmap_top20_genus_resolved.*, ",
    "heatmap_top30_genus_resolved.*\n", sep = "")

# ── Diagnostic: which genera differ between top-15 and top-20/30 panels? ──
# Helps identify which taxon/taxa likely carry the k=2 signal that the
# top-15 panel misses.
new_in_20 <- setdiff(res20$top_taxa, res15$top_taxa)
new_in_30 <- setdiff(res30$top_taxa, res15$top_taxa)

cat("\nGenera present in top-20 but NOT in top-15 (candidates for the\n",
    "discriminating taxa driving the k=2 split at wider panel sizes):\n", sep = "")
print(new_in_20)
cat("\nGenera present in top-30 but NOT in top-15:\n")
print(new_in_30)

write.csv(
  data.frame(Panel  = c("top15", "top20", "top30"),
             Best_k = c(res15$best_k, res20$best_k, res30$best_k)),
  out_path("dmm_bestk_by_panel_size.csv"), row.names = FALSE
)
write.csv(
  data.frame(Comparison = c("top15_vs_top20", "top15_vs_top30", "top20_vs_top30"),
             ARI        = c(ari_15_20, ari_15_30, ari_20_30)),
  out_path("dmm_sensitivity_ARI.csv"), row.names = FALSE
)

# ══════════════════════════════════════════════════════════════════════════
# 8. POST-HOC ONLY: association between DMM_1/DMM_2 (TOP-30 model — the
#    variant that consistently resolves 2 groups) and age at sampling
#    (14d / 24d / 30d). Age was NOT used to build DMM — this is purely a
#    test of whether the groups found by DMM happen to line up with
#    sampling day.
# ══════════════════════════════════════════════════════════════════════════
# FIX: the original section read "P/metadata.csv" (columns Sample/Age), which does not
# exist, so the test was silently skipped. It now uses metadata_pigs.csv via load_metadata().
metadata <- load_metadata()

age_df <- res30$assign %>%
  left_join(metadata[, c("sample", "Age")], by = c("Sample" = "sample"))
stopifnot(!anyNA(age_df$Age))

age_mat <- table(Age = age_df$Age, DMM_group = age_df$DMM_group)
cat("\n=== DMM group (top-30 model) distribution across sampling ages ===\n")
print(age_mat)

if (ncol(age_mat) < 2) {
  cat("Only one DMM group present: no DMM_1 vs DMM_2 split to test against age.\n")
} else {
  # Fisher's exact test (expected counts are small with n = 34)
  age_test <- fisher.test(age_mat)
  print(age_test)
  cat(sprintf("Fisher's exact test, age vs DMM group: p = %.4g\n", age_test$p.value))
  write.csv(as.data.frame.matrix(age_mat), out_path("dmm_by_age_contingency_top30.csv"))
}

cat("\n✓ Genus-resolved DMM re-analysis complete.\n")