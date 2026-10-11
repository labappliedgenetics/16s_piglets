# ══════════════════════════════════════════════════════════════════════════
# 02 — Alpha diversity (rarefied), beta diversity, PERMANOVA/PERMDISP,
#      genus barplot and top-50 heatmap.   Run from the repository root.
# Input : data/<variant>/feature-table.tsv, level-6.csv, data/metadata_pigs.csv
# Output: results/<variant>/ (alpha_diversity.csv, beta_*.csv, permanova_bray.csv, figures)
# ══════════════════════════════════════════════════════════════════════════
source("config.R")
# ── Load libraries ────────────────────────────────────────────────────────────
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)

# ── 1. Read feature table ─────────────────────────────────────────────────────
ft_raw <- read_feature_table(FT_PATH)
# ── Read metadata ─────────────────────────────────────────────────────────────
metadata <- load_metadata()
print(table(metadata$group))

# Transpose → rows = samples, cols = ASVs (required by vegan)
otu <- t(ft_raw)               # matrix: samples × ASVs
otu <- round(otu)              # ensure integers

# FIX: the original script called setdiff(rownames(otu), ...) before `otu` existed
missing_samples <- setdiff(rownames(otu), metadata$sample)
if (length(missing_samples) > 0)
  stop("Samples in feature table but not in metadata: ", paste(missing_samples, collapse = ", "))

# align metadata to the table: this order MUST be identical for adonis2()
metadata <- metadata %>% filter(sample %in% rownames(otu)) %>%
  arrange(match(sample, rownames(otu)))
stopifnot(identical(metadata$sample, rownames(otu)))

cat(sprintf("Samples: %d | ASVs: %d | min depth %d | median depth %d | max depth %d\n",
            nrow(otu), ncol(otu), min(rowSums(otu)), median(rowSums(otu)), max(rowSums(otu))))
# ── 2. ALPHA DIVERSITY ────────────────────────────────────────────────────────
set.seed(SEED_ALPHA)  # for reproducibility

# the depth is set in config.R; with organelle filtering the minimum depth can change
rare_depth <- RARE_DEPTH
if (min(rowSums(otu)) < rare_depth) {
  rare_depth <- min(rowSums(otu))
  warning(sprintf("Rarefaction depth lowered to the minimum sample depth: %d", rare_depth))
}
write(rare_depth, out_path("rarefaction_depth_used.txt"))
n_iter     <- 100

# Storage for each iteration's metrics
rare_results <- vector("list", n_iter)

for (i in seq_len(n_iter)) {
  
  otu_rare <- rrarefy(otu, sample = rare_depth)
  
  rare_results[[i]] <- data.frame(
    Sample     = rownames(otu_rare),
    Observed   = specnumber(otu_rare),
    Shannon    = diversity(otu_rare, index = "shannon"),
    Simpson    = diversity(otu_rare, index = "simpson"),
    InvSimpson = diversity(otu_rare, index = "invsimpson"),
    Chao1      = apply(otu_rare, 1, function(x) estimateR(x)["S.chao1"]),
    iter       = i,
    row.names  = NULL
  )
}

rare_all <- bind_rows(rare_results)

# Median across the 100 iterations, per sample
alpha <- rare_all %>%
  group_by(Sample) %>%
  summarise(
    Observed   = median(Observed),
    Shannon    = median(Shannon),
    Simpson    = median(Simpson),
    InvSimpson = median(InvSimpson),
    Chao1      = median(Chao1),
    .groups = "drop"
  )

# Merge with metadata
alpha <- left_join(alpha,
                   metadata,
                   by = c("Sample" = "sample"))

write.csv(alpha, out_path("alpha_diversity.csv"), row.names = FALSE)
head(alpha)
# ── Alpha diversity boxplots ──────────────────────────────────────────────────
alpha_long <- alpha %>%
  pivot_longer(
    cols = c(Observed, Shannon, Simpson, InvSimpson, Chao1),
    names_to = "Metric",
    values_to = "Value"
  )

# ── P-values for each alpha-diversity metric ───────────────────────────────────
p_values <- alpha_long %>%
  group_by(Metric) %>%
  summarise(
    p_value = kruskal.test(Value ~ group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("p = ", format.pval(p_value, digits = 3, eps = 0.001))
  )

print(p_values)


# ── Alpha diversity boxplots ──────────────────────────────────────────────────
p_alpha <- ggplot(alpha_long,
                  aes(x = group,
                      y = Value,
                      fill = group)) +
  geom_boxplot(
    alpha = 0.8,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.5,
    alpha = 0.7
  ) +
  
  facet_wrap(~Metric, scales = "free_y") +
  
  # p-value above each panel
  geom_text(
    data = p_values,
    aes(
      x = -Inf,
      y = Inf,
      label = p_label
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.5,
    size = 4
  ) +
  
  scale_fill_brewer(palette = "Set2") +
  
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "#4A90A4"),
    strip.text = element_text(color = "white", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  labs(
    title = "Alpha Diversity by Group",
    x = "Group",
    y = "Value"
  )


ggsave(
  out_path("alpha_diversity_boxplot.pdf"),
  p_alpha,
  width = 10,
  height = 6
)

ggsave(
  out_path("alpha_diversity_boxplot.png"),
  p_alpha,
  width = 10,
  height = 6,
  dpi = 300
)


# ── 3. BETA DIVERSITY ─────────────────────────────────────────────────────────
# Bray-Curtis dissimilarity on raw counts
bray   <- vegdist(otu, method = "bray")

# Jaccard (presence/absence)
jaccard <- vegdist(otu, method = "jaccard", binary = TRUE)

# Aitchison = Euclidean distance on CLR-transformed counts; zeros replaced by 0.5
# (simple pseudocount; this is NOT robust Aitchison/RPCA — see 02b for the latter)
clr_transform <- function(x) {
  x[x == 0] <- 0.5
  log(x) - mean(log(x))
}
otu_clr   <- t(apply(otu, 1, clr_transform))
aitchison <- dist(otu_clr, method = "euclidean")

# Save dissimilarity matrices
write.csv(as.matrix(bray),     out_path("beta_braycurtis.csv"))
write.csv(as.matrix(jaccard),  out_path("beta_jaccard.csv"))
write.csv(as.matrix(aitchison),out_path("beta_aitchison.csv"))
cat("Beta diversity matrices saved → beta_*.csv\n")

# ── PCoA ordination (Bray-Curtis) ─────────────────────────────────────────────
pcoa_res <- cmdscale(bray, k = 2, eig = TRUE)

eig_pct <- round(
  pcoa_res$eig[1:2] / sum(abs(pcoa_res$eig)) * 100,
  1
)

pcoa_df <- data.frame(
  Sample = rownames(pcoa_res$points),
  PC1 = pcoa_res$points[,1],
  PC2 = pcoa_res$points[,2]
)

# Add metadata
pcoa_df <- left_join(
  pcoa_df,
  metadata,
  by = c("Sample" = "sample")
)

p_pcoa <- ggplot(
  pcoa_df,
  aes(PC1, PC2, color = group)
) +
  geom_point(size = 4, alpha = 0.85) +
  stat_ellipse(aes(fill = group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  theme_bw(base_size = 13) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "PCoA (Bray-Curtis)",
    color = "Group",
    fill = "Group",
    x = paste0("PC1 (", eig_pct[1], "%)"),
    y = paste0("PC2 (", eig_pct[2], "%)")
  )

ggsave(out_path("beta_pcoa_braycurtis.pdf"), p_pcoa, width = 8, height = 6)
ggsave(out_path("beta_pcoa_braycurtis.png"), p_pcoa, width = 8, height = 6, dpi = 300)
# ── NMDS (Bray-Curtis) ────────────────────────────────────────────────────────
set.seed(42)

nmds_res <- metaMDS(
  otu,
  distance = "bray",
  k = 2,
  trymax = 100,
  trace = FALSE
)

cat(sprintf("NMDS stress = %.4f\n", nmds_res$stress))

nmds_df <- data.frame(
  Sample = rownames(nmds_res$points),
  NMDS1 = nmds_res$points[,1],
  NMDS2 = nmds_res$points[,2]
)

nmds_df <- left_join(
  nmds_df,
  metadata,
  by = c("Sample" = "sample")
)

p_nmds <- ggplot(
  nmds_df,
  aes(NMDS1, NMDS2, color = group)
) +
  geom_point(size = 4, alpha = 0.85) +
  stat_ellipse(aes(fill = group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = paste("Stress =", round(nmds_res$stress,4)),
    hjust = 1.1,
    vjust = -0.5
  ) +
  theme_bw(base_size = 13) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "NMDS (Bray-Curtis)",
    color = "Group",
    fill = "Group"
  )

ggsave(out_path("beta_nmds_braycurtis.pdf"), p_nmds, width = 8, height = 6)
ggsave(out_path("beta_nmds_braycurtis.png"), p_nmds, width = 8, height = 6, dpi = 300)


permanova <- adonis2(
  bray ~ group,
  data = metadata,
  permutations = 999
)

print(permanova)

write.csv(
  as.data.frame(permanova),
  out_path("permanova_bray.csv")
)

# ── PERMDISP: homogeneity of multivariate dispersion (needed to interpret PERMANOVA) ──
# NOTE: added to the repository; reports distance to group spatial median.
for (dn in c("bray", "jaccard", "aitchison")) {
  d  <- get(dn)
  bd <- betadisper(d, metadata$group)
  pt <- permutest(bd, permutations = 999)
  cat(sprintf("\nPERMDISP (%s): F = %.3f, p = %.4f\n", dn,
              pt$tab$F[1], pt$tab$`Pr(>F)`[1]))
  print(tapply(bd$distances, metadata$group, mean))
  write.csv(data.frame(Sample = names(bd$distances), Group = metadata$group,
                       DistToMedian = bd$distances),
            out_path(sprintf("betadisper_distances_%s.csv", dn)), row.names = FALSE)
  capture.output(print(pt), print(tapply(bd$distances, metadata$group, mean)),
                 file = out_path(sprintf("betadisper_summary_%s.txt", dn)))
  # PERMANOVA on every matrix (same seed-free 999 permutations as above)
  if (dn != "bray") {
    pm <- adonis2(d ~ group, data = metadata, permutations = 999)
    write.csv(as.data.frame(pm), out_path(sprintf("permanova_%s.csv", dn)))
  }
}



# ── 4. LEVEL-6 TAXA BARPLOT (ALL SAMPLES GROUPED BY CLUSTER) ─────────────────

if (file.exists(L6_PATH)) {
  
  # Read genus table (rows = samples) and drop manifest columns
  l6 <- read_level6()
  
  # ---------------- Genus-level resolution filter ----------------
  # Collapse each taxonomy string to its g__ genus name. Anything with no
  # resolved genus (e.g. f__Lachnospiraceae;g__ or ;__) or a generic
  # placeholder (uncultured/unclassified/metagenome) is routed into an
  # "Unassigned" bucket instead of being plotted as if it were a real genus.
  genus_names_all <- vapply(colnames(l6), extract_genus, character(1))
  genus_resolved  <- !is.na(genus_names_all)
  
  cat(sprintf("Excluded %d / %d taxa with no resolved genus-level assignment ",
              sum(!genus_resolved), length(genus_resolved)))
  cat("(reads routed into 'Unassigned').\n")
  
  # collapse genus-resolved columns by genus name (sum duplicate lineages)
  l6_genus <- as.data.frame(
    t(rowsum(t(l6[, genus_resolved, drop = FALSE]),
             group = genus_names_all[genus_resolved]))
  )
  
  # everything without a resolved genus goes into one "Unassigned" column
  if (any(!genus_resolved)) {
    l6_genus$Unassigned <- rowSums(l6[, !genus_resolved, drop = FALSE])
  }
  
  # Convert to relative abundance (genus-resolved taxa + Unassigned bucket)
  l6_rel <- sweep(l6_genus, 1, rowSums(l6_genus), "/") * 100
  
  # ---------------- Top taxa ----------------
  # Rank only genus-resolved taxa — "Unassigned" is never eligible for a
  # top-N slot, it's always folded into "Other" below regardless of size.
  top_n <- 20
  
  rankable <- setdiff(colnames(l6_rel), "Unassigned")
  mean_ab  <- colMeans(l6_rel[, rankable, drop = FALSE])
  
  top_taxa <- names(
    sort(mean_ab, decreasing = TRUE)[1:min(top_n, length(mean_ab))]
  )
  
  l6_top <- l6_rel[, top_taxa, drop = FALSE]
  
  other_cols <- setdiff(colnames(l6_rel), top_taxa)
  if (length(other_cols) > 0) {
    l6_top$Other <- rowSums(l6_rel[, other_cols, drop = FALSE])
  }
  
  # ---------------- Long format ----------------
  
  l6_top$Sample <- rownames(l6_top)
  
  l6_long <- pivot_longer(
    l6_top,
    cols = -Sample,
    names_to = "Taxon",
    values_to = "RelAbund"
  )
  
  # ---------------- Merge metadata ----------------
  
  l6_long <- left_join(
    l6_long,
    metadata,
    by = c("Sample" = "sample")
  )
  
  # ---------------- Order samples by cluster ----------------
  
  metadata_order <- metadata %>%
    arrange(group, sample)
  
  l6_long$Sample <- factor(
    l6_long$Sample,
    levels = metadata_order$sample
  )
  
  # Taxon names are already clean genus names (or "Other") — no further
  # shortening needed since collapsing already happened above.
  
  # Keep taxa ordered by abundance
  tax_order <- l6_long %>%
    group_by(Taxon) %>%
    summarise(mean_ab = mean(RelAbund), .groups = "drop") %>%
    arrange(desc(mean_ab)) %>%
    pull(Taxon)
  
  l6_long$Taxon <- factor(
    l6_long$Taxon,
    levels = rev(tax_order)
  )
  
  # ---------------- Colors ----------------
  
  n_taxa <- length(unique(l6_long$Taxon))
  
  pal <- colorRampPalette(
    brewer.pal(12, "Paired")
  )(n_taxa)
  
  # ---------------- Plot ----------------
  
  p_bar <- ggplot(
    l6_long,
    aes(
      x = Sample,
      y = RelAbund,
      fill = Taxon
    )
  ) +
    geom_bar(stat = "identity", width = 0.9) +
    
    facet_grid(
      ~group,
      scales = "free_x",
      space = "free_x"
    ) +
    
    scale_fill_manual(values = pal) +
    
    scale_y_continuous(
      expand = c(0,0)
    ) +
    
    theme_bw(base_size = 13) +
    
    theme(
      
      strip.background = element_rect(fill = "#4A90A4"),
      
      strip.text = element_text(
        colour = "white",
        face = "bold"
      ),
      
      axis.text.x = element_blank(),
      
      axis.ticks.x = element_blank(),
      
      panel.spacing.x = unit(0.3,"lines"),
      
      panel.grid = element_blank(),
      
      legend.text = element_text(size = 8),
      
      legend.key.size = unit(0.45,"cm")
      
    ) +
    
    labs(
      
      title = "Genus Composition of Individual Samples Grouped by Age",
      
      x = "Samples",
      
      y = "Relative abundance (%)",
      
      fill = "Genus"
      
    )
  
  ggsave(
    out_path("barplot_level6_clustered.pdf"),
    p_bar,
    width = 14,
    height = 7
  )
  
  ggsave(
    out_path("barplot_level6_clustered.png"),
    p_bar,
    width = 14,
    height = 7,
    dpi = 300
  )
  
  cat("✓ Clustered genus barplot saved.\n")
  
} else {
  
  cat("level-6.csv not found.\n")
  
}



# ── 5. HEATMAP OF TOP 50 GENERA (ALL SAMPLES GROUPED BY CLUSTER) ─────────────

if (file.exists(L6_PATH)) {
  
  l6 <- read_level6()
  
  # Relative abundance
  l6_rel <- sweep(l6, 1, rowSums(l6), "/") * 100
  
  # ---------------- Top 50 genera ----------------
  
  top_n <- 50
  
  mean_ab <- colMeans(l6_rel)
  
  top_taxa <- names(
    sort(mean_ab, decreasing = TRUE)[
      1:min(top_n, length(mean_ab))
    ]
  )
  
  heat <- l6_rel[, top_taxa, drop = FALSE]
  
  # ---------------- Order samples by cluster ----------------
  
  metadata_order <- metadata %>%
    arrange(group, sample)
  
  heat <- heat[metadata_order$sample, , drop = FALSE]
  
  # ---------------- Shorten taxonomy names ----------------
  
  colnames(heat) <- sapply(colnames(heat), shorten_taxa)
  
  # ---------------- Long format ----------------
  
  heat$Sample <- rownames(heat)
  
  heat_long <- pivot_longer(
    heat,
    cols = -Sample,
    names_to = "Genus",
    values_to = "Abundance"
  )
  
  heat_long <- left_join(
    heat_long,
    metadata,
    by = c("Sample"="sample")
  )
  
  # Sample order
  heat_long$Sample <- factor(
    heat_long$Sample,
    levels = metadata_order$sample
  )
  
  # Genus order
  genus_order <- heat_long %>%
    group_by(Genus) %>%
    summarise(mean_ab = mean(Abundance)) %>%
    arrange(mean_ab)
  
  heat_long$Genus <- factor(
    heat_long$Genus,
    levels = genus_order$Genus
  )
  
  # ---------------- Plot ----------------
  
  p_heat <- ggplot(
    heat_long,
    aes(
      x = Sample,
      y = Genus,
      fill = Abundance
    )
  ) +
    
    geom_tile() +
    
    facet_grid(
      ~group,
      scales = "free_x",
      space = "free_x"
    ) +
    
    scale_fill_gradientn(
      colours = c(
        "white",
        "#B3CDE3",
        "#6497B1",
        "#005B96",
        "#03396C"
      ),
      name = "Relative abundance (%)"
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      
      axis.text.x = element_blank(),
      
      axis.ticks.x = element_blank(),
      
      axis.text.y = element_text(size = 8),
      
      strip.background = element_rect(fill = "#4A90A4"),
      
      strip.text = element_text(
        colour = "white",
        face = "bold"
      ),
      
      panel.spacing.x = unit(0.2,"lines"),
      
      panel.grid = element_blank()
      
    ) +
    
    labs(
      
      title = "Top 50 Genera Across Samples Grouped by Age",
      
      x = "Samples",
      
      y = "Genus"
      
    )
  
  ggsave(
    out_path("heatmap_top50_clustered.pdf"),
    p_heat,
    width = 15,
    height = 10
  )
  
  ggsave(
    out_path("heatmap_top50_clustered.png"),
    p_heat,
    width = 15,
    height = 10,
    dpi = 300
  )
  
  cat("✓ Heatmap saved.\n")
  
}
