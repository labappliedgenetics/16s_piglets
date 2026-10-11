source("config.R")   # run from repository root
# ── ANCOM-BC2: сравнение 14 d vs 24 d vs 30 d на уровне рода ─────────────────
library(ANCOMBC)
library(phyloseq)
library(dplyr)
library(tidyr)

# ── 1. Чтение genus-level (level-6) таблицы counts ───────────────────────────
l6 <- read.csv(L6_PATH, row.names = 1, check.names = FALSE)
l6$`forward-absolute-filepath` <- NULL
l6$`reverse-absolute-filepath` <- NULL

# ── 2. Объединение одинаковых genus assignments ──────────────────────────────
# Колонки level-6.csv — это полные таксономические строки до уровня рода;
# разные строки могут сводиться к одному и тому же роду (например, если
# отличаются только вышестоящим неопределённым таксоном) — такие колонки нужно
# просуммировать до тестирования.
# shorten_taxa() is defined in config.R
genus_names <- sapply(colnames(l6), shorten_taxa)

l6_t <- t(l6)                                       # taxa x samples
rownames(l6_t) <- genus_names
l6_genus <- rowsum(l6_t, group = rownames(l6_t))    # суммирует дублирующиеся роды
l6_genus <- round(l6_genus)                         # целочисленные counts для ANCOM-BC2

cat(sprintf("Родов после объединения: %d (было колонок: %d)\n",
            nrow(l6_genus), ncol(l6)))

# ── 2b. Убираем неопределённые роды (g__ без названия) ───────────────────────
# после strip префикса "g__" у таких таксонов остаётся пустая строка —
# оставляем только строки с реальным названием рода
known <- nzchar(rownames(l6_genus))
cat(sprintf("Известных родов: %d (отброшено неклассифицированных g__: %d)\n",
            sum(known), sum(!known)))
l6_genus <- l6_genus[known, , drop = FALSE]

# ── 3. Метаданные (14_d / 24_d / 30_d) ────────────────────────────────────────
metadata <- load_metadata()
rownames(metadata) <- metadata$sample

# синхронизация выборок между таблицей counts и метаданными
common <- intersect(colnames(l6_genus), rownames(metadata))
l6_genus <- l6_genus[, common, drop = FALSE]
metadata  <- metadata[common, , drop = FALSE]
stopifnot(identical(colnames(l6_genus), rownames(metadata)))

# ── 4. Фильтр по относительной представленности ──────────────────────────────
# оставляем роды с относительной представленностью > 0.1% как минимум в 20% образцов
rel  <- sweep(l6_genus, 2, colSums(l6_genus), "/")
keep <- rowMeans(rel > 0.001) >= 0.20
cat(sprintf("Родов прошло фильтр (>0.1%% в >=20%% образцов): %d из %d\n",
            sum(keep), nrow(l6_genus)))
l6_filt <- l6_genus[keep, , drop = FALSE]

# ── 5. phyloseq-объект ────────────────────────────────────────────────────────
otu_tab <- otu_table(as.matrix(l6_filt), taxa_are_rows = TRUE)
samp    <- sample_data(metadata)
ps      <- phyloseq(otu_tab, samp)

# ── 6. ANCOM-BC2 ───────────────────────────────────────────────────────────
set.seed(SEED_ANCOM)
out <- ancombc2(
  data          = ps,
  assay_name    = "counts",
  fix_formula   = "group",
  rand_formula  = NULL,
  p_adj_method  = "BH",       # BH/FDR коррекция
  prv_cut       = 0,          # фильтрация по представленности уже сделана в п.4
  lib_cut       = 0,
  group         = "group",
  struc_zero    = TRUE,
  neg_lb        = TRUE,
  alpha         = 0.05,
  global        = TRUE,       # омнибус-тест across 14d/24d/30d
  pairwise      = TRUE,       # попарные сравнения (14d-24d, 14d-30d, 24d-30d)
  dunnet        = FALSE,
  trend         = FALSE,
  iter_control  = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
  em_control    = list(tol = 1e-5, max_iter = 100),
  lme_control   = NULL,
  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100)
)

# ── 7. Извлечение и сохранение результатов ───────────────────────────────────
res_primary <- out$res         # LFC каждой группы относительно референса (14_d)
res_global  <- out$res_global  # омнибус-тест по всем трём группам
res_pair    <- out$res_pair    # все попарные контрасты

write.csv(res_primary, out_path("ancombc2_genus_vsref14d.csv"), row.names = FALSE)
write.csv(res_global,  out_path("ancombc2_genus_global.csv"),   row.names = FALSE)
write.csv(res_pair,    out_path("ancombc2_genus_pairwise.csv"), row.names = FALSE)

cat("✓ Итоговые таблицы сохранены → results/<variant>/ancombc2_genus_*.csv\n")

# ── 8. Роды с q < 0.05 ─────────────────────────────────────────────────────
sig_global <- res_global %>% filter(q_val < 0.05)
cat(sprintf("\nРодов, значимых в глобальном (омнибус) тесте при q < 0.05: %d\n",
            nrow(sig_global)))
if (nrow(sig_global) > 0) print(sig_global$taxon)

q_pair_cols <- grep("^q_", colnames(res_pair), value = TRUE)
sig_pair <- res_pair %>%
  filter(if_any(all_of(q_pair_cols), ~ . < 0.05))
cat(sprintf("Родов, значимых хотя бы в одном попарном сравнении при q < 0.05: %d\n",
            nrow(sig_pair)))
if (nrow(sig_pair) > 0) print(sig_pair$taxon)

if (nrow(sig_global) == 0 && nrow(sig_pair) == 0) {
  cat("\n⚠ Ни один род не достиг q < 0.05 ни в глобальном, ни в попарных тестах.\n")
}