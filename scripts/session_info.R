# Record package versions for the Methods / supplementary: Rscript scripts/session_info.R
source("config.R")
for (p in c("vegan","ggplot2","dplyr","tidyr","reshape2","RColorBrewer","DirichletMultinomial",
            "mclust","ANCOMBC","phyloseq")) suppressWarnings(suppressMessages(library(p, character.only = TRUE)))
dir.create(OUT_DIR, showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
print(sessionInfo())
