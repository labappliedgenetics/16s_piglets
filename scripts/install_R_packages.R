cran <- c("vegan","ggplot2","dplyr","tidyr","reshape2","RColorBrewer","mclust","BiocManager")
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(c("DirichletMultinomial","ANCOMBC","phyloseq"), update = FALSE, ask = FALSE)
