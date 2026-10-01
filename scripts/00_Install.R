#Installation file : GSM multiomic analyses

#installation function
installPkg <- function(pkg) {
  if (!pkg %in% rownames(installed.packages())) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

#install biocmanager if not already present
if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager")
}

#packages to install
pkgs <- c("dada2", "phyloseq", "compositions", #microbiome processing
          "Biostrings", #save nucleotide sequences
          "maaslin3", #microbiome abundance and prevalence modelling
          "lme4", "lmerTest", #linear mixed models
          "variancePartition", #CCA analysis
          "ggraph", "igraph", "tidygraph", #for network manipulation
          "vegan", #for mantel test
          "ggpubr", #plot manipulation
          "tidyverse") # general data wrangling/manipulation

# Set a different CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

for (p in pkgs) {
  installPkg(p)
}

#check
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
print(missing)
