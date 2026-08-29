#### Install and library ####

# install.packages("biotools")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("illuminaHumanv4.db")
  BiocManager::install("clusterProfiler")
  BiocManager::install("org.Hs.eg.db")
}

library("illuminaHumanv4.db")
library("clusterProfiler")
library("org.Hs.eg.db")

#### read genes ####
exp_path <- '../../1_PreprocessData/data/Blood/exp_array.tsv'
genes <- read.csv(exp_path, sep='\t')[['gene_id']]

#### Preprocess the genes ####

src <- "ENSEMBL"
dst <- "SYMBOL"

if (substr(genes[1], 0, 3) == "ENS") {  # ENSEMBL + version
  genes <- sub("[.][0-9]*", "", genes)
} else if (substr(genes[1], 0, 3) == "NM_" | substr(genes[1], 0, 3) == "NR_") {  # RefSeq
  src <- "REFSEQ"
}

#### Convert probe Id to gene symbol ####

df <- data.frame(
  bitr(
    genes,
    fromType = src,
    toType = dst,
    OrgDb = org.Hs.eg.db,
    drop = TRUE
  )
)
colnames(df) <- c('initial_id', 'Gene')

# For Illumina (Such as Pancreas)
# df <- data.frame(Gene=unlist(mget(x = genes, envir = illuminaHumanv4SYMBOL)))

write.csv(x = df, file = './converted_genes.csv')
