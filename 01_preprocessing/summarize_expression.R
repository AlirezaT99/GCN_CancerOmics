cancer_types <- c('Brain', 'Breast', 'Colorectal', 'Lung', 'Nervous System', 'Pancreas', 'Uterus')
trouble_maker <- c('Blood', 'Pancreas')

prepare_expression <- function(data_folder) {
  exp_path <- sprintf('%s/exp_array.tsv', data_folder)
  exp_array <- read.csv(exp_path, sep='\t')[c('icgc_donor_id', 'gene_id', 'normalized_expression_value')]
  if (data_folder %in% trouble_maker) {
    mapper <- read.csv(sprintf("%s/%s_gene_mapper.csv", data_folder, data_folder))[c('gene_id', 'gene_symbol')]
    exp_array$gene_id <- sub("[.][0-9]*", "", exp_array$gene_id)
    exp_array1 <- merge(x = exp_array, y = mapper, by = 'gene_id', all.x = TRUE)
    exp_array1 <- na.omit(exp_array1)
    exp_array <- exp_array1[c('icgc_donor_id', 'gene_symbol', 'normalized_expression_value')]
    colnames(exp_array)[colnames(exp_array) == 'gene_symbol'] <- 'gene_id'
  }
  write.csv(x = exp_array, file = sprintf("%s/expression_data.tsv", data_folder), sep = "\t")
}

setwd("../../1_PreprocessData/data")

for (c_type in cancer_types) {
  print(sprintf("Working on %s", c_type))
  prepare_expression(data_folder = c_type)
}
