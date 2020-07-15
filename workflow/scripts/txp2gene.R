message("Started")

#
# Manage script inputs
#
gtf_file <- snakemake@input[['gtf']]
tgmap_file <- snakemake@output[[1]]
renv <- snakemake@params[['renv']]

#
# Manage R packages
#
if (renv) {
	renv::activate()
}
library(rtracklayer)
library(tidyverse)

#
# Main script
#
gtf_data <- import.gff2(gtf_file, feature.type = 'transcript')
stopifnot(all(gtf_data$type == 'transcript'))

tg_map <- mcols(gtf_data)[, c('transcript_id', 'transcript_version', 'gene_id', 'gene_version')]
tg_map <- tg_map %>% as_tibble() %>%
  unite('transcript_id', c('transcript_id', 'transcript_version'), sep = ".") %>% 
  unite('gene_id', c('gene_id', 'gene_version'), sep = ".")
tg_map <- unique(tg_map)

write.table(tg_map, tgmap_file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

sessioninfo::session_info()
message("Completed")
