library(tidyverse)
source('src/util.R')

acc_start = 7061266
meta_df <- read.csv('fastq/sra_metadata.csv')
meta_df['Run'] <- sapply(as.numeric(rownames(meta_df)), function(x) paste0('SRR', x+acc_start-1))
meta_df <- meta_df %>% filter(str_detect(Library.Name, 'nt.*minAtc'))
meta_df$Library.Name <- gsub(" ", "_", meta_df$Library.Name)



acc_files <- sapply(meta_df['Run'], function(x) paste0('fastq/', x, '-counts.txt'))
ct_df <- read.delim(acc_files[1], skip=1) %>% dplyr::rename(!!(meta_df$Library.Name[1]) := STDIN)
for (i in seq_along(acc_files[-1])){
	run_df <- read.delim(acc_files[i + 1], skip=1) %>% dplyr::select(STDIN) %>% dplyr::rename(!! meta_df$Library.Name[i + 1] := STDIN) 
	ct_df <- cbind(ct_df, run_df)		
}

ct_df <- ct_df[,-(2:6)]

outstem <- 'data/SRP142541_H37Rv'
write.table(ct_df, file=paste0(outstem, '_raw_count_std.txt'), quote=FALSE, row.names=FALSE, sep='\t')


