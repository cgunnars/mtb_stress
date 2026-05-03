source('src/util.R')
library(tidyverse)
library(readxl)
args <- io_sheet()
df <- read_excel(args$i, sheet = args$sheet)
colnames(df) <- fix_colnames(c(df[1,1:2], colnames(df)[3:4], df[1,5:length(df)])) 
df <- df[-1, ] %>% as.data.frame()

tbl <- load_rv_table()
row <- sapply(1:nrow(df),
	      function(x) {ifelse(is.na(df$Rv_number[x]) || is.na(rv_to_gene(df$Rv_number[x])), df$Gene[x], rv_to_gene(df$Rv_number[x], tbl))
		          }
	      )
df$Gene_name <- row
cpm <- df[,3:4]
write.table(cpm, file=paste0(args$o, '_filt-tpm.txt'), sep='\t', quote=FALSE, row.names=FALSE)

sig_df <- get_sig(df, 'p_value', 'Fold_Change')
head(sig_df)
sig_df_file <- file(paste0(args$o, '_pH-DE_5.7-vs-7.txt'))
writeLines(sig_df$Gene_name, sig_df_file)
close(sig_df_file)

