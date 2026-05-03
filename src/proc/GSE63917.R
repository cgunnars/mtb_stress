source('src/util.R')
library(tidyverse)
library(argparse)
parser <- ArgumentParser()
parser$add_argument('-i', type='character', nargs=2)
parser$add_argument('-o', type='character', nargs=1)
args <- parser$parse_args()

df_up <- read.delim(args$i[1], skip=3, check.names=FALSE)
df_down <- read.delim(args$i[2], skip=3, check.names=FALSE)

df <- rbind(df_up, df_down)
colnames(df) <- fix_colnames(colnames(df))
colnames(df)[1] <- 'MT_Number'

rv_to_gene((df$Rv_number[4]))
df$Gene_name <- sapply(as.numeric(rownames(df)), 
		       function(x) {ifelse(df$Rv_number[x] == "",
				    as.vector(df$MT_Number[x]), 
				    rv_to_gene(df$Rv_number[x]))})

cpm <- df[, 2:5]
rownames(cpm) <- df$Gene_name
write.table(cpm, file=paste0(args$o, '_filt-tpm.txt'), sep='\t', quote=FALSE)

sig_df <- get_sig(df, 'Adjusted_p_value', 'log2_Fold_change')
sig_df_file <- file(paste0(args$o, '_pH-DE_DMSO5.7-vs-DMSO7.txt'))
writeLines(sig_df$Gene_name, sig_df_file)
close(sig_df_file)



