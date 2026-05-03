source('src/util.R')
library(tidyverse)
args <- io_sheet() 
df <- read_sheet(args$i, args$sheet, header_row = 3) 
nrow(df)
cpm <- df[, 2:5]
head(cpm)
#cpm has MT_Number, RV_number, Gene_name, but some MTs lack gene names, and Gene_name has non-unique PEs
mtb_genes <- read.delim('data/Mycobacterium_tuberculosis_H37Rv_txt_v2.txt') %>% filter(Feature == 'CDS')
rownames(mtb_genes) <- mtb_genes$Locus
#convert, where available: gene name > rv number > mt number
rownames(cpm) <- sapply(as.numeric(rownames(cpm)), 
	     function(x) {ifelse( is.na(as.vector(mtb_genes[df$RV_number[x], 'Name'])), 
				  as.vector(df$MT_Number[x]), 
				  as.vector(mtb_genes[df$RV_number[x], 'Name']))
			 }
	     )
write.table(cpm, file=paste0(args$o, '_filt-tpm.txt'), sep='\t', quote=FALSE)

df$Gene_name <- rownames(cpm)

sig_df <- get_sig(df, 'Adjusted_p_value', 'Log2_Fold_change')
sig_df_file <- file(paste0(args$o, '_pH-DE_DMSO5-vs-DMSO7.txt'))
writeLines(sig_df$Gene_name, sig_df_file)
close(sig_df_file)
