io_base <- function() {
	library(argparse)
	parser <- ArgumentParser()
	parser$add_argument('-i', type='character', nargs=1)
	parser$add_argument('-o', type='character', nargs = 1)
	return(parser)
}

#wrappers
io <- function(){
	parser <- io_base()
	return(parser$parse_args())
}

io_sheet <- function(){
	parser <- io_base()
	parser$add_argument('-sheet', type='integer', nargs='+')
	return(parser$parse_args())
}

# read in df from file + specified sheet, using integer header_row for column names
read_sheet <- function(file, sheet, header_row = 1){
	library(readxl)
	library(tidyr)
	library(stringr)
	library(tidyverse)
	df <- read_excel(file, sheet = sheet) %>% as.data.frame()
	colnames(df) <- fix_colnames(df[header_row,])
	df <- df[-(1:header_row),]
	rownames(df) <- 1:nrow(df)
	return(df)
}

fix_colnames <- function(names) {
	fixed_names <- names %>% str_replace_all(c(' ' = '_', '-' = '_'))
	print(fixed_names)
	return(fixed_names)
}

get_sig <- function(df, pval, fc) {
	library(tidyverse)
	return(dplyr::filter(df, as.numeric(!!as.name(pval)) < 0.05 & abs(as.numeric(!!as.name(fc))) > 1.5))
}

volcano_helper <- function(plot){
	source('../mcsf-gmcsf/src/plotutils.R')
	plot <- plot + 
		geom_point(size = 0.5, stroke = 0, shape = 16) +
		theme_bw() + 
		theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
		get_default_theme() + 
		xlab(expression('Average log'[2]*' fold change')) +
		ylab(expression('-log'[10]*'BH-corrected p value')) +
		geom_vline(xintercept = 1.5, linetype='dotted', lwd=0.1) + geom_vline(xintercept=-1.5, linetype='dotted', lwd=0.1)
	return(plot)
}

#bad fix for slow
load_rv_table <- function() {
	mtb_genes <- read.delim('data/Mycobacterium_tuberculosis_H37Rv_txt_v2.txt') %>% filter(Feature == 'CDS')
	rownames(mtb_genes) <- mtb_genes$Locus
	return(mtb_genes)
}
#slow because reading a full txt file
rv_to_gene <- function(rv, mtb_genes=load_rv_table()){
	return(as.vector(mtb_genes[as.vector(rv), 'Name']))
}

