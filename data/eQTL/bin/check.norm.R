#!/usr/bin/Rscript

##  Check normalization
##  Diego Garrido Martín
##  14/09/2018

##  Load libraries 

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenABEL))
suppressPackageStartupMessages(library(preprocessCore))

##  Parse arguments

option_list = list(
  make_option(c("-i", "--input"), type="character",
              help="Gene expression file .gz (genes X samples)", metavar="character"),
  make_option(c("-o", "--output"), type="character", 
              help="Output PDF", metavar="character")
) 

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) 

if ( is.null(opt$input) ||is.null(opt$output) ){
  print_help(opt_parser)
  stop("Arguments 'input' and 'output' are required\n", call.=FALSE)
} 

##  Run

# Input gene expression
ge_file <- opt$input
if( grepl("\\.gz$", opt$input) ){
    ge_file <- paste0("zcat < '", ge_file, "'")
}

ge <- as.data.frame(suppressMessages(fread(ge_file, header = TRUE, skip = "#chr")))
colnames(ge) <- gsub("#", "", colnames(ge))

set.seed(123)
pdf(opt$output,  paper = 'a4r', width = 9, height = 6)
 par(mfrow=c(1,2))
 boxplot(ge[,c(7:17)], las = 2, arg.names = colnames(ge)[7:17], main = "Gene expression in samples 1-10")
 qqnorm(ge[2,-c(1:6)]); qqline(ge[1,-c(1:6)], col = 2)
dev <- dev.off()
