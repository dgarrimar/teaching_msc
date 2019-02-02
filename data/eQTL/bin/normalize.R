#!/usr/bin/Rscript

##  Normalize gene expression
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
              help="Gene expression file (genes X samples)", metavar="character"),
  make_option(c("-t","--threshold"), type="numeric", default = 0.1,
              help="Expression threshold [default %default]", metavar="numeric"),
  make_option(c("-o", "--output"), type="character", 
              help="Output file, normalized expression", metavar="character")
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

ge <- as.data.frame(fread(ge_file, header = TRUE, skip = "#chr"))
colnames(ge) <- gsub("#", "", colnames(ge))

# Filter
pass <- apply(ge[,-c(1:6)], 1, function(x) {median(x) > opt$threshold})
ge <- ge[pass, ]

# Normalize
X <- normalize.quantiles(as.matrix(ge[, -c(1:6)]))     # Quantile normalization (samples -> same distribution) check: boxplot(log10(X[,1:10]+1))
X <- t(apply(X, 1, rntransform))                       # Normalization (genes -> normal distribution)

colnames(X) <- colnames(ge)[-c(1:6)]
ge <- data.frame(ge[, c(1:6)], X)
colnames(ge)[1] <- paste("#", colnames(ge)[1], collapse = "", sep = "")

# Write output: normalized expression
write.table(ge, file = opt$output, quote = F, sep = "\t", col.names = T, row.names = F)
