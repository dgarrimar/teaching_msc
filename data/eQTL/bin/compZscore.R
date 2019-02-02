#!/usr/bin/Rscript

## Compute Z score file given a nominal pass file and a gene name
## Diego Garrido-Martín
## 29/09/2018

## Load libraries
library(optparse)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

## Parse arguments
option_list = list(
  make_option(c("-g", "--gene"), type = "character",
              help = "Gene name", metavar = "character"),
  make_option(c("-k", "--k"), type = "numeric",
              help = "Number of top variants to consider [default %default]", 
              default = 100, metavar = "numeric"),
  make_option(c("-n", "--nominals"), type = "character",
              help = "Output of QTLtools nominal pass", metavar = "character"),
  make_option(c("-o", "--output"), type="character",
              metavar = "character", help = "Output Z score file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gene) || is.null(opt$nominal) || is.null (opt$output)){
  print_help(opt_parser) 
  stop("Arguments 'gene', 'nominal' and 'output' are required\n", call. = FALSE)
}

## Run
set.seed(123)
tb <- as.data.frame(fread(opt$nominal, header = T))[,c(1,8,12,13)]
colnames(tb) <- c("pheno_id", "var_id", "pv", "slope")
tb <- subset(tb, pheno_id == opt$gene) 
tb <- tb[order(tb$pv), ]
tb <- tb[1:opt$k,]

tb$Z <- qnorm(tb$pv, lower.tail = FALSE)*sign(tb$slope)
tb <- subset(tb, select = c("var_id", "Z"))
write.table(tb, file = opt$output, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
