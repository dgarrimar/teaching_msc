#!/usr/bin/Rscript

##  Compute PEER factors from gene expression
##  Diego Garrido Martín
##  15/09/2018

##  Load libraries

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(peer))

option_list = list(
  make_option(c("-i", "--input"), type="character",
              help="Gene expression matrix", metavar="FILE"),
  make_option(c("-o", "--output"), type="character",
              help="PEER factors", metavar="FILE"),
  make_option(c("-p", "--peer"), type="numeric",
              help="Number of PEER factors",
              metavar="INTEGER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input
output_file <- opt$output
p <- opt$peer

if ( is.null(input_file) || is.null(output_file) || is.null(p)){
  print_help(opt_parser)
  stop("Arguments 'input', 'output' and 'peer' are required\n", call.=FALSE)
}

## Set seed
set.seed(1)

## Get PEER factors
dat <- read.table(gzfile(input_file), header = FALSE)
header <- unlist(strsplit(readLines(gzfile(input_file))[[1]], split = "\t"))

dat <- dat[,-c(1:6)]
samples <- header[-c(1:6)]
colnames(dat)<-samples

## Build PEER model
model <- PEER()
PEER_setPhenoMean(model, as.matrix(t(dat)))
PEER_setNk(model, p)

## Iterate
PEER_update(model)

## Results
factors <- PEER_getX(model)
# weights <- PEER_getW(model)
# precision <- PEER_getAlpha(model)
# residuals <- PEER_getResiduals(model)

colnames(factors) <- paste0("PEER",1:p)
factors <- data.frame(sampleID = samples, factors)
factors <- factors[order(factors$sampleID), ]

write.table(factors, file = output_file, sep ="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
