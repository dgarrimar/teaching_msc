#!/usr/bin/Rscript

##  Variant partition of expression given a set of covariates
##  Diego Garrido Martín
##  15/09/2018

##  Load libraries 

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(doParallel))

##  Parse arguments

option_list = list(
  make_option(c("-i", "--input"), type="character",
              help="Gene expression file", metavar="character"),
  make_option(c("-m","--metadata"), type="character",
              help="Metadata file", metavar="character"),
  make_option(c("-f","--formula"), type="character",
              help="formula for the linear model", metavar="character"),
  make_option(c("-o", "--output"), type="character", 
              help="Output pdf", metavar="character")
) 

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) 

if ( is.null(opt$input) || is.null(opt$metadata) || is.null(opt$formula) || is.null(opt$output) ){
  print_help(opt_parser)
  stop("Arguments 'input', 'metadata', 'formula' and 'output' are required\n", call.=FALSE)
} 

## Define funcions

OpenRead <- function(arg) {

  if (arg %in% c("-", "/dev/stdin")) {
    con <- file("stdin", open = "r", blocking = T)
  } else if (grepl("^/dev/fd/", arg)) {
    con <- fifo(arg, open = "r", blocking = T)
  } else {
    con <-file(arg, open = "r", blocking = T)
  }
  return (con)
}

## Run

# Input gene expression
ge_file <- opt$input
if( grepl("\\.gz$", opt$input) ){
    ge_file <- paste0("zcat < '", ge_file, "'")
}

ge <- as.data.frame(suppressMessages(fread(ge_file, header = TRUE, skip = "#chr")))
colnames(ge) <- gsub("#", "", colnames(ge))
rownames(ge) <- ge$gene
ge[, c(1:6)] <- NULL

# Input metadata
con <- OpenRead(opt$metadata)
metadata <- read.table(con, header = TRUE)
close(con)
rownames(metadata) <- metadata[,1]

# Parallel execution
cl <- makeCluster(4)
registerDoParallel(cl)

# Variance partition
form <- as.formula(opt$formula) 
varPart <- fitExtractVarPartModel( ge, form, metadata)
vp <- sortCols( varPart )

# Plot
pdf(opt$output,  paper = 'a4r', width = 9, height = 9)
 plotVarPart(vp)
dev <- dev.off()
