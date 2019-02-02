#!/usr/bin/Rscript

## Output pdf with QQplot and hist of p-values
## Diego Garrido-Martín
## 28/08/2018

## Load libraries
library(optparse)
library(data.table)

## Parse arguments
option_list = list(
  make_option(c("-i", "--input"), type = "character",
              help = "Input table", metavar = "character"),
  make_option(c("-c", "--col"), type="numeric", metavar = "numeric",
              help = "Index of the column that contains the p-values"),
  make_option(c("-H", "--header"), action = "store_true", default = FALSE, 
	      help = "Does the input file have a header? [default %default]"),
  make_option(c("-o", "--output"), type="character",
              metavar = "character", default = "pvdist.pdf",
              help = "Output pdf file [default %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$col)){
  print_help(opt_parser)
  stop("Arguments 'input' and 'col' are required\n", call. = FALSE)
}

## Define functions
QQplot <- function(pvector, main = NULL, ...) {
  ## See http://gettinggeneticsdone.blogspot.com/p/copyright.html
  o <- -log10(sort(pvector, decreasing = FALSE))
  e <- -log10(1:length(o)/length(o))
  plot(e, o, pch = 19, cex = 1, main = main, ...,
       xlab = expression(Expected~~-log[10](italic(p))),
       ylab = expression(Observed~~-log[10](italic(p))),
       xlim = c(0, max(e)), ylim = c(0,max(o)))
  lines(e, e, col = "red")
}

## Run
set.seed(123)

tb <- as.data.frame(fread(opt$input, header = opt$header))
pv <- tb[, opt$col]

if (is.numeric(pv) && all(pv <= 1) && all(pv > 0)){
  pdf(opt$output,  paper = 'a4r', width = 9, height = 6)
  par(mfrow=c(1,2))
    QQplot(sample(pv, 10000)) # Sample to avoid many points in the plot and large pdf size
    hist(pv)
  dev <- dev.off()
} else{
  stop(sprintf("Column '%s' of file '%s' is not a numeric value in (0,1]!"), opt$col, opt$input)
}

