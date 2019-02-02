#!/usr/bin/Rscript

## Perform multiple testing correction for QTLtools output
## Diego Garrido-Martín
## 28/08/2018

## Load libraries
library(optparse)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

## Parse arguments
option_list = list(
  make_option(c("-n", "--nominal"), type = "character",
              help = "Result of QTLtools nominal pass", metavar = "character"),
  make_option(c("-p", "--permutation"), type = "character",
              help = "Result of QTLtools permuted pass [default %default]. Only required for method 'perm-fdr'",
	      default = NULL, metavar = "character"),
  make_option(c("-m", "--method"), type="character", metavar = "character",
              help = "One of 'bonferroni', 'fdr', 'perm-fdr'"),
  make_option(c("-a", "--alpha"), type="numeric", metavar = "numeric",
              help = "Significance level [default %default]", default = 0.05),
  make_option(c("-o", "--output"), type="character",
              metavar = "character", help = "Output pdf file [default %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$nominal) || is.null(opt$method) || is.null(opt$output)){
  print_help(opt_parser) 
  stop("Arguments 'nominal', 'method' and 'output' are required\n", call. = FALSE)
}

## Define functions
traceBack <- function(nominals.gene, permuted.df, p_t){
  permuted.gene <- subset(permuted.df, pheno_id == nominals.gene$pheno_id[1])
  p_tn <- qbeta(p = p_t, shape1 = permuted.gene$shape1, shape2 = permuted.gene$shape2)
  nominals.gene$p_tn <- p_tn
  return(subset(nominals.gene, pv <= p_tn))
}


## Run
set.seed(123)

nominals <- as.data.frame(fread(opt$nominal, header = F))
colnames(nominals) <- c("pheno_id", "pheno_chr", "pheno_start", "pheno_end",
                        "pheno_strand", "cis_vars", "dist", "var_id", "var_chr",
                        "var_start", "var_end", "pv", "slope", "isTop")
 
if (opt$method == "perm-fdr" ){
  if (is.null(opt$permutation)){
    stop("the method 'perm-fdr' requires a permutation pass output as input file")
  } 
  permutations <- as.data.frame(fread(opt$permutation, header = F))
  colnames(permutations) <- c("pheno_id", "pheno_chr", "pheno_start", "pheno_end",
                             "pheno_strand", "cis_vars", "dist", "topvar_id", 
                             "topvar_chr", "topvar_start", "topvar_end", "df", 
                             "dummy", "shape1", "shape2", "topvar_pv", "topvar_slope", "epv_p", "epv_b")               

  permutations$fdr <- p.adjust(permutations$epv_b, method = "fdr")
  set0 <- subset(permutations, fdr <= opt$alpha)
  set1 <- subset(permutations, fdr > opt$alpha)
  p_t <- (sort(set1$epv_b)[1] - sort(-1.0 * set0$epv_b)[1]) / 2
  fdr_pt <- (sort(set1$fdr)[1] - sort(-1.0 * set0$fdr)[1]) / 2
  permutations <- set0
  nominals <- subset(nominals, pheno_id %in% permutations$pheno_id)
  err <- abs(fdr_pt - opt$alpha)
  if (err > 0.1 * opt$alpha) {
	warning(sprintf("|closest observed FDR - FDR threshold set at %0.3f| = %0.2e.", opt$alpha, err))
  }
  message(sprintf("Global empirical P-value threshold = %0.2e", p_t))
  res <- as.data.frame(dplyr::do(dplyr::group_by(nominals, pheno_id), traceBack(., permutations, p_t)))
  res <- res[order(res$pv),]  
  write.table(res, file = opt$output, sep = "\t" ,col.names = TRUE, row.names = FALSE, quote = FALSE)
              
} else {
  nominals$pv_adj <- p.adjust(nominals$pv, method = opt$method)
  nominals <- subset(nominals, pv_adj <= opt$alpha)
  nominals <- nominals[order(nominals$pv_adj),]
  colnames(nominals)[colnames(nominals)=="pv_adj"] <- opt$method
  write.table(nominals, file = opt$output, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}



