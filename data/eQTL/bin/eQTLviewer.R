#!/usr/bin/Rscript

##  Plot eQTL 
##  Diego Garrido Martín
##  28/08/2018

##  Load libraries 

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))

##  Parse arguments

option_list = list(
  make_option(c("-i", "--input"), type="character",
              help="eQTLs to plot (mtc.R output format)", metavar="character"),
  make_option(c("-g", "--genotype"), type="character",
              help="Indexed genotype VCF file", metavar="character"),
  make_option(c("-e","--expression"), type="character",
              help="Normalized gene expression BED file", metavar="character"),
  make_option(c("-G","--Gene"), type="character", default = NULL,
              help="Gene of interest", metavar="character"),
  make_option(c("-s","--snp"), type="character", default = NULL,
              help="SNP of interest", metavar="character"),
  make_option(c("-o", "--output"), type="character", 
              help="Output PDF", metavar="character"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print progress [default %default]")
) 

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) 

if ( is.null(opt$input) || is.null(opt$genotype) || is.null(opt$expression) || is.null(opt$output) ){
  print_help(opt_parser)
  stop("Arguments 'input', 'genotype', 'expression' and 'output' are required\n", call.=FALSE)
} 

## Define functions
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

read.chunk <- function(gr, file, header = TRUE) {
  bed <- tryCatch(unlist(Rsamtools::scanTabix(file, param = GenomicRanges::reduce(gr)), 
                         use.names = FALSE, recursive = FALSE), error = function(e) c())
  if (length(bed) == 0) {
    return(NULL)
  }
  ncol <- length(strsplit(bed[1], "\t")[[1]])
  bed <- matrix(unlist(strsplit(bed, "\t", fixed = TRUE), 
                       use.names = FALSE, recursive = FALSE), length(bed), 
                ncol, byrow = TRUE)
  bed <- data.table::data.table(bed)
  if (header) {
    data.table::setnames(bed, gsub("#", "", colnames(suppressMessages(data.table::fread(paste("zcat <", file), 
                                                             nrows = 1, skip = "#CHROM")))))
  }
  bed <- bed[, lapply(.SD, function(ee) utils::type.convert(as.character(ee), 
                                                            as.is = TRUE))]
  bed <- as.data.frame(bed)
  return(bed)
}

##  Run

con <- OpenRead(opt$input)
tb <- read.table(con, header = TRUE)
close(con)

pdf(opt$output, paper = "a4r", width = 12, height = 8 )

if(!is.null(opt$Gene)){
  tb <- subset(tb, pheno_id == opt$Gene)
}

if(!is.null(opt$snp)){
 tb <- subset(tb, var_id == opt$snp)
}

GE <- as.data.frame(suppressMessages(fread(paste("zcat <", opt$expression), header = TRUE, skip = "#chr")))
colnames(GE) <- gsub("#", "", colnames(GE))

for (x in 1:nrow(tb)){
 
  # Extract fields of interest
  chr <- tb[x, "var_chr"]
  start <- tb[x, "var_start"]
  end <- tb[x, "var_end"]
  gene <- tb[x, "pheno_id"]
  snp <- tb[x, "var_id"]
  slope <- tb[x, "slope"]
  pv <- tb[x, "pv"]

  if(opt$verbose){
    cat(sprintf("Plot %s/%s: gene %s and snp %s\n", x, nrow(tb), gene, snp))
  } 
 
  # Variant
  subset.snp <- read.chunk(gr = GenomicRanges::GRanges(chr, IRanges::IRanges(start,end)),
                           file = opt$genotype)
  if(is.null(subset.snp)){
    warning(sprintf("snp %s is not in the genotype file.", snp))
  }

  ref <- gsub(TRUE, "T", subset.snp$REF) # Why it reads T as TRUE? Temporal fix: gsub
  alt <- gsub(TRUE, "T", subset.snp$ALT)

  # Expression
  ge <- GE

  # Generate df to plot
  subset.gene <- ge[which(ge$gene == gene), ]
  comm.inds <- intersect(colnames(subset.gene), colnames(subset.snp))             # Get common individuals
  subset.snp <- subset.snp[, comm.inds]
  subset.gene <- subset.gene[, comm.inds]
 
  df <- data.frame(t(subset.snp), t(subset.gene), stringsAsFactors = F)
  colnames(df) <- c("genotype", "norm_exp")
  df <- df[!grepl("\\.", df$genotype), ]
  df$genotype <- gsub("0[/|]1", sprintf("Heterozygous %s/%s", ref, alt) , df$genotype)
  df$genotype <- gsub("1[/|]0", sprintf("Heterozygous %s/%s", ref, alt) , df$genotype)
  df$genotype <- gsub("0[/|]0", sprintf("Dominant homozigous %s/%s", ref, ref) , df$genotype) 
  df$genotype <- gsub("1[/|]1", sprintf("Recessive homozigous %s/%s", alt, alt) , df$genotype)
   
  # Plot
  p <- ggplot(df, aes(x = factor(genotype), y = norm_exp, fill = "whatever")) +
       geom_boxplot() +
       geom_jitter(width = 0.2) + 
       theme_bw() +
       theme(text = element_text(size=18)) +
       xlab("Genotype") +
       ylab("Normalized gene expression") +
       ggtitle(sprintf("Variant %s is an eQTL for gene %s (pv = %0.4e, slope = %0.4f)", snp, gene, pv, slope)) +
       theme(plot.title = element_text(hjust = 0.5, size = 15, face="bold")) +
       guides(fill = FALSE)
  print(p)
}

null <- dev.off()

