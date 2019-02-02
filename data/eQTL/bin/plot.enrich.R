#!/usr/bin/Rscript

##  Plot enrichment in functional annotations
##  Diego Garrido Martín
##  23/09/2018

##  Load libraries 

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))

##  Parse arguments

option_list = list(
  make_option(c("-i", "--input"), type="character",
              help="Input file", metavar="character"),
  make_option(c("-o", "--output"), type="character", 
              help="Output enrichment plot in pdf", metavar="character")
) 

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) 

if ( is.null(opt$input) || is.null(opt$output) ){
  print_help(opt_parser)
  stop("Arguments 'input' and 'output' are required\n", call.=FALSE)
} 

##  Run

# Input enrichment

D <- read.table(opt$input)

# Perform Fisher test
df <- c()
for (i in 1:nrow(D)){
     f <- fisher.test(matrix(c(D[i, 1], D[i, 2], round(D[i, 3]), D[i, 2]), ncol=2))
     df <- rbind(df, data.frame(feat = D[i, 9], or = f$estimate, ci.lower = f$conf.int[1], ci.upper = f$conf.int[2], pv = f$p.value))
}

# Perform FDR
df$fdr <- p.adjust(df$pv, method = "fdr")
row.names(df) <- NULL
df <- subset(df, fdr < 0.1)

# Plot
kmax <- max(df$ci.upper)
kmin <- min(df$ci.lower)
x <- with(df, reorder(feat, or , identity))


p <- ggplot(data = df, aes(label = feat)) +
  geom_point(aes(x = x, y = or, size = -log10(pv)))+
  geom_errorbar(aes(x = x, ymin = ci.lower, ymax = ci.upper), width = 0.5) +
  ylim(kmin, kmax) +
  coord_flip() +
  xlab("") +
  ylab("Odds Ratio") +
  geom_abline(slope = 0, intercept = 1, col = "red") +
  theme_bw() +
  labs(size = expression(-log["10"](p-value))) +
  theme(text = element_text(size=20))

pdf(opt$output,  paper = 'a4r', width = 9, height = 9)
p
null <- dev.off()
