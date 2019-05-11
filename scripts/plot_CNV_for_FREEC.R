#!/usr/bin/env Rscript

##############################################################
#  script: plot_CNV_for_FREEC.R
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.08.06
#  description: plot CNV for FREEC 
#  example: Rscript --vanilla --slave plot_CNV_for_FREEC.R --input input.bam_ratio.txt --ploidy ploidy --genome_fai genome.fa.fai --output output.cnv_plot.pdf
##############################################################

library("ggplot2")
library("optparse")
library(scales)

option_list <- list(
  make_option(c("--input"), type = "character", default = NULL, 
              help = "input bam_ratio.txt file name", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL, 
              help = "the output file name", metavar = "character"),
  make_option(c("--ploidy"), type = "integer", default = NULL, 
              help = "the normal ploidy number", metavar = "character"),
  make_option(c("--genome_fai"), type = "character", default = NULL, 
              help = "path to the genome.fai file", metavar = "character")
); 

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
ploidy <- opt$ploidy

data <- read.table(opt$input, header = TRUE, sep = "\t")
data$Ratio[which(data$Ratio == -1)] <- NA
data$color[data$CopyNumber < ploidy] <- colors()[461]
data$color[data$CopyNumber == ploidy] <- colors()[88]
data$color[data$CopyNumber > ploidy] <- colors()[136]

genome_fai <- read.table(opt$genome_fai, header = FALSE, sep = "\t")
colnames(genome_fai) <- c("chr", "length", "offset", "linebases", "linewidth")

# re-order the chromosomes

chr_list <- as.vector(unique(genome_fai$chr))
genome_fai$chr <- factor(chr_list, levels = chr_list)
chr_list_intersect <- intersect(genome_fai$chr, data$Chromosome)
data$Chromosome <- factor(data$Chromosome, levels = chr_list_intersect)


# plot
ggplot(data, aes(x = (Start + End)/2, y = Ratio * ploidy)) + 
  geom_point(color = "grey", size = 0.05) +
  geom_segment(aes(x = Start, xend = End, y = CopyNumber, yend = CopyNumber, color = as.factor(color))) +	
  scale_color_manual(name = "", breaks = sort(unique(data$color)), values = sort(unique(data$color)), guide = FALSE) +
  scale_x_continuous(name = "Chromosome (bp)", labels = comma_format()) +
  scale_y_continuous(name = "Copy number profile") + 
  coord_cartesian(ylim = c(0, ploidy * 4)) +
  ggtitle(opt$input) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  facet_wrap(~ Chromosome, ncol = 1, scales = "free")
ggsave(filename = opt$output, device = "pdf", width = 8, height = length(chr_list), limitsize = FALSE)

