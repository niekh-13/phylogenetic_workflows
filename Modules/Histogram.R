# Title     : Histogram visualisation
# Objective : visualisation with ggplot2
# Created by: Niek_Huijsmans
# Created on: 9-4-2020
library("optparse")
##options to Rscript
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="input dataset file name", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

library("tools")
file2 <- file_path_as_absolute(opt$file)
df <- read.delim(file = file2, header = FALSE)
library("ggplot2")
ggplot(data=df, aes(x=V3)) +
  geom_histogram(breaks=seq(75,100, by = 0.5))+
  labs(title="Histogram for ANI", x="ANI(%)", y="Count")
ggsave("histogram.png")
