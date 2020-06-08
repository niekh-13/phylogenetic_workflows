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
file2 <- opt$file
df <- read.delim(file = file2, header = TRUE)
library("ggplot2")
library("plyr")

#the_theme <- theme(
#  axis.title.x = element_text(size = 15, family = "Verdana"),
#  axis.text.x = element_text(size = 13, family = "Verdana"),
#  axis.title.y = element_text(size = 15, family = "Verdana"),
#  axis.text.y = element_text(size=13, family = "Verdana"))

c_meds <- ddply(df, .(assembler), summarise, med = median(complete_features))
complete_features <- ggplot(data=df, aes(x=assembler, y=complete_features)) +
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  xlab("Assemblers") +
  ylab("Complete features (%)") +
  geom_text(data = c_meds, aes(x = assembler, y = med, label = med), size = 5, vjust = -1) +
  theme_bw() +
  theme(text=element_text(size=16, colour = "black"))

pdf("complete_features.pdf", width = 6, height = 4)
complete_features
dev.off()
#ggsave(plot=complete_features, filename="complete_features.pdf", complete_features, width=3.5, height=3.5, dpi=320)
#ggsave(plot=complete_features, filename="complete_features.png", complete_features, width=3.5, height=3.5, dpi=320)
#unlink("complete_features.png")
#unlink("complete_features.pdf")

t_meds <- ddply(df, .(assembler), summarise, med = median(total_features))
total_features <- ggplot(data=df, aes(x=assembler, y=total_features)) +
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  xlab("Assembelers") +
  ylab("Total features (%)") +
  geom_text(data = t_meds, aes(x = assembler, y = med, label = med), size = 5, vjust = -1) +
  theme_bw() +
  theme(text=element_text(size=16, colour = "black"))

pdf("total_features.pdf", width = 6, height = 4)
total_features
dev.off()
#ggsave(plot=total_features, filename="total_features.pdf", width=3.5, height=3.5, dpi=320)
#ggsave(plot=total_features, filename="total_features.png", width=3.5, height=3.5, dpi=320)
#unlink("total_features.png")
#unlink("total_features.pdf")

N_meds <- ddply(df, .(assembler), summarise, med = median(NGA50))
NGA <- ggplot(data=df, aes(x=assembler, y=NGA50)) +
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  xlab("Assemblers") +
  ylab("NGA50") +
  scale_y_continuous() +
  geom_text(data = N_meds, aes(x = assembler, y = med, label = med), size = 5, vjust = -1) +
  theme_bw() +
  theme(text=element_text(size=16, colour = "black"))

pdf("NGA50.pdf", width = 6, height = 4)
NGA
dev.off()
#ggsave(plot = NGA, filename="NGA50.pdf", width=3.5, height=3.5, dpi=320)
#ggsave(plot = NGA, filename="NGA50.png", width=3.5, height=3.5, dpi=320)
#unlink("NGA50.pdf")
#unlink("NGA50.png")

contigs_meds <- ddply(df, .(assembler), summarise, med = median(contigs))
contigs <- ggplot(data=df, aes(x=assembler, y=contigs)) +
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  xlab("Assembelers") +
  ylab("#contigs") +
  geom_text(data = contigs_meds, aes(x = assembler, y = med, label = med), size = 5, vjust = -1) +
  theme_bw() +
  theme(text=element_text(size=16, colour = "black"))

pdf("contigs.pdf", width = 6, height = 4)
contigs
dev.off()
#ggsave(plot = contigs, filename="contigs.pdf", width=3.5, height=3.5, dpi=320)
#ggsave(plots = contigs, filename="contigs.png", width=3.5, height=3.5, dpi=320)
#unlink("contigs.pdf")
#unlink("contigs.png")

fraction_meds <- ddply(df, .(assembler), summarise, med = median(genome_fraction))
fraction <- ggplot(data=df, aes(x=assembler, y=genome_fraction)) +
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  xlab("Assembelers") +
  ylab("genome fraction (%)") +
  geom_text(data = fraction_meds, aes(x = assembler, y = med, label = med), size = 5, vjust = -1) +
  theme_bw() +
  theme(text=element_text(size=16, colour = "black"))

pdf("genome_fraction.pdf", width = 6, height = 4)
fraction
dev.off()
#ggsave(plot = fraction, file="genome_fraction.pdf", width=3.5, height=3.5, dpi=320)
#ggsave(plot = fraction, file="genome_fraction.png",  width=3.5, height=3.5, dpi=320)
#unlink("genome_fraction.pdf")
#unlink("genome_fraction.png")

fraction_meds <- ddply(df, .(assembler), summarise, med = median(misassemblies))
misassemblies <- ggplot(data=df, aes(x=assembler, y=misassemblies)) +
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  xlab("Assembelers") +
  ylab("#misassemblies") +
  geom_text(data = fraction_meds, aes(x = assembler, y = med, label = med), size = 5, vjust = -0.5) +
  theme_bw() +
  theme(text=element_text(size=16, colour = "black"))

pdf("misassemblies.pdf", width = 6, height = 4)
misassemblies
dev.off()
#ggsave(filename="misassemblies.pdf", width=3.5, height=3.5, dpi=320)
#ggsave(filename="misassemblies.png", width=3.5, height=3.5, dpi=320)
#unlink("missassemblies.pdf")
#unlink("missassemblies.png")