# Title     : plot_tree_differences
# Objective : visualisation of distance
# Created by: Niek_Huijsmans
# Created on: 9-6-2020
library("optparse")
## Load treespace, phytools ad adegenet
library("treespace") # for plottreediff function
#library("phytools") # to read in phylo tree
#library("adegenet") # for coloring of leaves
##options to Rscript
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="input dataset file name", metavar="character")
  make_option(c("-r", "--real"), type="character", default=NULL,
              help="input dataset of real tree", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$file) || is.null(opt$real)) {
  print_help(opt_parser)
  stop("At least two argument must be supplied (input file and real tree)", call.=FALSE)
}### Plotting tree differences in R using treespace

## Load in the trees you want to compare
## In this case, a tree of 723 strains consturcted using prokka, Roary, RAxML with normal GTRGAMMA (tree.normal) and ASC_GTRGAMMA (tree.asc)
tree.normal <- read.tree(opt$file)
tree.asc <- read.tree(opt$real)

### Compute differences between the trees beforehand
### This will save you a lot of time when trying to get graphical settings right
#wmTipDiff <- tipDiff(tree.normal,tree.asc, sizeOfDifferences=TRUE)
#
### Plot the difference between trees
#plotTreeDiff(tree.normal,tree.asc,tipDiff = wmTipDiff)
#
### Plot the difference with trees facing each other
#plotTreeDiff(tree.normal,tree.asc,tipDiff = wmTipDiff, treesFacing = TRUE)
#
### Plot the difference between phylograms (type), with:
### smaller font sizes (cex)
### matched tips without lines drawn between 'em (align.tip.label = 0)
### Color according to the palette spectral (see ?num2col for a list)
#plotTreeDiff(tree.normal, tree.asc,
#             tipDiff = wmTipDiff,
#             cex=0.05,
#             font=1,
#             tipMatch = TRUE,
#             type  = "phylogram",
#             #align.tip.label = 0,
#             treesFacing = TRUE,
#             colourMethod = "palette",
#             palette = wasp)

## compute the KC metric between two trees
treeDist(tree.normal, tree.asc)