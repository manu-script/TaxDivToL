library(iNEXT)
library(argparse)

parser <- ArgumentParser(description='Estimate asymptotic diversity per sample per taxon')
parser$add_argument("--taxon", help="taxon name")
parser$add_argument("--inc", help="input incidences")
parser$add_argument("--inext", help="output intrapolation and extrapolation of diversity")
parser$add_argument("--asyest", help="output asymptotic diversity estimates")

args <- parser$parse_args()

df <- read.table(args$inc, header=FALSE, sep='\t', row.names=1)

inc <- list()
inc[[args$taxon]] <- df$V2
inext <- iNEXT(inc, q=0, datatype="incidence_freq", se=TRUE, conf=0.95, nboot=100, size=seq(1, 250, by=1))

write.table(inext$iNextEst$size_based, file=args$inext, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(inext$AsyEst, file=args$asyest, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)