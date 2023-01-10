library(argparse)

parser <- ArgumentParser(description='Compute pairwise jaccard/sorensen similarity')
parser$add_argument("--inc", help="input incidence freq matrix")
parser$add_argument("--jaccard", help="output jaccard matrix")
args <- parser$parse_args()

library(SpadeR)

inc <- read.csv(args$inc, header=TRUE, sep='\t', row.names=1)

jaccard <- data.frame(matrix(nrow=length(colnames(inc)), ncol=length(colnames(inc))))
colnames(jaccard) <- colnames(inc)
rownames(jaccard) <- colnames(inc)

for (i in colnames(inc)){
  for (j in colnames(inc)){
    if (i == j){
      jaccard[i,j] <- 1
    } else {
      spader <- SimilarityPair(inc[,c(i,j)], datatype="incidence_freq", nboot=100)
      jaccard[i,j] <- spader$Empirical_richness["U02(q=0,Jaccard)", "Estimate"]
      jaccard[j,i] <- spader$Empirical_richness["U02(q=0,Jaccard)", "Estimate"]
    }
  }
}

write.table(jaccard, file=args$jaccard, sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)