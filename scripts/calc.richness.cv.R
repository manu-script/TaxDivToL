library(argparse)

parser <- ArgumentParser(description='Calculate coefficient of variation of richness')
parser$add_argument("--S27DEC19", help="Input stats table")
parser$add_argument("--S28DEC19", help="Input stats table")
parser$add_argument("--S29DEC19", help="Input stats table")
parser$add_argument("--cv", help="Output cv table")

args <- parser$parse_args()

S27DEC19 <- read.csv("diversity/richness/S27DEC19.richness.tsv", header=TRUE, sep="\t")
S28DEC19 <- read.csv("diversity/richness/S28DEC19.richness.tsv", header=TRUE, sep="\t")
S29DEC19 <- read.csv("diversity/richness/S29DEC19.richness.tsv", header=TRUE, sep="\t")

df <- rbind(rbind(S27DEC19, S28DEC19), S29DEC19)
df <- subset(df, depth==500)
mean_df <- aggregate(df$richness, by=list(taxon=df$taxon), mean)
colnames(mean_df)[colnames(mean_df) == "x"] <- "richness_mean"
sd_df <- aggregate(df$richness, by=list(taxon=df$taxon), sd)
colnames(sd_df)[colnames(sd_df) == "x"] <- "richness_sd"
df <- merge(mean_df, sd_df, by="taxon")
df$richness_cv <- df$richness_sd / df$richness_mean

write.table(df, file=args$cv, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
