library(argparse)

parser <- ArgumentParser(description='Phyloseq Analysis')
parser$add_argument("--otu", help="Path to count matrix")
parser$add_argument("--tax", help="Path to taxanomy table")
parser$add_argument("--meta", help="Path to metadata")
parser$add_argument("--domain", help="Select domian")
parser$add_argument("--nmds", help="Output NMDS plot")
parser$add_argument("--dist", help="Output distance matrix")
args <- parser$parse_args()

library(tidyverse)
library(phyloseq)
library(ggplot2)

otu <- read.table(args$otu, sep = "\t", header = TRUE, row.names=1)

tax <- read.table(args$tax, sep = "\t", header = TRUE, row.names=1)
tax <- tax %>% mutate_if(is.character, as.factor)

meta <- read.table(args$meta, sep = "\t", header = TRUE, row.names=1)
meta <- meta %>% mutate_if(is.character, as.factor)
meta$season <- factor(meta$season, levels = c("Summer", "Monsoon", "Winter"))
meta$sector <- factor(meta$sector, levels = c("Northern", "Central", "Southern"))
meta$month <- factor(meta$month, levels = c("DEC", "MAR", "JUL", "NOV"))

otu_mat <- as.matrix(otu)
tax_mat <- as.matrix(tax)

phylo_OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX <- tax_table(tax_mat)
phylo_samples <- sample_data(meta)

chilika <- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)

chilika <- subset_samples(chilika, year==2020)

domain <- paste0("d__", args$domain)

if (args$domain != "ToL"){
  chilika <- subset_taxa(chilika, Domain==domain)
}

chilika.rel <- transform_sample_counts(chilika, function(x) x / sum(x))

chilika.dist <- distance(chilika.rel, method="bray")

chilika.ord <- ordinate(chilika.rel, "NMDS", distance=chilika.dist)

chilika.nmds <- plot_ordination(chilika.rel, chilika.ord, type="samples", color="season")

chilika.nmds <- chilika.nmds + theme_bw(base_size=10) +
                geom_point(size=2) +
                geom_text(mapping=aes(label=station), size=3, show.legend=FALSE, vjust=1.75) +
                scale_colour_manual(labels=c("Summer", "Monsoon", "Winter"), values=c('darksalmon', 'darkturquoise', 'darkseagreen')) +
                theme(legend.position="bottom", legend.title=element_blank(), legend.box="horizontal")

ggsave(args$nmds, plot=chilika.nmds, width=4.33, height=4.33, units='in', dpi=300)
write.table(as.data.frame(as.matrix(chilika.dist)), file=args$dist, sep='\t', quote=FALSE, col.names=NA)
