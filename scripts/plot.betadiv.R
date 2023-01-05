library(argparse)

parser <- ArgumentParser(description='Plot spatio-temporal betadiversity')
parser$add_argument("--dist", help="Input distance table")
parser$add_argument("--plot", help="Output asyest plot")
args <- parser$parse_args()

library(ggplot2)

df <- read.table(args$dist, header=TRUE, sep='\t')
levels <- c("Tree of Life","Viruses","Archaea","Bacteria","Eukaryota")
plt <- ggplot(df, aes(x=factor(domain, levels=levels), y=bray, fill=type)) + theme_bw(base_size=12) +
              geom_boxplot(color="black", outlier.shape=NA, position=position_dodge(1))+
              #geom_jitter(position=position_jitter(0.2), color="darkgray") +
              scale_y_continuous(expand=c(0,0), limits=c(0,0.6)) +
              labs(y='Bray-Curtis dissimilarity') +
              theme(legend.position="bottom", legend.title=element_blank(), legend.box="horizontal", axis.title.x=element_blank(), axis.title.y=element_text(size=12, face='plain', margin=margin(r=10)))

ggsave(args$plot, plot=plt, width=7.02, height=4.33, units='in', dpi=300)