library(argparse)

parser <- ArgumentParser(description='Plot Asymptotic estimates')
parser$add_argument("--asyest", help="Input asyest table")
parser$add_argument("--plot", help="Output asyest plot")
args <- parser$parse_args()

library(ggplot2)

df <- read.table(args$asyest, header=TRUE, sep='\t')
df$Estimator <- round(df$Estimator)
levels <- c("Viruses","Archaea","Bacteria","Eukaryota","Metazoa","Fungi","Viridiplantae","Arthropoda","Chordata","Actinopteri","Mammalia","Aves","Amphibia")

plt <- ggplot(df, aes(x=factor(Assemblage, level=levels), y=Estimator)) + theme_bw(base_size=8) +
              geom_bar(stat="identity") + geom_errorbar(aes(ymin=Estimator-s.e., ymax=Estimator+s.e.), width=0.2) +
              geom_text(aes(y=Estimator-s.e., label=Estimator), vjust=1.5, color="white", size=2) +
              scale_y_continuous(expand=c(0,0), breaks=c(1,10,100,1000), limits=c(1,1000), trans="log10") +
              annotation_logticks(sides="l", outside = TRUE) + coord_cartesian(clip = "off") +
              labs(y='Asymptotic family richness') +
              theme(legend.position="none", axis.text.x=element_text(size=8, angle=45, hjust=1), axis.title.x=element_blank(),
                    axis.title.y=element_text(size=8, face='plain'), panel.grid.minor = element_blank(),
                    axis.text.y=element_text(margin=margin(r=10)))

ggsave(args$plot, plot=plt, width=4.33, height=2.165, units='in', dpi=300)