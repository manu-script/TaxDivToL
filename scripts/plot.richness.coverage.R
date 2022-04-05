library(argparse)

parser <- ArgumentParser(description='Plot Richness accumulation and coverage')
parser$add_argument("--stats", help="Input stats table")
parser$add_argument("--rich", help="Output richness plot")
parser$add_argument("--cov", help="Output coverage plot")

args <- parser$parse_args()

library(ggplot2)
library(scales)
options(bitmapType = 'cairo', device = 'png')

stats <- read.csv(args$stats, header = TRUE, sep = '\t')

rich <- subset(stats, taxon!="Tree of Life")

rich_plt <- ggplot(rich, aes(x=depth, y=richness)) + theme_bw(base_size=12) +
                    geom_smooth(method = "loess", span=1.4, alpha=0.2, color="black") +
                    facet_wrap(~ taxon, ncol=3, scales="free_y") +
                    scale_x_continuous(expand=c(0,0),breaks=c(0,500,1000), limits=c(00,1050)) +
                    scale_y_continuous(breaks=pretty_breaks()) +
                    labs(x='Sequencing depth (in millions of unique fragments)', y='Family richness', axis.title.x=element_text(size=12, face='plain'), axis.title.y=element_text(size=12, face='plain'))

ggsave(args$rich, plot=rich_plt, width=169, height=240, units='mm', dpi=300)

cov <- subset(stats, taxon=="Tree of Life")

cov_plt <- ggplot(cov, aes(x=depth, y=coverage)) + theme_bw(base_size=12) +
                  geom_smooth(method = "loess", span=1.4, color="black", alpha=0.2) +
                  scale_x_continuous(expand=c(0,0),breaks=seq(0,1000,100), limits=c(0,1025)) +
                  scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10), limits=c(-10,100)) +
                  coord_cartesian(xlim=c(0,1025), ylim=c(0,100)) +
                  labs(x="Sequencing depth (in millions of unique fragments)", y="Fraction of family richness") +
                  theme(legend.title=element_blank(), legend.position="none", axis.title.x=element_text(size=12, face="plain"),
                                      axis.title.y=element_text(size=12, face="plain"))

ggsave(args$cov, plot=cov_plt, width=169, height=80, units="mm", dpi=300)
