library(argparse)

parser <- ArgumentParser(description='Plot Domains Richness accumulation curve')
parser$add_argument("--inext", help="Input inext table")
parser$add_argument("--plot", help="Output accumulation plot")
args <- parser$parse_args()

library(ggplot2)

df <- read.table(args$inext, header=TRUE, sep='\t')
df$t <- (df$t * 1e8)/1e9
df.point <- df[which(df$Method=="Observed"),]
df.line <- df[which(df$Method!="Observed"),]
df.line$Method <- factor(df.line$Method,
                         c("Rarefaction", "Extrapolation"),
                         c("Rarefaction", "Extrapolation"))

plt <- ggplot(df, aes(x=t, y=qD)) + theme_bw(base_size=12) + facet_wrap(~Assemblage, ncol=2, scales = "free_y") +
              geom_point(color="steelblue4", size=5, data=df.point) +
              geom_line(aes(linetype=Method), color="steelblue4", lwd=1.5, data=df.line) +
              geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL), fill="steelblue4", alpha=0.2) +
              labs(x='Sequencing depth (billions)', y='Family richness') +
              theme(legend.position="bottom", legend.title=element_blank(), legend.box="horizontal", legend.key.width=unit(1,"cm"), plot.title=element_text(size=10, hjust=1),
                       axis.title.x=element_text(size=12, face='plain', margin=margin(t=10)), axis.title.y=element_text(size=12, face='plain', margin=margin(r=10)))

ggsave(args$plot, plot=plt, width=7.02, height=7.02, units='in', dpi=300)