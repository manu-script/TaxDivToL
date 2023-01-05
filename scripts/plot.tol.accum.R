library(argparse)

parser <- ArgumentParser(description='Plot ToL Richness accumulation curve')
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

plt <- ggplot(df, aes(x=t, y=qD)) + theme_bw(base_size=8) +
              geom_point(color="black", size=3, data=df.point) +
              geom_line(aes(linetype=Method), color="black", lwd=1, data=df.line) +
              geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL), color="black", alpha=0.2) +
              scale_x_continuous(expand=c(0,0), limits=c(0,25)) +
              scale_y_continuous(expand=c(0,0), limits=c(800,1100)) +
              labs(x='Sequencing depth (billions)', y='Family richness') +
              theme(legend.position="top", legend.title=element_blank(), legend.box="horizontal", legend.key.width=unit(1,"cm"),
                    legend.margin=margin(0,0,-5,0),
                    axis.title.x=element_text(size=8, face='plain'), axis.title.y=element_text(size=8, face='plain'))

ggsave(args$plot, plot=plt, width=4.33, height=2.16, units='in', dpi=300)