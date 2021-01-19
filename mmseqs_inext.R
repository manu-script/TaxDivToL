library(iNEXT)
library(ggplot2)
library(ggpubr)
library(ggthemes)
options(bitmapType = 'cairo', device = 'png')

Viruses_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Viruses_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		    		S28=as.matrix(read.csv('mmseqs/S28DEC19/Viruses_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		    		S29=as.matrix(read.csv('mmseqs/S29DEC19/Viruses_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Viruses_inext <- iNEXT(Viruses_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Viruses_type1 <- ggiNEXT(Viruses_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
						scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) +
						scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
						theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
						scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
						scale_y_continuous(breaks=c(40, 42, 44, 46, 48), limits=c(40, 48)) +
						coord_cartesian(xlim=c(0,30), ylim=c(40,48))
#------------------------------------------------------------------------------------------------------------------------------------

Archaea_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Archaea_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
                    S28=as.matrix(read.csv('mmseqs/S28DEC19/Archaea_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		    		S29=as.matrix(read.csv('mmseqs/S29DEC19/Archaea_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Archaea_inext <- iNEXT(Archaea_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Archaea_type1 <- ggiNEXT(Archaea_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
						scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
		        		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
			 			scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
						scale_y_continuous(breaks=c(48, 52, 56, 60, 64), limits=c(48, 66)) +
						coord_cartesian(xlim=c(0,30), ylim=c(48,66))
#------------------------------------------------------------------------------------------------------------------------------------

Bacteria_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Bacteria_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		     S28=as.matrix(read.csv('mmseqs/S28DEC19/Bacteria_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		     S29=as.matrix(read.csv('mmseqs/S29DEC19/Bacteria_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Bacteria_inext <- iNEXT(Bacteria_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Bacteria_type1 <- ggiNEXT(Bacteria_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
						scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
						theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
						scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
			  			scale_y_continuous(breaks=c(620, 640, 660, 680), limits=c(620, 681)) +
						coord_cartesian(xlim=c(0,30), ylim=c(620,681))

#------------------------------------------------------------------------------------------------------------------------------------

Eukaryota_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Eukaryota_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		      S28=as.matrix(read.csv('mmseqs/S28DEC19/Eukaryota_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		      S29=as.matrix(read.csv('mmseqs/S29DEC19/Eukaryota_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Eukaryota_inext <- iNEXT(Eukaryota_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Eukaryota_type1 <- ggiNEXT(Eukaryota_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
						scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
			   			theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
			    		scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
						scale_y_continuous(breaks=c(600, 700, 800, 900, 1000), limits=c(600, 1050)) +
						coord_cartesian(xlim=c(0,30), ylim=c(600,1050))

#------------------------------------------------------------------------------------------------------------------------------------

Protists_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Protists_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		     S28=as.matrix(read.csv('mmseqs/S28DEC19/Protists_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		     S29=as.matrix(read.csv('mmseqs/S29DEC19/Protists_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Protists_inext <- iNEXT(Protists_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Protists_type1 <- ggiNEXT(Protists_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
							scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
							theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
			 				scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
							scale_y_continuous(breaks=c(150, 200, 250, 300, 350), limits=c(150, 350)) +
							coord_cartesian(xlim=c(0,30), ylim=c(150,350))
#------------------------------------------------------------------------------------------------------------------------------------
Fungi_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Fungi_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
	          S28=as.matrix(read.csv('mmseqs/S28DEC19/Fungi_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
	          S29=as.matrix(read.csv('mmseqs/S29DEC19/Fungi_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Fungi_inext <- iNEXT(Fungi_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Fungi_type1 <- ggiNEXT(Fungi_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
						scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
						theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
			 			scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
						scale_y_continuous(breaks=c(140, 160, 180, 200), limits=c(140, 202)) +
						coord_cartesian(xlim=c(0,30), ylim=c(140,202))
#------------------------------------------------------------------------------------------------------------------------------------

Viridiplantae_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Viridiplantae_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		          S28=as.matrix(read.csv('mmseqs/S28DEC19/Viridiplantae_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		          S29=as.matrix(read.csv('mmseqs/S29DEC19/Viridiplantae_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Viridiplantae_inext <- iNEXT(Viridiplantae_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Viridiplantae_type1 <- ggiNEXT(Viridiplantae_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
							scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
							theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
							scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
							scale_y_continuous(breaks=c(100, 125, 150, 175, 200), limits=c(100, 200)) +
							coord_cartesian(xlim=c(0,30), ylim=c(100,200))
#------------------------------------------------------------------------------------------------------------------------------------

Metazoa_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Metazoa_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		    S28=as.matrix(read.csv('mmseqs/S28DEC19/Metazoa_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		    S29=as.matrix(read.csv('mmseqs/S29DEC19/Metazoa_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Metazoa_inext <- iNEXT(Metazoa_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Metazoa_type1 <- ggiNEXT(Metazoa_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
						scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
			 			theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
			  			scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
			 			scale_y_continuous(breaks=c(150, 200, 250, 300), limits=c(150, 305)) +
						coord_cartesian(xlim=c(0,30), ylim=c(150,305))
#------------------------------------------------------------------------------------------------------------------------------------

Arthropoda_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Arthropoda_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		       S28=as.matrix(read.csv('mmseqs/S28DEC19/Arthropoda_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		       S29=as.matrix(read.csv('mmseqs/S29DEC19/Arthropoda_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Arthropoda_inext <- iNEXT(Arthropoda_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Arthropoda_type1 <- ggiNEXT(Arthropoda_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
							scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
                            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
			     			scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
			    			scale_y_continuous(breaks=c(60, 70, 80, 90, 100), limits=c(60, 102)) +
							coord_cartesian(xlim=c(0,30), ylim=c(60,102))
#------------------------------------------------------------------------------------------------------------------------------------

Actinopteri_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Actinopteri_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		        S28=as.matrix(read.csv('mmseqs/S28DEC19/Actinopteri_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		        S29=as.matrix(read.csv('mmseqs/S29DEC19/Actinopteri_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Actinopteri_inext <- iNEXT(Actinopteri_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Actinopteri_type1 <- ggiNEXT(Actinopteri_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
							 scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
							 theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
							 scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
							 scale_y_continuous(breaks=c(25, 35, 45, 55, 65), limits=c(25, 69)) +
							 coord_cartesian(xlim=c(0,30), ylim=c(25,69))
#------------------------------------------------------------------------------------------------------------------------------------

Mammalia_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Mammalia_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		     S28=as.matrix(read.csv('mmseqs/S28DEC19/Mammalia_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
		     S29=as.matrix(read.csv('mmseqs/S29DEC19/Mammalia_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Mammalia_inext <- iNEXT(Mammalia_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Mammalia_type1 <- ggiNEXT(Mammalia_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
						scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) + scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
			  			theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
			  			scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
			  			scale_y_continuous(breaks=c(15, 20, 25, 30, 35), limits=c(15, 37.5)) +
						coord_cartesian(xlim=c(0,30), ylim=c(15,37.5))
#------------------------------------------------------------------------------------------------------------------------------------

Aves_inc <- list(S27=as.matrix(read.csv('mmseqs/S27DEC19/Aves_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
	         S28=as.matrix(read.csv('mmseqs/S28DEC19/Aves_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)),
	         S29=as.matrix(read.csv('mmseqs/S29DEC19/Aves_genus_level.tsv', header = TRUE, sep = '\t', row.names=1)))

Aves_inext <- iNEXT(Aves_inc, q=0, datatype="incidence_raw", endpoint=30, knots=30, se=TRUE, conf=0.95, nboot=100)

Aves_type1 <- ggiNEXT(Aves_inext, type=1, color.var='site') + theme_bw(base_size = 14) +
					scale_colour_manual(values=c('#0072B2','#009E73','#D55E00')) +
					scale_fill_manual(values=c('#0072B2','#009E73','#D55E00')) +
		      		theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
		      		scale_x_continuous(breaks=c(0,10,20,30), labels=c(0, 1, 2, 3)) +
		      		scale_y_continuous(breaks=c(4,8,12,16), limits=c(4,16)) +
					coord_cartesian(xlim=c(0,30), ylim=c(4,16))

#------------------------------------------------------------------------------------------------------------------------------------
type1 <- annotate_figure(ggarrange(Viruses_type1, Archaea_type1, Bacteria_type1, Eukaryota_type1, Protists_type1, Fungi_type1, Viridiplantae_type1,
			   Metazoa_type1, Arthropoda_type1, Actinopteri_type1, Mammalia_type1, Aves_type1, ncol=3, nrow=4, common.legend=TRUE,
			   legend="top", labels=c('a. Viruses', 'b. Archaea', 'c. Bacteria', 'd. Eukaryota', 'e. Protists', 'f. Fungi',
						     'g. Viridiplantae', 'h. Metazoa', 'i. Arthropoda', 'j. Actinopteri', 'k. Mammalia', 'l. Aves'),
			   font.label=list(size=12, face="plain", color="black"), align='hv', widths=c(1,1,1), heights=c(1,1,1,1),
			   label.x=c(0.15,0.15,0.15, 0.1,0.15,0.2, 0.05,0.15,0.1, 0.1,0.1,0.2),
			   label.y=0.3,
			   legend.grob=get_legend(Aves_type1 + theme(legend.box="horizontal", legend.title=element_blank(),
								     legend.position="top", legend.margin=margin()))),
						 left=text_grob('Genus richness', face='plain', size=14, rot = 90),
						 bottom=text_grob('Sequencing depth (in billions of unique PE reads)', face='plain', size=14))

ggsave('Type1.png', plot=type1, path='mmseqs', width=169, height=230, units='mm', dpi=300)
#------------------------------------------------------------------------------------------------------------------------------------
coverage <- function(taxon, inext){
	S27 <- (inext$iNextEst$S27$qD / subset(inext$AsyEst, Site=='S27' & Diversity=='Species richness')$Estimator) * 100
	S27.upr <- (inext$iNextEst$S27$qD.UCL / subset(inext$AsyEst, Site=='S27' & Diversity=='Species richness')$Estimator) * 100
	S27.lwr <- (inext$iNextEst$S27$qD.LCL / subset(inext$AsyEst, Site=='S27' & Diversity=='Species richness')$Estimator) * 100
	S27.df <- data.frame(taxon=taxon, depth=1:30, sample='S27', cov=S27, cov.lwr=S27.lwr, cov.upr=S27.upr)
	S28 <- (inext$iNextEst$S28$qD / subset(inext$AsyEst, Site=='S28' & Diversity=='Species richness')$Estimator) * 100
	S28.upr <- (inext$iNextEst$S28$qD.UCL / subset(inext$AsyEst, Site=='S28' & Diversity=='Species richness')$Estimator) * 100
	S28.lwr <- (inext$iNextEst$S28$qD.LCL / subset(inext$AsyEst, Site=='S28' & Diversity=='Species richness')$Estimator) * 100
	S28.df <- data.frame(taxon=taxon, depth=1:30, sample='S28', cov=S28, cov.lwr=S28.lwr, cov.upr=S28.upr)
	S29 <- (inext$iNextEst$S29$qD / subset(inext$AsyEst, Site=='S29' & Diversity=='Species richness')$Estimator) * 100
	S29.upr <- (inext$iNextEst$S29$qD.UCL / subset(inext$AsyEst, Site=='S29' & Diversity=='Species richness')$Estimator) * 100
	S29.lwr <- (inext$iNextEst$S29$qD.LCL / subset(inext$AsyEst, Site=='S29' & Diversity=='Species richness')$Estimator) * 100
	S29.df <- data.frame(taxon=taxon, depth=1:30, sample='S29', cov=S29, cov.lwr=S29.lwr, cov.upr=S29.upr)
	return(rbind(rbind(S27.df, S28.df),S29.df))
}

dat_type2 <- coverage('Viruses', Viruses_inext)
dat_type2 <- rbind(dat_type2, coverage('Archaea', Archaea_inext))
dat_type2 <- rbind(dat_type2, coverage('Bacteria', Bacteria_inext))
dat_type2 <- rbind(dat_type2, coverage('Eukaryota', Eukaryota_inext))
dat_type2 <- rbind(dat_type2, coverage('Protists', Protists_inext))
dat_type2 <- rbind(dat_type2, coverage('Fungi', Fungi_inext))
dat_type2 <- rbind(dat_type2, coverage('Viridiplantae', Viridiplantae_inext))
dat_type2 <- rbind(dat_type2, coverage('Metazoa', Metazoa_inext))
dat_type2 <- rbind(dat_type2, coverage('Arthropoda', Arthropoda_inext))
dat_type2 <- rbind(dat_type2, coverage('Actinopteri', Actinopteri_inext))
dat_type2 <- rbind(dat_type2, coverage('Mammalia', Mammalia_inext))
dat_type2 <- rbind(dat_type2, coverage('Aves', Aves_inext))

cbpalette <- c("#24ff24","#009292","#ff6db6","#490092","#000000","#006ddb","#b66dff","#b6dbff", "#920000","#924900","#db6d00","#ffff6d")
type2 <- ggplot(dat_type2, aes(x=depth, y=cov, color=taxon)) + theme_bw(base_size = 12) + geom_smooth(aes(color=taxon, fill=taxon), alpha=0.2) +
				scale_colour_manual(values=cbpalette) + scale_fill_manual(values=cbpalette) +
				scale_x_continuous(expand=c(0,0),breaks=c(5,10,15,20,25), labels=c(0.5, 1, 1.5, 2, 2.5), limits=c(0,35)) +
				scale_y_continuous(expand=c(0,0),breaks=c(65,70,75,80,85,90,95,100), limits=c(65,105)) +
				coord_cartesian(xlim=c(1,26), ylim=c(65,100)) +
				labs(x='Sequencing depth (in billions of unique PE reads)', y='Sample coverage (at the Genus level)') +
				theme(legend.title=element_blank(), legend.position="bottom", axis.title.x=element_text(size=12, face='plain'),
				      axis.title.y=element_text(size=12, face='plain'))


ggsave('Type2.png', plot=type2, path='mmseqs', width=169, height=120, units='mm', dpi=300)





