library(SpadeR)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggthemes)
options(bitmapType = 'cairo', device = 'png')

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Viruses_S27 <- rename(read_tsv('mmseqs/S27DEC19/Viruses_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Viruses_S28 <- rename(read_tsv('mmseqs/S28DEC19/Viruses_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Viruses_S29 <- rename(read_tsv('mmseqs/S29DEC19/Viruses_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Viruses_df <- Viruses_S27 %>% full_join(Viruses_S28, by="lineage") %>% full_join(Viruses_S29, by="lineage") %>% replace(is.na(.),0)
Viruses_inc <- Viruses_df %>% select(-lineage) %>% as.matrix
rownames(Viruses_inc) <- Viruses_df$lineage
Viruses_spader <- SimilarityMult(Viruses_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Archaea_S27 <- rename(read_tsv('mmseqs/S27DEC19/Archaea_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Archaea_S28 <- rename(read_tsv('mmseqs/S28DEC19/Archaea_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Archaea_S29 <- rename(read_tsv('mmseqs/S29DEC19/Archaea_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Archaea_df <- Archaea_S27 %>% full_join(Archaea_S28, by="lineage") %>% full_join(Archaea_S29, by="lineage") %>% replace(is.na(.),0)
Archaea_inc <- Archaea_df %>% select(-lineage) %>% as.matrix
rownames(Archaea_inc) <- Archaea_df$lineage
Archaea_spader <- SimilarityMult(Archaea_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Bacteria_S27 <- rename(read_tsv('mmseqs/S27DEC19/Bacteria_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Bacteria_S28 <- rename(read_tsv('mmseqs/S28DEC19/Bacteria_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Bacteria_S29 <- rename(read_tsv('mmseqs/S29DEC19/Bacteria_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Bacteria_df <- Bacteria_S27 %>% full_join(Bacteria_S28, by="lineage") %>% full_join(Bacteria_S29, by="lineage") %>% replace(is.na(.),0)
Bacteria_inc <- Bacteria_df %>% select(-lineage) %>% as.matrix
rownames(Bacteria_inc) <- Bacteria_df$lineage
Bacteria_spader <- SimilarityMult(Bacteria_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Eukaryota_S27 <- rename(read_tsv('mmseqs/S27DEC19/Eukaryota_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Eukaryota_S28 <- rename(read_tsv('mmseqs/S28DEC19/Eukaryota_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Eukaryota_S29 <- rename(read_tsv('mmseqs/S29DEC19/Eukaryota_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Eukaryota_df <- Eukaryota_S27 %>% full_join(Eukaryota_S28, by="lineage") %>% full_join(Eukaryota_S29, by="lineage") %>% replace(is.na(.),0)
Eukaryota_inc <- Eukaryota_df %>% select(-lineage) %>% as.matrix
rownames(Eukaryota_inc) <- Eukaryota_df$lineage
Eukaryota_spader <- SimilarityMult(Eukaryota_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Protists_S27 <- rename(read_tsv('mmseqs/S27DEC19/Protists_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Protists_S28 <- rename(read_tsv('mmseqs/S28DEC19/Protists_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Protists_S29 <- rename(read_tsv('mmseqs/S29DEC19/Protists_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Protists_df <- Protists_S27 %>% full_join(Protists_S28, by="lineage") %>% full_join(Protists_S29, by="lineage") %>% replace(is.na(.),0)
Protists_inc <- Protists_df %>% select(-lineage) %>% as.matrix
rownames(Protists_inc) <- Protists_df$lineage
Protists_spader <- SimilarityMult(Protists_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Fungi_S27 <- rename(read_tsv('mmseqs/S27DEC19/Fungi_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Fungi_S28 <- rename(read_tsv('mmseqs/S28DEC19/Fungi_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Fungi_S29 <- rename(read_tsv('mmseqs/S29DEC19/Fungi_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Fungi_df <- Fungi_S27 %>% full_join(Fungi_S28, by="lineage") %>% full_join(Fungi_S29, by="lineage") %>% replace(is.na(.),0)
Fungi_inc <- Fungi_df %>% select(-lineage) %>% as.matrix
rownames(Fungi_inc) <- Fungi_df$lineage
Fungi_spader <- SimilarityMult(Fungi_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Viridiplantae_S27 <- rename(read_tsv('mmseqs/S27DEC19/Viridiplantae_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Viridiplantae_S28 <- rename(read_tsv('mmseqs/S28DEC19/Viridiplantae_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Viridiplantae_S29 <- rename(read_tsv('mmseqs/S29DEC19/Viridiplantae_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Viridiplantae_df <- Viridiplantae_S27 %>% full_join(Viridiplantae_S28, by="lineage") %>% full_join(Viridiplantae_S29, by="lineage") %>% replace(is.na(.),0)
Viridiplantae_inc <- Viridiplantae_df %>% select(-lineage) %>% as.matrix
rownames(Viridiplantae_inc) <- Viridiplantae_df$lineage
Viridiplantae_spader <- SimilarityMult(Viridiplantae_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Metazoa_S27 <- rename(read_tsv('mmseqs/S27DEC19/Metazoa_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Metazoa_S28 <- rename(read_tsv('mmseqs/S28DEC19/Metazoa_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Metazoa_S29 <- rename(read_tsv('mmseqs/S29DEC19/Metazoa_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Metazoa_df <- Metazoa_S27 %>% full_join(Metazoa_S28, by="lineage") %>% full_join(Metazoa_S29, by="lineage") %>% replace(is.na(.),0)
Metazoa_inc <- Metazoa_df %>% select(-lineage) %>% as.matrix
rownames(Metazoa_inc) <- Metazoa_df$lineage
Metazoa_spader <- SimilarityMult(Metazoa_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Arthropoda_S27 <- rename(read_tsv('mmseqs/S27DEC19/Arthropoda_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Arthropoda_S28 <- rename(read_tsv('mmseqs/S28DEC19/Arthropoda_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Arthropoda_S29 <- rename(read_tsv('mmseqs/S29DEC19/Arthropoda_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Arthropoda_df <- Arthropoda_S27 %>% full_join(Arthropoda_S28, by="lineage") %>% full_join(Arthropoda_S29, by="lineage") %>% replace(is.na(.),0)
Arthropoda_inc <- Arthropoda_df %>% select(-lineage) %>% as.matrix
rownames(Arthropoda_inc) <- Arthropoda_df$lineage
Arthropoda_spader <- SimilarityMult(Arthropoda_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Actinopteri_S27 <- rename(read_tsv('mmseqs/S27DEC19/Actinopteri_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Actinopteri_S28 <- rename(read_tsv('mmseqs/S28DEC19/Actinopteri_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Actinopteri_S29 <- rename(read_tsv('mmseqs/S29DEC19/Actinopteri_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Actinopteri_df <- Actinopteri_S27 %>% full_join(Actinopteri_S28, by="lineage") %>% full_join(Actinopteri_S29, by="lineage") %>% replace(is.na(.),0)
Actinopteri_inc <- Actinopteri_df %>% select(-lineage) %>% as.matrix
rownames(Actinopteri_inc) <- Actinopteri_df$lineage
Actinopteri_spader <- SimilarityMult(Actinopteri_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Mammalia_S27 <- rename(read_tsv('mmseqs/S27DEC19/Mammalia_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Mammalia_S28 <- rename(read_tsv('mmseqs/S28DEC19/Mammalia_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Mammalia_S29 <- rename(read_tsv('mmseqs/S29DEC19/Mammalia_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Mammalia_df <- Mammalia_S27 %>% full_join(Mammalia_S28, by="lineage") %>% full_join(Mammalia_S29, by="lineage") %>% replace(is.na(.),0)
Mammalia_inc <- Mammalia_df %>% select(-lineage) %>% as.matrix
rownames(Mammalia_inc) <- Mammalia_df$lineage
Mammalia_spader <- SimilarityMult(Mammalia_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

Aves_S27 <- rename(read_tsv('mmseqs/S27DEC19/Aves_genus_level.tsv') %>% rename_all(function(x) paste0("S27_", x)), lineage = S27_lineage)
Aves_S28 <- rename(read_tsv('mmseqs/S28DEC19/Aves_genus_level.tsv') %>% rename_all(function(x) paste0("S28_", x)), lineage = S28_lineage)
Aves_S29 <- rename(read_tsv('mmseqs/S29DEC19/Aves_genus_level.tsv') %>% rename_all(function(x) paste0("S29_", x)), lineage = S29_lineage)

Aves_df <- Aves_S27 %>% full_join(Aves_S28, by="lineage") %>% full_join(Aves_S29, by="lineage") %>% replace(is.na(.),0)
Aves_inc <- Aves_df %>% select(-lineage) %>% as.matrix
rownames(Aves_inc) <- Aves_df$lineage
Aves_spader <- SimilarityMult(Aves_inc,"incidence_raw", units = c(14,14,26), q = 0, nboot = 100, "relative")

#-------------------------------------------------------------------------------------------------------------------------------------------------------

similarity <- function(taxon, spader){
  sorenson <- spader$estimated_richness['C0N(q=0,Sorensen)','Estimate']
  s.ucl <- spader$estimated_richness['C0N(q=0,Sorensen)','95%.UCL']
  s.lcl <- spader$estimated_richness['C0N(q=0,Sorensen)','95%.LCL']
  jaccard <- spader$estimated_richness['U0N(q=0,Jaccard)','Estimate']
  j.ucl <- spader$estimated_richness['U0N(q=0,Jaccard)','95%.UCL']
  j.lcl <- spader$estimated_richness['U0N(q=0,Jaccard)','95%.LCL']
  s.df <- data.frame(taxon=taxon, Index='Sorenson', similarity=sorenson, upr=s.ucl, lwr=s.lcl)
  j.df <- data.frame(taxon=taxon, Index='Jaccard', similarity=jaccard, upr=j.ucl, lwr=j.lcl)
  return(rbind(s.df, j.df))
}

dat <- similarity('Viruses', Viruses_spader)
dat <- rbind(dat, similarity('Archaea', Archaea_spader))
dat <- rbind(dat, similarity('Bacteria', Bacteria_spader))
dat <- rbind(dat, similarity('Eukaryota', Eukaryota_spader))
dat <- rbind(dat, similarity('Protists', Protists_spader))
dat <- rbind(dat, similarity('Fungi', Fungi_spader))
dat <- rbind(dat, similarity('Viridiplantae', Viridiplantae_spader))
dat <- rbind(dat, similarity('Metazoa', Metazoa_spader))
dat <- rbind(dat, similarity('Arthropoda', Arthropoda_spader))
dat <- rbind(dat, similarity('Actinopteri', Actinopteri_spader))
dat <- rbind(dat, similarity('Mammalia', Mammalia_spader))
dat <- rbind(dat, similarity('Aves', Aves_spader))

plt <- ggplot(dat, aes(x=as.factor(taxon), y=similarity, fill=Index)) + theme_bw(base_size = 12) +
       scale_fill_manual(values=c("#006ddb","#b66dff")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
       geom_bar(position=position_dodge(), stat="identity", colour='black', alpha=0.5) +
       geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2, position=position_dodge(.9)) +
       labs(y='Community similarity') +
       theme(legend.position="bottom", axis.title.x=element_blank(),
				      axis.title.y=element_text(size=12, face='plain'))

ggsave('Similarity.png', plot=plt, path='mmseqs', width=169, height=120, units='mm', dpi=300)

