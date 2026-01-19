library(pairwiseAdonis)
FUNGI_SOIL_100  <- readRDS("C:/R/NGS/Maria_Phyloseq_Objects_new new new/Microcosms Experiment/ITS/B-DIVERSITY/Fungi_SOIL_100.RDS")
levels(as.factor(sample_data(FUNGI_SOIL_100)$Pesticide))
levels(as.factor(sample_data(FUNGI_SOIL_100)$Soil_Type))
levels(as.factor(sample_data(FUNGI_SOIL_100)$Time_Rate))
levels(as.factor(sample_data(FUNGI_SOIL_100)$Dose))
##Relative abundance transformation 100%
FUNGI_SOIL_Pella <- prune_samples(!(sample_names(FUNGI_SOIL_Pella) %in% c("MK_001")), FUNGI_SOIL_Pella)

FUNGI_SOIL_Pella <- prune_taxa(taxa_sums(FUNGI_SOIL_Pella)>0,FUNGI_SOIL_Pella)


FUNGI_SOIL_Pella_100 <- transform_sample_counts(FUNGI_SOIL_Pella, function(OTU) 100*OTU/sum(OTU))
saveRDS(FUNGI_SOIL_Pella_100, file = "FUNGI_SOIL_Pella_100.RDS")

##Koutsoura ANALYSIS##
FUNGI_SOIL_Koutsoura_100 <- transform_sample_counts(FUNGI_SOIL_Koutsoura, function(OTU) 100*OTU/sum(OTU))
saveRDS(FUNGI_SOIL_Koutsoura_100, file = "FUNGI_SOIL_Koutsoura_100.RDS")

##Farsala ANALYSIS##
FUNGI_SOIL_Farsala <- prune_samples(setdiff(sample_names(FUNGI_SOIL_Farsala), c("MK_104", "MK_084")), FUNGI_SOIL_Farsala)


FUNGI_SOIL_Farsala <- prune_taxa(taxa_sums(FUNGI_SOIL_Farsala)>0,FUNGI_SOIL_Farsala)

FUNGI_SOIL_Farsala_100 <- transform_sample_counts(FUNGI_SOIL_Farsala, function(OTU) 100*OTU/sum(OTU))

saveRDS(FUNGI_SOIL_Farsala_100, file = "FUNGI_SOIL_Farsala_100.RDS")



##PERMANOVA FOR EACH SOIL TYPE##
mypermanova_FUNGI_SOIL_Pella <- adonis2(FUNGI_SOIL_Pella_100@otu_table ~ Dose + Time_Rate, method = "bray", data = data.frame(FUNGI_SOIL_Pella_100@sam_data), by = "terms")
write.table(data.frame(mypermanova_FUNGI_SOIL_Pella), file="mypermanova_FUNGI_SOIL_Pella.txt", quote = F,col.names = NA, sep="\t")

mypermanova_FUNGI_SOIL_Koutsoura <- adonis2(FUNGI_SOIL_Koutsoura_100@otu_table ~Dose + Time_Rate, method = "bray", data = data.frame(FUNGI_SOIL_Koutsoura_100@sam_data), by = "terms")
write.table(data.frame(mypermanova_FUNGI_SOIL_Koutsoura), file="mypermanova_FUNGI_SOIL_Koutsoura.txt", quote = F,col.names = NA, sep="\t")

mypermanova_FUNGI_SOIL_Farsala <- adonis2(FUNGI_SOIL_Farsala_100@otu_table ~Dose + Time_Rate, method = "bray", data = data.frame(FUNGI_SOIL_Farsala_100@sam_data), by = "terms")
write.table(data.frame(mypermanova_FUNGI_SOIL_Farsala), file="mypermanova_FUNGI_SOIL_Farsala.txt", quote = F,col.names = NA, sep="\t")


#-#-pairwise permanova (sto RA% object)-#-#
library(pairwiseAdonis)

mycmpfactor_FUNGI_100 <- interaction(data.frame(FUNGI_SOIL_100@sam_data)$Soil_Type, data.frame(FUNGI_SOIL_100@sam_data)$Dose) 

mympairwiseperm_FUNGI_100 <- pairwise.adonis(FUNGI_SOIL_100@otu_table, mycmpfactor_FUNGI_100, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")
write.table(mympairwiseperm_FUNGI_100,
            file = "mympairwiseperm_FUNGI_100_soil_types.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

####PELLA###
mycmpfactor_FUNGI_SOIL_Pella<- interaction(data.frame(FUNGI_SOIL_Pella_100@sam_data)$Dose) #data.frame(FUNGI_SOIL_Pella_100@sam_data)$Dose) 

mympairwiseperm_FUNGI_SOIL_Pella<- pairwise.adonis(FUNGI_SOIL_Pella_100@otu_table, mycmpfactor_FUNGI_SOIL_Pella, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")##dose vs dose

write.table(mympairwiseperm_FUNGI_SOIL_Pella,
            file = "mympairwiseperm_FUNGI_SOIL_Pella_dosevsdose.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

##FARSALA###
mycmpfactor_FUNGI_SOIL_Farsala<- interaction(data.frame(FUNGI_SOIL_Farsala_100@sam_data)$Dose)
#data.frame(FUNGI_SOIL_Farsala_100@sam_data)$Dose) 

mympairwiseperm_FUNGI_SOIL_Farsala <- pairwise.adonis(FUNGI_SOIL_Farsala_100@otu_table, mycmpfactor_FUNGI_SOIL_Farsala, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")##dose vs dose
write.table(mympairwiseperm_FUNGI_SOIL_Farsala,
            file = "mympairwiseperm_FUNGI_SOIL_Farsala_dosevsdose.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

##KOUTSOURA##
mycmpfactor_FUNGI_SOIL_Koutsoura<- interaction(data.frame(FUNGI_SOIL_Koutsoura_100@sam_data)$Dose)
#data.frame(FUNGI_SOIL_Koutsoura_100@sam_data)$Dose) 

mympairwiseperm_FUNGI_SOIL_Koutsoura <- pairwise.adonis(FUNGI_SOIL_Koutsoura_100@otu_table, mycmpfactor_FUNGI_SOIL_Koutsoura, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")##dose vs dose
write.table(mympairwiseperm_FUNGI_SOIL_Koutsoura,
            file = "mympairwiseperm_FUNGI_SOIL_Koutsoura_dosevsdose.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

##DONE
