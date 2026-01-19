FUNGI_SOIL  <- readRDS("C:/R/NGS/Maria_Phyloseq_Objects_new new new/Microcosms Experiment/ITS/B-DIVERSITY/Fungi_SSDS_MK.RDS")
FUNGI_SOIL.Cl <- prune_taxa(taxa_sums(FUNGI_SOIL)>0,FUNGI_SOIL)

##Remove Kingdom NA & Phylum NA##
table(tax_table(FUNGI_SOIL)[, "Kingdom"], exclude = NULL)
table(tax_table(FUNGI_SOIL)[, "Phylum"], exclude = NULL)

FUNGI_SOIL.Cl.2 <- subset_taxa(FUNGI_SOIL.Cl, !is.na(Kingdom) & !is.na(Phylum))

table(tax_table(FUNGI_SOIL.Cl.2)[, "Kingdom"], exclude = NULL)

FUNGI_SOIL.Cl.2 

table(tax_table(FUNGI_SOIL.Cl.2)[, "Phylum"], exclude = NULL)

FUNGI_SOIL.Cl.2 <- prune_taxa(taxa_sums(FUNGI_SOIL.Cl.2)>0,FUNGI_SOIL.Cl.2)

##Remove Mitochondria##
FUNGI_SOIL.Cl.3 <- subset_taxa(FUNGI_SOIL.Cl.2, !Family %in% c("Mitochondria"))

FUNGI_SOIL.Cl.3 <- prune_taxa(taxa_sums(FUNGI_SOIL.Cl.3)>0,FUNGI_SOIL.Cl.3)

##Remove Chloroplast##
FUNGI_SOIL.Cl.4 <- subset_taxa(FUNGI_SOIL.Cl.3, !Order %in% c("Chloroplast"))

FUNGI_SOIL.Cl.4 <- prune_taxa(taxa_sums(FUNGI_SOIL.Cl.4)>0,FUNGI_SOIL.Cl.4)

FUNGI_SOIL.Cl.4

##Annotation##
FUNGI_SOIL.Cl.2.An <- FUNGI_SOIL.Cl.2

for(i in 1:nrow(tax_table(FUNGI_SOIL.Cl.2.An))){
  for(j in 2:ncol(tax_table(FUNGI_SOIL.Cl.2.An))){
    if(is.na(tax_table(FUNGI_SOIL.Cl.2.An)[i,j])){
      tax_table(FUNGI_SOIL.Cl.2.An)[i,j] <- tax_table(FUNGI_SOIL.Cl.2.An)[i,j-1]
    }
  }
}   
#-#-#-#-#-#-Seperate projects per Soil Type-#-#-#-#-#-#-#
##Pella##
FUNGI_SOIL_Pella <- subset_samples(FUNGI_SOIL.Cl.2.An, Soil_Type %in% (c("Pella")))

FUNGI_SOIL_Pella <- prune_taxa(taxa_sums(FUNGI_SOIL_Pella)>0,FUNGI_SOIL_Pella)

View(data.frame(sample_data(FUNGI_SOIL_Pella)))


FUNGI_SOIL_Pella
saveRDS(FUNGI_SOIL_Pella, file = "FUNGI_SOIL_Pella.RDS")
##Koutrsoura##
FUNGI_SOIL_Koutsoura <- subset_samples(FUNGI_SOIL.Cl.2.An, Soil_Type %in% (c("Koutsoura")))

FUNGI_SOIL_Koutsoura <- prune_taxa(taxa_sums(FUNGI_SOIL_Koutsoura)>0,FUNGI_SOIL_Koutsoura)

View(data.frame(sample_data(FUNGI_SOIL_Koutsoura)))


FUNGI_SOIL_Koutsoura
saveRDS(FUNGI_SOIL_Koutsoura, file = "FUNGI_SOIL_Koutsoura.RDS")
##Farsala##
FUNGI_SOIL_Farsala <- subset_samples(FUNGI_SOIL.Cl.2.An, Soil_Type %in% (c("Farsala")))

FUNGI_SOIL_Farsala <- prune_taxa(taxa_sums(FUNGI_SOIL_Farsala )>0,FUNGI_SOIL_Farsala )

View(data.frame(sample_data(FUNGI_SOIL_Farsala )))


FUNGI_SOIL_Pella
saveRDS(FUNGI_SOIL_Pella, file = "FUNGI_SOIL_Pella.RDS")

#-#-#-#-start the analysis per soil seperately as the main factor is the dose effect-#-#-#-#
##NMDS#

levels(as.factor(sample_data(FUNGI_SOIL_100)$Pesticide))
levels(as.factor(sample_data(FUNGI_SOIL_100)$Soil_Type))
levels(as.factor(sample_data(FUNGI_SOIL_100)$Time_Rate))
levels(as.factor(sample_data(FUNGI_SOIL_100)$Dose))
##Relative abundance transformation 100%
FUNGI_SOIL_100 <- transform_sample_counts(FUNGI_SOIL, function(OTU) 100*OTU/sum(OTU))
View(data.frame(sample_data(FUNGI_SOIL_100)))
saveRDS(FUNGI_SOIL_100, file = "FUNGI_SOIL_100.RDS")

##Pella ANALYSIS##

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

###NMDS EACH SOIL SEPERATELY###


####PELLA#####


ord.nmds.bray_ALL <- ordinate(FUNGI_SOIL_Pella_100, method = "NMDS", distance = "bray")


sample_data(FUNGI_SOIL_Pella_100)$Dose <- factor(sample_data(FUNGI_SOIL_Pella_100)$Dose)


plot_pella <- plot_ordination(FUNGI_SOIL_Pella_100, ord.nmds.bray_ALL, 
                              color = "Dose", 
                              title = paste("NMDS (stress ", round(ord.nmds.bray_ALL$stress, 2), ")", sep = "")) + 
  geom_point(aes(shape = Time_Rate), size = 3) +
  stat_ellipse(aes(group = Dose, linetype = factor(Dose)), 
               type = "t", linewidth = 0.5) +
  scale_color_manual(values = c("#3288BD","#D3D3D3", "#5CC465" , "#D612E0", "#FDAE61","#D53E4F")) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

print(plot_pella)


# Print the plot
ggsave("FUNGI_plot_pella.pdf", plot = plot_pella, device = "pdf", width = 6, height = 6, dpi = 300)

####Koutsoura####

sample_data(FUNGI_SOIL_Koutsoura_100)$Dose <- factor(sample_data(FUNGI_SOIL_Koutsoura_100)$Dose)
ord.nmds.bray_ALL <- ordinate(FUNGI_SOIL_Koutsoura_100, method = "NMDS", distance = "bray")

# Create NMDS plot
plot_Koutsoura<- plot_ordination(FUNGI_SOIL_Koutsoura_100, ord.nmds.bray_ALL, 
                                 color = "Dose", 
                                 title = paste("NMDS (stress ", round(ord.nmds.bray_ALL$stress, 2), ")", sep = "")) + 
  geom_point(aes(shape = Time_Rate), size = 3) +
  stat_ellipse(aes(group = Dose, linetype = factor(Dose)), 
               type = "t", linewidth = 0.5) +  # Ellipses by Soil_Type with line types
  scale_color_manual(values = c("#3288BD","#D3D3D3", "#5CC465" , "#D612E0", "#FDAE61","#D53E4F")) +
  #facet_wrap(~Soil_Type) +  # Create separate panels for each soil type (optional)
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

print(plot_Koutsoura)
# Print the plot
ggsave("FUNGI_plot_Koutoura.pdf", plot = plot_Koutsoura, device = "pdf", width = 6, height = 6, dpi = 300)

####farsala####

sample_data(FUNGI_SOIL_Farsala_100)$Dose <- factor(sample_data(FUNGI_SOIL_Farsala_100)$Dose)
ord.nmds.bray_ALL <- ordinate(FUNGI_SOIL_Farsala_100, method = "NMDS", distance = "bray")

# Create NMDS plot
plot_Farsala<- plot_ordination(FUNGI_SOIL_Farsala_100, ord.nmds.bray_ALL, 
                               color = "Dose", 
                               title = paste("NMDS (stress ", round(ord.nmds.bray_ALL$stress, 2), ")", sep = "")) + 
  geom_point(aes(shape = Time_Rate), size = 3) +
  stat_ellipse(aes(group = Dose, linetype = factor(Dose)), 
               type = "t", linewidth = 0.5) +  # Ellipses by Soil_Type with line types
  scale_color_manual(values = c("#3288BD","#D3D3D3", "#5CC465" , "#D612E0", "#FDAE61","#D53E4F")) +
  #facet_wrap(~Soil_Type) +  # Create separate panels for each soil type (optional)
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank()
  )

print(plot_Farsala)
# Print the plot
ggsave("FUNGI_plot_Farsala.pdf", plot = plot_Farsala, device = "pdf", width = 6, height = 6, dpi = 300)



