############################################################
## ASV Modelling per Soil and Pesticide
############################################################

## Libraries
library(phyloseq)
library(microbiome)
library(dplyr)
library(tibble)
library(drc)
library(qpcR)
library(writexl)
library(ssdtools)
library(ggplot2)
library(gridExtra)
library(future)
library(doFuture)

############################################################
## Load phyloseq object
############################################################
ps <- readRDS("Fungi_F.RDS")

############################################################
## Sample data preparation
############################################################
write.table(
  data.frame(sample_data(ps)),
  file = "Sample_Data.txt",
  quote = FALSE,
  col.names = NA,
  sep = "\t"
)
 #need to edit the phyloseq template according to the factors you need to analysis (pesticide or soil type )
sample_data(ps) <- read.table(
  "Sample_Data.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

############################################################
## CLR transformation & taxonomic aggregation
############################################################
ps_tr <- transform(ps, "clr")
ps_tr <- tax_glom(ps_tr, taxrank = "Species")
ps_tr <- transform(ps_tr, "clr")

############################################################
## Clean taxonomy table
############################################################
remove_na_columns <- function(ps){
  tt <- as.data.frame(tax_table(ps))
  tt <- tt[, colSums(is.na(tt)) < nrow(tt)]
  tax_table(ps) <- tax_table(as.matrix(tt))
  ps
}

ps_tr <- remove_na_columns(ps_tr)

tax_df <- as.data.frame(tax_table(ps_tr))
tax_df$ASV_ID <- rownames(tax_df)
write_xlsx(tax_df, "taxonomy_table_ITS.xlsx")

############################################################
## Dose–response modelling
############################################################
# Get unique soil types and pesticide treatments (excluding NTC)
mysoils <- unique(sample_data(ps_tr)$Soil_Type)[which(unique(sample_data(ps)$Soil_Type) != "NTC")]
mypesticides <- unique(sample_data(ps_tr)$pesticide)[which(unique(sample_data(ps)$pesticide) != "NTC")]

# Open PDF for saving plots
pdf(file = "modelling_plots_tr_all_arhea_ed20.pdf", onefile = TRUE, width = 14, height = 18)

# Set plotting parameters
par(bty = "n", mar = c(4, 5, 4, 4), mfrow = c(4, 3))

# Initialize final results table
tab_mod_fin <- data.frame(
  Soil_Pest = character(),
  MyFeat = character(),
  ModRank = numeric(),
  mod = character(),
  ED = character(),
  Rsq = numeric(),
  IC = numeric(),
  Estimate = numeric(),
  `Std. Error` = numeric(),
  Lower = numeric(),
  Upper = numeric(),
  MaxNom = numeric(),
  check.names = FALSE
)

# Iterate over each soil type
for (mysoil in mysoils) {
  # Filter the phyloseq object for the current soil type
  ps_soil <- subset_samples(physeq = ps_tr, Soil_Type == mysoil)
  
  # Iterate over each pesticide
  for (mypesticide in mypesticides) {
    # Subset further for the pesticide type
    ps_sel <- subset_samples(physeq = ps_soil, pesticide == mypesticide | pesticide == "Metmethyl")
    ps_sel <- prune_taxa(taxa = taxa_sums(ps_sel) > 0, x = ps_sel)
    
    # Create model matrix
    mymodmat <- data.frame(
      pest = sample_data(ps_sel)$Dose,
      data.frame(otu_table(ps_sel))[row.names(sample_data(ps_sel)), ]
    )
    View(mymodmat)
    # Load required libraries
    library(drc)
    
    # Iterate over each feature (ASV)
    for (myfeat in colnames(mymodmat)[grep("ASV", colnames(mymodmat))]) {
      tryCatch({
        # Fit initial model
        dat_a.mLL.4 <- drm(mymodmat[, myfeat] ~ pest, data = mymodmat, fct = BC.5(names = c("Slope", "Lower Limit", "Upper Limit", "ED50", "Emax")))
        
        # Model selection based on AIC
        mod_comp <- mselect(dat_a.mLL.4, c(getMeanFunctions(3), getMeanFunctions(4), getMeanFunctions(5)), icfct = AIC)
        mod_comp <- mod_comp[complete.cases(mod_comp), ]
        sel_mod <- row.names(mod_comp)[1:5]
        
        # Fit top three models
        mod1 <- drm(mymodmat[, myfeat] ~ pest, data = mymodmat, fct = get(sel_mod[1])())
        mod2 <- drm(mymodmat[, myfeat] ~ pest, data = mymodmat, fct = get(sel_mod[2])())
        mod3 <- drm(mymodmat[, myfeat] ~ pest, data = mymodmat, fct = get(sel_mod[3])())
        mod4 <- drm(mymodmat[, myfeat] ~ pest, data = mymodmat, fct = get(sel_mod[4])())
        mod5 <- drm(mymodmat[, myfeat] ~ pest, data = mymodmat, fct = get(sel_mod[5])())
        # Plot results
        plot(mod1, type = "all", main = paste(mysoil, "\n", myfeat), 
             xlab = bquote(.(mypesticide) * " mg kg"^-1), 
             ylab = "RA (%)", cex.main = 2, lwd = 2, 
             cex = 1.2, cex.axis = 1.2, cex.lab = 1.2, 
             ylim = c(0, 1.2 * max(c(mymodmat[, myfeat], PR(mod1, seq(0, max(mymodmat$pest), 0.1))))))
        plot(mod2, add = TRUE, col = "red", lty = 2, lwd = 2)
        plot(mod3, add = TRUE, col = "blue", lty = 3, lwd = 2)
        plot(mod4, add = TRUE, col = "yellow", lty = 3, lwd = 2)
        plot(mod5, add = TRUE, col = "green", lty = 3, lwd = 2)
        legend("topleft", legend = c(getMeanFunctions(fname = sel_mod[1])[[1]]$text,
                                     getMeanFunctions(fname = sel_mod[2])[[1]]$text,
                                     getMeanFunctions(fname = sel_mod[3])[[1]]$text,
                                     getMeanFunctions(fname = sel_mod[4])[[1]]$text,
                                     getMeanFunctions(fname = sel_mod[5])[[1]]$text), 
               cex = 0.8, lty = 1:3, bty = "n", lwd = 2, col = c("black", "red", "blue","yellow","green"))
        
        # Collect R-squared values
        Rsq1 <- qpcR::Rsq(mod1)
        Rsq2 <- qpcR::Rsq(mod2)
        Rsq3 <- qpcR::Rsq(mod3)
        Rsq4 <- qpcR::Rsq(mod4)
        Rsq5 <- qpcR::Rsq(mod5)
        # Compile model information into a table
        for (model_num in 1:5) {
          tryCatch({
            Soil_Pest <- rep(paste(mysoil, mypesticide, sep = "_"), 3)
            MyFeat <- rep(myfeat, 3)
            ModRank <- rep(model_num, 3)
            mod <- rep(sel_mod[model_num], 3)
            ED <- c("ED20", "ED50", "ED90")
            Rsq <- rep(get(paste("Rsq", model_num, sep = "")), 3)
            IC <- rep(round(mod_comp[model_num, 2], 3), 3)
            EDxmod <- ED(get(paste("mod", model_num, sep = "")), respLev = c(20, 50, 90), interval = "delta")
            row.names(EDxmod) <- paste(mod, c("ED20", "ED50", "ED90"), sep = "_")
            tab_mod <- as.data.frame(cbind(Soil_Pest, MyFeat, ModRank, mod, ED, Rsq, IC, round(EDxmod, 3)))
            
            tab_mod_fin <- rbind(tab_mod_fin, tab_mod)
          }, error = function(e) {})
        }
      }, error = function(e) { print(paste("Model for", myfeat, "is not possible")) })
    }
  }
}

# Close the PDF file
dev.off()
write.table(tab_mod_fin,file=paste("report_tr_all_4_ITS.txt",sep=""),quote=F,col.names=NA,sep="\t")


############################################################
## SSD analysis
############################################################
mySoilPests <- unique(sample_data(myreport)$Soil_Pest)
myEDs <- unique(sample_data(myreport)$ED) 


# Initialize results storage
myhc_table <- data.frame(Soil_Pest = character(), ED = character(), hc5 = numeric(), hc20 = numeric())
myplot_list <- list()

for(mySoilPest in mySoilPests){
  for(myED in myEDs){
    # Filter data for the current Soil_Pest and ED combination
    myreport_fin_1 <- myreport[myreport$Soil_Pest == mySoilPest & 
                                 myreport$ModRank == 1 & 
                                 myreport$ED == myED, ]
    
    # Debugging
    print(paste("Processing:", mySoilPest, "-", myED))
    print(head(myreport_fin_1))
    
    # Clean and ensure numeric conversion
    myreport_fin_1$Estimate <- as.numeric(trimws(myreport_fin_1$Estimate))
    myreport_fin_1 <- myreport_fin_1[!is.na(myreport_fin_1$Estimate), ]
    
    # Remove rows with invalid MyFeat values
    myreport_fin_1 <- myreport_fin_1[!is.na(myreport_fin_1$MyFeat), ]
    
    # Filter Estimate values within desired range
    myreport_fin <- myreport_fin_1[myreport_fin_1$Estimate <= 4000 & myreport_fin_1$Estimate > 0, ]
    
    # Check if myreport_fin is empty before proceeding
    if (nrow(myreport_fin) == 0) next
    
    # Set row names to MyFeat
    row.names(myreport_fin) <- myreport_fin$MyFeat
    
    # Fit distributions and create CDF plots
    fits <- ssd_fit_dists(myreport_fin, left = "Estimate")
    myplot_list[[paste(mySoilPest, myED)]] <- ssd_plot_cdf(fits, label = "MyFeat")
    
    # Calculate hazard concentrations
    hc <- ssd_hc(fits, percent = c(5, 20))
    myhc_tableNew <- data.frame(Soil_Pest = mySoilPest, ED = myED, hc5 = hc$est[1], hc20 = hc$est[2])
    
    # Append results to the summary table
    myhc_table <- rbind(myhc_table, myhc_tableNew)
  }
}


write.table(myhc_table, file = "HC5_values_ITS_P.txt", sep = "\t", quote = FALSE, row.names = FALSE)

library(gridExtra)
ggsave(paste("SSD_plots_ITS_P.pdf",sep=""), marrangeGrob(grobs = myplot_list, layout_matrix = matrix(1:12, nrow = 4, ncol=3, byrow=TRUE)), device = "pdf", width = 14, height = 18)
dev.off()
str(myreport_fin)
colnames(myreport_fin) <- c("x","ID", "Microbe","ModRank","Mod","ED","Rsq", "IC", "Conc", "Std. Error", "Lower", "Upper","Phylum", "Class",	"Order",	"Family",	"Genus",	"Species")
ompleterecords <- na.omit(myreport_fin)

CNW_data <- ompleterecords

# Step 1: Fit SSD Models to Data
CNW_fits <- ssd_fit_dists(CNW_data, dists = ssd_dists_all())
print(CNW_fits)

# Step 2: Tidy Output for Easy Inspection
# If tidy() fails, you can summarize manually:
if ("tidy" %in% ls("package:ssdtools")) {
  tidy_fits <- tidy(CNW_fits) 
  print(tidy_fits)
} else {
  warning("tidy() is unavailable. Use summary or manual inspection.")
}

# Step 3: Assess Goodness of Fit
gof_results <- ssd_gof(CNW_fits)
print(gof_results)
write.table(gof_results, file = "gof_results_ITS_P.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Step 4: Register Multisession for Parallel Computing
doFuture::registerDoFuture()
future::plan(future::multisession)

# Step 5: Generate Predictions with Confidence Intervals
CNW_preds <- predict(CNW_fits, ci = TRUE)
print(CNW_preds)

# Step 6: Optional Data for Manual Labeling (Uncomment if required)
# x <- CNW_ED50_df$Conc  # Data for manual x-coordinates
# y <- CNW_preds$percent # Predicted percentiles for manual labeling

# Step 7: Create SSD Plot and Save to PDF
pdf("PELLA_ITS.pdf", width = 10, height = 6) # Customize dimensions per sample
ssd_plot(
  CNW_data, CNW_preds, 
  color = "Phylum", 
  label = "Genus", 
  xlab = "Concentration (mg/kg soil)", 
  ribbon = TRUE
) +
  expand_limits(x = 5000) + # Extend x-axis for better visualization
  ggtitle("Fungi Sensitivity for Metsulfuron Methyl") +
  # Uncomment next line for manual text labels
  # geom_text(aes(x = x, y = y, label = Microbe), color = "black", hjust = -0.2) +
  theme_bw()
dev.off()


############################################################
## End of script
############################################################
