############################################
#### 1. Package installation (if needed) ####
############################################

# Install required Bioconductor packages (run only once)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2")


################################
#### 2. Package loading ####
################################

# Load DADA2 and check version
library(dada2)
packageVersion("dada2")

# Set number of threads for multithreading
mymultthread <- 56


##########################################
#### 3. Working directory and input paths ####
##########################################

# Set working directory for all analysis outputs
setwd("/.............../")

# Define paths containing demultiplexed FASTQ files
path1 <- "~/ex.Fungi"
path2 <- "~/ex.Fungi"

# List FASTQ files in input directories
list.files(c(path1, path2))


##############################################
#### 4. File listing and sample name parsing ####
##############################################

# Identify forward and reverse read files
fnFs <- sort(list.files(c(path1, path2),
                        pattern = "_R1.fastq",
                        full.names = TRUE))
fnRs <- sort(list.files(c(path1, path2),
                        pattern = "_R2.fastq",
                        full.names = TRUE))

# Extract sample names from filenames
sample.names <- gsub("_lib_.+", "", basename(fnFs))


##############################################
#### 5. Raw read quality assessment ####
##############################################

# Plot per-base quality profiles for all samples
pdf("Initial Quality Files.pdf", onefile = TRUE)
for (i in seq_along(fnFs)) {
  plot1 <- plotQualityProfile(c(fnFs[i], fnRs[i]))
  print(plot1)
}
dev.off()


##############################################################
#### 6. Quality filtering and trimming of sequencing reads ####
##############################################################

# Define output paths for filtered reads
filtFs <- file.path("filtered2",
                    paste(sample.names, "_F_filt.fastq.gz", sep = ""))
filtRs <- file.path("filtered2",
                    paste(sample.names, "_R_filt.fastq.gz", sep = ""))

# Filter and trim reads
out <- filterAndTrim(fnFs, filtFs,
                     fnRs, filtRs,
                     trimLeft = 11,
                     trimRight = 50,
                     maxN = 0,
                     maxEE = c(2, 2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE,
                     matchIDs = TRUE)

# Inspect filtering summary
head(out)


############################################
#### 7. Error rate learning ####
############################################

# Learn error rates for forward and reverse reads
errF <- learnErrors(filtFs, multithread = mymultthread)
errR <- learnErrors(filtRs, multithread = mymultthread)

# Visualize error models
pdf("Estimated Error Rates.pdf", onefile = TRUE)
print(plotErrors(errF, nominalQ = TRUE))
print(plotErrors(errR, nominalQ = TRUE))
dev.off()


############################################
#### 8. Dereplication and ASV inference ####
############################################

# Dereplicate reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Assign sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer ASVs
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


############################################
#### 9. Paired-end read merging ####
############################################

# Merge forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs,
                      dadaRs, derepRs,
                      verbose = TRUE)


##################################################
#### 10. Sequence table construction and chimera removal ####
##################################################

# Construct ASV table
seqtab <- makeSequenceTable(mergers)

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method = "consensus",
                                    multithread = mymultthread,
                                    verbose = TRUE)

# Inspect table dimensions and retained reads
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)


############################################
#### 11. Read tracking through the pipeline ####
############################################

# Function to count reads
getN <- function(x) sum(getUniques(x))

# Track read counts
track <- cbind(out,
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered",
                     "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample.names

# Save read tracking table
write.table(track,
            file = "readQCFungi.txt",
            col.names = NA,
            sep = "\t",
            quote = FALSE)


############################################
#### 12. Taxonomic assignment ####
############################################

# Assign taxonomy using UNITE database
taxa <- assignTaxonomy(seqtab.nochim,
                       "/mnt/4abcb8af-e2eb-47e2-86d6-eed71f6d5304/00_home_guests_do_not_move/fotisbs/wshp_foldrs/dbs/UNITE_general_release_s_all_04.02.2020/sh_general_release_dynamic_s_all_04.02.2020.fasta",
                       minBoot = 80,
                       multithread = TRUE,
                       tryRC = TRUE)

# Save taxonomy table
write.table(taxa,
            file = "taxa.txt",
            col.names = NA,
            sep = "\t",
            quote = FALSE)


############################################
#### 13. Phyloseq object construction ####
############################################

library(phyloseq)
packageVersion("phyloseq")

# Load sample metadata
samdf <- read.table("SamdfFungi.txt",
                    header = TRUE,
                    row.names = 1,
                    sep = "\t",
                    check.names = FALSE)

# Construct phyloseq object
Fungi_Maria_2024 <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa)
)


##################################################
#### 14. ASV renaming and sequence storage ####
##################################################

library(Biostrings)
library(stringr)

# Store DNA sequences in phyloseq object
sequences <- DNAStringSet(taxa_names(Fungi_Maria_2024))
names(sequences) <- taxa_names(Fungi_Maria_2024)

ps <- merge_phyloseq(Fungi_Maria_2024, sequences)

# Rename ASVs
taxa_names(ps) <- paste0("ASV",
                         str_pad(seq_along(taxa_names(ps)),
                                 width = 4,
                                 pad = "0"))

Fungi_Maria_2024_Final <- ps


############################################
#### 15. Taxonomic summaries ####
############################################

table(tax_table(Fungi_Maria_2024_Final)[, "Kingdom"], exclude = NULL)
table(tax_table(Fungi_Maria_2024_Final)[, "Phylum"],  exclude = NULL)
table(tax_table(Fungi_Maria_2024_Final)[, "Class"],   exclude = NULL)
table(tax_table(Fungi_Maria_2024_Final)[, "Order"],   exclude = NULL)
table(tax_table(Fungi_Maria_2024_Final)[, "Family"],  exclude = NULL)
table(tax_table(Fungi_Maria_2024_Final)[, "Genus"],   exclude = NULL)


############################################
#### 16. Save phyloseq object ####
############################################

saveRDS(Fungi_Maria_2024_Final,
        file = "Fungi_Maria_2024_Final.RDS")


############################################
#### 17. Taxonomic filtering and final cleaning ####
############################################

# Reload phyloseq object
Fungi_Maria_2024_Final <- readRDS(".........../Fungi_Maria_2024_Final.RDS")

# Remove taxa without Kingdom assignment
Fungi_Maria_2024_Final.Cl <- subset_taxa(Fungi_Maria_2024_Final,
                                         !is.na(Kingdom))

# Retain Ascomycota and Basidiomycota only
Fungi_Maria_2024_Final.Cl2 <- subset_taxa(
  Fungi_Maria_2024_Final.Cl,
  !is.na(Phylum) &
    Phylum %in% c("p__Ascomycota", "p__Basidiomycota")
)

Fungi_Maria_2024_Final.Cl2 <- prune_taxa(
  taxa_sums(Fungi_Maria_2024_Final.Cl2) > 0,
  Fungi_Maria_2024_Final.Cl2
)


##################################################
#### 18. Propagate taxonomy to lower ranks ####
##################################################

# Assign missing taxonomic ranks from previous level
Fungi_Maria_2024_Final.Cl.2.An <- Fungi_Maria_2024_Final.Cl2

for (i in seq_len(nrow(tax_table(Fungi_Maria_2024_Final.Cl.2.An)))) {
  for (j in 2:ncol(tax_table(Fungi_Maria_2024_Final.Cl.2.An))) {
    if (is.na(tax_table(Fungi_Maria_2024_Final.Cl.2.An)[i, j])) {
      tax_table(Fungi_Maria_2024_Final.Cl.2.An)[i, j] <-
        tax_table(Fungi_Maria_2024_Final.Cl.2.An)[i, j - 1]
    }
  }
}

# Final phyloseq object
Fungi_Maria_2024 <- Fungi_Maria_2024_Final.Cl.2.An

Fungi_Maria_2024

# Continue with each file script separately to reproduce the analysis figures