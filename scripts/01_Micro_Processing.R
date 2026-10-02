#Changed into script 12th December
#Original Jan 2024
#Lauren Mee Git: DrLMee

# Required Input #
#assumed that trimmed 16S read data has been added to a folder
#within the current workspace termed input/microbiome/trimmed/
#training sets from SILVA SSU databases are put in input/microbiome/taxaDB/

##### set up #####
#set seed
set.seed(1049)

#libraries
libs <- c("dada2", "phyloseq", #microbiome processing
          "Biostrings", #save nucleotide sequences
          "ggpubr", "tidyverse") #general plots etc

for (pkg in libs) {
  library(pkg, character.only = T)
}

#output directory architecture
dirs <- c("processed/", "processed/microbiome/", 
          "output/", "output/microbiome/", "output/microbiome/fastqQC/",
          "output/microbiome/community/", "output/figures/")

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d)
  }
}

#visualisations
#set colour palette
vanPal <- c("#5d83df", "#df5d83", "#40B0A6",
            "#00ffff",   "black", "#cdb126",
            "#8dbab1", "#F585F5", "#65348C")

#ggplot themes
theme_set(theme_bw(base_size = 13))
theme_update(
  strip.background = element_rect(fill = "lightgrey", colour = "white"),
  strip.text = element_text(colour = "black", face = "bold")
)

##### Process Reads Using dada2 ####
#list files, splitting into forward and reverse per sample
fqs <- sort(list.files("input/microbiome/trimmed/", full.names = T))
#store for and rev fastq files
fqF <- fqs[grepl("R1", fqs)]
#pull out sample IDs
microS <- sapply(strsplit(basename(fqF), "_"), '[', 1)

#plot sample quality, 9 at a time
batches <- ceiling(length(fqF) / 9)
for (i in 1:batches) {
  idx <- ((i-1) * 9 + 1):min(i * 9, length(fqF))
  p <- plotQualityProfile(fqF[idx]) 
  ggsave(paste0("output/microbiome/fastqQC/FastqQC_", i, ".pdf"))
}

#prepare forward read filepaths
filtFs <- file.path(paste0("processed/microbiome/filteredFastq/"), 
                    paste0(microS, "_F_filt.fastq.gz"))
#add sample names
names(filtFs) <- microS

#using standard filtering parameters
out <- suppressMessages(filterAndTrim(fqF, filtFs,
                                      maxN = 0, #dada2 will not except Ns
                                      maxEE = 2, #set the maximum number of expected errors allowed in a read
                                      truncQ = 2, #truncate reads at the first instance of a quality score less than 2
                                      rm.phix = T, #discard reads that match against the phiX genome
                                      compress = T, 
                                      multithread = T #not sure if this will work
))

#estimate error rate
errF <- learnErrors(filtFs, multithread = T)

#sample inference
dadaFs <- dada(filtFs, err = errF, multithread = T)

#construct the sequence table
seqtab <- makeSequenceTable(dadaFs)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = T, verbose = T)

#track everything that's happened thus far
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               rowSums(seqtab.nochim))
#build dataframe
colnames(track) <- c("input", "filtered", "denoisedF", 
                     "nonchim")
rownames(track) <- microS

#convert to dataframe
trackDF <- track %>%
  as.data.frame() %>%
  mutate(PropPassed = nonchim/input)

#how many reads remain
rownames(seqtab.nochim) <- microS
print(paste0(round(sum(seqtab.nochim)/sum(seqtab)*100, 2), "% reads remain"))

#read in metadata
met <- read.csv("input/SampleMetadata.csv") %>%
  # factorise group
  mutate(Group = factor(Group, levels = c("Control",
                                          "Recruitment",
                                          "Follow up")))

##### Assign taxonomy #####
#bayesian approach employed by dada2
#standard taxonomy (97% identity)
taxa <- assignTaxonomy(seqtab.nochim, 
                       "input/microbiome/taxaDB/silva_nr99_v138.1_train_set.fa",
                       multithread = T)

#using species level taxonomy using silva database
taxa <- addSpecies(taxa, "input/microbiome/taxaDB/silva_species_assignment_v138.1.fa")

#assess (manual)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

#as a good few of later analyses are genus-level, I want to make
#sure that there is a genus-level entry for every ASV
#this is very difficult to do in a character matrix, so going to 
#convert to a dataframe to apply the changes and then revert
tax <- data.frame(taxa)
#replace NA "genus"
#starting with the most to least fine tuned phylum classification
for (i in 1:nrow(tax)) {
  if (is.na(tax$Genus[i])) {
    if (!is.na(tax$Family[i])) {
      tax$Genus[i] <- paste0(tax$Family[i], "_unclassified")
    } else {
      if (!is.na(tax$Order[i])) {
        tax$Genus[i] <- paste0(tax$Order[i], "_unclassified")
      } else {
        if (!is.na(tax$Class[i])) {
          tax$Genus[i] <- paste0(tax$Class[i], "_unclassified")
        } else {
          if (!is.na(tax$Phylum[i])) {
            tax$Genus[i] <- paste0(tax$Phylum[i], "_unclassified")
          } else {
            if(!is.na(tax$Kingdom[i])) {
              tax$Genus[i] <- paste0(tax$Kingdom[i], "_unclassified")
            }
          }
        }
      }
    }
  }
}

#convert back to a character matrix so that phyloseq can read it
tax <- suppressWarnings(tax_table(tax))
#return matrix column/row names
dimnames(tax) <- dimnames(taxa)

#check
tail(tax, n = 3)

##### Contamination Checks ####
#we do not have negative controls so this is done by considering
#certain assumptions about how kit / reagent contaminants would behave
#1. Prevalence = high. However, vaginal microbiomes are so dominated by
#single taxa (Lactobacillus) that moderate prevalence may also be
#evidence of a contamination, in combination with other factors.
#2. Abundance = low. Though detected across these samples, they
#should not be at particularly high abundance except in the case of
#point 3.
#3. Contaminants have higher relative abundance at lower library 
#sizes - perhaps because there is less biological sample mixed with
#reagents to be amplified up during amplicon sequencing. When plotted
#relative abundance (y) versus libary size, there is a characteristic 
#high peak at the lower library sizes followed by a sharp decline and 
#flat plateau for the rest of the distribution
#note this is not bullet proof but we must do something to try
#and control for this

#prepare phyloseq amenable sample data
samdf <- met %>%
  column_to_rownames(var = "Sample") 

#make object
phylo <- phyloseq(tax_table(tax),
                  sample_data(samdf),
                  otu_table(seqtab.nochim, taxa_are_rows = FALSE))

#Compute total library size per sample
libSizes <- sample_sums(phylo)
#Extract count matrix and taxonomy
cntMat <- data.frame(otu_table(phylo))
taxDF <- data.frame(tax_table(phylo))
#Compute relative abundance matrix
relAMat <- cntMat / rowSums(cntMat)

#For each ASV compute:
#prevalence (as proportion of samples where present)
#Spearman correlation of relative abundance vs library size
contamDF <- map_dfr(colnames(relAMat), function(x) {
  rel   <- relAMat[, x]
  cnt   <- cntMat[, x]
  prev  <- mean(cnt > 0)        
  genus <- taxDF[x, "Genus"]
  species <- taxDF[x, "Species"]
  # Only test ASVs present in at least 2 samples
  # (correlation requires variance)
  if (sum(cnt > 0) < 2) {
    return(data.frame(
      Seq = x,
      Genus = genus,
      Species = species,
      Prevalence = prev,
      rho = NA,
      pval = NA,
      adjP = NA,
      Flag = "Singleton — untestable"
    ))
  }
  #otherwise, test for negative correlation between relA and library size
  ct <- cor.test(rel, libSizes, method = "spearman", exact = FALSE)
  #store
  out <- data.frame(Seq = x, Genus = genus, Species = species,
                    Prevalence = prev, rho = unname(ct$estimate),
                    pval = ct$p.value, adjP = NA,   # filled after map
                    Flag = NA)
}) %>%
  mutate(adjP = p.adjust(pval, method = "BH"),
         #Contaminants expected to be:
         # - moderately to highly prevalent 
         # - negatively correlated with library size
         Flag = case_when(
           Flag == "Singleton — untestable" ~ "Singleton — untestable",
           # Strong likelihood of contaminant: high prevalence + significant negative correlation
           Prevalence >= 0.40 & rho < 0 & adjP < 0.05
           ~ "Likely contaminant",
           # Possible contaminant: moderate prevalence + suggestive negative correlation
           Prevalence >= 0.15 & rho < 0 & adjP < 0.20
           ~ "Possible contaminant",
           # Prevalent but no frequency signal — worth manual review
           # could be genuine low-abundance commensal or contaminant 
           # whose signal is masked by biological noise
           Prevalence >= 0.75 & (is.na(rho) | rho >= 0)
           ~ "High prevalence — review",
           TRUE ~ "Retain"
         )
  )

#assess flagged ASVs visually
#store flags
posContam <- contamDF %>%
  filter(!Flag == "Retain",
         !is.na(rho)) %>%
  select(Seq) %>%
  pull()

#prepare IDs
posContamDF <- contamDF %>%
  filter(Seq %in% posContam) %>%
  mutate(conID = paste0("Contam", 1:length(posContam)))

#taxDF edit
taxDF <- taxDF %>%
  rownames_to_column(var = "Seq")

#visualise
relA <- relAMat %>%
  data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, values_to = "RelA", names_to = "Seq") %>%
  filter(Seq %in% posContam) 
conPlot <- libSizes  %>%
  data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  inner_join(., relA, by = "Sample") %>%
  inner_join(., taxDF, by = "Seq") %>%
  inner_join(., posContamDF[,c(1,8:9)]) %>%
  mutate(Presence = ifelse(RelA == 0, "Not present", "Present"),
         Presence = factor(Presence, levels = c("Present", "Not present")),
         conID = factor(conID, levels = paste0("Contam", 1:length(posContam))))
#plot
ggplot(conPlot, aes(x = ., y = RelA)) + 
  geom_point(aes(colour = Presence, shape = Flag)) + 
  facet_wrap(~conID, scales = "free") +
  scale_colour_manual(values = c("Present" = "black",
                                 "Not present" = "lightgrey")) + 
  labs(x = "Library size",
       y = "Relative abundance") + 
  theme(legend.position = "bottom")
ggsave("processed/microbiome/PossibleContaminants.png",
       width = 20, height = 20, unit = "cm")

#from looking at this data - and looking at the literature I will be removing
#the following:
#Contam1 and 2 - both Acinetobacter, both notorious reagent contaminants
#Contam3 - Delftia - despite there being evidence of Delftia and the vaginal
#microbiome on closer inspection none of these top cited articles try
#to control for potential contaminants and there is not enough evidence
#to overcome the signal in this data suggesting it is a reagent microbe
#Delftia are largely water-based bacteria found in lakes, gas effluent, etc
#can be opportunistic pathogens but again - the distribution of relative
#abundance to library size is textbook contamination patterns and it
#is abundant in the environment and has been found in sequencing kits multiple times
#across the literature
#contam6 - Rhizobium complex. These are plant microbes and also known kit
#contaminants.
#contam7 - Comamonas. Again graphically look like contaminants though could
#be skin microbes that have migrated. 
conToRem <- contamDF %>%
  filter(Genus %in% c("Acinetobacter", "Delftia", 
                      "Comamonas", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium")) %>%
  select(Seq) %>%
  pull()
#remove
seqtab.nocont <- seqtab.nochim[, !colnames(seqtab.nochim) %in% conToRem]

#samples that were more than 50% suspected contamination will be removed
mostCon <- (rowSums(seqtab.nocont) / rowSums(seqtab.nochim)) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  filter(as.numeric(.) < 0.5) %>%
  select(Sample) %>%
  pull()

#remove
seqtab.nocont <- seqtab.nocont[! rownames(seqtab.nocont) %in% mostCon, ]

##### PhyloSeq: Final Processing ####
#Remove samples that 1) were removed from this analysis due to 
#sequencing QC and 2) do not have the second time point 
met <- met %>%
  filter(Sample %in% rownames(seqtab.nocont))

#now determine which atrophy patients do not have
#both visits
oneTime <- met %>%
  filter(Status == "Atrophy") %>%
  group_by(Patient.ID) %>%
  tally() %>%
  filter(n == 1) %>%
  select(Patient.ID) %>%
  pull()
#remove
met <- met %>%
  filter(!Patient.ID %in% oneTime)

#generate counts
cnts <- data.frame(seqtab.nocont)
#remove samples without > 1 timepoint
cnts <- cnts[rownames(cnts) %in% met$Sample, ]

#make new phyloseq object with the removed samples (necessary for contamination
#exploration, but cannot be used in downstream analyses)
samdf <- met %>%
  column_to_rownames(var = "Sample") 
phylo <- phyloseq(tax_table(tax),
                  sample_data(samdf),
                  otu_table(cnts, taxa_are_rows = FALSE))

#sequence names are annoying so
#store the sequences for use if needed
dna <- DNAStringSet(taxa_names(phylo))
names(dna) <- taxa_names(phylo)
phylo <- merge_phyloseq(phylo, dna)
taxa_names(phylo) <- paste0("ASV", seq(ntaxa(phylo)))

#also apply to reads
names(dna) <- taxa_names(phylo)

#filtering
#taxonomic filtering
#ensure mitochondrial sequences aren't retained
if (any(grepl("Mitochondria", tax))) {
  phylo <- subset_taxa(phylo, Family != "Mitochondria")
}

#keep only ASVs that have been identified as bacteria
phylo <- subset_taxa(phylo, Kingdom == "Bacteria")

#remove very rare taxa
phylo <- prune_taxa(taxa_sums(phylo) > 10, phylo)

#check that samples have survived contaminant removal and other filtering
lowBio <- rownames(cnts)[rowSums(data.frame(otu_table(phylo))) < 1000]
#add corresponding paired samples
lowPat <- met$Patient.ID[met$Sample %in% lowBio]
toKeep <- met$Sample[!met$Patient.ID %in% lowPat]

#remove samples with low biological signal
phylo <- prune_samples(toKeep, phylo)
#remove from metadata
met <- met %>%
  filter(Sample %in% toKeep)

#collapsing to other levels
#collapse counts to genus level; summing species within
phyGen <- tax_glom(phylo, taxrank = "Genus")

#remove any genera that appear in < 5% of all samples
phyGen <- filter_taxa(phyGen, function(x) sum(x > 0) > (0.05 * nsamples(phylo)), TRUE)

#collapse counts to phylum level
phyPhy <- tax_glom(phylo, taxrank = "Phylum")

##### Community Composition #####
#extract counts from phylum object
#rename ASV ids with taxa names
#add metadata
phyCnts <- data.frame(otu_table(phyPhy))
colnames(phyCnts) <- data.frame(tax_table(phyPhy))$Phylum

phyCnts <- phyCnts %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample) %>%
  inner_join(., met, by = "Sample") %>%
  dplyr::rename("Phylum" = "name", 
                "counts" = "value") %>%
  #add relative abundance per sample
  group_by(Sample) %>%
  mutate(relA = counts / sum(counts), .after = counts)

#order the patient IDs correctly
pIds <-  unique(phyCnts$Patient.ID[phyCnts$Status == "Atrophy"])
ordNo <- pIds[order(formatC(pIds, width = nchar(max(pIds)),
                            flag = "0"))]

#phylum level colour palette
phyPal <- vanPal
names(phyPal)[1] <- "Firmicutes"
names(phyPal)[2] <- "Bacteroidota"
names(phyPal)[4] <- "Campylobacterota"
names(phyPal)[5] <- "Actinobacteriota"
names(phyPal)[7] <- "Fusobacteriota"
names(phyPal)[8] <- "Proteobacteria"

#produce the community plots by group
pBA <- phyCnts %>%
  filter(Status == "Atrophy") %>%
  mutate(Patient.ID = droplevels(factor(Patient.ID)),
         Patient.ID = factor(Patient.ID, 
                             levels = ordNo),
         Group = factor(Group,
                            levels = c("Recruitment", "Follow up"))) %>%
  ggplot(aes(x = (Patient.ID), y = relA)) + 
  theme_bw(base_size = 13) +
  geom_bar(stat = "identity", aes(fill = Phylum),
           colour = "black", linewidth = 0.25) +
  theme(legend.position = "bottom",
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank()) +
  labs(y = "Relative abundance",
       x = "Patient") +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = phyPal) +
  facet_grid(Group~.) +
  guides(fill = guide_legend(nrow = 2))
#control sample plot
pC <- phyCnts %>%
  filter(!Status == "Atrophy") %>%
  mutate(Patient.ID = droplevels(factor(Patient.ID))) %>%
  ggplot(aes(x = (Patient.ID), y = relA)) + 
  theme_bw(base_size = 13) +
  geom_bar(stat = "identity", aes(fill = Phylum),
           colour = "black", linewidth = 0.25) +
  theme(legend.position = "bottom",
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank()) +
  labs(y = "",
       x = "") +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = phyPal) +
  facet_grid(Group~.) +
  guides(fill = guide_legend(nrow = 2))
#arrange
phyCom <- ggarrange(pC, pBA, ncol = 1,
                    common.legend = T, heights = c(0.85, 1.85),
                    legend = "bottom")
#save
ggsave("output/microbiome/community/CommunityComposition_Phy.png")

#genus-level plots
#extract counts from phylum object
#rename ASV ids with taxa names
#add metadata
genCnts <- data.frame(otu_table(phyGen))
colnames(genCnts) <- data.frame(tax_table(phyGen))$Genus
genCnts <- genCnts %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample) %>%
  inner_join(., met, by = "Sample") %>%
  dplyr::rename("Genus" = "name", 
                "counts" = "value") %>%
  #add relative abundance per sample
  group_by(Sample) %>%
  mutate(relA = counts / sum(counts), .after = counts)

#determine dominant genera
#average relative abundances ordered by size (largest to smallest)
#and take the 15 most commonly abundant genera
topbyRel <- genCnts %>%
  mutate(Genus = gsub("_unclassified", " (unclassified)", Genus)) %>%
  group_by(Genus) %>%
  mutate(AvgRel = mean(relA), .after = Sample) %>%
  ungroup() %>%
  select(AvgRel, Genus) %>%
  unique() %>%
  arrange(-AvgRel) %>%
  select(Genus) %>%
  head(n = 15) %>%
  pull()
#prepare palettes
palFun <- colorRampPalette(vanPal)
plotPal <- palFun(length(topbyRel))
plotPal <- c(plotPal, "white")
names(plotPal) <- c(topbyRel, "Other")
#this produces two indistinguishable colours so will manually fix
plotPal["Prevotella"] <- "#007474"
plotPal["Bifidobacterium"] <- "#ffc862"
plotPal["Corynebacterium"] <- "#ff4500"
plotPal["Gardnerella"] <- "#0000d9"
plotPal["Streptococcus"] <- "#32cd32"
plotPal["Escherichia-Shigella"] <- "#9A1814"
plotPal["Prevotellaceae (unclassified)"] <- "#ffff00"
plotPal["Dialister"] <- "#800080"
plotPal["Atopobium"] <- "#df5d83"

#prepare plot metrics
plotDF <- genCnts %>%
  group_by(Sample, Genus) %>%
  mutate(relA = sum(relA)) %>%
  ungroup() %>%
  select(Sample, Patient.ID, relA, Genus, Group) %>%
  unique() %>%
  mutate(Genus = gsub("_unclassified", " (unclassified)", Genus),
         Genus = ifelse(Genus %in% topbyRel, Genus, "Other")) %>%
  group_by(Sample, Genus) %>%
  mutate(relA = sum(relA),
         Group = factor(Group,
                            levels = c("Control", "Recruitment", "Follow up")),
         Genus = factor(Genus, levels = c(topbyRel, "Other"))) %>%
  ungroup() %>%
  select(Sample, Patient.ID, relA, Genus, Group) %>%
  unique() 
#plot cohorts
pBA <- plotDF %>%
  filter(!Group == "Control") %>%
  mutate(Patient.ID = droplevels(factor(Patient.ID)),
         Patient.ID = factor(Patient.ID, 
                             levels = ordNo)) %>%
  ggplot(aes(x = (Patient.ID), y = relA)) + 
  theme_bw(base_size = 13) +
  geom_bar(stat = "identity", aes(fill = Genus),
           colour = "black", linewidth = 0.25) +
  theme(legend.position = "bottom",
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(face = "italic")) +
  labs(y = "Relative abundance",
       x = "Patient") +
  guides(fill = guide_legend(ncol = 7)) +
  scale_fill_manual(values = plotPal) +
  facet_grid(Group~.) + 
  guides(fill = guide_legend(nrow = 4))
#plot control subjects
pC <- plotDF %>%
  filter(Group == "Control") %>%
  mutate(Patient.ID = droplevels(factor(Patient.ID))) %>%
  ggplot(aes(x = (Patient.ID), y = relA)) + 
  theme_bw(base_size = 13) +
  geom_bar(stat = "identity", aes(fill = Genus),
           colour = "black", linewidth = 0.25) +
  theme(legend.position = "bottom",
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        legend.text = element_text(face = "italic")) +
  labs(y = "",
       x = "") +
  guides(fill = guide_legend(ncol = 7)) +
  scale_fill_manual(values = plotPal) +
  facet_grid(Group~.) + 
  guides(fill = guide_legend(nrow = 4))
#arrange
genCom <- ggarrange(pC, pBA, ncol = 1,
                    common.legend = T, heights = c(0.85, 1.85),
                    legend = "bottom")
#save
ggsave("output/microbiome/community/CommunityComposition_Gen.png")

#store prevalence
prevDF <- genCnts %>%
  ungroup() %>%
  mutate(Presence = ifelse(counts > 0, 1, 0),
         .after = Sample) %>%
  select(Genus, Presence, Group) %>%
  group_by(Genus) %>%
  mutate(Prevalence = sum(Presence) / length(unique(genCnts$Sample)), #tot samples
         .after = Genus) %>%
  ungroup() %>%
  select(Genus, Prevalence) %>%
  unique() %>%
  arrange(-Prevalence)

##### Outputs #####
#manuscript figure
ggarrange(plotlist = list(phyCom, genCom),
          ncol = 1, labels = "AUTO")
ggsave("output/figures/Fig-community.png",
       height = 15, width = 10, units = "in")
ggsave("output/figures/Fig-community.pdf",
       height = 15, width = 10, units = "in")
ggsave("output/figures/Fig-community.jpg",
       height = 15, width = 10, units = "in")

#contamination investigations
contamDF %>%
  #add contamination IDs in for those that were more closely investigated
  full_join(., posContamDF[c("Seq", "conID")], by = "Seq") %>%
  rename("ContaminantID" = "conID") %>%
  write.csv(., "processed/microbiome/Contam_All.csv",
          row.names = F, quote = F)
write.csv(posContamDF, "output/microbiome/PotentialContam.csv",
          row.names = F, quote = F)

#genus colour scheme
genPal <- plotPal
save(genPal, file = "processed/microbiome/AbundantGeneraPalette.rds")

#output of preprocessing - metadata + dada2 stats
trackDF %>%
  rownames_to_column(var = "Sample") %>%
  inner_join(met, .) %>% 
  write.csv(., "output/microbiome/MicroPreProcessing_Data.csv",
          quote = F, row.names = F)

#dada2 error rates
png("processed/microbiome/dada2_QC_ErrorRates.png")
plotErrors(errF, nominalQ = T)
dev.off()

#dada2 descriptive statistics
write.csv(trackDF, "processed/microbiome/ReadDada2_tracked.csv",
          row.names = T, quote = F)

#phyloseq objects
save(phylo, phyGen, phyPhy, file = "processed/microbiome/PhyloSeqObjs.rds")

#prevalence
write.csv(prevDF, "output/microbiome/community/Prevalence.csv",
          quote = F, row.names = F)

#record samples used in analysis
met %>%
  filter(Sample %in% rownames(cnts)) %>%
  mutate(Microbiome = "Present") %>%
  select(Sample, Patient.ID, Microbiome) %>%
  write.csv(., "processed/Sample_AnalysisRegister.csv",
            quote = F, row.names = F)

#session information
writeLines(capture.output(sessionInfo()), "output/microbiome/PreProcessing_sessionInfo.txt")
