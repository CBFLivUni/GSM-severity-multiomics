#Changed into script 13th December
#Original mid 2024
#Lauren Mee

##### Set up #####
#set seed
set.seed(1049)

#libraries
libs <- c("variancePartition", #CCA analysis
          "lme4", "lmerTest",
          "ggpubr", "tidyverse") #general plots etc

for (pkg in libs) {
  library(pkg, character.only = T)
}

#load data
#metadata
met <- read.csv("input/SampleMetadata.csv") %>%
  mutate(Group = factor(Group, levels = c("Control", "Recruitment", "Follow up")))

#nmr data
nmr <- read.csv("input/nmr/VANS_tampons_master_data_matrix_updated_Sep2024_missing_values_replaced.csv")

#output directory architecture
dirs <- c("processed/", "processed/nmr/",
          "output/", "output/nmr/")

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d)
  }
}

#misc
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

##### Functions ####
#extract ID from NMR data that is in a different format to the microbiome 
#data I've been using so far
getID <- function(sampleID){
  #extract sample ID no
  dig <- str_extract_all(gsub("_0+", "", sampleID), pattern = "\\d+")
  #determine if there's a B or not
  treat <- ifelse(str_detect(sampleID, "B"), "B", "")
  id <- paste("VAN", treat, dig, sep = "")
  return(id)
}

#fix metabolite ID
#unknowns retain their numeric part of the ID
fixMet <- function(metabolite) {
  tmp <- metabolite
  tmp <- gsub("_(?=[Uu]nknown)", "-", tmp, perl = TRUE)
  tmp <- gsub("(?<=[Uu]nknown)_", "-", tmp, perl = TRUE)
  tmp <- gsub("_(?=[0-9])", "-", tmp, perl = TRUE)
  tmp <- gsub("\\.", "-", tmp)
  tmp <- gsub("-\\(", " (", tmp)
  tmp <- gsub("O_", "O-", tmp)
  tmp <- gsub("tinephos", "tine phos", tmp)
  tmp <- gsub("_", "\n", tmp)
  return(tmp)
}

##### Preprocessing ####
#convert the nmr ids to those that match the materials used thus far
nmr <- nmr %>%
  mutate(Sample = unlist(map(ID, getID)),
         .after = ID) %>%
  select(-ID, -Group)

#add group variable
nmr <- nmr %>%
  inner_join(met[, c("Sample", "Group")], ., by = "Sample")

#remove any atrophy samples that don't have both timepoints
#id these patients
oneTime <- met %>%
  filter(Sample %in% nmr$Sample,
         Status == "Atrophy") %>%
  group_by(Patient.ID) %>%
  tally() %>%
  filter(n == 1) %>%
  inner_join(., met, by = "Patient.ID") %>%
  select(Sample) %>%
  pull()
#remove from data
nmr <- nmr %>%
  filter(!Sample %in% oneTime)

##### Normalise and scale ####
#visualise spectra
nmr %>%
  select(-Group) %>%
  pivot_longer(-Sample) %>%
  inner_join(., met, by = "Sample") %>%
  ggplot(aes(x = Sample, y = value, group = name)) +
  theme_bw(base_size = 13) +
  geom_line(aes(colour = Sample), alpha = 0.75) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom") +
  #scale_colour_manual(values = vanPal) +
  labs(x = "Sample",
       y = "Intensity") +
  guides(colour = "none")
#save
ggsave("processed/nmr/RawSpectra.png")

#by group
nmr %>%
  select(-Group) %>%
  pivot_longer(-Sample) %>%
  inner_join(., met, by = "Sample") %>%
  mutate(Group = factor(Group,
                            levels = c("Control",
                                       "Recruitment",
                                       "Follow up"))) %>%
  ggplot(aes(x = Sample, y = value)) +
  theme_bw(base_size = 13) +
  geom_boxplot(aes(colour = Group), alpha = 0.75) +
  guides(colour = "none") +
  labs(x = "Sample",
       y = "Intensity") +
  theme(axis.text.x = element_blank()) +
  scale_colour_manual(values = vanPal)
#save
ggsave("processed/nmr/RawSpectra_Boxplot_ByGroup.png")

#convert nmr df to matrix
#metabolites end up rows, samples = columns
nmrMat <- nmr %>%
  column_to_rownames(var = "Sample") %>%
  select(-Group) %>%
  as.matrix() %>%
  t()

#compute reference spectrum
ref <- apply(nmrMat, 1, median)
q <- nmrMat / ref
qMed <- apply(q, 2, median)
pqnMat <- t(nmrMat)/qMed

#log2 transform (zeros have already been imputed by collaborator)
logMatLim <- log2(pqnMat)
#visualise
logMatLim %>%
  data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample) %>%
  inner_join(., met, by = "Sample") %>%
  mutate(Group = factor(Group,
                            levels = c("Control",
                                       "Recruitment",
                                       "Follow up"))) %>%
  ggplot(aes(x = Sample, y = value, group = Sample)) +
  theme_bw(base_size = 13) +
  geom_boxplot(aes(colour = Group), alpha = 0.75) +
  scale_colour_manual(values = vanPal) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom") +
  labs(y = "log2(PQN-Intensity)")
#save
ggsave("processed/nmr/log2PQN_Spectra_Boxplot_ByGroup.png")

#plot PCA
pca <- prcomp(logMatLim, scale = F)
exVars <- round(summary(pca)$importance[2,]*100, 2)
exVars <- exVars[exVars >= 1]
plotDF <- pca$x %>%
  data.frame() %>%
  select(PC1, PC2, PC3, PC4) %>%
  rownames_to_column(var = "Sample") %>%
  inner_join(., met, by = "Sample") 
cen <- plotDF %>%
  group_by(Group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2)) %>%
  as.data.frame()
plotDF %>%
  mutate(Group = factor(Group,
                            levels = c("Control",
                                       "Recruitment",
                                       "Follow up"))) %>%
  ggplot(aes(x = PC1, y = PC2, colour = Group, fill = Group)) +
  theme_bw(base_size = 13) +
  geom_point(aes(colour = Group)) +
  scale_colour_manual(values = vanPal) +
  scale_fill_manual(values = vanPal) +
  labs(x = paste0("PC1 (", exVars[1], "%)"),
       y = paste0("PC2 (", exVars[2], "%)")) +
  geom_point(data = cen, size = 5, 
             shape = 21, colour = "black", 
             aes(fill = Group),
             show.legend = FALSE) +
  stat_ellipse(aes(colour = Group)) 
#save
ggsave("output/nmr/PCA_AllSpectra.png")


##### Metabolites and NGAT / DIVA: LMM ####
#collinearity checks from previous need to be redone as this is 
#a slightly different subset of the full cohort
#subset met to include only those that are kept in the nmr analysis
met <- met[met$Sample %in% rownames(logMatLim), ]

#now check what covariates have multiple levels - factor categorical
metGSM <- met %>%
  filter(Status == "Atrophy") %>%
  mutate(Patient.ID = as.factor(Patient.ID),
         Group = factor(Group, levels = c("Recruitment",
                                             "Follow up")),
         Ethnicity = as.factor(Ethnicity),
         Smoking.status = factor(Smoking.status), 
         Menopausal.status = factor(Menopausal.status))
#will only be running LMM on useable factors
summary(metGSM)

#n = 42 per GSM timepoint cohort. Changes - we now have 2 smokers, and a 
#single pre-menopausal patient. Unfortunately this is really too small
#numbers to stably estimate variance against 40/41 others so these will not be
#considered in the modelling
#check other variables for collinearity
form <- ~ Group + Age + BMI + DIVA + NGAT
corMat <- canCorPairs(form, metGSM)
#plot
corMat %>%
  data.frame() %>%
  rownames_to_column(var = "one") %>%
  pivot_longer(-one, names_to = "two", values_to = "Correlation") %>%
  ggplot(aes(x = one, y = two)) +
  theme_bw(base_size = 13) +
  geom_tile(aes(fill = Correlation)) +
  geom_label(aes(label = round(Correlation,2))) +
  scale_fill_gradient2(low = vanPal[[6]],
                       mid = vanPal[[7]],
                       high = "#050808",
                       midpoint = 0.5) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank()) 
#save
ggsave("output/nmr/clinicalMetricCollinear_CCA.png")

#again, nothing really to worry about save Age & NGAT (not as
#correlated as they are in the microbiome cohort)

#prepare to run LMMs
#there will be two models: full and reduced
#full: met ~ DIVA + NGAT + Group + BMI + Age + (1|Patient.ID)
#reduced: met ~ DIVA + NGAT + Group + (1|Patient.ID)
modDat <- logMatLim %>%
  data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  inner_join(., metGSM[c("Sample", "NGAT", "DIVA",
                         "Group", "Patient.ID",
                         "BMI", "Age")],
             by = "Sample")
#store metabolites
metVec <- colnames(logMatLim)

#full model
fullForm <- "~ DIVA + NGAT + Group + BMI + Age + (1|Patient.ID)"
fullOutDF <- map_dfr(metVec, function(m) {
  #run model
  tmod <- lmer(paste0(m,  fullForm),
               data = modDat)
  #extract summary
  tSum <- summary(tmod)
  #extract coefficients
  est <- tSum$coefficients
  #extract variance 
  var <- as.data.frame(VarCorr(tmod))
  #store vars
  piVar <- var$vcov[var$grp == "Patient.ID"]
  resVar <- var$vcov[var$grp == "Residual"]
  #collate
  out <- data.frame(
    Metabolite = m,
    Model = "Full",
    Formula = fullForm,
    Patient.IDVariance = piVar,
    ResidualVariance = resVar,
    #variance of patient
    Patient_ICC = piVar / (piVar + resVar),
    #collect estimate (coefficient), pvalue and SE per term
    DIVA_estimate = est["DIVA", "Estimate"],
    DIVA_se = est["DIVA", "Std. Error"],
    DIVA_pvalue = est["DIVA", "Pr(>|t|)"],
    NGAT_estimate = est["NGAT", "Estimate"],
    NGAT_se = est["NGAT", "Std. Error"],
    NGAT_pvalue = est["NGAT", "Pr(>|t|)"],
    Group_estimate = est["GroupFollow up", "Estimate"],
    Group_se = est["GroupFollow up", "Std. Error"],   
    Group_pvalue = est["GroupFollow up", "Pr(>|t|)"],
    BMI_estimate = est["BMI", "Estimate"],
    BMI_se   = est["BMI", "Std. Error"],
    BMI_pvalue = est["BMI", "Pr(>|t|)"],
    Age_estimate = est["Age", "Estimate"],
    Age_se   = est["Age", "Std. Error"],
    Age_pvalue = est["Age", "Pr(>|t|)"],
    #singular fit warning
    Singular = isSingular(tmod)
  )
})
#apply corrections
fullOutDF <- fullOutDF %>%
  mutate(DIVA_adjP = p.adjust(DIVA_pvalue, method = "BH"),
         NGAT_adjP = p.adjust(NGAT_pvalue, method = "BH"),
         Age_adjP = p.adjust(Age_pvalue, method = "BH"),
         Group_adjP = p.adjust(Group_pvalue, method = "BH"))

#hits for DIVA
if(any(fullOutDF$DIVA_adjP < 0.05)) {
  DIVAFullHits <- fullOutDF %>%
    filter(Singular == FALSE,
           DIVA_adjP < 0.05) 
  print(nrow(DIVAFullHits))
  print(head(DIVAFullHits))
} else {
  print("No metabolite significantly associates with DIVA score")
}


#NGAT
if(any(fullOutDF$NGAT_adjP < 0.05)) {
  NGATFullHits <- fullOutDF %>%
    filter(Singular == FALSE,
           NGAT_adjP < 0.05) %>%
    arrange(NGAT_adjP)
  print(nrow(NGATFullHits))
  #28 metabolites associated when controlling for patient variability
  print(head(NGATFullHits))
} else {
  print("No metabolite significantly associates with NGAT score")
}

#reduced model
redForm <- "~ DIVA + NGAT + Group + (1|Patient.ID)"
redOutDF <- map_dfr(metVec, function(m) {
  #run model
  tmod <- lmer(paste0(m,  redForm),
               data = modDat)
  #extract summary
  tSum <- summary(tmod)
  #extract coefficients
  est <- tSum$coefficients
  #extract variance 
  var <- as.data.frame(VarCorr(tmod))
  #store vars
  piVar <- var$vcov[var$grp == "Patient.ID"]
  resVar <- var$vcov[var$grp == "Residual"]
  #collate
  out <- data.frame(
    Metabolite = m,
    Model = "Reduced",
    Formula = redForm,
    Patient.IDVariance = piVar,
    ResidualVariance = resVar,
    Patient_ICC = piVar / (piVar + resVar),
    #collect estimate (coefficient), pvalue and SE per term
    DIVA_estimate = est["DIVA", "Estimate"],
    DIVA_se = est["DIVA", "Std. Error"],
    DIVA_pvalue = est["DIVA", "Pr(>|t|)"],
    NGAT_estimate = est["NGAT", "Estimate"],
    NGAT_se = est["NGAT", "Std. Error"],
    NGAT_pvalue = est["NGAT", "Pr(>|t|)"],
    Group_estimate = est["GroupFollow up", "Estimate"],
    Group_se = est["GroupFollow up", "Std. Error"],   
    Group_pvalue = est["GroupFollow up", "Pr(>|t|)"],
    #singular fit warning
    Singular = isSingular(tmod)
  )
})
#apply corrections
redOutDF <- redOutDF %>%
  mutate(DIVA_adjP = p.adjust(DIVA_pvalue, method = "BH"),
         NGAT_adjP = p.adjust(NGAT_pvalue, method = "BH"),
         Group_adjP = p.adjust(Group_pvalue, method = "BH"))

#hits for DIVA
if(any(redOutDF$DIVA_adjP < 0.05)) {
  DIVARedHits <- redOutDF %>%
    filter(Singular == FALSE,
           DIVA_adjP < 0.05) 
  print(nrow(DIVARedHits))
  print(head(DIVARedHits))
} else {
  print("No metabolite significantly associates with DIVA score")
}

#NGAT
if(any(redOutDF$NGAT_adjP < 0.05)) {
  NGATRedHits <- redOutDF %>%
    filter(Singular == FALSE,
           NGAT_adjP < 0.05) %>%
    arrange(NGAT_adjP)
  print(nrow(NGATRedHits))
  #81 metabolites associated when controlling for patient variability
  print(head(NGATRedHits))
} else {
  print("No metabolite significantly associates with NGAT score")
}

#visualise top hits
#Standardise post hoc rather than refitting: for a linear model the
#standardised coefficient is exactly Estimate * sd(predictor) / sd(response),
#so no model needs to be run again. This makes estimates comparable BOTH
#across metabolites (which differ in dynamic range) and across predictors
#(only interested in visualising NGAT, DIVA, Age (as NGAT confounder) and Group.
#get metabolite standard deviations
sdMat  <- apply(logMatLim, 2, sd)
#get predictor sds
sdPred <- c(NGAT  = sd(modDat$NGAT),
            DIVA  = sd(modDat$DIVA),
            Age   = sd(modDat$Age),
            #binary term: SD of the 0/1 indicator, so it is on the same
            #footing as the continuous predictors (Gelman 2008)
            Group = sd(as.numeric(modDat$Group) - 1))

#need SE, pvalue, estimate per term, per model
allModDF <- bind_rows(fullOutDF, redOutDF)
#prepare to plot
plotDF <- allModDF %>%
  #remove singular fits
  filter(Singular == "FALSE") %>%
  select(Metabolite, Model,
         matches("^(NGAT|DIVA|Age|Group)(_estimate|_pvalue|_se|_adjP)")) %>%
  pivot_longer(-c(Metabolite, Model), 
               names_to = c("Term", ".value"),
               names_pattern = "^(NGAT|DIVA|Age|Group)_(estimate|pvalue|se|adjP)$") %>%
  #remove Age + Reduced pairing as there are no accompanying stats
  filter(!(Term == "Age" & Model == "Reduced")) %>%
  #add SD metrics
  mutate(scaleFac = sdPred[Term] / sdMat[Metabolite],
         StdEst = estimate * scaleFac,
         Lower = (estimate - 1.96 * se) * scaleFac,
         Upper = (estimate + 1.96 * se) * scaleFac,
         #order required variables
         Model = factor(Model, levels = c("Reduced", "Full")),
         Term = factor(Term, levels = c("NGAT", "DIVA", "Group", "Age")),
         #assign significance
         Sig = factor(ifelse(adjP < 0.05, "BH p < 0.05", "n.s."),
                           levels = c("BH p < 0.05", "n.s."))) %>%
  data.frame()
#two versions - top significant metabolites under reduced model (n = 30,
#to reduce if too busy), and those significant in both models

#top 30 most sig
top30Red <- NGATRedHits %>%
  head(n = 30) %>%
  select(Metabolite) %>%
  pull()
NGATRedHits <- NGATRedHits %>%
  mutate(MetLab = map_vec(Metabolite, fixMet))
pTop30LMM <- plotDF %>%
  filter(Metabolite %in% top30Red) %>%
  mutate(#add nice labels
         MetLab = map_vec(Metabolite, fixMet),
         #order alphabetical so metabolite bins are together
         MetLab = factor(MetLab, levels = rev(sort(NGATRedHits$MetLab))))%>%
  ggplot(aes(x = StdEst, y = MetLab, colour = Model, shape = Sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_linerange(aes(xmin = Lower, xmax = Upper),
                 position = position_dodge(width = 0.55), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.55), 
             size = 2.8, fill = "white") +
  facet_grid(~Term) +
  scale_colour_manual(values = c(vanPal[5], vanPal[7])) +
  scale_shape_manual(values = c("BH p < 0.05" = 16, "n.s." = 21)) +
  labs(y = NULL,
       x = "Standardised effect\n(SD units)",
       shape = NULL, colour = NULL)
ggsave("output/nmr/LMM_CombinedEffectPlot_top30RedBins.png")

#28 from full model
NGATFullHits <- NGATFullHits %>%
  mutate(MetLab = map_vec(Metabolite, fixMet))
pFullLMM <- plotDF %>%
  filter(Metabolite %in% NGATFullHits$Metabolite) %>%
  #add nice labels
  mutate(MetLab = map_vec(Metabolite, fixMet),
    #order alphabetical so metabolite bins are together
    MetLab = factor(MetLab, levels = rev(sort(NGATFullHits$MetLab)))) %>%
  ggplot(aes(x = StdEst, y = MetLab, colour = Model, shape = Sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_linerange(aes(xmin = Lower, xmax = Upper),
                 position = position_dodge(width = 0.55), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.55), 
             size = 2.8, fill = "white") +
  facet_grid(~Term) +
  scale_colour_manual(values = c(vanPal[5], vanPal[7])) +
  scale_shape_manual(values = c("BH p < 0.05" = 16, "n.s." = 21)) +
  labs(y = NULL,
       x = "Standardised effect\n(SD units)",
       shape = NULL, colour = NULL) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 7))
ggsave("output/nmr/LMM_CombinedEffectPlot_All28FullModBins.png")

#slim version (just severity metrics)
pFullLMMslim <- plotDF %>%
  filter(Metabolite %in% NGATFullHits$Metabolite) %>%
  #add nice labels
  mutate(MetLab = map_vec(Metabolite, fixMet),
         #order alphabetical so metabolite bins are together
         MetLab = factor(MetLab, levels = rev(sort(NGATFullHits$MetLab)))) %>%
  #keep only severity metrics
  filter(Term %in% c("NGAT", "DIVA")) %>%
  ggplot(aes(x = StdEst, y = MetLab, colour = Model, shape = Sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_linerange(aes(xmin = Lower, xmax = Upper),
                 position = position_dodge(width = 0.55), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.55), 
             size = 2.8, fill = "white") +
  facet_grid(~Term) +
  scale_colour_manual(values = c(vanPal[5], vanPal[7])) +
  scale_shape_manual(values = c("BH p < 0.05" = 16, "n.s." = 21)) +
  labs(y = NULL,
       x = "Standardised effect\n(SD units)",
       shape = NULL, colour = NULL) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 7)) +
  guides(colour = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2, ncol = 1))
#save output later - manuscript figure of choice

##### Metabolites and NGAT / DIVA: Correlation ####
#assess correlation directly between NGAT and significant hits 
#LMM results are primary - the following is just to assess if the
#bins significant above are correlated with NGAT at both visits

#run correlations
#NGAT is ordinal and ties will exist - will use
#Kendall to try and compensate
#looking at recruitment first
metas <- NGATFullHits$Metabolite
resRec <- map_dfr(metas, function(x) {
  df <- modDat %>%
    filter(Group == "Recruitment")
  test <- cor.test(df[[x]], df$NGAT,
                   method = "kendall")
  out <- data.frame(Cohort = "Recruitment",
                    MetaboliteBin = x, 
                    Statistic = names(test$statistic),
                    StatValue = unname(test$statistic),
                    tau = unname(test$estimate),
                    Exact = (names(test$statistic) == "T"),
                    pval = test$p.value)
  return(out)
})
#correct
resRec$adjP <- p.adjust(resRec$pval, method = 'BH')
#are they all still significant? 
all(resRec$adjP < 0.05)
#nope
resRec %>%
  filter(adjP < 0.05)

#follow up
resFu <- map_dfr(metas, function(x) {
  df <- modDat %>%
    filter(Group == "Follow up")
  test <- cor.test(df[[x]], df$NGAT,
                   method = "kendall")
  out <- data.frame(Cohort = "Follow up",
                    MetaboliteBin = x, 
                    Statistic = names(test$statistic),
                    StatValue = unname(test$statistic),
                    tau = unname(test$estimate),
                    Exact = (names(test$statistic) == "T"),
                    pval = test$p.value)
  return(out)
})
#correct
resFu$adjP <- p.adjust(resFu$pval, method = 'BH')
#are they all still significant? 
all(resFu$adjP < 0.05)
#nope
resFu %>%
  filter(adjP < 0.05)

#combine
corOut <- bind_rows(resRec, resFu)

#produce a table of the metabolites of interest 
#stats reporting will come from reduced model
sumDF <- redOutDF %>%
  #order by effect size
  #arrange(-abs(NGAT_estimate)) %>%
  arrange(NGAT_adjP) %>%
  #get metabolites of each bin
  mutate(BinLab = map_vec(Metabolite, fixMet),
         Met = map_vec(BinLab, function(m) {
           ifelse(grepl("^Unknown", m), m, gsub("-[0-9]+$", "", m, perl = TRUE))
           })) %>%
  #get metabolite level metrics
  group_by(Met) %>%
  mutate(NBins = n(),
         NSig = sum(NGAT_adjP < 0.05),
         #estimates from significant bins only
         LowerEst = min(NGAT_estimate[NGAT_adjP < 0.05]),
         UpperEst = max(NGAT_estimate[NGAT_adjP < 0.05]), 
         Estimates = ifelse(LowerEst == UpperEst, 
                            as.character(round(LowerEst, 3)),
                            paste(formatC(round(LowerEst, 3), digits = 3, 
                                          format = "f", flag = "#"), 
                                  formatC(round(UpperEst, 3), digits = 3, 
                                          format = "f", flag = "#"), 
                                  sep = " to ")),
         Direction = ifelse(LowerEst < 0, "Negative", "Positive"),
         SigBins = paste(NSig, NBins, sep = " of ")) %>%
  #keep only the metabolites that are significant in both LMM models
  filter(Metabolite %in% intersect(NGATRedHits$Metabolite, 
                                   NGATFullHits$Metabolite)) %>%
  #slim down
  select(Met, SigBins, Estimates, Direction) %>%
  unique() %>%
  data.frame()

#get the correlation results to share
sumDF <- corOut %>%
  mutate(BinLab = map_vec(MetaboliteBin, fixMet),
         Met = map_vec(BinLab, function(m) {
           ifelse(grepl("^Unknown", m), m, gsub("-[0-9]+$", "", m, perl = TRUE))
         })) %>%
  #keep only significant correlations
  filter(adjP < 0.05) %>%
  #collapse by metabolite
  group_by(Met) %>%
  mutate(Concordance = paste(unique(Cohort), collapse = "; "),
         Concordance = ifelse(grepl(";", Concordance), "Both visits", Concordance)) %>%
  data.frame() %>%
  select(Met, Concordance) %>%
  unique() %>%
  inner_join(sumDF, ., by = "Met")

#visualise distributions of "interesting" MOI
#looking at top 12 most significant by reduced model
top12 <- NGATRedHits$Metabolite[1:12]

#prepare to visualise
logMatLim %>%
  data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "MetaboliteBin", values_to = "logPQN") %>%
  #only stick to metabolites of interest
  filter(MetaboliteBin %in% top12) %>%
  #add cohort in
  inner_join(., met, by = "Sample") %>%
  mutate(Group = factor(Group, levels = c("Control",
                                          "Recruitment",
                                          "Follow up"))) %>%
  mutate(BinLab = map_vec(MetaboliteBin, fixMet),
         Met = map_vec(BinLab, function(m) {
           ifelse(grepl("^Unknown", m), m, gsub("-[0-9]+$", "", m, 
                                                perl = TRUE))})
    ) %>%
  ggplot(aes(x = Group, y = logPQN)) +
  geom_boxplot(aes(fill = Group)) +
  scale_fill_manual(values = vanPal) +
  facet_wrap(~Met, scales = "free_y") + 
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("output/nmr/top12SigredModOut_MetaboliteInt_FullCohort.png")
  
#visualisation per metabolite > bin
#this will be metabolites that are significant in both models AND
#significantly correlated at both timepoints
consistent <- corOut %>%
  group_by(MetaboliteBin) %>%
  mutate(Concordance = paste(unique(Cohort), collapse = ", ")) %>%
  select(MetaboliteBin, Concordance) %>%
  filter(grepl(",", Concordance)) %>%
  unique()
#go back to metabolite intensities
logMatLim %>%
  data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "MetaboliteBin", values_to = "logPQN") %>%
  #keep only consistent metabolite bins
  filter(MetaboliteBin %in% consistent$MetaboliteBin) %>%
  #add metabolite groups
  mutate(BinLab = map_vec(MetaboliteBin, fixMet),
         Met = map_vec(BinLab, function(m) {
           ifelse(grepl("^Unknown", m), m, gsub("-[0-9]+$", "", m, perl = TRUE))
         })) %>%
  #bring in grouping metadata
  inner_join(., met[,c("Sample", "Group")]) %>%
  #plot
  ggplot(aes(x = Group, y = logPQN)) +
  geom_boxplot(aes(fill = Group), alpha = 0.75) + 
  #geom_point(aes(colour = Group), position = "jitter", size = 1) +
  scale_colour_manual(values = vanPal) +
  scale_fill_manual(values = vanPal) +
  facet_wrap(~Met, scales = "free_y", ncol = 4) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")
#save
ggsave("output/nmr/mostConsistent_MetaboliteInt_FullCohort.png",
       height = 20, unit = "cm")

#there's a further break down of metabolites that seem particularly interesting
#the below consistently appear in the 3 different analyses when ordering by adjP
#or absolute effect - and all are consistently associated with NGAT at both
#timepoints
moi2 <- c("Formate_5", "Unknown_153", "Unknown_23",
          "O_phosphocholine_120", "Lactate_48", 
          "Glutamine_154", "Glucose_116", "Glucose_101",
          "Glucose_102", "Glucose_103", "Glucose_104", "Acetate_197")
#plot
moiBox <- logMatLim %>%
  data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "MetaboliteBin", values_to = "logPQN") %>%
  #keep only consistent metabolite bins
  filter(MetaboliteBin %in% moi2) %>%
  #add metabolite groups
  mutate(BinLab = map_vec(MetaboliteBin, fixMet),
         Met = map_vec(BinLab, function(m) {
           ifelse(grepl("^Unknown", m), m, gsub("-[0-9]+$", "", m, perl = TRUE))
         })) %>%
  #bring in grouping metadata
  inner_join(., met[,c("Sample", "Group")]) %>%
  #plot
  ggplot(aes(x = Group, y = logPQN)) +
  geom_boxplot(aes(fill = Group), alpha = 0.75) + 
  #geom_point(aes(colour = Group), position = "jitter", size = 1) +
  scale_colour_manual(values = vanPal) +
  scale_fill_manual(values = vanPal) +
  facet_wrap(~Met, scales = "free_y", nrow = 2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")
#save
ggsave("output/nmr/mostInteresting_MetaboliteInt_FullCohort.png")

##### Output ####
#log(PQN) intensities
write.csv(logMatLim, "processed/nmr/logPQN_NMR.csv",
            quote = F, row.names = T)

#LMM results
bind_rows(fullOutDF, redOutDF) %>%
  write.csv(., "output/nmr/lme4_allresults.csv",
            row.names = F, quote = F)

#correlation results
write.csv(corOut, "output/nmr/NGAT_Kendall_PerVisit.csv",
          row.names = F, quote = F)

#summary of metabolites (summing bins)
write.csv(sumDF, "output/nmr/Metabolite-NGAT_SigMets_Summary.csv",
          row.names = F, quote = F)

#options for manuscript
#option 1
pFullLMMslim
ggsave("output/figures/lmm-moi-slim.png",
       width = 12, height = 15, unit = "cm")
ggsave("output/figures/lmm-moi-slim.pdf",
       width = 12, height = 15, unit = "cm")
ggsave("output/figures/lmm-moi-slim.jpg",
       width = 12, height = 15, unit = "cm")
save(pFullLMMslim, file = "output/figures/lmm-moi-slim.rds")

#record samples used in analysis
reg <- read.csv("processed/Sample_AnalysisRegister.csv")
met %>%
  mutate(Metabolome = "Present") %>%
  select(Sample, Patient.ID, Metabolome) %>%
  full_join(reg, ., by = c("Sample", "Patient.ID")) %>%
  mutate(Microbiome = ifelse(is.na(Microbiome), "Absent", Microbiome),
         Metabolome = ifelse(is.na(Metabolome), "Absent", Metabolome)) %>%
  write.csv(., "processed/Sample_AnalysisRegister.csv",
            row.names = F, quote = F)

#Session information
writeLines(capture.output(sessionInfo()),
           "output/nmr/NMRAnalysis_sessionInfo.txt")
