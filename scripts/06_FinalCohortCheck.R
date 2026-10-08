#14th September 2026
#Lauren Mee
#Final bits

##### Set up #####
#set seed
set.seed(1248)

#libraries
libs <- c("tidyverse") #general data manipulation

for (pkg in libs) {
  library(pkg, character.only = T)
}

#output directory architecture
dirs <- c("output/")

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d)
  }
}

##### Functions ####
#convert pvalues to a significance notation
getSig <- function(pval) {
  if (pval < 0.05) {
    dec <- "*"
    if (pval < 0.01) {
      dec <- "**"
      if (pval < 0.001) {
        dec <- "***"
      }
    }
  } else {
    dec <- "NS"
  }
  return(dec)
}

##### Load data #####
#full metadata, factorised as elsewhere in the pipeline
met <- read.csv("input/SampleMetadata.csv") %>%
  mutate(Group = factor(Group, levels = c("Control",
                                          "Recruitment",
                                          "Follow up")),
         Patient.ID = factor(Patient.ID))

#cohort "registry"
#note - only samples that are included in at least one arm of the analyses
#are in here - no need for further filtering
reg <- read.csv("processed/Sample_AnalysisRegister.csv") %>%
  mutate(Patient.ID = factor(Patient.ID))

##### Build the cohort #####
#any sample used in any of the three arms of analysis (106 samples originally)
metCo <- met %>%
  filter(Sample %in% reg$Sample)
  #only 2 samples were lost entirely from the analysis

#restrict to GSM patients: controls are single-visit by design and cannot
#contribute to a between-visit comparison
metGSM <- metCo %>%
  filter(Status == "Atrophy") %>%
  droplevels()

#Sort both subsets by Patient.ID and check they match before testing.
rec <- metGSM %>% 
  filter(Group == "Recruitment") %>% 
  arrange(Patient.ID)
fup <- metGSM %>% 
  filter(Group == "Follow up") %>% 
  arrange(Patient.ID)

#order
paired <- intersect(rec$Patient.ID, fup$Patient.ID)
rec <- rec %>% 
  filter(Patient.ID %in% paired)
fup <- fup %>% 
  filter(Patient.ID %in% paired)
#check
identical(as.character(rec$Patient.ID), as.character(fup$Patient.ID))

##### Tests #####
#prepare to store output
stats <- vector(mode = "list")

#NGAT
#ties exist - include continuity correction with wilcox, as in script 02
t0 <- wilcox.test(rec$NGAT, fup$NGAT,
                  exact = FALSE, continuity = TRUE, paired = TRUE)
stats[[1]] <- data.frame(group1 = "Recruitment",
                         group2 = "Follow up",
                         .y. = "NGAT",
                         pairs = length(paired),
                         V = t0$statistic,
                         p = t0$p.value,
                         y.position = 16) %>%
  mutate(p.label = unlist(map(p, getSig)))
names(stats)[1] <- "NGAT_RecruitmentVFollowUp_WilcoxonSignedRankTest_Paired_ContinuityCorrection"

#DIVA
#ties exist
t0 <- wilcox.test(rec$DIVA, fup$DIVA,
                  exact = FALSE, continuity = TRUE, paired = TRUE)
stats[[2]] <- data.frame(group1 = "Recruitment",
                         group2 = "Follow up",
                         .y. = "DIVA",
                         pairs = length(paired),
                         V = t0$statistic,
                         p = t0$p.value,
                         y.position = 95) %>%
  mutate(p.label = unlist(map(p, getSig)))
names(stats)[2] <- "DIVA_RecruitmentVFollowUp_WilcoxonSignedRankTest_Paired_ContinuityCorrection"

##### DIVA v NGAT correlation ####
#as NGAT and DIVA both measure atrophy severity, they should correlate. Tested 
#per timepoint as each patient contributes two non-independent samples.
#ties exist in both scores so kendall is used, as has been case throughout pipeline.
metaCors <- data.frame(Cohort = character(),
                   Var1 = character(),
                   Var2 = character(),
                   Statistic = character(),
                   StatValue = numeric(),
                   tau = numeric(),
                   Exact = logical(),
                   p = numeric())

#Recruitment
test <- cor.test(rec$DIVA, rec$NGAT,
                 method = "kendall")
metaCors[1,] <- c("Recruitment", "DIVA", "NGAT", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))
#Follow up
test <- cor.test(fup$DIVA, fup$NGAT,
                 method = "kendall")
metaCors[2,] <- c("Follow up", "DIVA", "NGAT", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))

#correct
metaCors <- metaCors %>%
  mutate(p = as.numeric(p)) %>%
  mutate(p.adj = p.adjust(p, method = "fdr"),
         p.label = map_vec(p.adj, getSig))

stats[[3]] <- metaCors
names(stats)[3] <- "DIVAvNGAT_PerTimepoint_KendallRankCorrelation"

##### Age v severity scores correlation ####
#in microbiome/metabolome subcohorts age and NGAT associate, DIVA and age do not.
ageCors <- data.frame(Cohort = character(),
                   Var1 = character(),
                   Var2 = character(),
                   Statistic = character(),
                   StatValue = numeric(),
                   tau = numeric(),
                   Exact = logical(),
                   p = numeric())

#Recruitment
test <- cor.test(rec$Age, rec$NGAT,
                 method = "kendall")
ageCors[1,] <- c("Recruitment", "Age", "NGAT", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", 
                  unname(test$p.value))
#Follow up
test <- cor.test(fup$Age, fup$NGAT,
                 method = "kendall")
ageCors[2,] <- c("Follow up", "Age", "NGAT", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", 
                  unname(test$p.value))
#Recruitment
test <- cor.test(rec$Age, rec$DIVA,
                 method = "kendall")
ageCors[3,] <- c("Recruitment", "Age", "DIVA", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))
#Follow up
test <- cor.test(fup$Age, fup$DIVA,
                 method = "kendall")
ageCors[4,] <- c("Follow up", "Age", "DIVA", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))

#correct within age-severity correlation tests
ageCors <- ageCors %>%
  mutate(p = as.numeric(p)) %>%
  mutate(p.adj = p.adjust(p, method = "fdr"),
         p.label = map_vec(p.adj, getSig))

stats[[4]] <- ageCors
names(stats)[4] <- "AgevSeverityScores_PerTimepoint_KendallRankCorrelation"

##### Demography #####
#cohort demographies
healthy <- metCo %>%
  filter(Group == "Control") %>%
  select(Sample) %>% pull()

#add registry to GSM cohort
metGSM <- metGSM %>%
  inner_join(., reg, by = c("Sample", "Patient.ID"))

#microbiome
microGSM <- metGSM %>%
  filter(Microbiome == "Present") %>%
  select(Sample) %>% pull()

#metabolome
metaGSM <- metGSM %>%
  filter(Metabolome == "Present") %>%
  select(Sample) %>% pull()

#network
netGSM <- metGSM %>%
  filter(Network == "Present") %>%
  select(Sample) %>% pull()

demoDF <- data.frame(
  Cohort = c("Healthy", "GSM: Microbiome Subset", "GSM: NMR Subset",
             "GSM: Network Subset"),
  sampleN = c(length(healthy), length(microGSM), length(metaGSM), length(netGSM)),
  PatientN = c(length(unique(met$Patient.ID[met$Sample %in% healthy])),
               length(unique(met$Patient.ID[met$Sample %in% microGSM])),
               length(unique(met$Patient.ID[met$Sample %in% metaGSM])),
               length(unique(met$Patient.ID[met$Sample %in% netGSM]))),
  Age = c(round(mean(met$Age[met$Sample %in% healthy]), 2),
          round(mean(met$Age[met$Sample %in% microGSM &
                               met$Group == "Recruitment"]), 2),
          round(mean(met$Age[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"]), 2),
          round(mean(met$Age[met$Sample %in% netGSM &
                               met$Group == "Recruitment"]), 2)),
  SDAge = c(round(sd(met$Age[met$Sample %in% healthy]), 2),
            round(sd(met$Age[met$Sample %in% microGSM &
                               met$Group == "Recruitment"]), 2),
            round(sd(met$Age[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"]), 2),
            round(sd(met$Age[met$Sample %in% netGSM &
                               met$Group == "Recruitment"]), 2)),
  BMI = c(round(mean(met$BMI[met$Sample %in% healthy], na.rm = T), 2),
          round(mean(met$BMI[met$Sample %in% microGSM &
                               met$Group == "Recruitment"]), 2),
          round(mean(met$BMI[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"]), 2),
          round(mean(met$BMI[met$Sample %in% netGSM &
                               met$Group == "Recruitment"]), 2)),
  SDBMI = c(round(sd(met$BMI[met$Sample %in% healthy], na.rm = T), 2),
            round(sd(met$BMI[met$Sample %in% microGSM &
                               met$Group == "Recruitment"]), 2),
            round(sd(met$BMI[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"]), 2),
            round(sd(met$BMI[met$Sample %in% netGSM &
                               met$Group == "Recruitment"]), 2)),
  NGAT_Recruitment = c(round(mean(met$NGAT[met$Sample %in% healthy]), 2),
          round(mean(met$NGAT[met$Sample %in% microGSM &
                               met$Group == "Recruitment"]), 2),
          round(mean(met$NGAT[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"]), 2),
          round(mean(met$NGAT[met$Sample %in% netGSM &
                               met$Group == "Recruitment"]), 2)),
  SDNGAT_Recruitment = c(round(sd(met$NGAT[met$Sample %in% healthy]), 2),
            round(sd(met$NGAT[met$Sample %in% microGSM &
                               met$Group == "Recruitment"]), 2),
            round(sd(met$NGAT[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"]), 2),
            round(sd(met$NGAT[met$Sample %in% netGSM &
                               met$Group == "Recruitment"]), 2)),
  NGAT_FollowUp = c(NA,
          round(mean(met$NGAT[met$Sample %in% microGSM &
                               met$Group == "Follow up"]), 2),
          round(mean(met$NGAT[met$Sample %in% metaGSM &
                               met$Group == "Follow up"]), 2),
          round(mean(met$NGAT[met$Sample %in% netGSM &
                               met$Group == "Follow up"]), 2)),
  SDNGAT_FollowUp = c(NA,
            round(sd(met$NGAT[met$Sample %in% microGSM &
                               met$Group == "Follow up"]), 2),
            round(sd(met$NGAT[met$Sample %in% metaGSM &
                               met$Group == "Follow up"]), 2),
            round(sd(met$NGAT[met$Sample %in% netGSM &
                               met$Group == "Follow up"]), 2)),
  DIVA_Recruitment = c(round(mean(met$DIVA[met$Sample %in% healthy]), 2),
                       round(mean(met$DIVA[met$Sample %in% microGSM &
                                             met$Group == "Recruitment"]), 2),
                       round(mean(met$DIVA[met$Sample %in% metaGSM &
                                             met$Group == "Recruitment"]), 2),
                       round(mean(met$DIVA[met$Sample %in% netGSM &
                                             met$Group == "Recruitment"]), 2)),
  SDDIVA_Recruitment = c(round(sd(met$DIVA[met$Sample %in% healthy]), 2),
                         round(sd(met$DIVA[met$Sample %in% microGSM &
                                             met$Group == "Recruitment"]), 2),
                         round(sd(met$DIVA[met$Sample %in% metaGSM &
                                             met$Group == "Recruitment"]), 2),
                         round(sd(met$DIVA[met$Sample %in% netGSM &
                                             met$Group == "Recruitment"]), 2)),
  DIVA_FollowUp = c(NA,
                    round(mean(met$DIVA[met$Sample %in% microGSM &
                                         met$Group == "Follow up"]), 2),
                    round(mean(met$DIVA[met$Sample %in% metaGSM &
                                         met$Group == "Follow up"]), 2),
                    round(mean(met$DIVA[met$Sample %in% netGSM &
                                         met$Group == "Follow up"]), 2)),
  SDDIVA_FollowUp = c(NA,
                      round(sd(met$DIVA[met$Sample %in% microGSM &
                                         met$Group == "Follow up"]), 2),
                      round(sd(met$DIVA[met$Sample %in% metaGSM &
                                         met$Group == "Follow up"]), 2),
                      round(sd(met$DIVA[met$Sample %in% netGSM &
                                         met$Group == "Follow up"]), 2)),
  Ancestry_WhiteBritish = c(
    round(sum((met$Ethnicity[met$Sample %in% healthy] == "White British") /
                                  length(healthy)) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% microGSM &
                             met$Group == "Recruitment"] == "White British") /
          length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% metaGSM &
                           met$Group == "Recruitment"] == "White British") /
          length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% netGSM &
                             met$Group == "Recruitment"] == "White British") /
          length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)),
  Ancestry_AsianPakistani = c(
    round(sum((met$Ethnicity[met$Sample %in% healthy] == "Asian Pakistani") /
                length(healthy)) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% microGSM &
                               met$Group == "Recruitment"] == "Asian Pakistani") /
                length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"] == "Asian Pakistani") /
                length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% netGSM &
                               met$Group == "Recruitment"] == "Asian Pakistani") /
                length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)),
  Ancestry_OtherWhite = c(
    round(sum((met$Ethnicity[met$Sample %in% healthy] == "Other white") /
                length(healthy)) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% microGSM &
                               met$Group == "Recruitment"] == "Other white") /
                length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"] == "Other white") /
                length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% netGSM &
                               met$Group == "Recruitment"] == "Other white") /
                length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)),
  Ancestry_MixedWhiteHispanic = c(
    round(sum((met$Ethnicity[met$Sample %in% healthy] == "Mixed White/Hispanic") /
                length(healthy)) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% microGSM &
                               met$Group == "Recruitment"] == "Mixed White/Hispanic") /
                length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"] == "Mixed White/Hispanic") /
                length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% netGSM &
                               met$Group == "Recruitment"] == "Mixed White/Hispanic") /
                length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)),
  Ancestry_BritishAsian = c(
    round(sum((met$Ethnicity[met$Sample %in% healthy] == "British Asian") /
                length(healthy)) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% microGSM &
                               met$Group == "Recruitment"] == "British Asian") /
                length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"] == "British Asian") /
                length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% netGSM &
                               met$Group == "Recruitment"] == "British Asian") /
                length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)),
  Ancestry_MixedWhiteBlack = c(
    round(sum((met$Ethnicity[met$Sample %in% healthy] == "Mixed White/Black") /
                length(healthy)) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% microGSM &
                               met$Group == "Recruitment"] == "Mixed White/Black") /
                length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"] == "Mixed White/Black") /
                length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% netGSM &
                               met$Group == "Recruitment"] == "Mixed White/Black") /
                length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)),
  Ancestry_African = c(
    round(sum((met$Ethnicity[met$Sample %in% healthy] == "African") /
                length(healthy)) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% microGSM &
                               met$Group == "Recruitment"] == "African") /
                length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"] == "African") /
                length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% netGSM &
                               met$Group == "Recruitment"] == "African") /
                length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)),
  Ancestry_MixedWhiteChinese = c(
    round(sum((met$Ethnicity[met$Sample %in% healthy] == "Mixed White/Chinese") /
                length(healthy)) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% microGSM &
                               met$Group == "Recruitment"] == "Mixed White/Chinese") /
                length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"] == "Mixed White/Chinese") /
                length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Ethnicity[met$Sample %in% netGSM &
                               met$Group == "Recruitment"] == "Mixed White/Chinese") /
                length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)),
  SmokingStatus_NonSmoker = c(
    round(sum((met$Smoking.status[met$Sample %in% healthy] == "Non-smoker") /
                length(healthy)) * 100, 2),
    round(sum((met$Smoking.status[met$Sample %in% microGSM &
                               met$Group == "Recruitment"] == "Non-smoker") /
                length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Smoking.status[met$Sample %in% metaGSM &
                               met$Group == "Recruitment"] == "Non-smoker") /
                length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Smoking.status[met$Sample %in% netGSM &
                               met$Group == "Recruitment"] == "Non-smoker") /
                length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)),
  MenopausalStatus_Menopausal = c(
    round(sum((met$Menopausal.status[met$Sample %in% healthy] == "Menopausal") /
                length(healthy)) * 100, 2),
    round(sum((met$Menopausal.status[met$Sample %in% microGSM &
                                    met$Group == "Recruitment"] == "Menopausal") /
                length(unique(met$Patient.ID[met$Sample %in% microGSM]))) * 100, 2),
    round(sum((met$Menopausal.status[met$Sample %in% metaGSM &
                                    met$Group == "Recruitment"] == "Menopausal") /
                length(unique(met$Patient.ID[met$Sample %in% metaGSM]))) * 100, 2),
    round(sum((met$Menopausal.status[met$Sample %in% netGSM &
                                    met$Group == "Recruitment"] == "Menopausal") /
                length(unique(met$Patient.ID[met$Sample %in% netGSM]))) * 100, 2)))
#omg. worst code I've written in years. I cannot be bothered to streamline, its being left. 
#sorry future, embarrassed me.


##### Output ####
#stats
writeLines(capture.output(print(stats)),
           "output/BetweenVisits_UsedCohort_Summary.txt")
#demography summary
write.csv(demoDF, "output/VAN_DemographyBreakdown.csv",
          row.names = F, quote = F)
#registration matched with input metadata
met %>%
  left_join(., reg, by = c("Sample", "Patient.ID")) %>%
  mutate(Microbiome = ifelse(is.na(Microbiome), "Absent", Microbiome),
         Metabolome = ifelse(is.na(Metabolome), "Absent", Metabolome),
         Network = ifelse(is.na(Network), "Absent", Network)) %>%
  write.csv(., "output/EndofAnalysis_Metadata.csv",
            quote = F, row.names = F)