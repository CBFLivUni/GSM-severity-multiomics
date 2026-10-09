#Changed into script 12th December
#Original Early 2024
#Lauren Mee

##### Set up #####
#set seed
set.seed(1049)

#libraries
libs <- c("vegan", "phyloseq", #microbiome processing
          "maaslin3", "lme4", "lmerTest", #microbiome abundance and prevalence modelling
          "variancePartition", #collinearity plot (cca)
          "ggpubr", "tidyverse") #general plots etc

for (pkg in libs) {
  library(pkg, character.only = T)
}

#load data
#metadata and factorise appropriately
met <- read.csv("output/microbiome/MicroPreProcessing_Data.csv") %>%
  mutate(Group = factor(Group, levels = c("Control",
                                                  "Recruitment",
                                                  "Follow up")),
         Patient.ID = factor(Patient.ID))
#phyloseq objects
load("processed/microbiome/PhyloSeqObjs.rds")
#genus-level based dataframe
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

#output directory architecture
dirs <- c("processed/", "processed/microbiome/",
          "output/", "output/microbiome/", "output/microbiome/community/",
          "output/microbiome/alpha/", "output/microbiome/beta/",
          "output/microbiome/clinicalMetric/", "output/microbiome/maaslin3/")

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
#genus palette
load("processed/microbiome/AbundantGeneraPalette.rds")

#ggplot themes
theme_set(theme_bw(base_size = 13))
theme_update(
  strip.background = element_rect(fill = "lightgrey", colour = "white"),
  strip.text = element_text(colour = "black", face = "bold")
)

#for plots require more than default # of colours
palFun <- colorRampPalette(vanPal)

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

#determine best k (number of dimensions) for nmds plot with 
#input distance matrix (dist)
chooseK <- function(dist, 
                    type = "k", #other options: table, figure
                    method = "bray" 
) {
  strs <- vector(length = 10)
  for (i in 1:10) {
    mds <- metaMDS(dist, distance = method, autotransform = F, k = i)
    strs[i] <- mds$stress
    out <- data.frame(NoDimensions = 1:10, Stress = strs)
  }
  #if stress isn't great, abandon
  if (min(out$Stress) > 0.2) {
    stop("Stress is quite high. Perhaps reassess data.")
  }
  #determine best k using elbow
  #normalise k (x) and stress (y)
  xNorm <- (out$NoDimensions - min(out$NoDimensions)) / 
    (max(out$NoDimensions) - min(out$NoDimensions))
  yNorm <- (out$Stress - min(out$Stress)) / (max(out$Stress) - min(out$Stress))
  #produce line from first and last points
  p1 <- c(xNorm[1], yNorm[1])
  p2 <- c(xNorm[length(xNorm)], yNorm[length(yNorm)])
  #compute how far each point sits from said line
  distances <- map2_vec(xNorm, yNorm, function(x, y) {
    p3 <- c(x, y)
    abs((p2[2]-p1[2])*p3[1] - (p2[1]-p1[1])*p3[2] + p2[1]*p1[2] - p2[2]*p1[1]) /
      sqrt((p2[2]-p1[2])^2 + (p2[1]-p1[1])^2)
  })
  #store largest point-to-line distance as elbow point
  kDim <- which.max(distances)
  if (type == "k") {
    return(kDim)
  }
  if (type == "table") {
    return(out)
  }
  if (type == "plot") {
    p <- out %>%
      mutate(PointCol = ifelse(Stress > 0.2, "Bad", "OK")) %>%
        ggplot(aes(y = Stress, x = NoDimensions)) +
      geom_hline(yintercept = out$Stress[kDim], linetype = "dashed",
                 alpha = 0.5) +
      geom_vline(xintercept = out$NoDimensions[kDim], linetype = "dashed",
                 alpha = 0.5) +
      geom_point(aes(colour = PointCol)) +
      scale_colour_manual(values = c("Bad" = "red3",
                                     "OK" = "black")) + 
      labs(y = "Stress", 
           x = "Number of Dimensions",
           colour = "",
           subtitle = "Stress versus Dimensionality") +
      theme_bw(base_size = 13) +
      guides(colour = "none") +
      geom_point(data = out[kDim,], aes(x = NoDimensions, y = Stress),
                 shape = 21, colour = "darkgrey", fill = NA, 
                 size = 6, stroke = 1.5) +
      scale_x_continuous(breaks = 1:10, labels = 1:10)
    return(p)
  }
}

##### Alpha diversity ####
#rarefaction
#extract ASV counts from genus collapsed object
cnts <- data.frame(otu_table(phyGen))

#set number of iterations
n <- 10000 

#initialise empty matrix for shannon calculations
rMat <- matrix(0, ncol = 1, nrow = nrow(cnts))
#initialise empty vector fo rrichness calculations
rich <- rep(0, nrow(cnts))

#prepare rarefaction level
nLmin <- min(rowSums(cnts))

#shannon for richness and evenness, 
#species richness for ... richness.
for (i in 1:n){
  if (i %% 500 == 0) {
    print(i)
  }
  #rarefy count table by the lowest total reads in one sample
  #I am using raw counts but vegan keeps warning that 
  #there are no singletons - I will suppress this for now
  rCnts <- suppressWarnings(rrarefy(cnts, nLmin))
  #compute alpha diversity metric (shannons)
  iMat <- as.matrix(diversity(rCnts, index = "shannon")) 
  #accumulate counts
  rMat <- rMat + iMat
  #record richness
  rich <- rich + specnumber(rCnts)
}

#produce dataframe
alpha <- data.frame(Sample  = rownames(cnts),
                    Shannon = rMat / n,
                    Richness = rich / n) %>%
  #add metadata
  inner_join(., met, by = "Sample")

#prepare space to store stats throughout analysis
stats <- vector(mode = "list")
#add test statistics
test <- wilcox.test(alpha$Shannon[alpha$Group == "Recruitment"],
                    alpha$Shannon[alpha$Group == "Follow up"], 
                    paired = TRUE)
stats[[1]] <- data.frame(group1 = "Recruitment",
                    group2 = "Follow up",
                    .y. = "Shannon",
                    y.position = 3,
                    p = round(test$p.value, 3),
                    V = test$statistic) %>%
  mutate(p.label = unlist(map(p, getSig)))

#store stats
names(stats)[1] <- "Alpha_BeforeVAfter_WilcoxonSignedRankSum_Paired"

#is shannon's index associated with either clinical score of atrophy? 
#(linear modelling) independently per score
alphaGSM <- alpha %>%
  filter(!Group == "Control") %>%
  mutate(Group = factor(Group, levels = c("Recruitment", "Follow up")))
#full  model
fam <- lmer(Shannon ~ NGAT + DIVA + Age + Group + BMI + (1 | Patient.ID), 
            data = alphaGSM)
famOut <- summary(fam)
#save
capture.output(famOut,
               file = "output/microbiome/alpha/Shannon_LMM_FullModel.txt")
#reduced NGAT
ram <- lmer(Shannon ~ NGAT + DIVA + Group + (1 | Patient.ID), 
           data = alphaGSM)
ramOut <- summary(ram)
#save
capture.output(ramOut,
               file = "output/microbiome/alpha/Shannon_LMM_ReducedModel.txt")
#richness
#add test statistics
test <- wilcox.test(alpha$Richness[alpha$Group == "Recruitment"],
                    alpha$Richness[alpha$Group == "Follow up"], 
                    paired = TRUE)
stats[[2]] <- data.frame(group1 = "Recruitment",
                         group2 = "Follow up",
                         .y. = "Richness",
                         y.position = 57,
                         p = round(test$p.value, 3),
                         V = test$statistic) %>%
  mutate(p.label = unlist(map(p, getSig)))

#store stats
names(stats)[2] <- "Richness_BeforeVAfter_WilcoxonSignedRankSum_Paired"

#is richness associated with either clinical score of atrophy? 
#(linear modelling) independently per score
#full  model
famR <- lmer(Richness ~ NGAT + DIVA + Age + Group + BMI + (1 | Patient.ID), 
            data = alphaGSM)
famOutR <- summary(famR)

#save
capture.output(famOutR,
               file = "output/microbiome/alpha/Richness_LMM_FullModel.txt")
#reduced NGAT
ramR <- lmer(Richness ~ NGAT + DIVA + Group + (1 | Patient.ID), 
            data = alphaGSM)
ramOutR <- summary(ramR)
#save
capture.output(ramOutR,
               file = "output/microbiome/alpha/Richness_LMM_ReducedModel.txt")

#visualise
#Standardise predictors. Raw coefficients cannot share an axis and cannot be
#not comparable "per point". Standardising makes every estimate "SDs of the 
#metric per SD of the predictor". Factors are left alone so Group stays 
#interpretable as a between-visit difference.
mods <- list(Full = "~ NGAT + DIVA + Age + Group + BMI + (1 | Patient.ID)",
             Reduced = "~ NGAT + DIVA + Group + (1 | Patient.ID)")

coefTab <- map_dfr(c("Shannon", "Richness"), function(m) {
  map_dfr(names(mods), function(mn) {
    #subset and rescale data
    d <- alphaGSM %>% 
      #convert continuous variables to z scores 
      mutate(y = as.numeric(scale(.data[[m]])), 
             NGAT = as.numeric(scale(NGAT)),
             DIVA = as.numeric(scale(DIVA)), 
             Age = as.numeric(scale(Age)), 
             BMI = as.numeric(scale(BMI)))
    #run model
    fit <- lmer(as.formula(paste("y", mods[[mn]])), data = d)
    #store beta coefficients
    co  <- summary(fit)$coefficients
    #calculate confidence intervals
    ci  <- confint(fit, method = "Wald")     
    #output 
    tibble(Metric = m, Model = mn, Term = rownames(co),
           Estimate = co[, 1], P = co[, ncol(co)],
           Lower = ci[rownames(co), 1], Upper = ci[rownames(co), 2])
  })
})

#BMI stays in the model but is not plotted: BMI is associated with 
#neither score nor outcome.
keep <- c(NGAT = "NGAT", DIVA = "DIVA", Age = "Age",
          `GroupFollow up` = "Follow up vs\nRecruitment")

#plot
alphaPlotDF <- coefTab %>%
  filter(Term %in% names(keep)) %>%
  mutate(Term = factor(keep[Term], levels = rev(unname(keep))),
         Metric = factor(Metric, levels = c("Shannon", "Richness")),
         Model = factor(Model, levels = c("Reduced", "Full")),
         Sig = ifelse(P < 0.05, "p < 0.05", "n.s."),
         Sig = factor(Sig, levels = c("p < 0.05", "n.s.")))

alphaPlot <- ggplot(alphaPlotDF, 
                    aes(x = Estimate, y = Term, colour = Model, shape = Sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_linerange(aes(xmin = Lower, xmax = Upper),
                 position = position_dodge(width = 0.55), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.55), 
             size = 2.8, fill = "white") +
  facet_wrap(~ Metric) +
  scale_colour_manual(values = c(vanPal[5], vanPal[7])) +
  scale_shape_manual(values = c("p < 0.05" = 16, "n.s." = 21)) +
  labs(x = "Standardised effect\n(SD units)",
       y = NULL, colour = NULL, shape = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2, ncol = 1))
ggsave("output/microbiome/alpha/LMM_AlphaDiversityMetrics_vs_ClinicalMetrics.png")


##### Beta diversity: Bray-Curtis #####
#also requires rarefaction
#considering two patient timepoints in separation, often
#pre-existing or unchanged treatment, "Follow up" samples cannot be
#treated as a clean post-intervention timepoint, and pooling both
#timepoints per patient would pseudoreplicate the beta diversity model
grpz <- c("Recruitment", "Follow up")
samps <- map(grpz, function(timeP) {
  return(met$Sample[met$Group == timeP])
})

#prepare cnts to convert into distance matrices
betaCnts <- map(samps, function(s) {
  return(cnts[rownames(cnts) %in% s,])
})

#note - using bray-curtis here, but will also be using robust.aitchison
#later. Bray-Curtis is sensitive to single taxa dominance - which
#is exactly the case with vaginal microbiome 16S data (Lactobacillus)
#same limit as alpha rarefaction
distsBC <- map(betaCnts, function(b) {
  return(avgdist(b, dmethod = "bray", iterations = 10000,
                 sample = nLmin))
})

#nmds
#plot stress v dimensionality plots
map2(distsBC, grpz, function(d, g) {
  #plot
  chooseK(d, type = "plot") +
    labs(title = paste0(g, ": Bray-Curtis"))
  #save
  lab <- gsub(" ", "-", g)
  ggsave(paste0("output/microbiome/beta/NMDS_Bray-Curtis_", lab,
                "_StressVsDimensions.png"))
})

#decide on k parameter
kListBC <- map(distsBC, chooseK)

#nmds
nmdsListBC <- map2(distsBC, kListBC, function(d, k){
  metaMDS(d, k = k, trymax = 100, trace = F, distance = "bray")
}) 

#plot stressplot
map2(nmdsListBC, grpz, function(n, g) {
  lab <- gsub(" ", "-", g)
  png(paste0("output/microbiome/beta/NMDS_Bray-Curtis_", lab, "_Stressplot.png"))
    stressplot(n)
  dev.off()
})

#produce multiple plots
nmdsBCDFs <- map(nmdsListBC, function(n) {
  out <- scores(n) %>%
    data.frame() %>%
    select(NMDS1, NMDS2) %>%
    rownames_to_column(var = "Sample") %>%
    #add sample metadata 
    inner_join(., met, by = "Sample") %>%
    rename("NMDS1.BC" = "NMDS1",
           "NMDS2.BC" = "NMDS2")
})

#store variables of interest for simple plots
vois <- c("NGAT", "DIVA",
          "Age", "BMI",
          "Ethnicity")

nmdsPlotListBC <- map2(nmdsBCDFs, grpz, function(n, g) {
  #prepare to store simple explorations
  plots <- vector(mode = "list", length = length(vois))
  #figure out how many categorical levels I have to plot
  values <- vector()
  for (i in 1:length(vois)) {
    if (!is.numeric(met[[vois[i]]])) {
      values <- c(values, unique(met[[vois[i]]]))
    }
  }
  #iterate through variables and plot and save
  for (i in 1:length(vois)) {
    plotDF <- n %>%
      mutate(Group = n[[vois[i]]])
    lab <- gsub("\\.", " ", vois[[i]])
    #store base plot
    plots[[i]] <- plotDF %>%
      ggplot(aes(x = NMDS1.BC, y = NMDS2.BC)) +
      geom_point(aes(colour = Group),
                 size = 3, alpha = 0.75) +
      labs(subtitle = lab,
           colour = "",
           x = "NMDS1",
           y = "NMDS2") +
      theme(plot.subtitle = element_text(face = "bold"))
    #figure out colour schemes
    if (!is.numeric(plotDF$Group)) {
      plots[[i]] <- plots[[i]] + 
        scale_colour_manual(values = vanPal[c(4,9,1,8,7)]) 
    } else {
      midP <- ((max(plotDF$Group) - min(plotDF$Group)) / 2) + min(plotDF$Group)
      plots[[i]] <- plots[[i]] +
        scale_colour_gradient2(low = vanPal[[6]],
                               mid = vanPal[[7]],
                               high = "#050808",
                               midpoint = midP)
    }
  }
  p <- ggarrange(plotlist = plots, ncol = 2, nrow = 3)
  p <- annotate_figure(p, top = text_grob(g, size = 15, face = "bold"))
  #save
  title <- gsub(" ", "-", g)
  ggsave(paste0("output/microbiome/beta/NMDS_Bray-Curtis_", title,
                "_MetaVariables.pdf"), scale = 2)
  ggsave(paste0("output/microbiome/beta/NMDS_Bray-Curtis_", title,
                "_MetaVariables.png"), scale = 2)
  return(p)
})

#more indepth plots
#its likely lactobacillus being so ubiquitous and at times
#close to 100% abundant has a large impact on the NMDS
#also = likely to have a considerable affect on bray-curtis distance
#as mentioned above
lacto <- genCnts %>%
  data.frame() %>%
  #keep only Lactobacillus
  filter(Genus == "Lactobacillus")
#factorise lactobacillus level
lacto$LactoRelABin[lacto$relA < 0.2] <- "0-19%"
lacto$LactoRelABin[lacto$relA >= 0.8] <- ">= 80%"
lacto$LactoRelABin[lacto$relA >= 0.2 &
                     lacto$relA < 0.4] <- "20-39%"
lacto$LactoRelABin[lacto$relA >= 0.4 &
                     lacto$relA < 0.6] <- "40-59%"
lacto$LactoRelABin[lacto$relA >= 0.6 &
                     lacto$relA < 0.8] <- "60-79%"
#order
lacto$LactoRelABin <- factor(lacto$LactoRelABin, 
                             levels = c("0-19%",
                                        "20-39%",
                                        "40-59%",
                                        "60-79%",
                                        "> 80%"))
nmdsPlotList2 <- map2(nmdsBCDFs, grpz, function(n, g) {
  #plot
  p <- lacto %>%
    dplyr::select(Sample, relA) %>%
    inner_join(., n, by = "Sample") %>%
    ggplot(aes(x = NMDS1.BC, y = NMDS2.BC, colour = relA)) +
    geom_point(size = 3, alpha = 0.75) +
    labs(colour = "Relative\nabundance",
         title = g,
         subtitle = "Lactobacillus ",
         x = "NMDS1",
         y = "NMDS2") +
    theme_bw(base_size = 13) +
    theme(legend.position = "right", 
          plot.subtitle = element_text(face = "italic")) +
    scale_colour_gradient2(low = "#aec1ef",
                           mid = "#5D83DF", #lactobacillus colour from other plots
                           high = "#0a1633",
                           midpoint = 0.5)
  #save
  title <- gsub(" ", "-", g)
  ggsave(paste0("output/microbiome/beta/Bray-Curtis_", title, "_LactoRelA_NMDS.png"))
  return(p)
})

#what's dominant when lacto isn't?
#for a taxa to be dominant it needs to have relative abundance >= 30% 
#following precedent set by Labeer et al 2022, Brooks et al 2016
topDF <- genCnts %>%
  filter(Status == "Atrophy") %>%
  group_by(Sample) %>%
  filter(relA == max(relA)) %>%
  ungroup() %>%
  #if there is no taxa at >=30%, sample is considered to have no dominance
  mutate(Genus = ifelse(relA < 0.3, "No dominance", Genus)) %>%
  dplyr::select(Sample, Patient.ID, Genus, Group) %>%
  dplyr::rename("DominantGenus" = "Genus") %>%
  mutate(Group = factor(Group, 
                            levels = c("Control",
                                       "Recruitment",
                                       "Follow up"))) %>%
  inner_join(lacto[,c("Sample", "LactoRelABin")], .,  by = "Sample") %>%
  ungroup() 
  
#slim down
topSlim <- topDF %>% 
  select(LactoRelABin, DominantGenus, Group) %>%
  unique() %>%
  group_by(LactoRelABin, Group) %>%
  mutate(DominantGenus = paste(DominantGenus, collapse = ", ")) %>%
  unique()
#get numbers per lactobacillus bin and save
genCnts %>%
  filter(Status == "Atrophy") %>%
  inner_join(., lacto[, c("Sample", "LactoRelABin")], by = "Sample") %>%
  select(Sample, Group, LactoRelABin) %>%
  unique() %>%
  group_by(Group, LactoRelABin) %>%
  tally() %>%
  inner_join(topSlim, by = c("Group", "LactoRelABin")) %>%
  unique() %>%
  write.csv("output/microbiome/community/DominantGenus_Per_LactoGroup.csv",
            quote = F, row.names = F)

#combine data
nmdsBCDFs <- map(nmdsBCDFs, function(n) {
  n <- n %>%
    #add sample metadata 
    inner_join(., topDF[, c("Sample", "DominantGenus")], by = "Sample") %>%
    ungroup()
  #add in binned NGAT severities
  #add in binned NGAT severities
  n$NGATBin[n$NGAT > 10] <- "Severe"
  n$NGATBin[n$NGAT <= 10 & 
              n$NGAT >= 5] <- "Moderate"
  n$NGATBin[n$NGAT < 5] <- "Mild"
  return(n)
})
nmdsPlot <- bind_rows(nmdsBCDFs) %>%
  mutate(Group = factor(Group, levels = c("Recruitment",
                                          "Follow up")),
         DominantGenus = 
           factor(DominantGenus, levels = c(setdiff(unique(DominantGenus), 
                                                    "No dominance"),
                                            "No dominance")))

#add no dominance to genPal
genPal["No dominance"] <- "black"

#plot
ggplot(nmdsPlot, aes(x = NMDS1.BC, y = NMDS2.BC, colour = DominantGenus)) +
  labs(subtitle = "Bray-Curtis",
       col = "Dominant Genus",
       shape = "Atrophy Severity",
       x = "NMDS1",
       y = "NMDS2") +
  geom_point(size = 4, aes(shape = NGATBin)) +
  theme_bw(base_size = 13) +
  scale_colour_manual(values = genPal,
                      guide = guide_legend(label.theme = element_text(face = "italic"))) +
  facet_grid(~Group) +
  theme(legend.position = "right") 
#save
ggsave("output/microbiome/beta/DominantGenus_AtrophySeverity_NMDS_Bray-Curtis.png")

#There are no variables that can be checked for dispersion with
#betadisper as there are no suitable (categorical) variables going into 
#PERMANOVA. Ethnicity is 5 levels, but 4 of those are one individuals
#each which makes modelling too unstable.
table(met$Ethnicity, met$Group)
#log testable variables (non-continuous)
#"Smoking.status" also excluded, as there's no smokers in the cohort
#permanova

fullMod <- map2(distsBC, grpz, function(d, g) {
  fm <- adonis2(d ~ NGAT + DIVA + Age + BMI,
          data = nmdsPlot[nmdsPlot$Group == g, ], 
          by = "margin",
          permutations = 10000, method = "bray")
  lab <- gsub(" ", "-", g)
  capture.output(fm,
                 file = paste0("output/microbiome/beta/PERMANOVA_Bray-Curtis_",
                               lab, "_FullModel.txt"))
  return(fm)
  
})  

#check other variables for collinearity
form <- ~ Age + BMI + Ethnicity + DIVA + NGAT
corMat <- canCorPairs(form, nmdsPlot)
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
ggsave("output/microbiome/clinicalMetric/clinicalMetricCollinear_CCA.png")

#deeper look at NGAT and age = atrophy severity and age
#there's no real way of unlinking these - atrophy is age related disease
#Age and DIVA correlations with NGAT are computed separately per group
#(Recruitment vs Follow up) rather than pooled across both timepoints,
#since pooling would pseudoreplicate the test (each patient contributes
#two, non-independent, rows)
#there will be a number of NGAT correlations so I'm going to store and
#correct together
metaCors <- data.frame(Cohort = character(),
                   Var1 = character(),
                   Var2 = character(),
                   Statistic = character(),
                   StatValue = numeric(),
                   tau = numeric(),
                   Exact = logical(),
                   p = numeric())
recDF <- nmdsPlot %>%
  filter(Group == "Recruitment")
fuDF <- nmdsPlot %>%
  filter(Group == "Follow up")

#Recruitment
test <- cor.test(recDF$Age, 
                 recDF$NGAT,
                 method = "kendall")
metaCors[1,] <- c("Recruitment", "Age", "NGAT", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))
#Follow up
test <- cor.test(fuDF$Age,
                 fuDF$NGAT,
                 method = "kendall")
metaCors[2,] <- c("Follow up", "Age", "NGAT", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))

#for completeness, more closely assess DIVA and age
#Recruitment
test <- cor.test(recDF$Age, 
                 recDF$DIVA,
                 method = "kendall")
metaCors[3,] <- c("Recruitment", "Age", "DIVA", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))
#Follow up
test <- cor.test(fuDF$Age,
                 fuDF$DIVA,
                 method = "kendall")
metaCors[4,] <- c("Follow up", "Age", "DIVA", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))

#deeper look at NGAT and DIVA
#as NGAT and DIVA both measure atrophy severity, one would hope
#they correlate
#Recruitment
test <- cor.test(recDF$DIVA, recDF$NGAT,
                 method = "kendall")
metaCors[5,] <- c("Recruitment", "DIVA", "NGAT", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))
#Follow up
test <- cor.test(fuDF$DIVA,
                 fuDF$NGAT,
                 method = "kendall")
metaCors[6,] <- c("Follow up", "DIVA", "NGAT", names(test$statistic),
                  unname(test$statistic), unname(test$estimate),
                  names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                  unname(test$p.value))

#correct
metaCors <- metaCors %>%
  mutate(p = as.numeric(p)) %>%
  mutate(p.adj = p.adjust(p, method = "fdr"))

#NGAT and age are correlated - but no way to say which is driving the 
#other. Just looking at metrics (variables of interest)
redMod <- map2(distsBC, grpz, function(d, g) {
  rm <- adonis2(d ~ NGAT + DIVA,
                data = nmdsPlot[nmdsPlot$Group == g, ], 
                by = "margin",
                permutations = 10000, method = "bray")
  lab <- gsub(" ", "-", g)
  capture.output(rm,
                 file = paste0("output/microbiome/beta/PERMANOVA_Bray-Curtis_",
                               lab, "_ReducedModel.txt"))
  return(rm)
  
})  

#compare the two models
fullModAll <- map2(distsBC, grpz, function(d, g) {
  fma <- adonis2(d ~ NGAT + DIVA + Age + BMI,
                 data = nmdsPlot[nmdsPlot$Group == g, ], 
                 permutations = 10000, method = "bray")
  lab <- gsub(" ", "-", g)
  capture.output(fma,
                 file = paste0("output/microbiome/beta/PERMANOVA_Bray-Curtis_",
                               lab, "_FullModel_byModel.txt"))
  return(fma)
})
redModAll <- map2(distsBC, grpz, function(d, g) {
  rma <- adonis2(d ~ NGAT + DIVA,
                 data = nmdsPlot[nmdsPlot$Group == g, ], 
                 permutations = 10000, method = "bray")
  lab <- gsub(" ", "-", g)
  capture.output(rma,
                 file = paste0("output/microbiome/beta/PERMANOVA_Bray-Curtis_",
                               lab, "_ReducedModel_byModel.txt"))
  return(rma)
})

#check DIVA and NGAT
plotsBC <- vector(mode = "list")
#prepare labels for facet plot
labBCNGAT <- data.frame(Group = factor(c("Recruitment", 
                                  "Follow up")),
                        x = c(-0.25, -0.25),
                        y = c(-0.35, -0.35),
                        redR2 = c(round(redMod[[1]]["NGAT", c("R2")], 3), 
                               round(redMod[[2]]["NGAT", c("R2")], 3)),
                        redpVal = c(redMod[[1]]["NGAT", c("Pr(>F)")],
                                    redMod[[2]]["NGAT", c("Pr(>F)")]),
                        fullR2 = c(round(fullMod[[1]]["NGAT", c("R2")], 3), 
                                   round(fullMod[[2]]["NGAT", c("R2")], 3)),
                        fullpVal = c(fullMod[[1]]["NGAT", c("Pr(>F)")],
                                     fullMod[[2]]["NGAT", c("Pr(>F)")])) %>%
  mutate(pValLabRed = ifelse(redpVal < 0.001, "< 0.001", paste("=", round(redpVal, 3))),
         pValLabFull = ifelse(fullpVal < 0.001, "< 0.001", paste("=", round(fullpVal, 3))),
         label = paste0("\nPERMANOVA\nReduced: R² = ", redR2, 
                        ", p ", pValLabRed, "\nFull: R² = ", fullR2, 
                        ", p ", pValLabFull))


#NGAT
midP <- ((max(nmdsPlot$NGAT) - min(nmdsPlot$NGAT)) / 2) + min(nmdsPlot$NGAT)
plotsBC[[1]] <- ggplot(nmdsPlot, aes(x = NMDS1.BC, y = NMDS2.BC, colour = NGAT)) +
  labs(col = "Score",
       subtitle = "NGAT") +
  geom_point(size = 3) +
  scale_colour_gradient2(low = vanPal[[6]],
                         mid = vanPal[[7]],
                         high = "#050808",
                         midpoint = midP) +
  facet_grid(~Group) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(legend.position = "right",
        plot.subtitle = element_text(face = "bold"))  + 
  geom_text(data = labBCNGAT, size = 3,
            aes(x = x, y = y, label = label, fontface = "italic"),
            inherit.aes = F, hjust = 0)
#DIVA
#prepare labels for facet plot
labBCDIVA <- data.frame(Group = factor(c("Recruitment", 
                                         "Follow up")),
                        x = c(-0.25, -0.25),
                        y = c(-0.35, -0.35),
                        redR2 = c(round(redMod[[1]]["DIVA", c("R2")], 3), 
                                  round(redMod[[2]]["DIVA", c("R2")], 3)),
                        redpVal = c(redMod[[1]]["DIVA", c("Pr(>F)")],
                                    redMod[[2]]["DIVA", c("Pr(>F)")]),
                        fullR2 = c(round(fullMod[[1]]["DIVA", c("R2")], 3), 
                                   round(fullMod[[2]]["DIVA", c("R2")], 3)),
                        fullpVal = c(fullMod[[1]]["DIVA", c("Pr(>F)")],
                                     fullMod[[2]]["DIVA", c("Pr(>F)")])) %>%
  mutate(pValLabRed = ifelse(redpVal < 0.001, "< 0.001", paste("=", round(redpVal, 3))),
         pValLabFull = ifelse(fullpVal < 0.001, "< 0.001", paste("=", round(fullpVal, 3))),
         label = paste0("\nPERMANOVA\nReduced: R² = ", redR2, 
                        ", p ", pValLabRed, "\nFull: R² = ", fullR2, 
                        ", p ", pValLabFull))

midP <- ((max(nmdsPlot$DIVA) - min(nmdsPlot$DIVA)) / 2) + min(nmdsPlot$DIVA)
plotsBC[[2]] <- ggplot(nmdsPlot, aes(x = NMDS1.BC, y = NMDS2.BC, colour = DIVA)) +
  labs(col = "Score",
       subtitle = "DIVA") +
  geom_point(size = 3) +
  scale_colour_gradient2(low = vanPal[[6]],
                         mid = vanPal[[7]],
                         high = "#050808",
                         midpoint = midP) +
  facet_grid(~Group) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(legend.position = "right",
        plot.subtitle = element_text(face = "bold")) + 
  geom_text(data = labBCDIVA, size = 3,
            aes(x = x, y = y, label = label, fontface = "italic"),
            inherit.aes = F, hjust = 0)
pBetaClinical <- ggarrange(plotlist = plotsBC, ncol = 1)
#save
ggsave("output/microbiome/beta/AtrophySeverity_Bray-Curtis_DIVA_NGAT_NMDS.png")

##### Beta diversity: Robust Aitchinson #####
#note - robust aitchinson is said to be scale invariant (Aitchison J (1982))
#and as such it is said rarefaction is not needed - however, more recent
#publications dispute this (Schloss 2024 / 2026). For completeness,
#rarefaction again will be used here. However, Aitchison should handle
#the compositionality much better than Bray-Curtis, and be less 
#sensitive to Lactobacillus dominance allowing for more "community" level
#separations.
distsRA <- map(betaCnts, function(b) {
  return(avgdist(b, dmethod = "robust.aitchison", iterations = 10000,
                 sample = nLmin))
})

#nmds
#plot stress v dimensionality plots
map2(distsRA, grpz, function(d, g) {
  #plot
  out <- chooseK(d, type = "plot", method = "robust.aitchison") +
    labs(title = paste0(g, ": Robust Aitchison"))
  #save
  lab <- gsub(" ", "-", g)
  ggsave(paste0("output/microbiome/beta/RobustAichison", lab,
                "_StressVsDimensions.png"))
  return(out)
})

#decide on k parameter
kListRA <- map(distsRA, chooseK, method = "robust.aitchison")

#nmds
nmdsListRA <- map2(distsRA, kListRA, function(d, k){
  metaMDS(d, k = k, trymax = 100, trace = F, distance = "robust.aitchison")
}) 

#plot stressplot
map2(nmdsListRA, grpz, function(n, g) {
  lab <- gsub(" ", "-", g)
  png(paste0("output/microbiome/beta/NMDS_RobustAitchison_", lab, "_Stressplot.png"))
  stressplot(n)
  dev.off()
})

#produce multiple plots
nmdsRADFs <- map(nmdsListRA, function(n) {
  out <- scores(n) %>%
    data.frame() %>%
    select(NMDS1, NMDS2) %>%
    rownames_to_column(var = "Sample") %>%
    #add sample metadata 
    inner_join(., met, by = "Sample") %>%
    rename("NMDS1.RA" = "NMDS1",
           "NMDS2.RA" = "NMDS2")
})

nmdsPlotListRA <- map2(nmdsRADFs, grpz, function(n, g) {
  #prepare to store simple explorations
  plots <- vector(mode = "list", length = length(vois))
  #figure out how many categorical levels I have to plot
  values <- vector()
  for (i in 1:length(vois)) {
    if (!is.numeric(met[[vois[i]]])) {
      values <- c(values, unique(met[[vois[i]]]))
    }
  }
  plotPal <- palFun(length(values))
  #iterate through variables and plot and save
  for (i in 1:length(vois)) {
    plotDF <- n %>%
      mutate(Group = n[[vois[i]]])
    lab <- gsub("\\.", " ", vois[[i]])
    #store base plot
    plots[[i]] <- plotDF %>%
      ggplot(aes(x = NMDS1.RA, y = NMDS2.RA)) +
      geom_point(aes(colour = Group),
                 size = 3, alpha = 0.75) +
      labs(subtitle = lab,
           colour = "",
           x = "NMDS1",
           y = "NMDS2") +
      theme(plot.subtitle = element_text(face = "bold"))
    #figure out colour schemes
    if (!is.numeric(plotDF$Group)) {
      plots[[i]] <- plots[[i]] + 
        scale_colour_manual(values = vanPal[c(4,9,1,8,7)]) 
    } else {
      midP <- ((max(plotDF$Group) - min(plotDF$Group)) / 2) + min(plotDF$Group)
      plots[[i]] <- plots[[i]] +
        scale_colour_gradient2(low = vanPal[[6]],
                               mid = vanPal[[7]],
                               high = "#050808",
                               midpoint = midP)
    }
  }
  p <- ggarrange(plotlist = plots, ncol = 2, nrow = 3)
  p <- annotate_figure(p, top = text_grob(g, size = 15, face = "bold"))
  #save
  title <- gsub(" ", "-", g)
  ggsave(paste0("output/microbiome/beta/NMDS_RobustAitchison_", title,
                "_MetaVariables.pdf"), scale = 2)
  ggsave(paste0("output/microbiome/beta/NMDS_RobustAitchison_", title,
                "_MetaVariables.png"), scale = 2)
  return(p)
})


nmdsPlotList2 <- map2(nmdsRADFs, grpz, function(n, g) {
  #plot
  p <- lacto %>%
    dplyr::select(Sample, relA) %>%
    inner_join(., n, by = "Sample") %>%
    ggplot(aes(x = NMDS1.RA, y = NMDS2.RA, colour = relA)) +
    labs(colour = "Relative\nabundance",
         title = g,
         subtitle = "Lactobacillus ",
         x = "NMDS1",
         y = "NMDS2") +
    geom_point(size = 3, alpha = .75) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right", 
          plot.subtitle = element_text(face = "italic")) +
    scale_colour_gradient2(low = "#aec1ef",
                           mid = "#5D83DF", #lactobacillus colour from other plots
                           high = "#0a1633",
                           midpoint = 0.5)
  #save
  title <- gsub(" ", "-", g)
  ggsave(paste0("output/microbiome/beta/RobustAitchison_", title, "_LactoRelA_NMDS.png"))
  return(p)
})

#combine all data
#combine data
nmdsBCDFs <- map(nmdsBCDFs, function(n) {
  n <- n %>%
    #add sample metadata 
    inner_join(., topDF[, c("Sample", "DominantGenus")], by = "Sample") %>%
    ungroup()
  #add in binned NGAT severities
  #add in binned NGAT severities
  n$NGATBin[n$NGAT > 10] <- "Severe"
  n$NGATBin[n$NGAT <= 10 & 
              n$NGAT >= 5] <- "Moderate"
  n$NGATBin[n$NGAT < 5] <- "Mild"
  return(n)
})

nmdsPlot <- bind_rows(nmdsRADFs) %>%
  select(Sample, NMDS1.RA, NMDS2.RA) %>%
  inner_join(., nmdsPlot)

#plot
ggplot(nmdsPlot, aes(x = NMDS1.RA, y = NMDS2.RA, colour = DominantGenus)) +
  labs(col = "Dominant Genus",
       shape = "Atrophy Severity") +
  geom_point(size = 4, aes(shape = NGATBin)) +
  theme_bw(base_size = 13) +
  scale_colour_manual(values = genPal,
                      guide = guide_legend(label.theme = element_text(face = "italic"))) +
  facet_grid(~Group) +
  theme(legend.position = "right") 
#save
ggsave("output/microbiome/beta/DominantGenus_RobustAitchision_AtrophySeverity_NMDS.png")

#permanova
#Group dropped from all models below: beta diversity is restricted to
#Per visit samples, so Group is constant and cannot be modelled
fullMod <- map2(distsRA, grpz, function(d, g) {
  fm <- adonis2(d ~ NGAT + DIVA + Age + BMI,
                data = nmdsPlot[nmdsPlot$Group == g, ], 
                by = "margin",
                permutations = 10000, method = "robust.aitchison")
  lab <- gsub(" ", "-", g)
  capture.output(fm,
                 file = paste0("output/microbiome/beta/PERMANOVA_RobustAitchison_",
                               lab, "_FullModel.txt"))
  return(fm)
  
})  


#NGAT and age are correlated - but no way to say which is driving the 
#other. Just looking at metrics (variables of interest)
redMod <- map2(distsRA, grpz, function(d, g) {
  rm <- adonis2(d ~ NGAT + DIVA,
                data = nmdsPlot[nmdsPlot$Group == g, ], 
                by = "margin",
                permutations = 10000, method = "robust.aitchison")
  lab <- gsub(" ", "-", g)
  capture.output(rm,
                 file = paste0("output/microbiome/beta/PERMANOVA_RobustAitchison_",
                               lab, "_ReducedModel.txt"))
  return(rm)
  
})  

#compare the two models
fullModAll <- map2(distsRA, grpz, function(d, g) {
  fma <- adonis2(d ~ NGAT + DIVA + Age + BMI,
                 data = nmdsPlot[nmdsPlot$Group == g, ], 
                 permutations = 10000, method = "robust.aitchison")
  lab <- gsub(" ", "-", g)
  capture.output(fma,
                 file = paste0("output/microbiome/beta/PERMANOVA_RobustAitchison_",
                               lab, "_FullModel_byModel.txt"))
  return(fma)
})

redModAll <- map2(distsRA, grpz, function(d, g) {
  rma <- adonis2(d ~ NGAT + DIVA,
                 data = nmdsPlot[nmdsPlot$Group == g, ], 
                 permutations = 10000, method = "robust.aitchison")
  lab <- gsub(" ", "-", g)
  capture.output(rma,
                 file = paste0("output/microbiome/beta/PERMANOVA_RobustAitchison_",
                               lab, "_ReducedModel_byModel.txt"))
  return(rma)
})


#check DIVA and NGAT
plots <- vector(mode = "list")
#prepare labels for facet plot
labRANGAT <- data.frame(Group = factor(c("Recruitment", 
                                         "Follow up")),
                        x = c(-15, -15),
                        y = c(-6.75, -6.75),
                        redR2 = c(round(redMod[[1]]["NGAT", c("R2")], 3), 
                                  round(redMod[[2]]["NGAT", c("R2")], 3)),
                        redpVal = c(redMod[[1]]["NGAT", c("Pr(>F)")],
                                    redMod[[2]]["NGAT", c("Pr(>F)")]),
                        fullR2 = c(round(fullMod[[1]]["NGAT", c("R2")], 3), 
                                   round(fullMod[[2]]["NGAT", c("R2")], 3)),
                        fullpVal = c(fullMod[[1]]["NGAT", c("Pr(>F)")],
                                     fullMod[[2]]["NGAT", c("Pr(>F)")])) %>%
  mutate(pValLabRed = ifelse(redpVal < 0.001, "< 0.001", paste("=", round(redpVal, 3))),
         pValLabFull = ifelse(fullpVal < 0.001, "< 0.001", paste("=", round(fullpVal, 3))),
         label = paste0("\nPERMANOVA\nReduced: R² = ", redR2, 
                        ", p ", pValLabRed, "\nFull: R² = ", fullR2, 
                        ", p ", pValLabFull))


#NGAT
midP <- ((max(nmdsPlot$NGAT) - min(nmdsPlot$NGAT)) / 2) + min(nmdsPlot$NGAT)
plots[[1]] <- ggplot(nmdsPlot, aes(x = NMDS1.RA, y = NMDS2.RA, colour = NGAT)) +
  labs(col = "Score",
       subtitle = "NGAT") +
  geom_point(size = 3) +
  theme_bw(base_size = 13) +
  scale_colour_gradient2(low = vanPal[[6]],
                         mid = vanPal[[7]],
                         high = "#050808",
                         midpoint = midP) +
  facet_grid(~Group) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(legend.position = "right",
        plot.subtitle = element_text(face = "bold"))  + 
  geom_text(data = labRANGAT, size = 3,
            aes(x = x, y = y, label = label, fontface = "italic"),
            inherit.aes = F, hjust = 0)
#DIVA
#prepare labels for facet plot
labRADIVA <- data.frame(Group = factor(c("Recruitment", 
                                         "Follow up")),
                        x = c(-18.5, -18.5),
                        y = c(-6.75, -6.75),
                        redR2 = c(round(redMod[[1]]["DIVA", c("R2")], 3), 
                                  round(redMod[[2]]["DIVA", c("R2")], 3)),
                        redpVal = c(redMod[[1]]["DIVA", c("Pr(>F)")],
                                    redMod[[2]]["DIVA", c("Pr(>F)")]),
                        fullR2 = c(round(fullMod[[1]]["DIVA", c("R2")], 3), 
                                   round(fullMod[[2]]["DIVA", c("R2")], 3)),
                        fullpVal = c(fullMod[[1]]["DIVA", c("Pr(>F)")],
                                     fullMod[[2]]["DIVA", c("Pr(>F)")])) %>%
  mutate(pValLabRed = ifelse(redpVal < 0.001, "< 0.001", paste("=", round(redpVal, 3))),
         pValLabFull = ifelse(fullpVal < 0.001, "< 0.001", paste("=", round(fullpVal, 3))),
         label = paste0("\nPERMANOVA\nReduced: R² = ", redR2, 
                        ", p ", pValLabRed, "\nFull: R² = ", fullR2, 
                        ", p ", pValLabFull))

midP <- ((max(nmdsPlot$DIVA) - min(nmdsPlot$DIVA)) / 2) + min(nmdsPlot$DIVA)
plots[[2]] <- ggplot(nmdsPlot, aes(x = NMDS1.RA, y = NMDS2.RA, colour = DIVA)) +
  labs(col = "Score",
       subtitle = "DIVA") +
  geom_point(size = 3) +
  theme_bw(base_size = 13) +
  scale_colour_gradient2(low = vanPal[[6]],
                         mid = vanPal[[7]],
                         high = "#050808",
                         midpoint = midP) +
  facet_grid(~Group) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(legend.position = "right",
        plot.subtitle = element_text(face = "bold")) + 
  geom_text(data = labRADIVA, size = 3,
            aes(x = x, y = y, label = label, fontface = "italic"),
            inherit.aes = F, hjust = 0)
pBetaClinicalRA <- ggarrange(plotlist = plots, ncol = 1)
#save
ggsave("output/microbiome/beta/AtrophySeverity_RobustAitchison_DIVA_NGAT_NMDS.png",
       height = 18, width = 20, unit = "cm")

##### Lactobacillus ####
lactoRelA <- lacto %>% 
  mutate(Group = factor(Group, levels = c("Control", "Recruitment", "Follow up"))) %>%
  ggplot(aes(x = Group, y = relA)) +
  theme_bw(base_size = 13) +
  geom_violin(alpha = 0.5, width = 0.5, aes(fill = Group)) +
  geom_boxplot(width = 0.1, alpha = 0.5, aes(colour = Group)) +
  geom_point(aes(colour = Group), size = 2, position = "jitter") +
  scale_colour_manual(values = vanPal) +
  scale_fill_manual(values = vanPal) +
  labs(y = expression(italic("Lactobacillus")* " relative abundance")) +
  guides(fill = "none",
         colour = "none")
ggsave("output/microbiome/community/LactobacillusRelA_All.png")

#breakdown by species
taxDF2 <- tax_table(phylo) %>%
  data.frame() %>%
  rownames_to_column(var = "ASV")
lactoSpec <- otu_table(phylo) %>%
  data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample) %>%
  dplyr::rename("count"= "value",
                "ASV" = "name") %>%
  group_by(Sample) %>%
  mutate(relA = count / sum(count)) %>%
  #get sample metdata
  inner_join(., met, by = "Sample") %>%
  #get taxonomic information
  inner_join(., taxDF2, by = "ASV") %>%
  #keep only Lactobacillus
  filter(Genus == "Lactobacillus") %>%
  mutate(Label = paste(Genus, Species),
         Label = gsub("NA", "(unclassified)", Label))
#plot
lactoSpec %>%
  ggplot(aes(x = Group, y = relA)) +
  theme_bw(base_size = 13) +
  geom_point(aes(colour = Group), size = 2) +
  scale_colour_manual(values = vanPal) +
  facet_wrap(.~Label,
             labeller = label_wrap_gen(width = 10)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  labs(y = "Relative abundance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
#save
ggsave("output/microbiome/community/Lactobacillus_relA_BySpecies.png")

#prevalence
sampsPerTreat <- lactoSpec %>%
  select(Sample, Group) %>%
  unique() %>%
  group_by(Group) %>%
  tally()
lactoSpec %>%
  ungroup() %>%
  select(Sample, Label, Group, count) %>%
  group_by(Label, Sample) %>%
  mutate(totCount = sum(count)) %>%
  mutate(Presence = ifelse(totCount > 0, 1, 0),
         .after = Label) %>%
  ungroup() %>%
  select(Sample, Label, Presence, Group) %>%
  inner_join(., sampsPerTreat, by = "Group") %>%
  unique() %>%
  group_by(Label, Group) %>%
  mutate(TotalTreat = sum(Presence), .after = Presence) %>%
  select(Label, TotalTreat, n, Group) %>%
  unique() %>%
  mutate(Prevalence = TotalTreat / n) %>%
  ungroup() %>%
  select(Label, Group, Prevalence) %>%
  unique() %>%
  ggplot(aes(x = Group, y = Label)) +
  geom_tile(aes(fill = Prevalence)) +
  scale_fill_gradient2(low = "white",
                         mid = "#5D83DF", #lactobacillus colour from other plots
                         high = "#0a1633",
                         midpoint = 0.5) +
  theme(axis.text.y = element_text(face = "italic"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
#save
ggsave("output/microbiome/community/Lactobacillus_relA_BySpecies_Heatmap.jpg")

#what is the genus prevalence?
all(cnts$ASV1 > 0)
#100%

#causally, there is an already characterised mechanism to suggest
#Lactobacillus level will be associated with GSM severity
#with Lactobacillus also being 100% prevalent it is worth an attempt
#at correlating RELATIVE abundance with clinical metrics (NGAT/DIVA)
#this of course does not mean more / less Lactobacillus = association
#either way as we don't have an absolute quantification - rather
#its level relative to the whole community may associate.
#again, using kendall with NGAT/DIVA as ties are likely, nonparametric etc
lactoCors <- data.frame(Cohort = character(),
                        Var1 = character(),
                        Var2 = character(),
                        Statistic = character(),
                        StatValue = numeric(),
                        tau = numeric(),
                        Exact = logical(),
                        p = numeric())
#NGAT recruitment
test <- cor.test(lacto$relA[lacto$Group == "Recruitment"],
                 lacto$NGAT[lacto$Group == "Recruitment"], 
                 method = "kendall")
lactoCors[1,] <- c("Recruitment", "LactoRelA", "NGAT", names(test$statistic),
                   unname(test$statistic), unname(test$estimate),
                   names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                   unname(test$p.value))
#NGAT follow up
test <- cor.test(lacto$relA[lacto$Group == "Follow up"],
                 lacto$NGAT[lacto$Group == "Follow up"], 
                 method = "kendall")
lactoCors[2,] <- c("Follow up", "LactoRelA", "NGAT", names(test$statistic),
                   unname(test$statistic), unname(test$estimate),
                   names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                   unname(test$p.value))
#DIVA recruitment
test <- cor.test(lacto$relA[lacto$Group == "Recruitment"],
                 lacto$DIVA[lacto$Group == "Recruitment"], 
                 method = "kendall")
lactoCors[3,] <- c("Recruitment", "LactoRelA", "DIVA", names(test$statistic),
                   unname(test$statistic), unname(test$estimate),
                   names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                   unname(test$p.value))
#DIVA follow up
test <- cor.test(lacto$relA[lacto$Group == "Follow up"],
                 lacto$DIVA[lacto$Group == "Follow up"], 
                 method = "kendall")
lactoCors[4,] <- c("Follow up", "LactoRelA", "DIVA", names(test$statistic),
                   unname(test$statistic), unname(test$estimate),
                   names(test$statistic) == "T", #if exact pvalues are computed tau's test statistic is T
                   unname(test$p.value))
#correct
lactoCors <- lactoCors %>%
  mutate(p = as.numeric(p)) %>%
  mutate(p.adj = p.adjust(p, method = "fdr"))
#prepare to use in visualisation
NGATlab <- lactoCors %>%
  filter(Var2 == "NGAT") %>%
  mutate(y = c(.8, .7),
         x = c(9, 11),
         tau = as.numeric(tau),
         p.lab = ifelse(p.adj < 0.001, '< 0.001', round(p.adj, 3)), 
         label = paste0("tau: ", round(tau, 3), 
                        "\nadj.P: ", p.lab)) %>%
  rename("Group" = "Cohort") %>%
  mutate(Group = factor(Group, levels = c("Recruitment", "Follow up")))
DIVAlab <- lactoCors %>%
  filter(Var2 == "DIVA") %>%
  mutate(y = c(0.30, .2),
         x = c(20, 70),
         tau = as.numeric(tau),
         p.lab = ifelse(p.adj < 0.001, '< 0.001', round(p.adj, 3)), 
         label = paste0("tau: ", round(tau, 3), 
                        "\nadj.P: ", p.lab)) %>%
  rename("Group" = "Cohort") %>%
  mutate(Group = factor(Group, levels = c("Recruitment", "Follow up")))

#visualise
plotsLacto <- vector(mode = "list", length = 2)
plotsLacto[[1]] <- lacto %>%
  filter(!Group == "Control") %>%
  mutate(Group = factor(Group, levels = c("Recruitment", "Follow up"))) %>%
  ggplot(aes(x = NGAT, y = relA)) + 
  geom_hline(yintercept = 1, colour = "darkgrey", alpha = 0.75, 
             linetype = "dashed") + 
  geom_hline(yintercept = 0, colour = "darkgrey", alpha = 0.75, 
             linetype = "dashed") +
  stat_smooth(method = "lm", colour = "darkgrey") +
  geom_point() +
  facet_grid(~Group) +
  geom_text(data = NGATlab, size = 3, fontface = "italic",
            aes(x = x, y = y, label = label),
            hjust = 0) +
  labs(y = expression(atop(italic("Lactobacillus"), "relative abundance"))) +
  coord_cartesian(ylim = c(-0.05, 1.05), expand = FALSE) +
  theme(panel.spacing = unit(1.2, "lines"))
plotsLacto[[2]] <- lacto %>%
  filter(!Group == "Control") %>%
  mutate(Group = factor(Group, levels = c("Recruitment", "Follow up"))) %>%
  ggplot(aes(x = DIVA, y = relA)) + 
  geom_hline(yintercept = 1, colour = "darkgrey", alpha = 0.75, 
             linetype = "dashed") + 
  geom_hline(yintercept = 0, colour = "darkgrey", alpha = 0.75, 
             linetype = "dashed") +
  stat_smooth(method = "lm", colour = "darkgrey") +
  geom_point() +
  facet_grid(~Group) +
  geom_text(data = DIVAlab, size = 3, fontface = "italic",
            aes(x = x, y = y, label = label), hjust = 0) +
  labs(y = expression(atop(italic("Lactobacillus"), "relative abundance"))) 
corLacto <- ggarrange(plotlist = plotsLacto, nrow = 2, widths = 15)
ggsave("output/microbiome/clinicalMetric/LactobacillusRelA_vs_NGATDIVA.png")


##### Microbiome-Clinical Association: MaAslin3 ####
#a separate MaAsLin3 run with identical specification will be fitted for 
#each clinical score. NGAT/DIVA are not correlated at all, so adjusting
#for either does not affect the other. There is no benefit to a joint
#model here, but there would be a considerable multiple testing correction cost
#Recruitment is already reference factor level for group in NMDSplot
#can include both timepoints as patientID can be modelled as a random
#effect. 

#will be running NGAT/DIVA separately purely due the otherwise huge
#multiple testing cost (correction is applied across taxa tested * variable tested
#so gets very considerable very quickly)

#prepare count matrices - MaAslin3 requires samples as columns
#remove controls
gsm <- met$Sample[!met$Group == "Control"]
tMCnts <- t(cnts[rownames(cnts) %in% gsm,])
#add genus names > ASV identifiers
genLab <- taxDF2$Genus[match(rownames(tMCnts), taxDF2$ASV)]
rownames(tMCnts) <- genLab

#add read depth from total counts per patient
#note - these counts have been taken as sum of bacterial ASV counts after
#removal of reagent contaminants, mitochondrial and non-bacterial sequences, 
#and ASVs with ≤10 reads across the dataset.
modDF <- data.frame(colSums(tMCnts)) %>%
  rownames_to_column(var = "Sample") %>%
  inner_join(., nmdsPlot) %>%
  rename("Reads" = "colSums.tMCnts.") %>%
  column_to_rownames(var = "Sample")

#NGAT
#NGAT and Age are very closely correlated. Causally, this makes sense -
#if NGAT is a good measure of atrophy severity then
#age > menopause > eostrogen drop > Lacto reduction + GSM [NGAT]
#in this Age is an upstream mediator not a confounder
#for NGAT the model will be ran twice, with and without age
#and compared for discussion in the manuscript
NGATout1 <- maaslin3(input_data = tMCnts,
                    input_metadata = modDF,
                    output = 'output/microbiome/maaslin3/NGAT_noAge/',
                    formula = '~ NGAT + Group + 
                        Reads + (1 | Patient.ID)',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.05,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    plot_summary_plot = TRUE,
                    summary_plot_first_n = 25,
                    verbosity = "WARN")

NGATout2 <- maaslin3(input_data = tMCnts,
                    input_metadata = modDF,
                    output = 'output/microbiome/maaslin3/NGAT/',
                    formula = '~ NGAT + Age + Group + 
                        Reads + (1 | Patient.ID)',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.05,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    plot_summary_plot = TRUE,
                    summary_plot_first_n = 25,
                    verbosity = "WARN")

#DIVA
#DIVA and Age are not at all correlated, and so age will be included
DIVAout <- maaslin3(input_data = tMCnts,
                    input_metadata = modDF,
                    output = 'output/microbiome/maaslin3/DIVA/',
                    formula = '~ DIVA + Group + Age +
                        Reads + (1 | Patient.ID)',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.05,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    plot_summary_plot = TRUE,
                    summary_plot_first_n = 25,
                    verbosity = "WARN")

##### Combined figures ####
#manuscript figure: alpha / beta diversity and NGAT/DIVA
alphaBeta <- ggarrange(alphaPlot, pBetaClinical,
          widths = c(1.5, 4),
          labels = "AUTO")
ggsave("output/figures/alpha-beta.jpg", 
       width = 28, height = 19, unit = "cm")
ggsave("output/figures/alpha-beta.png", 
       width = 28, height = 19, unit = "cm")
ggsave("output/figures/alpha-beta.pdf", 
       width = 28, height = 19, unit = "cm")
save(alphaBeta, file = "output/figures/alpha-beta.rds")

#manuscript figure: lactobacillus correlations and levels
lactoFacet <- ggarrange(lactoRelA, corLacto,
          widths = c(1.5, 3),
          labels = "AUTO")
ggsave("output/figures/lacto-facet.jpg", 
       width = 21.5, height = 12, unit = "cm",
       scale = 1.2)
ggsave("output/figures/lacto-facet.png", 
       width = 21, height = 12, unit = "cm",
       scale = 1.2)
ggsave("output/figures/lacto-facet.pdf", 
       width = 21, height = 12, unit = "cm",
       scale = 1.2)
save(lactoFacet, file = "output/figures/lacto-facet.rds")

##### Outputs #####
#all data
vanData <- alpha %>%
  select(Sample, Shannon, Richness) %>%
  inner_join(., met, by = "Sample") %>%
  left_join(., nmdsPlot[, c("Sample", "NMDS1.BC", "NMDS2.BC",
                            "NMDS1.RA", "NMDS2.RA")], by = "Sample") %>%
  left_join(., lacto[, c("Sample", "relA", "LactoRelABin")], by = "Sample") %>%
  dplyr::rename("LactoRelA" = "relA")
write.csv(vanData,
          "output/microbiome/MicroAnalyses_Data.csv",
          row.names = F, quote = F)

#beta diversity distance matrices
save(distsBC,
     file = "processed/microbiome/BrayCurtisDistance.rds")
save(distsRA,
     file = "processed/microbiome/RobustAitchisionDistance.rds")

#rarefaction level
writeLines(capture.output(print(min(rowSums(cnts)))),
           "output/microbiome/alpha/RarefactionLevel.txt")

#save alpha statistics
write.csv(alpha, "output/microbiome/alpha/ShannonDiversity_Data.csv",
            row.names = F, quote = F)

#save tests
writeLines(capture.output(print(stats)),
           "output/microbiome/Stat_Tests.txt")
cors <- list(metaCors, lactoCors)
names(cors) <- c("Covariate correlation",
                 "Lactobacillus relative abundance correlations")
writeLines(capture.output(print(cors)),
           "output/microbiome/Correlation_Tests.txt")

#session information
writeLines(capture.output(sessionInfo()),
           "output/microbiome/Analyses_sessionInfo.txt")
