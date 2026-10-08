#November 2024, rewritten 21st Jan 2025
#Lauren Mee
#Produce and compare microbiome/metabolome networks
#built on Spearman's correlations

##### Set up #####
#set seed
set.seed(1152)

#libraries
libs <- c("phyloseq", #to access microbiome material (taxonomy lookups)
          "compositions", #clr transformation
          "ggraph", "igraph", "tidygraph", #for network manipulation
          "tidyverse") #general data manipulation, plots


for (pkg in libs) {
  library(pkg, character.only = T)
}

#load data
#metadata 
met <- read.csv("input/SampleMetadata.csv")

#microbiome data
load("processed/microbiome/PhyloSeqObjs.rds")
#extract taxonomy
tax <- tax_table(phyGen) %>%
  data.frame() %>%
  rownames_to_column(var = "ASV")

#nmr (already scaled and log-transformed)
#script 03 writes this with write.csv(): it is COMMA separated with a leading
#empty header field for the rownames. Read it as such.
y <- read.csv("processed/nmr/logPQN_NMR.csv", row.names = 1) %>%
  data.frame()

#output directory architecture
dirs <- c("output/", "output/network/", 
          "output/network/Recruitment", "output/network/FollowUp",
          "processed/", "processed/network/")

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

##### Functions #####
#plots a ggraph from an igraph object with
#all nodes labelled
plotGraph <- function(graphObj, 
                      labels = TRUE,
                      title = NULL) {
  p <- ggraph(graphObj, layout = "kk") +
    geom_edge_link(aes(colour = Relationship, 
                       width = Correlation)) +
    scale_edge_width(range = c(0.5, 3)) + 
    geom_node_point(aes(size = Degree, colour = Feature)) +
    scale_edge_colour_manual(values = vanPal[c(7,9)]) +
    scale_colour_manual(values = vanPal) +
    theme_void() +
    theme(legend.margin = margin(t = 5, r = 10, b = 5, l = 5))
  if (labels == TRUE) {
    p <- p + 
      geom_node_label(aes(label = name), repel = T)
  }
  if (!is.null(title)) {
    p <- p + 
      ggtitle(title) + 
      theme(plot.title = element_text(hjust = 0.5, 
                                      face = "bold",
                                      size = 16))
  }
  return(p)
}

##### Preprocessing #####
#metadata requirements
#NB NGAT change per patient (NGATDifference / NGATResponse) used to be
#derived here, but was only ever consumed by dead code paths in script 05.
#Clinical-response measures are now built properly in script 06.

#Keep only the atrophy patients and make sure there is a both recruitment and
#follow up sample patient.
atrophy <- met %>%
  filter(Status == "Atrophy") %>%
  select(Sample) %>%
  pull()
#remove samples from metadata that don't have both timepoints
toRem <- met %>%
  group_by(Patient.ID) %>%
  tally() %>%
  filter(n == 1) %>%
  select(Patient.ID) %>%
  pull()
#subset metadata
met <- met %>%
  filter(Sample %in% atrophy) %>%
  filter(!Patient.ID %in% toRem)

#genus-level counts from the phyloseq object built in script 01. Columns stay
#as ASV identifiers: the taxonomy lookup below joins on ASV and adds the genus
#name at that point.
cnts <- data.frame(otu_table(phyGen))

#first, I need to further filter out the microbiome data so that
#only microbes found in at least 30% of the samples are kept
#this should reduce the likelihood of spurious associations 
#occurring when "0" hits are CLR transformed into sightly different
#negative values, that can become associations that obviously
#then don't exist
lt30 <- apply(cnts, 2, function(x) sum(x > 0)) <= nrow(cnts) * 0.30
x <- cnts[, !colnames(cnts) %in% names(lt30)[lt30 == T]]

#apply clr transformations
x[x == 0] <- 0.5
x <- clr(x) %>%
  as.data.frame()

#Now to start limiting the dataframes to samples that are found in 
#both datasets.
both <- intersect(rownames(x), rownames(y))
#apply
y <- y[both, ]
x <- x[rownames(x) %in% both, ]

#check must be TRUE before continuing
all(rownames(x) == rownames(y))

##### Run correlation and permutation #####
#just GSM patients
grpz <- c("Recruitment", "Follow up")
samps <- map(grpz, function(run) {
  return(met$Sample[met$Group == run])
})

#subset dataframes
xBA <- map(samps, function(set) {
  return(x[rownames(x) %in% set,])
})
yBA <- map(samps, function(set) {
  return(y[rownames(y) %in% set,])
})

#generate correlations between microbiome and metabolome data
corsBA <- map2(xBA, yBA, 
    function(x2, y2) {
      #compute observed spearman correlation
      #cor = pairwise correlations between columns of the matrices
      corMat <- cor(x2, y2, method = "spearman")
      return(corMat)
    })

#run permutations to generate pvalues
pvalBA <- map2(xBA, yBA, function(x2, y2) {
  #prepare space for pvalues
  pMat <- matrix(NA, ncol(x2), ncol(y2))
  #loop through each pair of columns
  for (i in 1:ncol(x2)) {
    for (j in 1:ncol(y2)) {
      #get observed spearman correlation for the pair
      obsCor <- cor(x2[,i], y2[,j], method = "spearman")
      #perform permutation test
      perCor <- replicate(1000, {
        xPer <- sample(x2[,i]) #shuffle the ith column of x
        cor(xPer, y2[,j], method = "spearman")
      })
      #compute the empirical pvalue
      pval <- mean(abs(perCor) >= abs(obsCor))
      #store pvalue in the matrix
      pMat[i, j] <- pval
    }
  }
  return(pMat)
})

aPvalBA <- map(pvalBA, function(pMat) {
    #adjust all values together ...
    #correcting via a matrix introduces a correction per column,
    #meaning pairs that are less significant in the grand scheme of things
    #are potentially upweighted
    pVec <- c(pMat)
    #adjust
    pAVec <- p.adjust(pVec, method = "BH")
    #return to a matrix
    pAMat <- matrix(pAVec, nrow = nrow(pMat), ncol = ncol(pMat))
    return(pAMat)
})

#keep just significant correlations
sigCorDFBA <- map2(aPvalBA, corsBA, function(pAMat, corMat) { 
  #produce significant only matrix
  sig <- pAMat < 0.05
  #mask non-significant results
  sigCorMat <- ifelse(sig, corMat, NA)
  #add feature names
  rownames(sigCorMat) <- colnames(x)
  colnames(sigCorMat) <- colnames(y)
  sigInd <- which(!is.na(sigCorMat), arr.ind = T)
  #create a dataframe with variable pairs and correlation values
  corDFSig <- data.frame(Microbe = rownames(sigCorMat)[sigInd[,1]],
                         Metabolite = colnames(sigCorMat)[sigInd[,2]],
                         Correlation = sigCorMat[sigInd])
  #add genus names to microbes
  corDFSig <- corDFSig %>%
    inner_join(., tax[,c("ASV", "Genus")], by = c("Microbe" = "ASV")) %>%
    mutate(Relationship = ifelse(Correlation > 0, "Positive", "Negative"),
           Correlation = abs(Correlation)) %>%
    dplyr::select(Genus, Metabolite, Correlation, Relationship)
  return(corDFSig)
})

#produce graph objects
graphsBA <- map(sigCorDFBA, function(corDFSig) {
  #create igraph object
  g1 <- graph_from_data_frame(d = corDFSig, 
                              directed = F)
  #add microbe/metabolite annotation
  V(g1)$Feature <- ifelse(V(g1)$name %in% tax$Genus, "Microbe", "Metabolite")
  V(g1)$Degree <- degree(g1)
  return(g1)
})

#convert to tidygraphs objects
tidygraphsBA <- map(graphsBA, as_tbl_graph)

##### Output ####
#clr-transformed counts
write.csv(x, "processed/microbiome/CLRTransformedCounts_NetworkSubset.csv",
            row.names = T, quote = F)

#save data objects
save(met, pvalBA, aPvalBA, 
     corsBA, sigCorDFBA, 
     graphsBA, tidygraphsBA,
     file = "processed/network/NetworkBuild_Objs.rds")

#save network edges with correlations
scOut <- map2(sigCorDFBA, grpz, function(x,y) {
  out <- x %>%
    mutate(Group = y)
})
scOut <- bind_rows(scOut)
write.csv(scOut, "output/network/SigCorr_NetworkBuilds.csv",
          row.names = F, quote = F)

#combine all data to be saved as supplementary (nonsig included)
allDat <- map2(corsBA, pvalBA, function(cor, pval) {
  d1 <- cor %>%
    data.frame() %>%
    rownames_to_column(var = "ASV") %>%
    pivot_longer(-ASV, values_to = "SpearmanRho",
                 names_to = "Metabolite")
  d2 <- pval
  colnames(d2) <- colnames(cor)
  rownames(d2) <- rownames(cor)
  d2 <- d2 %>%
    data.frame() %>%
    rownames_to_column(var = "ASV") %>%
    pivot_longer(-ASV, values_to = "Pvalue",
                 names_to = "Metabolite")
  out <- inner_join(d1, d2, by = c("ASV", "Metabolite"))
})
allDat <- map2(allDat, aPvalBA, function(df, q) {
  d3 <- q
  colnames(d3) <- colnames(corsBA[[1]])
  rownames(d3) <- rownames(corsBA[[1]])
  d3 <- d3 %>%
    data.frame() %>%
    rownames_to_column(var = "ASV") %>%
    pivot_longer(-ASV, values_to = "AdjPval",
                 names_to = "Metabolite")
  out <- inner_join(df, d3, by = c("ASV", "Metabolite"))
}) 
#add genus labels in place of ASV
allDat <- map2(allDat, grpz, function(x, y) {
  x <- x %>%
    inner_join(., tax, by = "ASV") %>%
    select(Genus, Metabolite, SpearmanRho, Pvalue, AdjPval) %>%
    mutate(Network = y)
})
#flatten into dataframe
allDat <- bind_rows(allDat)
#save
write.csv(allDat, "output/network/NetworkBuild_AllData.csv",
          row.names = F, quote = F)

#visualisations
#groups
#filesystem-safe versions of grpz for directory names (matches dirs created above)
labs <- c("Recruitment", "FollowUp")
pmap(list(grpz, labs, tidygraphsBA), function(x, fileN, y) {
  plotGraph(y, labels = F, title = x)
  ggsave(paste0("output/network/", fileN, "/SigPairwiseSpearmanCorr_NetworkGraph_Unlabelled.png"))
  ggsave(paste0("output/network/", fileN, "/SigPairwiseSpearmanCorr_NetworkGraph_Unlabelled.pdf"))
  ggsave(paste0("output/network/", fileN, "/SigPairwiseSpearmanCorr_NetworkGraph_Unlabelled.jpg"))
})

#record samples used in analysis
reg <- read.csv("processed/Sample_AnalysisRegister.csv")
met %>%
  #keep only the samples that went into the network build
  filter(Sample %in% rownames(x)) %>%
  mutate(Network = "Present") %>%
  select(Sample, Patient.ID, Network) %>%
  full_join(reg, ., by = c("Sample", "Patient.ID")) %>%
  mutate(Network = ifelse(is.na(Network), "Absent", Network)) %>%
  write.csv(., "processed/Sample_AnalysisRegister.csv",
            row.names = F, quote = F)

#Session information
writeLines(capture.output(sessionInfo()),
           "output/network/BuildingNetworks_sessionInfo.txt")
