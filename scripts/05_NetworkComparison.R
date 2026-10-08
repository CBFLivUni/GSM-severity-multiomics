#21st Jan 2025
#Lauren Mee
#Produce and compare microbiome/metabolome networks
#built on Spearman's correlations

##### Set up #####
#set seed
set.seed(1152)

#libraries
libs <- c("ggraph", "igraph", "tidygraph", #for network manipulation
          "UpSetR", #for venn-like graphs
          "phyloseq", #for handling phyloseq objects
          "vegan", #for mantel test
          "ggrepel", "ggpubr", "tidyverse") #general data manipulation, plots


for (pkg in libs) {
  library(pkg, character.only = T)
}

#load data
#network build objects
load("processed/network/NetworkBuild_Objs.rds")
#microbiome phyloseq objects
load("processed/microbiome/PhyloSeqObjs.rds")

#microbiome data (clr-transformed)
x <- read.csv("processed/microbiome/CLRTransformedCounts_NetworkSubset.csv",
              row.names = 1)

#nmr (already scaled and log-transformed)
y <- read.csv("processed/nmr/logPQN_NMR.csv", row.names = 1) %>%
  data.frame() %>%
  #keep only samples in both datasets
  filter(rownames(.) %in% rownames(x))

#get top metabolite bins that correlated with NGAT
metNGAT <- read.csv("output/nmr/lme4_allresults.csv")

#taxonomy
tax <- tax_table(phyGen) %>%
  data.frame() %>%
  rownames_to_column(var = "ASV") 

#misc
#set colour palette
vanPal <- c("#5d83df", "#df5d83", "#40B0A6",
            "#00ffff",   "black", "#cdb126",
            "#8dbab1", "#F585F5", "#65348C")

##### Functions ####
#extracts degrees/betweenness from igraphs objects
#and compares the two group networks
#produces a violin plot and a dataframe with metrics
#for later perusal
compTopCraft <- function(graphs,
                         metric = "Degree" #Degree or Betweenness
) {
  lab1 <- "Recruitment"
  lab2 <- "Follow up"
  fileN <- "RecruitmentvFollowUp"
  if (metric == "Degree") {
    one <- data.frame(degree(graphs[[1]]))
    two <- data.frame(degree(graphs[[2]]))
  }
  if (metric == "Betweenness") {
    one <- data.frame(betweenness(graphs[[1]]))
    two <- data.frame(betweenness(graphs[[2]]))
  }
  bva <- map2(list(one, two), c(lab1, lab2), function(df, l) {
    df <- df %>%
      rownames_to_column(var = "Node") %>%
      mutate(Network = l)
    names(df)[2] <- "Metric"
    return(df)
  })
  bva <- bind_rows(bva) %>%
    mutate(Network = factor(Network, levels = c(lab1, lab2)))
  ggplot(bva, aes(x = Network, y = Metric)) +
    theme_bw(base_size = 13) +
    geom_violin(aes(fill = Network), alpha = 0.75) + 
    geom_boxplot( width = 0.1) + 
    scale_fill_manual(values = vanPal[c(2:3)]) + 
    guides(colour = 'none', fill = "none") +
    labs(y = metric) 
  ggsave(paste0("output/network/", fileN, "_Network", metric, ".png"))
  names(bva)[2] <- metric
  write.csv(bva, paste0("output/network/", fileN, "_Network", metric, ".csv"),
              row.names = F, quote = F)
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

##### Data preparation #####
#combine objects to work together
grpz <- c("Recruitment", "Follow up") 

#produce dataframes with pvals and correlations per pair 
statsDF <- map2(aPvalBA, corsBA, function(pvals, cors) {
  #apply metabolite/ASV IDS
  rownames(pvals) <- rownames(cors)
  colnames(pvals) <- colnames(cors)
  #convert pval matrix into long form dataframe
  pvals <- pvals %>%
    data.frame() %>%
    rownames_to_column(var = "ASV") %>%
    pivot_longer(-ASV,
                 names_to = "Metabolite",
                 values_to = "adjPval")
  #convert correlation matrix into long form dataframe
  cors <- cors %>%
    data.frame() %>%
    rownames_to_column(var = "ASV") %>%
    pivot_longer(-ASV,
                 names_to = "Metabolite",
                 values_to = "Correlation")
  #combine both, and add adjusted pvalues
  df <- cors %>%
    inner_join(., pvals, by = c("ASV", "Metabolite")) 
  return(df)
})


##### Hubs #####
#determine hubs in networks made from all pairwise connections
hubs <- map2(grpz, graphsBA, function(x, y) {
  #figure out threshold for top 5% most connected nodes
  hubThres <- quantile(V(y)$Degree, 0.95)
  #extract top 5% highest degree nodes
  hubs <- V(y)[V(y)$Degree >= hubThres] 
  #get just hub names
  hubs <- rownames(as.data.frame(hubs))
  #write text files with list of hubs
  write.table(as.data.frame(hubs),
              paste0("output/network/", gsub(" u", "U", x),
                     "/Hub_Nodes_FullGraph_",
                     round(hubThres, 3), ".txt"),
              col.names = F, quote = F, row.names = F)
  return(hubs)
})
names(hubs) <- grpz

##### Network Visualisations ####
#add NGAT correlation
#consider only significant metabolite bins in reduced model
ngatSig <- metNGAT[metNGAT$NGAT_adjP < 0.05 &
                     metNGAT$Model == "Reduced",]
moi <- ngatSig$Metabolite

#visualise networks themselves
#add extra parameters required for aesthetic reasons
tidyList <- map2(tidygraphsBA, hubs, function(tg, h) {
  tg <- tg %>%
    #add hubs
    mutate(Label = ifelse(name %in% h, name, ""),
           Type = ifelse(name %in% h, "Hub", "Other"),
           #add if DM metabolite or not 
           Label = ifelse(name %in% moi, name, Label),
           Type = ifelse(name %in% moi, "DM Met", Type),
           #add Lactobacillus
           Label = ifelse(name == "Lactobacillus", "Lactobacillus", Label),
           #add Streptococcus as antagonising microbe
           Label = ifelse(name == "Streptococcus", "Streptococcus", Label),
           #fix Metabolite labels
           Label = map_vec(Label, fixMet),
           #allow for italicisation of microbe taxa names
           Face = ifelse(Feature == "Microbe", "bold.italic", "bold"))
  return(tg)
})

#apply
tidyList <- map(tidyList, function(tg) {
  #filter to significant metabolites only
  #keeping all in makes network graphs way too busy
  #apply
  tg <- tg %>%
    #add NGAT label if name is in most significant NGAT metabolites
    mutate(Label = ifelse(name %in% moi, 
                          map_vec(name, fixMet), Label),
           Label = map_vec(Label, function(m) {
             ifelse(grepl("^Unknown", m) | grepl("^Glucose", m), 
                    m, gsub("-[0-9]+$", "", m,
                            perl = TRUE))
           }),
      #add NGAT correlation direction as Type
      Type = ifelse(name %in% ngatSig$Metabolite[ngatSig$NGAT_estimate > 0], 
                    "NGAT+",
        ifelse(name %in% ngatSig$Metabolite[ngatSig$NGAT_estimate < 0], 
               "NGAT-", Type)))
  return(tg)
})

#plot (all labels)
netListAll <- map2(tidyList, grpz, function(tg, l) {
  filename <- gsub(" u", "U", l)
  p <- ggraph(tg, layout = "kk") + 
    geom_edge_link(aes(colour = Relationship,
                       width = Correlation)) + 
    scale_edge_width(range = c(0.5, 3)) + 
    geom_node_point(aes(size = Degree, shape = Feature,
                        colour = Type)) + 
    scale_size_continuous(range = c(3, 11)) + 
    scale_edge_colour_manual(values = vanPal[c(9,7)], guide = "legend") +
    scale_colour_manual(values = c("Other" = "black",
                                   "Hub" =  "#008080",
                                   "NGAT+" = "#d8783d",
                                   "NGAT-" = "#0179c2"), 
                        guide = "legend") +
    scale_shape_manual(values = c(15, 16)) +
    theme_void() +
    theme(legend.margin = margin(t = 5, r = 10, b = 5, l = 5)) + 
    guides(size = "none", fill = "none") + 
    geom_node_label(aes(label = name,
                        fontface = Face,
                        fill = Type),
                    repel = T, box.padding = 1,
                    alpha = 0.8) +
    scale_fill_manual(values = c("Other" = "darkgrey",
                                 "Hub" =  "#00a7a7",
                                 "NGAT+" = "#f2d3c0",
                                 "NGAT-" = "#b0e1ff")) +
    ggtitle(l) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    face = "bold",
                                    size = 16))
  ggsave(paste0("output/network/", filename, 
                "/SigPairwiseSpearmanNetwork_AllLabels.png"))
  ggsave(paste0("output/network/", filename, 
                "/SigPairwiseSpearmanNetwork_AllLabels.pdf"))
  return(p)
})

#plot list just moi and hubs.
netList <- map2(tidyList, grpz, function(tg, l) {
  p <- ggraph(tg, layout = "kk") + 
    geom_edge_link(aes(colour = Relationship,
                       width = Correlation)) + 
    scale_edge_width(range = c(0.5, 3)) + 
    geom_node_point(aes(size = Degree, shape = Feature,
                        colour = Type)) + 
    scale_size_continuous(range = c(3, 11)) + 
    scale_edge_colour_manual(values = vanPal[c(9,7)], guide = "legend") +
    scale_colour_manual(values = c("Other" = "black",
                                   "Hub" =  "#008080",
                                   "NGAT+" = "#d8783d",
                                   "NGAT-" = "#0179c2"), 
                        guide = "legend") +
    scale_shape_manual(values = c(15, 16)) +
    theme_void() +
    theme(legend.margin = margin(t = 5, r = 10, b = 5, l = 5)) + 
    guides(size = "none", fill = "none") + 
    geom_node_label(aes(label = Label,
                        fontface = Face,
                        fill = Type),
                    repel = T, box.padding = 1.3,
                    alpha = 0.8, size = 3) +
    scale_fill_manual(values = c("Other" = "darkgrey",
                                 "Hub" =  "#00a7a7",
                                 "NGAT+" = "#f2d3c0",
                                 "NGAT-" = "#b0e1ff")) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    face = "bold",
                                    size = 16),
          legend.position = "left")
  ggsave(paste0("output/network/", gsub(" u", "U", l),
                "/SigPairwiseSpearmanNetwork_FocusedLabels.png"))
  ggsave(paste0("output/network/", gsub(" u", "U", l),
                "/SigPairwiseSpearmanNetwork_FocusedLabels.pdf"))
  return(p)
})

##### Comparing Topology ####
#general
genTop <- map(graphsBA, function(g) {
  nodes <- vcount(g)
  edges <- ecount(g)
  meandeg <- mean(degree(g))
  meddeg <- median(degree(g))
  sddeg <- sd(degree(g))
  meanbet <- mean(betweenness(g))
  medbet <- median(betweenness(g))
  sdbet <- sd(betweenness(g))
  out <- c(nodes, edges, meddeg, meandeg, sddeg, medbet, meanbet, sdbet)
  return(out)
})
genTop <- as.data.frame(do.call(rbind, genTop))
names(genTop) <- c("NoNodes", "NoEdges", "MedianDegrees",
                   "MeanDegrees", "SDDegrees", "MedianBetweenness",
                   "MeanBetweenness", "SDBetweenness")
rownames(genTop) <- grpz
#compute SE
genTop <- genTop %>%
  mutate(SEDegrees = SDDegrees / sqrt(NoNodes),
         SEBetweenness = SDBetweenness / sqrt(NoNodes))


#shared nodes and edges
nodes <- map(graphsBA, function(g) {
  return(V(g)$name)
})
edgez <- map(graphsBA, function(g) {
  e <- as_edgelist(g)
  return(paste(e[,1], e[,2], sep = "-"))
})

#groups
bVa.Node <- intersect(nodes[[1]], nodes[[2]])
bVa.SharedNodesNo <- length(bVa.Node)
if (bVa.SharedNodesNo == 0) {
  bVa.Node <- NA
}

bVa.Edge <- intersect(edgez[[1]], edgez[[2]])
bVa.SharedEdgesNo <- length(bVa.Edge)
if (bVa.SharedEdgesNo == 0) {
  bVa.Edge <- NA
}

#compare matrix
compTop <- data.frame(SharedNodeNo = bVa.SharedNodesNo,
           SharedNodes = paste(bVa.Node, collapse = ", "),
           SharedEdgeNo = bVa.SharedEdgesNo,
           SharedEdges = paste(bVa.Edge, collapse = ", "))
rownames(compTop) <- "RecruitmentVsFollowUp"

#degrees / betweenness
compTopCraft(graphsBA)
compTopCraft(graphsBA, "Betweenness")

##### Mantel Test #####
#test global congruence of each network
#prepare to subset to get matrices to test
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

#ensure rownames are in the same orders
yBA <- map2(xBA, yBA, function(x, y) {
  y <- y[rownames(x), ]
  return(y)
})

#double check
map2(xBA, yBA, function(x, y) {
  return(print(all(rownames(x) == rownames(y))))
})

#convert into distance matrices
#CLR is compositional so Euclidean is more apt
#than bray or the other distances I'd use with untransformed
#counts
dist.x <- map(xBA, ~dist(.x, method = "euclidean"))
dist.y <- map(yBA, ~dist(.x, method = "euclidean"))

#spearman to keep test non-parametric
manOut <- map2(dist.x, dist.y, ~mantel(.x, .y, method = "spearman", permutations = 10000))
names(manOut) <- grpz

##### Output ####
#potential figure
metMicroNetworks <- ggarrange(plotlist = netList, common.legend = T, 
          legend = "left", labels = "AUTO", nrow = 2) + 
  theme(plot.margin = margin(12, 12, 12, 12))
ggsave("output/figures/met-micro-networks.pdf",
       height = 11, width = 10, units = "in")
ggsave("output/figures/met-micro-networks.png",
       height = 11, width = 10, units = "in")
ggsave("output/figures/met-micro-networks.jpg",
       height = 11, width = 10, units = "in")
#save rds
save(metMicroNetworks, file = "output/figures/met-micro-networks.rds")

#topography comparison
write.csv(genTop,
            "output/network/GeneralTopologyStats.csv",
            row.names = T, quote = F)
write.csv(compTop,
            "output/network/SharedEdgesNodes.csv",
            row.names = T, quote = F)

#Mantel Tests output
writeLines(capture.output(print(manOut)),
           "output/network/MantelTestResults.txt")

#Session information
writeLines(capture.output(sessionInfo()),
           "output/network/Comparisons_sessionInfo.txt")
