library(tidyverse)
library(magrittr)
library(readxl)
library(fgsea)
library(data.table)
library(MASS)

inputPath <- "/home/rstudio/workspace/Data Collection/"
intermediatePath <- "/home/rstudio/workspace/Data Collection/intermediate/"
outPath <- "/home/rstudio/workspace/Data Collection/output/"

source("/home/rstudio/repos/Comp_Dataset_Call.R")

### Read in Data
MS_BEADdel_ProF <- read.csv("/home/rstudio/workspace/Data Collection/QCed/MS/BEADdel/ProDat_QCed.csv", row.names = 1)
DEL_ProMap <- read_excel("/home/rstudio/workspace/Data Collection/QCed/MS/BEADdel/BEADdepletions_ProteinMapping.xlsx") %>% as.data.frame()
colnames(MS_BEADdel_ProF) <- sapply(colnames(MS_BEADdel_ProF), function(x){DEL_ProMap[which(DEL_ProMap$Protein.Group==x),"Genes"]}) ### Using the Entrez Gene Naming scheme

MS_NONdel_ProF <- read.csv("/home/rstudio/workspace/Data Collection/QCed/MS/NONdel/ProDat_QCed.csv", row.names = 1)
NON_ProMap <- read_excel("/home/rstudio/workspace/Data Collection/QCed/MS/NONdel/NONDEPLETED_ProteinMapping.xlsx") %>% as.data.frame()
colnames(MS_NONdel_ProF) <- sapply(colnames(MS_NONdel_ProF), function(x){NON_ProMap[which(NON_ProMap$Protein.Group==x),"Genes"]})

Olink_ProF <- read.csv("/home/rstudio/workspace/Data Collection/QCed/Olink/CSF/CSF_Customised_ProDataClean.csv", row.names = 1)

Clinic_MetaF_All <- read.csv("/home/rstudio/workspace/Data Collection/DA76/da76-wb-fleeting-cabbage-4050/Data/20251025_als_clinical_with_ambrosia_data_pseudonymized_FINAL.csv") %>% filter(CSF_OLINK_MANIFEST != "")
rownames(Clinic_MetaF_All) <- Clinic_MetaF_All$CSF_OLINK_MANIFEST

# Only the base line samples will be considered in Clinic_MetaF
Clinic_MetaF_All %<>% 
  mutate(
    # C9 vs Non-C9 classification
    C9vsNonC9 = ifelse(
      PATHOGENIC_VARIANT %in% c("KIF5A", "MAPT", "SBMA", "SOD1", "SPG7", "TBK1"),
      "NEGATIVE",
      PATHOGENIC_VARIANT
    ),
    
    # ALS vs HC
    ALSvsHC = ifelse(
      GROUP %in% c("AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL"),
      GROUP,
      NA
    ),
    
    # ALS vs DC
    ALSvsDC = ifelse(
      GROUP %in% c("AMYOTROPHIC LATERAL SCLEROSIS", "DISEASE CONTROL"),
      GROUP,
      NA
    ),
    
    # BULBAR vs SPINAL classification
    BULBARvSpinal = case_when(
      WEAKNESS_SITE == "BULBAR" ~ "BULBAR",
      WEAKNESS_SITE %in% c("GENERALISED", "LOWER LIMB", "RESPIRATORY", "TRUNCAL", "UPPER LIMB") ~ "SPINAL",
      TRUE ~ ""
    ),
    
    # Transformations
    ALSFRSR_RATE = log(ALSFRSR_RATE),
    ECAS_SCORE = boxcoxTranform(ECAS_SCORE, 134),
    
    # Factor conversions
    GROUP = factor(GROUP),
    ALSvsHC = factor(ALSvsHC),
    GENDER = factor(GENDER),
    WEAKNESS_SITE = factor(WEAKNESS_SITE),
    PATHOGENIC_VARIANT = factor(PATHOGENIC_VARIANT),
    COHORT = factor(COHORT),
    BULBARvSpinal = factor(BULBARvSpinal),
    C9vsNonC9 = factor(C9vsNonC9)
  )

Clinic_MetaF <- Clinic_MetaF_All %>% filter(LONGITUDINAL_ENCOUNTER_NUMBER == 1) 

### Merge the clinical meta and proteomic data
commonID1 <- intersect(rownames(MS_BEADdel_ProF), rownames(Clinic_MetaF))
MS_BEADdel_merged <- cbind(MS_BEADdel_ProF[commonID1,], Clinic_MetaF[commonID1,]) 

commonID2 <- intersect(rownames(MS_NONdel_ProF), rownames(Clinic_MetaF))
MS_NONdel_merged <- cbind(MS_NONdel_ProF[commonID2,], Clinic_MetaF[commonID2,])

commonID3 <- intersect(rownames(Olink_ProF), rownames(Clinic_MetaF))
Olink_merged <- cbind(Olink_ProF[commonID3,], Clinic_MetaF[commonID3,])

########################
########################
# DE Between MS & OLINK
########################
########################
# Direct compare the correlation coefficient
individualCor_BEADdel_Olink_Pro <- getIndividualRelation(MS_BEADdel_ProF, Olink_ProF, "Protein", "spearman")
individualCor_NONdel_Olink_Pro <- getIndividualRelation(MS_NONdel_ProF, Olink_ProF, "Protein", "spearman")
individualCor_BEADdel_NONdel_Pro <- getIndividualRelation(MS_BEADdel_ProF, MS_NONdel_ProF, "Protein", "spearman")

corPro_BEADdel_Olink <- names(individualCor_BEADdel_Olink_Pro)[which(abs(individualCor_BEADdel_Olink_Pro)>0.5)]
corPro_NONdel_Olink <- names(individualCor_NONdel_Olink_Pro)[which(abs(individualCor_NONdel_Olink_Pro)>0.5)]
corPro_BEADdel_NONdel <- names(individualCor_BEADdel_NONdel_Pro)[which(abs(individualCor_BEADdel_NONdel_Pro)>0.5)]

hist(individualCor_BEADdel_Olink_Pro, breaks=100)
hist(individualCor_NONdel_Olink_Pro, breaks=100)
hist(individualCor_BEADdel_NONdel_Pro, breaks=100)

# Enrichment test
geneSymuniverse <- names(individualCor_BEADdel_Olink_Pro)
geneList <- corPro_BEADdel_Olink

geneSymuniverse <- names(individualCor_NONdel_Olink_Pro)
geneList <- corPro_NONdel_Olink

geneSymuniverse <- names(individualCor_BEADdel_NONdel_Pro)
geneList <- corPro_BEADdel_NONdel

GOBPGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c5.go.bp.v2023.2.Hs.symbols.gmt"))
GOMFGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c5.go.mf.v2023.2.Hs.symbols.gmt"))
GOCCGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c5.go.cc.v2023.2.Hs.symbols.gmt"))
REACTOMEGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.reactome.v2023.2.Hs.symbols.gmt"))
KEGGGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt"))
BiocartaGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.biocarta.v2023.2.Hs.symbols.gmt"))
WikiPathwaysGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt"))

### Within each module perform the enrichment test based on the selected database
pathBase <- c("GOBP","GOMF","GOCC","REACTOME","KEGG","Biocarta","WikiPathways") ### user define which databases to look into
SigPath.List <- lapply(pathBase, function(x) vector("list", length = 1))

GeneSetList <- sapply(paste0(pathBase,"GeneSets"),function(x){get(x)})

gSk=1
for(GeneSetsX in GeneSetList){
  EnrichPath <- fgsea::fora(GeneSetsX, geneList,geneSymuniverse, minSize = 5, maxSize = 1000) 
  
  fwrite(EnrichPath, file=paste0(intermediatePath,pathBase[gSk],".BEADdelvsNONdel.txt"), sep="\t", sep2=c("", " ", ""))
  
  SigPath.List[[gSk]] <-  EnrichPath[which(EnrichPath$padj<0.05),] ### a list store the significant items for further plot
  
  names(SigPath.List[[gSk]])[1] <- "BEADdelvsNONdel"
  
  gSk <- gSk+1 
}

pdf(paste0(outPath, "Comp.pdf"))

### Bubble plot to show the significantly enriched pathways
files <- paste0(pathBase,".rdata")

for (i in 1:length(pathBase)){
  makeBubble(intermediatePath, pathBase[i], 5, 5, files[i])
}

individualAdjP_BEADdel_Olink_Samp <- getIndividualRelation(MS_BEADdel_ProF, Olink_ProF, "Sample", "spearman")
individualAdjP_NONdel_Olink_Samp <- getIndividualRelation(MS_NONdel_ProF, Olink_ProF, "Sample", "spearman")
individualAdjP_BEADdel_NONdel_Samp <- getIndividualRelation(MS_BEADdel_ProF, MS_NONdel_ProF, "Sample", "spearman")

par(mfrow=c(1,3))
hist(individualAdjP_BEADdel_Olink_Samp, breaks=100)
hist(individualAdjP_NONdel_Olink_Samp, breaks=100)
hist(individualAdjP_BEADdel_NONdel_Samp, breaks=100)

corSamp_BEADdel_Olink <- names(individualAdjP_BEADdel_Olink_Samp)[which(abs(individualAdjP_BEADdel_Olink_Samp)>0.7)]
corSamp_NONdel_Olink <- names(individualAdjP_NONdel_Olink_Samp)[which(abs(individualAdjP_NONdel_Olink_Samp)>0.7)]
corSamp_BEADdel_NONdel <- names(individualAdjP_BEADdel_NONdel_Samp)[which(abs(individualAdjP_BEADdel_NONdel_Samp)>0.7)]

# Compare with global, whether there is something important for the sample characteristics
CorSampReport <- Clinic_MetaF_All[corSamp_BEADdel_NONdel,]

# Define variables
categorical_vars <- c("ALSvsHC", "ALSvsDC", "BULBARvSpinal", "C9vsNonC9", "GENDER", "COHORT")
numeric_vars     <- c("AGE_AT_SYMPTOM_ONSET", "AGE_AT_SAMPLING", "BODY_MASS_INDEX", "ALSFRSR_RATE", "ECAS_SCORE")

all_vars <- c(categorical_vars, numeric_vars)

# Initialize result vector
CorSamp_vs_AllSamp <- setNames(rep(NA, length(all_vars)), all_vars)

# Compare categorical variables using Fisher's Exact Test
for(var in categorical_vars){
  CorSampVar <- CorSampReport %>% 
    filter(!is.na(var)) %>% 
    pull(var)
  
  CorSampAll <- Clinic_MetaF %>%
    filter(!is.na(var)) %>% 
    pull(var)
  
  AllSampVar <- Clinic_MetaF[,var]
  testF <- fisher.test(table(CorSampVar), table(CorSampAll))
  CorSamp_vs_AllSamp[var] <- testF$p.value
}

# Compare numeric variables using Wilcoxon rank-sum test (non-parametric)
for(var in numeric_vars){
  CorSampVar <- CorSampReport %>% 
    filter(!is.na(var)) %>% 
    pull(var)
  
  CorSampAll <- Clinic_MetaF %>%
    filter(!is.na(var)) %>% 
    pull(var)
  
  testKS <- ks.test(table(CorSampVar), table(CorSampAll))
  CorSamp_vs_AllSamp[var] <- testKS$p.value
}

CorSamp_vs_AllSamp <- p.adjust(CorSamp_vs_AllSamp, method="BH")

# View results
CorSamp_vs_AllSamp

# --- Plot distributions for significant variables ---

sig_vars <- names(CorSamp_vs_AllSamp)[CorSamp_vs_AllSamp < 0.05]

for(var in sig_vars){
  if(var %in% categorical_vars){
    # Bar plot for categorical variables
    df <- rbind(
      data.frame(value = CorSampReport[[var]], group = "Subset"),
      data.frame(value = Clinic_MetaF[[var]], group = "Full")
    )
    
    p <- ggplot(df, aes(x = value, fill = group)) +
      geom_bar(position = "dodge") +
      labs(title = paste("Distribution of", var),
           x = var, y = "Count") +
      theme_minimal()
    print(p)
    
  } else if(var %in% numeric_vars){
    # Histogram / density for numeric variables
    df <- rbind(
      data.frame(value = CorSampReport[[var]], group = "Subset"),
      data.frame(value = Clinic_MetaF[[var]], group = "Full")
    )
    
    p <- ggplot(df, aes(x = value, fill = group)) +
      geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
      labs(title = paste("Distribution of", var),
           x = var, y = "Count") +
      theme_minimal()
    print(p)
  }
}

dev.off()

#######################
#######################
# DE Between Platforms
#######################
#######################
# --------------------------
# Generate DE P value Table
# --------------------------
varList1 = c("ALSvsHC", "ALSvsDC", "BULBARvSpinal", "C9vsNonC9", "BODY_MASS_INDEX", "ALSFRSR_FINE_MOTOR", "ALSFRSR_GROSS_MOTOR", "ALSFRSR_RATE", "ECAS_SCORE")
varList2 = c("AGE_AT_SYMPTOM_ONSET", "AGE_AT_SAMPLING","GENDER", "COHORT")
varList = c(varList1,varList2)

### check whether the categorical variables are set as factor
for(id in c(1,2,3,4, 12,13)){
  print(class(Clinic_MetaF[,varList[id]]))
} 

result_tbl_adj_BEADdel <- get_DE_Pvalue_Table(MS_BEADdel_ProF, MS_BEADdel_merged, varList2, "/home/rstudio/workspace/Data Collection/output/", "BEADdel", varList1)
result_tbl_adj_NONdel <- get_DE_Pvalue_Table(MS_NONdel_ProF, MS_NONdel_merged, varList2, "/home/rstudio/workspace/Data Collection/output/", "NONdel", varList1)
result_tbl_adj_Olink <- get_DE_Pvalue_Table(Olink_ProF, Olink_merged, varList2, "/home/rstudio/workspace/Data Collection/output/", "Olink", varList1)

# Extract the common significantly expressed proteins
# Significance threshold
alpha <- 0.05

sig_BEADdel <- sig_NONdel <- sig_Olink <- list()
names(sig_BEADdel) <- names(sig_NONdel) <- names(sig_Olink)
for(var in varList1){
  sig_BEADdel[[var]] <- rownames(result_tbl_adj_BEADdel)[which(result_tbl_adj_BEADdel[,var] < 0.05)]
  sig_NONdel[[var]] <- rownames(result_tbl_adj_NONdel)[which(result_tbl_adj_NONdel[,var] < 0.05)]
  sig_Olink[[var]]  <- rownames(result_tbl_adj_Olink)[which(result_tbl_adj_Olink[,var] < 0.05)]
}
save(sig_BEADdel, sig_NONdel, sig_Olink, file=paste0(outPath, "sig.rdat"))

# Construct a table, rows as clinical variables, columns as three platforms, entries are significant proteins
SigProSummary <- data.frame(
  ClinicalVar = varList1,
  BEADdel = sapply(varList1, function(v) paste(sig_BEADdel[[v]], collapse = "/")),
  NONdel  = sapply(varList1, function(v) paste(sig_NONdel [[v]], collapse = "/")),
  Olink   = sapply(varList1, function(v) paste(sig_Olink  [[v]], collapse = "/")),
  stringsAsFactors = FALSE
)

write.csv(SigProSummary, file = paste0(outPath, "SigProSummary.csv"), row.names = FALSE)

# Make scatter plot of protein expression between Olink and MS for ALS vs HC
ALSvsHC_SigPro_List <- unique(c(sig_BEADdel[["ALSvsHC"]], sig_NONdel[["ALSvsHC"]], sig_Olink[["ALSvsHC"]]))

ALSvsHC_SigPro_Exist_Platform <- matrix(NA, nrow = length(ALSvsHC_SigPro_List), ncol=3)
rownames(ALSvsHC_SigPro_Exist_Platform) <- ALSvsHC_SigPro_List
colnames(ALSvsHC_SigPro_Exist_Platform) <- c("Olink", "MS_BEADdel", "MS_NONdel")

for(ProS in ALSvsHC_SigPro_List){ 
  ALSvsHC_SigPro_Exist_Platform[ProS, "Olink"] <- ifelse(ProS %in% colnames(Olink_ProF), "Y", "N") 
  ALSvsHC_SigPro_Exist_Platform[ProS, "MS_BEADdel"] <- ifelse(ProS %in% colnames(MS_BEADdel_ProF), "Y", "N") 
  ALSvsHC_SigPro_Exist_Platform[ProS, "MS_NONdel"] <- ifelse(ProS %in% colnames(MS_NONdel_ProF), "Y", "N")
}

write.csv(ALSvsHC_SigPro_Exist_Platform, file = paste0(outPath, "ALSvsHC_SigPro_Exist_Platform.csv"))

pdf(paste0(outPath, "ALSvsHC_SigPro_In_Platform.pdf"))

# ---------- scatter plot between two platforms ----------
plot_scatter <- function(df1, df2, protein, label1, label2) {
  # Find common samples
  comSamp <- intersect(rownames(df1), rownames(df2))

  plot_df <- tibble(
    x = df1[comSamp, protein],
    y = df2[comSamp, protein]
  )
  
  ggplot(plot_df, aes(x, y)) +
    geom_point(alpha = 0.7) +
    labs(
      title = paste("Protein:", protein),
      x = label1,
      y = label2
    ) +
    theme_bw(base_size = 14)
}

for (ProS in ALSvsHC_SigPro_List) {
  
  has_Olink   <- ProS %in% colnames(Olink_ProF)
  has_BEADdel <- ProS %in% colnames(MS_BEADdel_ProF)
  has_NONdel  <- ProS %in% colnames(MS_NONdel_ProF)
  
  # Olink vs BEADdel
  if (has_Olink && has_BEADdel) {
    print(
      plot_scatter(
        Olink_ProF, MS_BEADdel_ProF,
        protein = ProS,
        label1 = "Olink",
        label2 = "MS_BEADdel"
      )
    )
  }
  
  # Olink vs NONdel
  if (has_Olink && has_NONdel) {
    print(
      plot_scatter(
        Olink_ProF, MS_NONdel_ProF,
        protein = ProS,
        label1 = "Olink",
        label2 = "MS_NONdel"
      )
    )
  }
}

dev.off()

# compute overlaps per var and build a summary table
summary_list <- lapply(varList1, function(var) {
  s_bead <- sig_BEADdel[[var]]
  s_non  <- sig_NONdel[[var]]
  s_ol   <- sig_Olink[[var]]
  
  overlap_BEADdel_NONdel <- intersect(s_bead, s_non)
  overlap_BEADdel_Olink  <- intersect(s_bead, s_ol)
  overlap_NONdel_Olink   <- intersect(s_non, s_ol)
  overlap_all3           <- Reduce(intersect, list(s_bead, s_non, s_ol))
  
  data.frame(
    var = var,
    Comparison = c("BEADdel vs NONdel", "BEADdel vs Olink", "NONdel vs Olink", "All three"),
    nProteins = c(
      length(overlap_BEADdel_NONdel),
      length(overlap_BEADdel_Olink),
      length(overlap_NONdel_Olink),
      length(overlap_all3)
    ),
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_list)

# Function to compute overlaps for a single var
compute_overlap <- function(var) {
  bead  <- sig_BEADdel[[var]]
  nond  <- sig_NONdel[[var]]
  olink <- sig_Olink[[var]]
  
  list(
    BEADdel_vs_NONdel   = intersect(bead, nond),
    BEADdel_vs_Olink    = intersect(bead, olink),
    NONdel_vs_Olink     = intersect(nond, olink),
    All_three           = Reduce(intersect, list(bead, nond, olink))
  )
}

# Compute overlaps for each var
overlap_list <- lapply(varList1, compute_overlap)
names(overlap_list) <- varList1

# Write summary
write.csv(summary_df, file = paste0(outPath, "overlap_summary_counts.csv"), 
          row.names = FALSE)

# Also save element names for full inspection
save(overlap_list, file = paste0(outPath, "overlap_elements.rdat"))

summary_df

# -------------------------------------------
# Significantly expressed protein expression 
# -------------------------------------------


# ------------
# Volcano Plot
# ------------

# -------------------
# Enrichment Testing
# -------------------
# fgseaPath <- fgsea::fgsea(pathways = all_gene_sets, stats=ranks, scoreType = whichType, eps = 0.0) ### we will not limit the size (minSize and maxSize) here, just output everything, but would select the proper size for further plotting. 

