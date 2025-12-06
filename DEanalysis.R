library(tidyverse)
library(magrittr)
library(readxl)
library(fgsea)
library(data.table)
library(MASS)
library(lme4)
library(nnet)
library(MASS)
library(ordinal)

inputPath <- "/home/rstudio/workspace/Data Collection/"
intermediatePath <- "/home/rstudio/workspace/Data Collection/intermediate/"
outPath <- "/home/rstudio/workspace/Data Collection/output_baseline/"

source("/home/rstudio/repos/DE_functions_call.R")

### Read in Data
MS_BEADdel_ProF <- read.csv("/home/rstudio/workspace/Data Collection/QCed/MS/BEADdel/ProDat_QCed.csv", row.names = 1)
# "/home/rstudio/workspace/Data Collection/oxf-qc2-als-mass-spectometry-wb-genial-hazelnut-7025/Data/output/BEADdel/SAMPLE/QCed_Frame.csv"
DEL_ProMap <- read_excel("/home/rstudio/workspace/Data Collection/QCed/MS/BEADdel/BEADdepletions_ProteinMapping.xlsx") %>% as.data.frame()
colnames(MS_BEADdel_ProF) <- sapply(colnames(MS_BEADdel_ProF), function(x){DEL_ProMap[which(DEL_ProMap$Protein.Group==x),"Genes"]}) ### Using the Entrez Gene Naming scheme

MS_NONdel_ProF <- read.csv("/home/rstudio/workspace/Data Collection/QCed/MS/NONdel/ProDat_QCed.csv", row.names = 1)
# "/home/rstudio/workspace/Data Collection/oxf-qc2-als-mass-spectometry-wb-genial-hazelnut-7025/Data/output/NONdel/SAMPLE/QCed_Frame.csv"
NON_ProMap <- read_excel("/home/rstudio/workspace/Data Collection/QCed/MS/NONdel/NONDEPLETED_ProteinMapping.xlsx") %>% as.data.frame()
colnames(MS_NONdel_ProF) <- sapply(colnames(MS_NONdel_ProF), function(x){NON_ProMap[which(NON_ProMap$Protein.Group==x),"Genes"]})

###  CV% < 0.8
# Olink_ProF <- read.csv("/home/rstudio/workspace/Data Collection/QCed/Olink/CSF/CSF_Customised_ProDataClean.csv", row.names = 1)
Olink_ProF <- read.csv("/home/rstudio/workspace/Data Collection/Data-CSF/oxf-da28-opdc-als-multiple-cohorts-olink-csf/CSF_Customised_ProDataClean.csv", sep= ",", row.names = 1)

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
    
    # DC vs HC
    DCvsHC = ifelse(
      GROUP %in% c("DISEASE CONTROL", "HEALTHY CONTROL"),
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
    
    # SubjectID for random effect regression
    SubjectID_Random = CSF_OLINK_MANIFEST,
    
    # Factor conversions
    GROUP = factor(GROUP),
    ALSvsHC = factor(ALSvsHC),
    ALSvsDC = factor(ALSvsDC),
    DCvsHC = factor(DCvsHC),
    GENDER = factor(GENDER),
    WEAKNESS_SITE = factor(WEAKNESS_SITE),
    PATHOGENIC_VARIANT = factor(PATHOGENIC_VARIANT),
    COHORT = factor(COHORT),
    BULBARvSpinal = factor(BULBARvSpinal),
    C9vsNonC9 = factor(C9vsNonC9),
    SubjectID_Random = factor(SubjectID_Random)
  )

Clinic_MetaF <- Clinic_MetaF_All %>% filter(LONGITUDINAL_ENCOUNTER_NUMBER == 1) 

### Merge the clinical meta and proteomic data
commonID1 <- intersect(rownames(MS_BEADdel_ProF), rownames(Clinic_MetaF))
MS_BEADdel_merged <- cbind(
  MS_BEADdel_ProF[commonID1, , drop=FALSE],
  Clinic_MetaF[commonID1, , drop=FALSE]
)
MS_BEADdel_ProF_matchClinic <- MS_BEADdel_ProF[commonID1, , drop = FALSE]

commonID2 <- intersect(rownames(MS_NONdel_ProF), rownames(Clinic_MetaF))
MS_NONdel_merged <- cbind(
  MS_NONdel_ProF[commonID2, , drop = FALSE],
  Clinic_MetaF[commonID2, , drop = FALSE]
)
MS_NONdel_ProF_matchClinic <- MS_NONdel_ProF[commonID2, , drop = FALSE]

commonID3 <- intersect(rownames(Olink_ProF), rownames(Clinic_MetaF))
Olink_merged <- cbind(
  Olink_ProF[commonID3, , drop = FALSE],
  Clinic_MetaF[commonID3, , drop = FALSE]
)
Olink_ProF_matchClinic <- Olink_ProF[commonID3, ,drop = FALSE]

########################
########################
# DE Between MS & OLINK
########################
########################
# Direct compare the correlation coefficient
# individualCor_BEADdel_Olink_Pro <- getIndividualRelation(MS_BEADdel_ProF, Olink_ProF, "Protein", "spearman")
# individualCor_NONdel_Olink_Pro <- getIndividualRelation(MS_NONdel_ProF, Olink_ProF, "Protein", "spearman")
# individualCor_BEADdel_NONdel_Pro <- getIndividualRelation(MS_BEADdel_ProF, MS_NONdel_ProF, "Protein", "spearman")
# 
# corPro_BEADdel_Olink <- names(individualCor_BEADdel_Olink_Pro)[which(abs(individualCor_BEADdel_Olink_Pro)>0.5)]
# corPro_NONdel_Olink <- names(individualCor_NONdel_Olink_Pro)[which(abs(individualCor_NONdel_Olink_Pro)>0.5)]
# corPro_BEADdel_NONdel <- names(individualCor_BEADdel_NONdel_Pro)[which(abs(individualCor_BEADdel_NONdel_Pro)>0.5)]
# 
# hist(individualCor_BEADdel_Olink_Pro, breaks=100)
# hist(individualCor_NONdel_Olink_Pro, breaks=100)
# hist(individualCor_BEADdel_NONdel_Pro, breaks=100)
# 
# # Enrichment test
# geneSymuniverse <- names(individualCor_BEADdel_Olink_Pro)
# geneList <- corPro_BEADdel_Olink
# 
# geneSymuniverse <- names(individualCor_NONdel_Olink_Pro)
# geneList <- corPro_NONdel_Olink
# 
# geneSymuniverse <- names(individualCor_BEADdel_NONdel_Pro)
# geneList <- corPro_BEADdel_NONdel
# 
# GOBPGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c5.go.bp.v2023.2.Hs.symbols.gmt"))
# GOMFGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c5.go.mf.v2023.2.Hs.symbols.gmt"))
# GOCCGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c5.go.cc.v2023.2.Hs.symbols.gmt"))
# REACTOMEGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.reactome.v2023.2.Hs.symbols.gmt"))
# KEGGGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt"))
# BiocartaGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.biocarta.v2023.2.Hs.symbols.gmt"))
# WikiPathwaysGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt"))
# 
# ### Within each module perform the enrichment test based on the selected database
# pathBase <- c("GOBP","GOMF","GOCC","REACTOME","KEGG","Biocarta","WikiPathways") ### user define which databases to look into
# SigPath.List <- lapply(pathBase, function(x) vector("list", length = 1))
# 
# GeneSetList <- sapply(paste0(pathBase,"GeneSets"),function(x){get(x)})
# 
# gSk=1
# for(GeneSetsX in GeneSetList){
#   EnrichPath <- fgsea::fora(GeneSetsX, geneList,geneSymuniverse, minSize = 5, maxSize = 1000) 
#   
#   fwrite(EnrichPath, file=paste0(intermediatePath,pathBase[gSk],".BEADdelvsNONdel.txt"), sep="\t", sep2=c("", " ", ""))
#   
#   SigPath.List[[gSk]] <-  EnrichPath[which(EnrichPath$padj<0.05),] ### a list store the significant items for further plot
#   
#   names(SigPath.List[[gSk]])[1] <- "BEADdelvsNONdel"
#   
#   gSk <- gSk+1 
# }
# 
# pdf(paste0(outPath, "Comp.pdf"))
# 
# ### Bubble plot to show the significantly enriched pathways
# files <- paste0(pathBase,".rdata")
# 
# for (i in 1:length(pathBase)){
#   makeBubble(intermediatePath, pathBase[i], 5, 5, files[i])
# }
# 
# individualAdjP_BEADdel_Olink_Samp <- getIndividualRelation(MS_BEADdel_ProF, Olink_ProF, "Sample", "spearman")
# individualAdjP_NONdel_Olink_Samp <- getIndividualRelation(MS_NONdel_ProF, Olink_ProF, "Sample", "spearman")
# individualAdjP_BEADdel_NONdel_Samp <- getIndividualRelation(MS_BEADdel_ProF, MS_NONdel_ProF, "Sample", "spearman")
# 
# par(mfrow=c(1,3))
# hist(individualAdjP_BEADdel_Olink_Samp, breaks=100)
# hist(individualAdjP_NONdel_Olink_Samp, breaks=100)
# hist(individualAdjP_BEADdel_NONdel_Samp, breaks=100)
# 
# corSamp_BEADdel_Olink <- names(individualAdjP_BEADdel_Olink_Samp)[which(abs(individualAdjP_BEADdel_Olink_Samp)>0.7)]
# corSamp_NONdel_Olink <- names(individualAdjP_NONdel_Olink_Samp)[which(abs(individualAdjP_NONdel_Olink_Samp)>0.7)]
# corSamp_BEADdel_NONdel <- names(individualAdjP_BEADdel_NONdel_Samp)[which(abs(individualAdjP_BEADdel_NONdel_Samp)>0.7)]
# 
# # Compare with global, whether there is something important for the sample characteristics
# CorSampReport <- Clinic_MetaF_All[corSamp_BEADdel_NONdel,]
# 
# # Define variables
# categorical_vars <- c("ALSvsHC", "ALSvsDC", "BULBARvSpinal", "C9vsNonC9", "GENDER", "COHORT")
# numeric_vars     <- c("AGE_AT_SYMPTOM_ONSET", "AGE_AT_SAMPLING", "BODY_MASS_INDEX", "ALSFRSR_RATE", "ECAS_SCORE")
# 
# all_vars <- c(categorical_vars, numeric_vars)
# 
# # Initialize result vector
# CorSamp_vs_AllSamp <- setNames(rep(NA, length(all_vars)), all_vars)
# 
# # Compare categorical variables using Fisher's Exact Test
# for(var in categorical_vars){
#   CorSampVar <- CorSampReport %>% 
#     filter(!is.na(var)) %>% 
#     pull(var)
#   
#   CorSampAll <- Clinic_MetaF %>%
#     filter(!is.na(var)) %>% 
#     pull(var)
#   
#   AllSampVar <- Clinic_MetaF[,var]
#   testF <- fisher.test(table(CorSampVar), table(CorSampAll))
#   CorSamp_vs_AllSamp[var] <- testF$p.value
# }
# 
# # Compare numeric variables using Wilcoxon rank-sum test (non-parametric)
# for(var in numeric_vars){
#   CorSampVar <- CorSampReport %>% 
#     filter(!is.na(var)) %>% 
#     pull(var)
#   
#   CorSampAll <- Clinic_MetaF %>%
#     filter(!is.na(var)) %>% 
#     pull(var)
#   
#   testKS <- ks.test(table(CorSampVar), table(CorSampAll))
#   CorSamp_vs_AllSamp[var] <- testKS$p.value
# }
# 
# CorSamp_vs_AllSamp <- p.adjust(CorSamp_vs_AllSamp, method="BH")
# 
# # View results
# CorSamp_vs_AllSamp
# 
# # --- Plot distributions for significant variables ---
# 
# sig_vars <- names(CorSamp_vs_AllSamp)[CorSamp_vs_AllSamp < 0.05]
# 
# for(var in sig_vars){
#   if(var %in% categorical_vars){
#     # Bar plot for categorical variables
#     df <- rbind(
#       data.frame(value = CorSampReport[[var]], group = "Subset"),
#       data.frame(value = Clinic_MetaF[[var]], group = "Full")
#     )
#     
#     p <- ggplot(df, aes(x = value, fill = group)) +
#       geom_bar(position = "dodge") +
#       labs(title = paste("Distribution of", var),
#            x = var, y = "Count") +
#       theme_minimal()
#     print(p)
#     
#   } else if(var %in% numeric_vars){
#     # Histogram / density for numeric variables
#     df <- rbind(
#       data.frame(value = CorSampReport[[var]], group = "Subset"),
#       data.frame(value = Clinic_MetaF[[var]], group = "Full")
#     )
#     
#     p <- ggplot(df, aes(x = value, fill = group)) +
#       geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
#       labs(title = paste("Distribution of", var),
#            x = var, y = "Count") +
#       theme_minimal()
#     print(p)
#   }
# }
# 
# dev.off()

#######################
#######################
# DE Between Platforms
#######################
#######################
# --------------------------
# Generate DE P value Table
# --------------------------
varList1 = c("ALSvsHC", "DCvsHC", "ALSvsDC", "BULBARvSpinal", "C9vsNonC9", "BODY_MASS_INDEX", "ALSFRSR_FINE_MOTOR", "ALSFRSR_GROSS_MOTOR", "ALSFRSR_RATE", "ECAS_SCORE")
varList2 = c("AGE_AT_SYMPTOM_ONSET", "AGE_AT_SAMPLING","GENDER", "COHORT")
varList = c(varList1,varList2)

### check whether the categorical variables are set as factor
for(id in c(1,2,3,4,5,13,14)){
  print(class(Clinic_MetaF[,varList[id]]))
} 

result_tbl_adj_BEADdel <- get_DE_Pvalue_Table(MS_BEADdel_ProF_matchClinic, MS_BEADdel_merged, varList2, "/home/rstudio/workspace/Data Collection/output_baseline/", "BEADdel", varList1)
result_tbl_adj_NONdel <- get_DE_Pvalue_Table(MS_NONdel_ProF_matchClinic, MS_NONdel_merged, varList2, "/home/rstudio/workspace/Data Collection/output_baseline/", "NONdel", varList1)
result_tbl_adj_Olink <- get_DE_Pvalue_Table(Olink_ProF_matchClinic, Olink_merged, varList2, "/home/rstudio/workspace/Data Collection/output_baseline/", "Olink", varList1)

# Extract the common significantly expressed proteins
# Significance threshold
alpha <- 0.05

sig_BEADdel <- sig_NONdel <- sig_Olink <- list()
names(sig_BEADdel) <- names(sig_NONdel) <- names(sig_Olink)
for(var in varList1){
  sig_BEADdel[[var]] <- rownames(result_tbl_adj_BEADdel)[order(result_tbl_adj_BEADdel[, var])][which(result_tbl_adj_BEADdel[order(result_tbl_adj_BEADdel[, var]), var] < 0.05)]
  
  sig_NONdel[[var]] <- rownames(result_tbl_adj_NONdel)[order(result_tbl_adj_NONdel[, var])][which(result_tbl_adj_NONdel[order(result_tbl_adj_NONdel[, var]), var] < 0.05)]
  
  sig_Olink[[var]] <- rownames(result_tbl_adj_Olink)[order(result_tbl_adj_Olink[, var])][which(result_tbl_adj_Olink[order(result_tbl_adj_Olink[, var]), var] < 0.05)]
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

#----------------------
#----------------------
# Special for ALSvsHC
#----------------------
#----------------------
# Make scatter plot of protein expression between Olink and MS for ALS vs HC, HC as the reference group
ALSvsHC_logFC_Olink_F <- extract_logFC(Olink_merged, Olink_ProF_matchClinic, "AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL")
ALSvsHC_logFC_MS_BEADdel_F <- extract_logFC(MS_BEADdel_merged, MS_BEADdel_ProF_matchClinic, "AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL")
ALSvsHC_logFC_MS_NONdel_F <- extract_logFC(MS_NONdel_merged, MS_NONdel_ProF_matchClinic, "AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL")

ALSvsHC_SigPro_List <- unique(c(sig_BEADdel[["ALSvsHC"]], sig_NONdel[["ALSvsHC"]], sig_Olink[["ALSvsHC"]]))

ALSvsHC_SigPro_Exist_Platform <- matrix(NA, nrow = length(ALSvsHC_SigPro_List), ncol=9)
rownames(ALSvsHC_SigPro_Exist_Platform) <- ALSvsHC_SigPro_List
colnames(ALSvsHC_SigPro_Exist_Platform) <- c("Olink", "MS_BEADdel", "MS_NONdel", "padj_Olink", "logFC_Olink", "padj_MS_BEADdel", "logFC_MS_BEADdel", "padj_MS_NONdel", "logFC_MS_NONdel")

for(ProS in ALSvsHC_SigPro_List){ 
  ALSvsHC_SigPro_Exist_Platform[ProS, "Olink"] <- ifelse(ProS %in% colnames(Olink_ProF), "YES", "NO") 
  ALSvsHC_SigPro_Exist_Platform[ProS, "MS_BEADdel"] <- ifelse(ProS %in% colnames(MS_BEADdel_ProF), "YES", "NO") 
  ALSvsHC_SigPro_Exist_Platform[ProS, "MS_NONdel"] <- ifelse(ProS %in% colnames(MS_NONdel_ProF), "YES", "NO")
}

for(ProS in ALSvsHC_SigPro_List){
  ALSvsHC_SigPro_Exist_Platform[ProS, "logFC_Olink"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "Olink"] == "YES", signif(ALSvsHC_logFC_Olink_F[ProS],4), NA)
  ALSvsHC_SigPro_Exist_Platform[ProS, "padj_Olink"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "Olink"] == "YES", signif(result_tbl_adj_Olink[ProS,"ALSvsHC"],4), NA)
  
  ALSvsHC_SigPro_Exist_Platform[ProS, "logFC_MS_BEADdel"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "MS_BEADdel"] == "YES", signif(ALSvsHC_logFC_MS_BEADdel_F[ProS],4), NA)
  ALSvsHC_SigPro_Exist_Platform[ProS, "padj_MS_BEADdel"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "MS_BEADdel"] == "YES", signif(result_tbl_adj_BEADdel[ProS,"ALSvsHC"],4), NA)
  
  ALSvsHC_SigPro_Exist_Platform[ProS, "logFC_MS_NONdel"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "MS_NONdel"] == "YES", signif(ALSvsHC_logFC_MS_NONdel_F[ProS],4), NA)
  ALSvsHC_SigPro_Exist_Platform[ProS, "padj_MS_NONdel"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "MS_NONdel"] == "YES", signif(result_tbl_adj_NONdel[ProS,"ALSvsHC"],4), NA)
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
  
  # BEADdel vs NONdel
  if (has_BEADdel && has_NONdel) {
    print(
      plot_scatter(
        MS_BEADdel_ProF, MS_NONdel_ProF,
        protein = ProS,
        label1 = "MS_BEADdel",
        label2 = "MS_NONdel"
      )
    )
  }
  
}

dev.off()

#----------------------
#----------------------
# Special for DCvsHC
#----------------------
#----------------------
# Make scatter plot of protein expression between Olink and MS for ALS vs HC, HC as the reference group
DCvsHC_logFC_Olink_F <- extract_logFC(Olink_merged, Olink_ProF_matchClinic, "DISEASE CONTROL", "HEALTHY CONTROL")
DCvsHC_logFC_MS_BEADdel_F <- extract_logFC(MS_BEADdel_merged, MS_BEADdel_ProF_matchClinic, "DISEASE CONTROL", "HEALTHY CONTROL")
DCvsHC_logFC_MS_NONdel_F <- extract_logFC(MS_NONdel_merged, MS_NONdel_ProF_matchClinic, "DISEASE CONTROL", "HEALTHY CONTROL")

DCvsHC_SigPro_List <- unique(c(sig_BEADdel[["DCvsHC"]], sig_NONdel[["DCvsHC"]], sig_Olink[["DCvsHC"]]))

###########################
###########################
# SANITY CHECK FOR ALSvsHC
###########################
###########################
#-------
# Limma 
#-------

library(limma)
library(edgeR)   

Olink_merged_test <- Olink_merged %>% filter(!(is.na(ALSvsHC) | is.na(AGE_AT_SAMPLING) | is.na(GENDER)))
Olink_ProF_test <- Olink_ProF[rownames(Olink_merged_test), ]

expr <- as.matrix(t(Olink_ProF_test))

# 1. Set reference level
Olink_merged_test$ALSvsHC <- relevel(Olink_merged_test$ALSvsHC, ref = "HEALTHY CONTROL")

# 2. Build design matrix
design <- model.matrix(~ ALSvsHC + AGE_AT_SAMPLING + GENDER, data = Olink_merged_test)
colnames(design) <- make.names(colnames(design))
colnames(design)

fit <- lmFit(expr, design)

# 3. Define contrast
contrast.matrix <- makeContrasts(
  ALS_vs_HC = ALSvsHCAMYOTROPHIC.LATERAL.SCLEROSIS,   # the coefficient comparing ALS to baseline (HC)
  levels = design
)

# 4. Fit contrasts and eBayes
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

colnames(fit2$coefficients)

results <- topTable(fit2, coef = "ALS_vs_HC", number = Inf)
head(results)

DEGs <- subset(results, adj.P.Val < 0.05)

ALSvsHC_limma_DE <- rownames(DEGs)

ALSvsHC_Olink_DE <- sig_Olink$ALSvsHC

library(VennDiagram)

venn.plot <- draw.pairwise.venn(
  area1 = length(limma_DE),
  area2 = length(Olink_DE),
  cross.area = length(common_DE),
  category = c("limma DE", "Olink DE"),
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  lty = 1,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20)
)

grid.newpage()
grid.draw(venn.plot)

Olink_merged_test <- Olink_merged %>% filter(!(is.na(DCvsHC) | is.na(AGE_AT_SAMPLING) | is.na(GENDER)))
Olink_ProF_test <- Olink_ProF[rownames(Olink_merged_test), ]

expr <- as.matrix(t(Olink_ProF_test))

# 1. Set reference level
Olink_merged_test$DCvsHC <- relevel(Olink_merged_test$DCvsHC, ref = "HEALTHY CONTROL")

# 2. Build design matrix
design <- model.matrix(~ DCvsHC + AGE_AT_SAMPLING + GENDER, data = Olink_merged_test)
colnames(design) <- make.names(colnames(design))
colnames(design)

fit <- lmFit(expr, design)

# 3. Define contrast
contrast.matrix <- makeContrasts(
  DC_vs_HC = DCvsHCDISEASE.CONTROL,   # the coefficient comparing ALS to baseline (HC)
  levels = design
)

# 4. Fit contrasts and eBayes
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

colnames(fit2$coefficients)

results <- topTable(fit2, coef = "DC_vs_HC", number = Inf)
head(results)

DEGs <- subset(results, adj.P.Val < 0.05)

limma_DE <- rownames(DEGs)

# ---------------------
# Logistic Regression
# ---------------------
### Merge the clinical meta and proteomic data
commonID1 <- intersect(rownames(MS_BEADdel_ProF), rownames(Clinic_MetaF_All))
MS_BEADdel_merged_All <- cbind(
  MS_BEADdel_ProF[commonID1, , drop=FALSE],
  Clinic_MetaF_All[commonID1, , drop=FALSE]
)
MS_BEADdel_ProF_matchClinic_All <- MS_BEADdel_ProF[commonID1, , drop = FALSE]

commonID2 <- intersect(rownames(MS_NONdel_ProF), rownames(Clinic_MetaF_All))
MS_NONdel_merged_All <- cbind(
  MS_NONdel_ProF[commonID2, , drop = FALSE],
  Clinic_MetaF_All[commonID2, , drop = FALSE]
)
MS_NONdel_ProF_matchClinic_All <- MS_NONdel_ProF[commonID2, , drop = FALSE]

commonID3 <- intersect(rownames(Olink_ProF), rownames(Clinic_MetaF_All))
Olink_merged_All <- cbind(
  Olink_ProF[commonID3, , drop = FALSE],
  Clinic_MetaF_All[commonID3, , drop = FALSE]
)
Olink_ProF_matchClinic_All <- Olink_ProF[commonID3, ,drop = FALSE]

### baseline sample only
ALSvsHC_Olink_result_tbl_adj_v1_GML <- get_DE_Pvalue_Table_GLM(Olink_ProF_matchClinic, Olink_merged, varList2, outPath, "Olink", "logistic", "ALSvsHC", FALSE)
ALSvsHC_Olink_logisticReg_v1_DE <- rownames(ALSvsHC_Olink_result_tbl_adj_v1_GML)[which(ALSvsHC_Olink_result_tbl_adj_v1_GML < 0.05)]

ALSvsHC_MS_BEADdel_result_tbl_adj_v1_GML <- get_DE_Pvalue_Table_GLM(MS_BEADdel_ProF_matchClinic, MS_BEADdel_merged, varList2, outPath, "MS_BEADdel", "logistic", "ALSvsHC", FALSE)
ALSvsHC_MS_BEADdel_logisticReg_v1_DE <- rownames(ALSvsHC_MS_BEADdel_result_tbl_adj_v1_GML)[which(ALSvsHC_MS_BEADdel_result_tbl_adj_v1_GML < 0.05)]

ALSvsHC_MS_NONdel_result_tbl_adj_v1_GML <- get_DE_Pvalue_Table_GLM(MS_NONdel_ProF_matchClinic, MS_NONdel_merged, varList2, outPath, "MS_NONdel", "logistic", "ALSvsHC", FALSE)
ALSvsHC_MS_NONdel_logisticReg_v1_DE <- rownames(ALSvsHC_MS_NONdel_result_tbl_adj_v1_GML)[which(ALSvsHC_MS_NONdel_result_tbl_adj_v1_GML < 0.05)]

ALSvsHC_lgReg_SigPro_List <- unique(c(ALSvsHC_Olink_logisticReg_v1_DE, ALSvsHC_MS_BEADdel_logisticReg_v1_DE, ALSvsHC_MS_NONdel_logisticReg_v1_DE))

ALSvsHC_lgReg_SigPro_Exist_Platform <- matrix(NA, nrow = length(ALSvsHC_lgReg_SigPro_List), ncol=9)
rownames(ALSvsHC_lgReg_SigPro_Exist_Platform) <- ALSvsHC_lgReg_SigPro_List
colnames(ALSvsHC_lgReg_SigPro_Exist_Platform) <- c("Olink", "MS_BEADdel", "MS_NONdel", "padj_Olink", "logFC_Olink", "padj_MS_BEADdel", "logFC_MS_BEADdel", "padj_MS_NONdel", "logFC_MS_NONdel")

for(ProS in ALSvsHC_lgReg_SigPro_List){ 
  ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "Olink"] <- ifelse(ProS %in% colnames(Olink_ProF), "YES", "NO") 
  ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "MS_BEADdel"] <- ifelse(ProS %in% colnames(MS_BEADdel_ProF), "YES", "NO") 
  ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "MS_NONdel"] <- ifelse(ProS %in% colnames(MS_NONdel_ProF), "YES", "NO")
}

for(ProS in ALSvsHC_lgReg_SigPro_List){
  ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "logFC_Olink"] <- ifelse(ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "Olink"] == "YES", signif(ALSvsHC_logFC_Olink_F[ProS],4), NA)
  ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "padj_Olink"] <- ifelse(ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "Olink"] == "YES", signif(ALSvsHC_Olink_result_tbl_adj_v1_GML[ProS,1],4), NA)
  
  ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "logFC_MS_BEADdel"] <- ifelse(ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "MS_BEADdel"] == "YES", signif(ALSvsHC_logFC_MS_BEADdel_F[ProS],4), NA)
  ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "padj_MS_BEADdel"] <- ifelse(ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "MS_BEADdel"] == "YES", signif(ALSvsHC_MS_BEADdel_result_tbl_adj_v1_GML[ProS,1],4), NA)
  
  ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "logFC_MS_NONdel"] <- ifelse(ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "MS_NONdel"] == "YES", signif(ALSvsHC_logFC_MS_NONdel_F[ProS],4), NA)
  ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "padj_MS_NONdel"] <- ifelse(ALSvsHC_lgReg_SigPro_Exist_Platform[ProS, "MS_NONdel"] == "YES", signif(ALSvsHC_MS_NONdel_result_tbl_adj_v1_GML[ProS,1],4), NA)
}

write.csv(ALSvsHC_lgReg_SigPro_Exist_Platform, file = paste0(outPath, "ALSvsHC_lgReg_SigPro_Exist_Platform.csv"))

pdf(paste0(outPath, "ALSvsHC_lgReg_SigPro_In_Platform.pdf"))

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

for (ProS in ALSvsHC_lgReg_SigPro_List) {
  
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
  
  # BEADdel vs NONdel
  if (has_BEADdel && has_NONdel) {
    print(
      plot_scatter(
        MS_BEADdel_ProF, MS_NONdel_ProF,
        protein = ProS,
        label1 = "MS_BEADdel",
        label2 = "MS_NONdel"
      )
    )
  }
  
}

dev.off()

#######################################################################################################################################################################################
DCvsHC_Olink_result_tbl_adj_v1_GML <- get_DE_Pvalue_Table_GLM(Olink_ProF_matchClinic, Olink_merged, varList2, outPath, "Olink", "logistic", "DCvsHC", FALSE)
DCvsHC_Olink_logisticReg_v1_DE <- rownames(DCvsHC_Olink_result_tbl_adj_v1_GML)[which(DCvsHC_Olink_result_tbl_adj_v1_GML<0.05)]

DCvsHC_MS_BEADdel_result_tbl_adj_v1_GML <- get_DE_Pvalue_Table_GLM(MS_BEADdel_ProF_matchClinic, MS_BEADdel_merged, varList2, outPath, "MS_BEADdel", "logistic", "DCvsHC", FALSE)
DCvsHC_MS_BEADdel_logisticReg_v1_DE <- rownames(DCvsHC_MS_BEADdel_result_tbl_adj_v1_GML)[which(DCvsHC_MS_BEADdel_result_tbl_adj_v1_GML < 0.05)]

DCvsHC_MS_NONdel_result_tbl_adj_v1_GML <- get_DE_Pvalue_Table_GLM(MS_NONdel_ProF_matchClinic, MS_NONdel_merged, varList2, outPath, "MS_NONdel", "logistic", "DCvsHC", FALSE)
DCvsHC_MS_NONdel_logisticReg_v1_DE <- rownames(DCvsHC_MS_NONdel_result_tbl_adj_v1_GML)[which(DCvsHC_MS_NONdel_result_tbl_adj_v1_GML < 0.05)]

DCvsHC_SigPro_List <- unique(c(DCvsHC_Olink_logisticReg_v1_DE, DCvsHC_MS_BEADdel_logisticReg_v1_DE, DCvsHC_MS_NONdel_logisticReg_v1_DE))

DCvsHC_logFC_Olink_F <- extract_logFC(Olink_merged, Olink_ProF_matchClinic, "DISEASE CONTROL", "HEALTHY CONTROL")
DCvsHC_logFC_MS_BEADdel_F <- extract_logFC(MS_BEADdel_merged, MS_BEADdel_ProF_matchClinic, "DISEASE CONTROL", "HEALTHY CONTROL")
DCvsHC_logFC_MS_NONdel_F <- extract_logFC(MS_NONdel_merged, MS_NONdel_ProF_matchClinic, "DISEASE CONTROL", "HEALTHY CONTROL")

DCvsHC_SigPro_Exist_Platform <- matrix(NA, nrow = length(DCvsHC_SigPro_List), ncol=9)
rownames(DCvsHC_SigPro_Exist_Platform) <- DCvsHC_SigPro_List
colnames(DCvsHC_SigPro_Exist_Platform) <- c("Olink", "MS_BEADdel", "MS_NONdel", "padj_Olink", "logFC_Olink", "padj_MS_BEADdel", "logFC_MS_BEADdel", "padj_MS_NONdel", "logFC_MS_NONdel")

for(ProS in DCvsHC_SigPro_List){ 
  DCvsHC_SigPro_Exist_Platform[ProS, "Olink"] <- ifelse(ProS %in% colnames(Olink_ProF), "YES", "NO") 
  DCvsHC_SigPro_Exist_Platform[ProS, "MS_BEADdel"] <- ifelse(ProS %in% colnames(MS_BEADdel_ProF), "YES", "NO") 
  DCvsHC_SigPro_Exist_Platform[ProS, "MS_NONdel"] <- ifelse(ProS %in% colnames(MS_NONdel_ProF), "YES", "NO")
}

for(ProS in DCvsHC_SigPro_List){
  DCvsHC_SigPro_Exist_Platform[ProS, "logFC_Olink"] <- ifelse(DCvsHC_SigPro_Exist_Platform[ProS, "Olink"] == "YES", signif(DCvsHC_logFC_Olink_F[ProS],4), NA)
  DCvsHC_SigPro_Exist_Platform[ProS, "padj_Olink"] <- ifelse(DCvsHC_SigPro_Exist_Platform[ProS, "Olink"] == "YES", signif(DCvsHC_Olink_result_tbl_adj_v1_GML[ProS,1],4), NA)
  
  DCvsHC_SigPro_Exist_Platform[ProS, "logFC_MS_BEADdel"] <- ifelse(DCvsHC_SigPro_Exist_Platform[ProS, "MS_BEADdel"] == "YES", signif(DCvsHC_logFC_MS_BEADdel_F[ProS],4), NA)
  DCvsHC_SigPro_Exist_Platform[ProS, "padj_MS_BEADdel"] <- ifelse(DCvsHC_SigPro_Exist_Platform[ProS, "MS_BEADdel"] == "YES", signif(DCvsHC_MS_BEADdel_result_tbl_adj_v1_GML[ProS,1],4), NA)
  
  DCvsHC_SigPro_Exist_Platform[ProS, "logFC_MS_NONdel"] <- ifelse(DCvsHC_SigPro_Exist_Platform[ProS, "MS_NONdel"] == "YES", signif(DCvsHC_logFC_MS_NONdel_F[ProS],4), NA)
  DCvsHC_SigPro_Exist_Platform[ProS, "padj_MS_NONdel"] <- ifelse(DCvsHC_SigPro_Exist_Platform[ProS, "MS_NONdel"] == "YES", signif(DCvsHC_MS_NONdel_result_tbl_adj_v1_GML[ProS,1],4), NA)
}

write.csv(DCvsHC_SigPro_Exist_Platform, file = paste0(outPath, "DCvsHC_SigPro_Exist_Platform.csv"))

pdf(paste0(outPath, "DCvsHC_SigPro_In_Platform.pdf"))

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

for (ProS in DCvsHC_SigPro_List) {
  
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
  
  # BEADdel vs NONdel
  if (has_BEADdel && has_NONdel) {
    print(
      plot_scatter(
        MS_BEADdel_ProF, MS_NONdel_ProF,
        protein = ProS,
        label1 = "MS_BEADdel",
        label2 = "MS_NONdel"
      )
    )
  }
  
}

dev.off()
# ------------------------------------------
# Different Regression Model Comparison
# ------------------------------------------
# Linear-regression DE list
### 1. Linear-regression DE list
ALSvsHC_SigPro_List <- unique(c(
  sig_BEADdel[["ALSvsHC"]],
  sig_NONdel[["ALSvsHC"]],
  sig_Olink[["ALSvsHC"]]
))

### 2. Logistic-regression DE list
ALSvsHC_lgReg_SigPro_List <- unique(c(
  ALSvsHC_Olink_logisticReg_v1_DE, 
  ALSvsHC_MS_BEADdel_logisticReg_v1_DE, 
  ALSvsHC_MS_NONdel_logisticReg_v1_DE
))

### 3. Platform-specific DE sets
# Linear
lin_olink  <- sig_Olink[["ALSvsHC"]]
lin_bead   <- sig_BEADdel[["ALSvsHC"]]
lin_nondel <- sig_NONdel[["ALSvsHC"]]

# Logistic
log_olink  <- ALSvsHC_Olink_logisticReg_v1_DE
log_bead   <- ALSvsHC_MS_BEADdel_logisticReg_v1_DE
log_nondel <- ALSvsHC_MS_NONdel_logisticReg_v1_DE

### 4. Load library
library(VennDiagram)

### 5. Function to draw and save platform-wise Venn diagrams
plot_lin_log_venn <- function(lin_set, log_set, title, filename){
  venn.plot <- venn.diagram(
    x = list(
      Linear   = lin_set,
      Logistic = log_set
    ),
    filename = NULL,
    fill = c("#1f78b4", "#33a02c"),
    alpha = 0.5,
    cex = 1.0,
    cat.cex = 1.0,
    cat.pos = 0,
    main = title,
    main.cex = 1.2
  )
  
  png(filename, width = 1800, height = 1800, res = 300)
  grid::grid.draw(venn.plot)
  dev.off()
}

### 6. Generate Venn diagrams for each platform

# Olink
plot_lin_log_venn(
  lin_set = lin_olink,
  log_set = log_olink,
  title = "Olink: Linear vs Logistic DE Proteins",
  filename = paste0(outPath,"venn_olink_lin_vs_log.png")
)

# MS_BEADdel
plot_lin_log_venn(
  lin_set = lin_bead,
  log_set = log_bead,
  title = "MS_BEADdel: Linear vs Logistic DE Proteins",
  filename = paste0(outPath,"venn_bead_lin_vs_log.png")
)

# MS_NONdel
plot_lin_log_venn(
  lin_set = lin_nondel,
  log_set = log_nondel,
  title = "MS_NONdel: Linear vs Logistic DE Proteins",
  filename = paste0(outPath,"venn_nondel_lin_vs_log.png")
)


# -------------------
# Enrichment Testing
# -------------------
# fgseaPath <- fgsea::fgsea(pathways = all_gene_sets, stats=ranks, scoreType = whichType, eps = 0.0) ### we will not limit the size (minSize and maxSize) here, just output everything, but would select the proper size for further plotting. 

