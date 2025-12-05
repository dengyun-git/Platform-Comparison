# This script performs differential annalysis using all samples not only baseline samples, where including sampleID as random effect
library(tidyverse)
library(magrittr)
library(readxl)
library(fgsea)
library(data.table)
library(MASS)

inputPath <- "/home/rstudio/workspace/Data Collection/"
intermediatePath <- "/home/rstudio/workspace/Data Collection/intermediate/"
outPath <- "/home/rstudio/workspace/Data Collection/output_AllSamp/"

source("/home/rstudio/repos/DE_functions_call.R")

### Read in Data
MS_BEADdel_ProF <- read.csv("/home/rstudio/workspace/Data Collection/QCed/MS/BEADdel/ProDat_QCed.csv", row.names = 1)
DEL_ProMap <- read_excel("/home/rstudio/workspace/Data Collection/QCed/MS/BEADdel/BEADdepletions_ProteinMapping.xlsx") %>% as.data.frame()
colnames(MS_BEADdel_ProF) <- sapply(colnames(MS_BEADdel_ProF), function(x){DEL_ProMap[which(DEL_ProMap$Protein.Group==x),"Genes"]}) ### Using the Entrez Gene Naming scheme

MS_NONdel_ProF <- read.csv("/home/rstudio/workspace/Data Collection/QCed/MS/NONdel/ProDat_QCed.csv", row.names = 1)
NON_ProMap <- read_excel("/home/rstudio/workspace/Data Collection/QCed/MS/NONdel/NONDEPLETED_ProteinMapping.xlsx") %>% as.data.frame()
colnames(MS_NONdel_ProF) <- sapply(colnames(MS_NONdel_ProF), function(x){NON_ProMap[which(NON_ProMap$Protein.Group==x),"Genes"]})

Olink_ProF <- read.csv("/home/rstudio/workspace/Data Collection/QCed/Olink/CSF/CSF_Customised_ProDataClean.csv", row.names = 1)

Clinic_MetaF <- read.csv("/home/rstudio/workspace/Data Collection/DA76/da76-wb-fleeting-cabbage-4050/Data/20251025_als_clinical_with_ambrosia_data_pseudonymized_FINAL.csv") %>% filter(CSF_OLINK_MANIFEST != "")
rownames(Clinic_MetaF) <- Clinic_MetaF$CSF_OLINK_MANIFEST

# Only the base line samples will be considered in Clinic_MetaF
Clinic_MetaF %<>% 
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

### Merge the clinical meta and proteomic data
commonID1 <- intersect(rownames(MS_BEADdel_ProF), rownames(Clinic_MetaF))
MS_BEADdel_merged <- cbind(MS_BEADdel_ProF[commonID1,], Clinic_MetaF[commonID1,]) 

commonID2 <- intersect(rownames(MS_NONdel_ProF), rownames(Clinic_MetaF))
MS_NONdel_merged <- cbind(MS_NONdel_ProF[commonID2,], Clinic_MetaF[commonID2,])

commonID3 <- intersect(rownames(Olink_ProF), rownames(Clinic_MetaF))
Olink_merged <- cbind(Olink_ProF[commonID3,], Clinic_MetaF[commonID3,])


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

result_tbl_adj_BEADdel <- get_DE_Pvalue_Table_AllSamp(MS_BEADdel_ProF, MS_BEADdel_merged, varList2, "/home/rstudio/workspace/Data Collection/output/", "BEADdel", varList1)
result_tbl_adj_NONdel <- get_DE_Pvalue_Table_AllSamp(MS_NONdel_ProF, MS_NONdel_merged, varList2, "/home/rstudio/workspace/Data Collection/output/", "NONdel", varList1)
result_tbl_adj_Olink <- get_DE_Pvalue_Table_AllSamp(Olink_ProF, Olink_merged, varList2, "/home/rstudio/workspace/Data Collection/output/", "Olink", varList1)

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

