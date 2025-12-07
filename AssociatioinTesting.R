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
library(patchwork)
library(ggrepel)
library(plotly)


inputPath <- "/home/rstudio/workspace/Data Collection/"
intermediatePath <- "/home/rstudio/workspace/Data Collection/intermediate_AssociationTesting/"
outPath <- "/home/rstudio/workspace/Data Collection/output_AssociationTesting/"

source("/home/rstudio/repos/DE_functions_call.R")

##################
################## 
# Load Data
################## 
################## 
BEADdel_ProF <- read.csv(paste0(inputPath, "oxf-qc2-als-mass-spectometry-wb-genial-hazelnut-7025/Data/output/BEADdel/SAMPLE/QCed_Frame.csv"), row.names = 1)[,-1]
DEL_ProMap <- read_excel(paste0(inputPath, "oxf-qc2-als-mass-spectometry-wb-genial-hazelnut-7025/Data/output/BEADdel/SAMPLE/BEADdepletions_ProteinMapping.xlsx"))
colnames(BEADdel_ProF) <- sapply(colnames(BEADdel_ProF), function(x){DEL_ProMap[which(DEL_ProMap$Protein.Group==x),"Genes"]}) ### Using the Entrez Gene Naming scheme

NONdel_ProF <- read.csv(paste0(inputPath, "oxf-qc2-als-mass-spectometry-wb-genial-hazelnut-7025/Data/output/NONdel/SAMPLE/QCed_Frame.csv"), row.names = 1)[,-1]
NON_ProMap <- read_excel(paste0(inputPath, "oxf-qc2-als-mass-spectometry-wb-genial-hazelnut-7025/Data/output/NONdel/SAMPLE/NONDEPLETED_ProteinMapping.xlsx"))
colnames(NONdel_ProF) <- sapply(colnames(NONdel_ProF), function(x){NON_ProMap[which(NON_ProMap$Protein.Group==x),"Genes"]})

Olink_ProF <- read.csv(paste0(inputPath, "Data-CSF/oxf-da28-opdc-als-multiple-cohorts-olink-csf/CSF_Customised_ProDataClean.csv"), sep= ",", row.names = 1)

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
    ALSFRSR_FINE_MOTOR = boxcoxTranform(ALSFRSR_FINE_MOTOR, 12),
    ALSFRSR_GROSS_MOTOR = boxcoxTranform(ALSFRSR_GROSS_MOTOR, 12),
    
    # SubjectID for random effect regression
    SubjectID_Random = CSF_OLINK_MANIFEST,
    
    # Factor conversions
    GROUP = factor(GROUP),
    ALSvsHC = relevel(factor(ALSvsHC), ref = "HEALTHY CONTROL"),
    ALSvsDC = relevel(factor(ALSvsDC), ref = "DISEASE CONTROL"),
    DCvsHC  = relevel(factor(DCvsHC),  ref = "HEALTHY CONTROL"),
    GENDER  = relevel(factor(GENDER),  ref = "FEMALE"),
    WEAKNESS_SITE = factor(WEAKNESS_SITE),
    PATHOGENIC_VARIANT = factor(PATHOGENIC_VARIANT),
    COHORT  = factor(COHORT),
    BULBARvSpinal = relevel(factor(BULBARvSpinal), ref = "BULBAR"),
    C9vsNonC9 = relevel(factor(C9vsNonC9), ref = "C9ORF72"),
    SubjectID_Random = factor(SubjectID_Random)
  )

### Normality transformation
Clinic_MetaF_All$ALSFRSR_RATE <- ifelse(Clinic_MetaF_All$ALSFRSR_RATE > 0,
                                        log(Clinic_MetaF_All$ALSFRSR_RATE),
                                        NA)

Clinic_MetaF_All$ECAS_SCORE <- boxcoxTranform(Clinic_MetaF_All$ECAS_SCORE, 134)
Clinic_MetaF_All$ALSFRSR_FINE_MOTOR <- boxcoxTranform(Clinic_MetaF_All$ALSFRSR_FINE_MOTOR, 12)
Clinic_MetaF_All$ALSFRSR_GROSS_MOTOR <- boxcoxTranform(Clinic_MetaF_All$ALSFRSR_GROSS_MOTOR, 12)

Clinic_MetaF <- Clinic_MetaF_All %>% filter(LONGITUDINAL_ENCOUNTER_NUMBER == 1) # Clinic_MetaF only contain the baseline samples

### Merge the clinical meta and proteomic data
commonID1 <- intersect(rownames(BEADdel_ProF), rownames(Clinic_MetaF))
BEADdel_merged <- cbind(
  BEADdel_ProF[commonID1, , drop=FALSE],
  Clinic_MetaF[commonID1, , drop=FALSE]
)
BEADdel_ProF_matchClinic <- BEADdel_ProF[commonID1, , drop = FALSE]

commonID2 <- intersect(rownames(NONdel_ProF), rownames(Clinic_MetaF))
NONdel_merged <- cbind(
  NONdel_ProF[commonID2, , drop = FALSE],
  Clinic_MetaF[commonID2, , drop = FALSE]
)
NONdel_ProF_matchClinic <- NONdel_ProF[commonID2, , drop = FALSE]

commonID3 <- intersect(rownames(Olink_ProF), rownames(Clinic_MetaF))
Olink_merged <- cbind(
  Olink_ProF[commonID3, , drop = FALSE],
  Clinic_MetaF[commonID3, , drop = FALSE]
)
Olink_ProF_matchClinic <- Olink_ProF[commonID3, ,drop = FALSE]

save(Clinic_MetaF_All, Clinic_MetaF, 
     Olink_merged, Olink_ProF, Olink_ProF_matchClinic, 
     BEADdel_merged, BEADdel_ProF, BEADdel_ProF_matchClinic,
     NONdel_merged, NONdel_ProF, NONdel_ProF_matchClinic, file = paste0(inputPath, "output_report/Load_Data_For_Report.rdata"))

##########################
##########################
# ClinicaL Variable View
##########################
##########################
# Define variable lists
cat_vars <- c(
  "GROUP", "ALSvsHC", "ALSvsDC", "DCvsHC", "GENDER", 
  "WEAKNESS_SITE", "PATHOGENIC_VARIANT", "COHORT",
  "BULBARvSpinal", "C9vsNonC9"
)

num_vars <- c(
  "AGE_AT_SAMPLING",
  "AGE_AT_SYMPTOM_ONSET",
  "BODY_MASS_INDEX",
  "ALSFRSR_FINE_MOTOR",
  "ALSFRSR_GROSS_MOTOR",
  "ALSFRSR_RATE",
  "ECAS_SCORE"
)

varList <- c(cat_vars, num_vars)

All_ClinicVar_View <- ViewClinincalVar(cat_vars, num_vars, Clinic_MetaF)

Olink_ClinicVar_View <- ViewClinincalVar(cat_vars, num_vars, Olink_merged)

BEADdel_ClinicVar_View <- ViewClinincalVar(cat_vars, num_vars, BEADdel_merged)

NONdel_ClinicVar_View <- ViewClinincalVar(cat_vars, num_vars, NONdel_merged)

save(All_ClinicVar_View, Olink_ClinicVar_View, BEADdel_ClinicVar_View, NONdel_ClinicVar_View, file = paste0(intermediatePath, "ClinicVar_View.Rdata"))
########################
########################
# Association Testing
######################## 
########################
# ---------------------
# Logistic Regression
# ---------------------
logi_vars <- c("ALSvsHC", "ALSvsDC", "DCvsHC", "BULBARvSpinal", "C9vsNonC9")
line_var <- c("ALSFRSR_FINE_MOTOR", "ALSFRSR_GROSS_MOTOR", "ALSFRSR_RATE", "ECAS_SCORE")
varList_Conf <- c("AGE_AT_SYMPTOM_ONSET", "AGE_AT_SAMPLING","GENDER", "COHORT")
varAll_List <- c(logi_vars, line_var)

result_tbl_adj_BEADdel_logi <- get_DE_Pvalue_Table_GLM(ProExpF = BEADdel_ProF_matchClinic, 
                                                       Merged = BEADdel_merged, 
                                                       covariateList = varList_Conf, 
                                                       outPath = intermediatePath, 
                                                       whichPlatform = "BEADdel", 
                                                       regMethod = "logistic",
                                                       varList1 = logi_vars,
                                                       use_random_YN = FALSE)

result_tbl_adj_NONdel_logi <- get_DE_Pvalue_Table_GLM(ProExpF = NONdel_ProF_matchClinic, 
                                                      Merged = NONdel_merged, 
                                                      covariateList = varList_Conf, 
                                                      outPath = intermediatePath, 
                                                      whichPlatform = "NONdel", 
                                                      regMethod = "logistic",
                                                      varList1 = logi_vars,
                                                      use_random_YN = FALSE) 

result_tbl_adj_Olink_logi <- get_DE_Pvalue_Table_GLM(ProExpF = Olink_ProF_matchClinic, 
                                                     Merged = Olink_merged, 
                                                     covariateList = varList_Conf, 
                                                     outPath = intermediatePath, 
                                                     whichPlatform = "Olink", 
                                                     regMethod = "logistic",
                                                     varList1 = logi_vars,
                                                     use_random_YN = FALSE)


# ---------------------
# Linear Regression
# ---------------------
result_tbl_adj_BEADdel_lin <- get_DE_Pvalue_Table_GLM(ProExpF = BEADdel_ProF_matchClinic, 
                                                      Merged = BEADdel_merged, 
                                                      covariateList = varList_Conf, 
                                                      outPath = intermediatePath, 
                                                      whichPlatform = "BEADdel", 
                                                      regMethod = "linear",
                                                      varList1 = line_var,
                                                      use_random_YN = FALSE)

result_tbl_adj_NONdel_lin <- get_DE_Pvalue_Table_GLM(ProExpF = NONdel_ProF_matchClinic, 
                                                     Merged = NONdel_merged, 
                                                     covariateList = varList_Conf, 
                                                     outPath = intermediatePath, 
                                                     whichPlatform = "NONdel", 
                                                     regMethod = "linear",
                                                     varList1 = line_var,
                                                     use_random_YN = FALSE) 

result_tbl_adj_Olink_lin <- get_DE_Pvalue_Table_GLM(ProExpF = Olink_ProF_matchClinic, 
                                                    Merged = Olink_merged, 
                                                    covariateList = varList_Conf, 
                                                    outPath = intermediatePath, 
                                                    whichPlatform = "Olink", 
                                                    regMethod = "linear",
                                                    varList1 = line_var,
                                                    use_random_YN = FALSE)

# -----------------------------------------------
# Combine All the Clinical Meta Testing Results
# -----------------------------------------------
result_tbl_adj_BEADdel <- cbind(result_tbl_adj_BEADdel_logi, result_tbl_adj_BEADdel_lin[rownames(result_tbl_adj_BEADdel_logi), ])
result_tbl_adj_NONdel <- cbind(result_tbl_adj_NONdel_logi, result_tbl_adj_NONdel_lin[rownames(result_tbl_adj_NONdel_logi), ])
result_tbl_adj_Olink <- cbind(result_tbl_adj_Olink_logi, result_tbl_adj_Olink_lin[rownames(result_tbl_adj_Olink_logi), ])

# Extract the common significantly expressed proteins
# Significance threshold
alpha <- 0.05

sig_BEADdel <- sig_NONdel <- sig_Olink <- list()
names(sig_BEADdel) <- names(sig_NONdel) <- names(sig_Olink)
for(var in varAll_List){
  sig_BEADdel[[var]] <- rownames(result_tbl_adj_BEADdel)[order(result_tbl_adj_BEADdel[, var])][which(result_tbl_adj_BEADdel[order(result_tbl_adj_BEADdel[, var]), var] < 0.05)]
  
  sig_NONdel[[var]] <- rownames(result_tbl_adj_NONdel)[order(result_tbl_adj_NONdel[, var])][which(result_tbl_adj_NONdel[order(result_tbl_adj_NONdel[, var]), var] < 0.05)]
  
  sig_Olink[[var]] <- rownames(result_tbl_adj_Olink)[order(result_tbl_adj_Olink[, var])][which(result_tbl_adj_Olink[order(result_tbl_adj_Olink[, var]), var] < 0.05)]
}

### Save all the significant protein names per platform
save(sig_BEADdel, sig_NONdel, sig_Olink, file=paste0(outPath, "sig.rdat"))

# Construct a table, rows as clinical variables, columns as three platforms, entries are significant proteins
SigProSummary <- data.frame(
  ClinicalVar = varAll_List,
  BEADdel = sapply(varAll_List, function(v) paste(sig_BEADdel[[v]], collapse = "/")),
  NONdel  = sapply(varAll_List, function(v) paste(sig_NONdel [[v]], collapse = "/")),
  Olink   = sapply(varAll_List, function(v) paste(sig_Olink  [[v]], collapse = "/")),
  stringsAsFactors = FALSE
)

write.csv(SigProSummary, file = paste0(outPath, "SigProSummary.csv"), row.names = FALSE)

# compute overlaps per var and build a summary table
summary_list <- lapply(varAll_List, function(var) {
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
overlap_list <- lapply(varAll_List, compute_overlap)
names(overlap_list) <- varAll_List

# Write summary
write.csv(summary_df, file = paste0(outPath, "overlap_summary_counts.csv"), 
          row.names = FALSE)

# Also save element names for full inspection
save(overlap_list, file = paste0(outPath, "overlap_elements.rdat"))

summary_df

# -----------------------------------------------
# Combine and Write Logistic and Linear Results
# -----------------------------------------------
# Load effect size files 
Olink_eff_logi  <- read_df(paste0(intermediatePath, "Olink_logistic_RandomEffectFALSE_result_tbl_EffectSize.csv"))
BEADdel_eff_logi <- read_df(paste0(intermediatePath, "BEADdel_logistic_RandomEffectFALSE_result_tbl_EffectSize.csv"))
NONdel_eff_logi <- read_df(paste0(intermediatePath,  "NONdel_logistic_RandomEffectFALSE_result_tbl_EffectSize.csv"))
Olink_eff_lin  <- read_df(paste0(intermediatePath, "Olink_linear_RandomEffectFALSE_result_tbl_EffectSize.csv"))
BEADdel_eff_lin <- read_df(paste0(intermediatePath, "BEADdel_linear_RandomEffectFALSE_result_tbl_EffectSize.csv"))
NONdel_eff_lin <- read_df(paste0(intermediatePath,  "NONdel_linear_RandomEffectFALSE_result_tbl_EffectSize.csv"))

# Load p value files 
Olink_p_logi  <- read_df(paste0(intermediatePath, "Olink_logistic_RandomEffectFALSE_result_tbl_p.csv"))
BEADdel_p_logi <- read_df(paste0(intermediatePath, "BEADdel_logistic_RandomEffectFALSE_result_tbl_p.csv"))
NONdel_p_logi <- read_df(paste0(intermediatePath, "NONdel_logistic_RandomEffectFALSE_result_tbl_p.csv"))
Olink_p_lin  <- read_df(paste0(intermediatePath, "Olink_linear_RandomEffectFALSE_result_tbl_p.csv"))
BEADdel_p_lin <- read_df(paste0(intermediatePath, "BEADdel_linear_RandomEffectFALSE_result_tbl_p.csv"))
NONdel_p_lin <- read_df(paste0(intermediatePath, "NONdel_linear_RandomEffectFALSE_result_tbl_p.csv"))

# Olink
common_proteins_olink <- intersect(rownames(Olink_eff_logi), rownames(Olink_eff_lin))
Olink_eff_logi <- Olink_eff_logi[common_proteins_olink, , drop = FALSE]
Olink_eff_lin  <- Olink_eff_lin[common_proteins_olink, , drop = FALSE]
Olink_p_logi   <- Olink_p_logi[common_proteins_olink, , drop = FALSE]
Olink_p_lin    <- Olink_p_lin[common_proteins_olink, , drop = FALSE]

# BEADdel
common_proteins_beaddel <- intersect(rownames(BEADdel_eff_logi), rownames(BEADdel_eff_lin))
BEADdel_eff_logi <- BEADdel_eff_logi[common_proteins_beaddel, , drop = FALSE]
BEADdel_eff_lin  <- BEADdel_eff_lin[common_proteins_beaddel, , drop = FALSE]
BEADdel_p_logi   <- BEADdel_p_logi[common_proteins_beaddel, , drop = FALSE]
BEADdel_p_lin    <- BEADdel_p_lin[common_proteins_beaddel, , drop = FALSE]

# NONdel
common_proteins_nondel <- intersect(rownames(NONdel_eff_logi), rownames(NONdel_eff_lin))
NONdel_eff_logi <- NONdel_eff_logi[common_proteins_nondel, , drop = FALSE]
NONdel_eff_lin  <- NONdel_eff_lin[common_proteins_nondel, , drop = FALSE]
NONdel_p_logi   <- NONdel_p_logi[common_proteins_nondel, , drop = FALSE]
NONdel_p_lin    <- NONdel_p_lin[common_proteins_nondel, , drop = FALSE]

# Effect sizes
Olink_eff_combined   <- cbind(Olink_eff_logi, Olink_eff_lin)
BEADdel_eff_combined <- cbind(BEADdel_eff_logi, BEADdel_eff_lin)
NONdel_eff_combined  <- cbind(NONdel_eff_logi, NONdel_eff_lin)

# P-values
Olink_p_combined   <- cbind(Olink_p_logi, Olink_p_lin)
BEADdel_p_combined <- cbind(BEADdel_p_logi, BEADdel_p_lin)
NONdel_p_combined  <- cbind(NONdel_p_logi, NONdel_p_lin)

write.csv(Olink_eff_combined,  file = paste0(outPath, "Olink_EffectSize.csv"))
write.csv(BEADdel_eff_combined, file = paste0(outPath, "BEADdel_EffectSize.csv"))
write.csv(NONdel_eff_combined,  file = paste0(outPath, "NONdel_EffectSize.csv"))

write.csv(Olink_p_combined,  file = paste0(outPath, "Olink_pvalue.csv"))
write.csv(BEADdel_p_combined, file = paste0(outPath, "BEADdel_pvalue.csv"))
write.csv(NONdel_p_combined,  file = paste0(outPath, "NONdel_pvalue.csv"))

write.csv(result_tbl_adj_Olink, file = paste0(outPath, "Olink_adjp.csv"))
write.csv(result_tbl_adj_BEADdel, file = paste0(outPath, "BEADdel_adjp.csv"))
write.csv(result_tbl_adj_NONdel, file = paste0(outPath, "NONdel_adjp.csv"))

########################
########################
# Vocano Plot
######################## 
########################
platformList <- c("Olink", "BEADdel", "NONdel")

Olink_eff <- read.csv(paste0(outPath, "Olink_EffectSize.csv"), row.names = 1)
BEADdel_eff <- read.csv(paste0(outPath, "BEADdel_EffectSize.csv"), row.names = 1)
NONdel_eff <- read.csv(paste0(outPath, "NONdel_EffectSize.csv"), row.names = 1)

effect_list <- list(
  Olink = Olink_eff,
  BEADdel = BEADdel_eff,
  NONdel = NONdel_eff
)

Olink_adjp <- read_df(paste0(outPath, "Olink_adjp.csv"))
BEADdel_adjp <- read_df(paste0(outPath, "BEADdel_adjp.csv"))
NONdel_adjp <- read_df(paste0(outPath, "NONdel_adjp.csv"))

adjp_list <- list(
  Olink = Olink_adjp,
  BEADdel = BEADdel_adjp,
  NONdel = NONdel_adjp
)

volcano_panel_list <- vector("list", length = length(varAll_List))
names(volcano_panel_list) <- varAll_List

 # Loop over traits to create and store volcano plots panels
for (trait in varAll_List) {
  volcano_panel_list[[trait]] <- make_volcano_panel(trait, effect_list, adjp_list)
  
  out_file <- paste0(outPath, trait, "_volcano_panel.pdf")
  
  ggsave(
    filename = out_file,
    plot = volcano_panel_list[[trait]],
    device = "pdf",
    width = 14,   # landscape width
    height = 6,   # landscape height
    units = "in"
  )
}

# Save the list of volcano panels
save(volcano_panel_list, file = paste0(intermediatePath, "volcano_panel_list.rdata"))

# Interactive volcano plot
plotly_volcano_list <- list()

# Loop through each trait and platform
for (trait in varAll_List) {
  plotly_volcano_list[[trait]] <- list()
  
  for (plat in platformList) {
    df_eff <- effect_list[[plat]]
    df_p   <- adjp_list[[plat]]
    
    plotly_volcano_list[[trait]][[plat]] <- make_volcano_plotly(df_eff, df_p, trait, plat)
  }
}

save(plotly_volcano_list, file = paste0(intermediatePath, "plotly_volcano_list.rdata"))

########################
########################
# Enrichment Testing
######################## 
########################
### Ranking matrix calculation
Ranks_Olink_AllVar <- Calculate_Ranking_Metric(paste0(outPath, "Olink_pvalue.csv"), 
                                               paste0(outPath, "Olink_EffectSize.csv"))

Ranks_BEADdel_AllVar <- Calculate_Ranking_Metric(paste0(outPath, "BEADdel_pvalue.csv"), 
                                                  paste0(outPath, "BEADdel_EffectSize.csv"))

Ranks_NONdel_AllVar <- Calculate_Ranking_Metric(paste0(outPath, "NONdel_pvalue.csv"), 
                                                  paste0(outPath, "NONdel_EffectSize.csv"))

### Prepare genesets used for our enrichment test
GOBPGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c5.go.bp.v2023.2.Hs.symbols.gmt"))
GOMFGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c5.go.mf.v2023.2.Hs.symbols.gmt"))
GOCCGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c5.go.cc.v2023.2.Hs.symbols.gmt"))
REACTOMEGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.reactome.v2023.2.Hs.symbols.gmt"))
KEGGGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt"))
BiocartaGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.biocarta.v2023.2.Hs.symbols.gmt"))
WikiPathwaysGeneSets <- fgsea::gmtPathways(paste0(inputPath,"OnlineResource/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt"))

pathBaseList <- c("GOBP","GOMF","GOCC","REACTOME","KEGG","Biocarta","WikiPathways") 

GeneSetList <- sapply(paste0(pathBaseList,"GeneSets"),function(x){get(x)})

SigPath_List <- lapply(varAll_List, function(v) {
  
  # For each variable, create platform layer
  platformContainer <- lapply(platformList, function(p) {
    
    # For each platform, create gene set layer
    geneSetContainer <- setNames(
      vector("list", length(pathBaseList)),
      pathBaseList
    )
    
    return(geneSetContainer)
  })
  
  # Name platforms inside this variable
  platformContainer <- setNames(platformContainer, platformList)
  
  return(platformContainer)
})

# Name the top layer with variable names
SigPath_List <- setNames(SigPath_List, varAll_List)

for (ClinicVar in varAll_List) {
  
  for (platform in platformList) {
    
    # Extract numeric ranking vector
    ranks_M <- get(paste0("Ranks_", platform, "_AllVar"))
    ranks <- ranks_M[,ClinicVar]
    names(ranks) <- rownames(ranks_M)
    ranks <- ranks[which(!is.na(ranks))]
    
    # Determine scoreType
    if (sign(min(ranks)) * sign(max(ranks)) == -1) {
      whichType <- "std"
    } else if (min(ranks) > 0) {
      whichType <- "pos"
    } else {
      whichType <- "neg"
    }
    
    for (GScont in seq_along(GeneSetList)) {
      
      GeneSetsX <- GeneSetList[[GScont]]
      
      fgseaPath <- fgsea::fgsea(
        pathways = GeneSetsX,
        stats = ranks,
        scoreType = whichType,
        eps = 0.0,
        minSize = 10,
        maxSize = 1000
      )
      
      # write result table
      fwrite(
        fgseaPath,
        file = paste0(
          intermediatePath,
          pathBaseList[GScont], "_", ClinicVar, "_", platform, ".txt"
        ),
        sep = "\t",
        sep2 = c("", " ", "")
      )
      
      sig_fgseaPath <- fgseaPath %>%
        filter(padj < 0.05) %>%
        arrange(padj)
      
      # save into SigPath_List (IMPORTANT FIX: moved inside loop!)
      SigPath_List[[ClinicVar]][[platform]][[pathBaseList[GScont]]] <- sig_fgseaPath
    }
  }
}

### Bubble plot to show the significantly enriched pathways
Enrich_Rdat_files <- outer(
  pathBaseList,
  varAll_List,
  Vectorize(function(p, v) paste0(p, "_", v, ".rdata"))
)

Enrich_Rdat_files <- as.vector(Enrich_Rdat_files)

pdf(paste0(outPath, "AssociationTesting_PathwayEnrichment.pdf"))

fileCt <- 1
for (pathBase in pathBaseList){
  for(ClinicVar in varAll_List){
    patternHere <- paste0(pathBase, "_", ClinicVar)
    makeBubble(intermediatePath, patternHere, 5, 5, Enrich_Rdat_files[fileCt])
    fileCt <- fileCt + 1
  }
}

dev.off()
