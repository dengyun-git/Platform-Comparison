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
intermediatePath <- "/home/rstudio/workspace/Data Collection/intermediate_DifferentialExpression/"
outPath <- "/home/rstudio/workspace/Data Collection/output_DifferentialExpression/"

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

##############################
##############################
# DE per Trait per Platform 
##############################
##############################
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

result_tbl_adj_BEADdel <- get_DE_Pvalue_Table(BEADdel_ProF_matchClinic, BEADdel_merged, varList2, outPath, "BEADdel", varList1)
result_tbl_adj_NONdel <- get_DE_Pvalue_Table(NONdel_ProF_matchClinic, NONdel_merged, varList2, outPath, "NONdel", varList1)
result_tbl_adj_Olink <- get_DE_Pvalue_Table(Olink_ProF_matchClinic, Olink_merged, varList2, outPath, "Olink", varList1)

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
ALSvsHC_logFC_BEADdel_F <- extract_logFC(BEADdel_merged, BEADdel_ProF_matchClinic, "AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL")
ALSvsHC_logFC_NONdel_F <- extract_logFC(NONdel_merged, NONdel_ProF_matchClinic, "AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL")

ALSvsHC_SigPro_List <- unique(c(sig_BEADdel[["ALSvsHC"]], sig_NONdel[["ALSvsHC"]], sig_Olink[["ALSvsHC"]]))

ALSvsHC_SigPro_Exist_Platform <- matrix(NA, nrow = length(ALSvsHC_SigPro_List), ncol=9)
rownames(ALSvsHC_SigPro_Exist_Platform) <- ALSvsHC_SigPro_List
colnames(ALSvsHC_SigPro_Exist_Platform) <- c("Olink", "BEADdel", "NONdel", "padj_Olink", "logFC_Olink", "padj_BEADdel", "logFC_BEADdel", "padj_NONdel", "logFC_NONdel")

for(ProS in ALSvsHC_SigPro_List){ 
  ALSvsHC_SigPro_Exist_Platform[ProS, "Olink"] <- ifelse(ProS %in% colnames(Olink_ProF), "YES", "NO") 
  ALSvsHC_SigPro_Exist_Platform[ProS, "BEADdel"] <- ifelse(ProS %in% colnames(BEADdel_ProF), "YES", "NO") 
  ALSvsHC_SigPro_Exist_Platform[ProS, "NONdel"] <- ifelse(ProS %in% colnames(NONdel_ProF), "YES", "NO")
}

for(ProS in ALSvsHC_SigPro_List){
  ALSvsHC_SigPro_Exist_Platform[ProS, "logFC_Olink"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "Olink"] == "YES", signif(ALSvsHC_logFC_Olink_F[ProS],4), NA)
  ALSvsHC_SigPro_Exist_Platform[ProS, "padj_Olink"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "Olink"] == "YES", signif(result_tbl_adj_Olink[ProS,"ALSvsHC"],4), NA)
  
  ALSvsHC_SigPro_Exist_Platform[ProS, "logFC_BEADdel"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "BEADdel"] == "YES", signif(ALSvsHC_logFC_BEADdel_F[ProS],4), NA)
  ALSvsHC_SigPro_Exist_Platform[ProS, "padj_BEADdel"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "BEADdel"] == "YES", signif(result_tbl_adj_BEADdel[ProS,"ALSvsHC"],4), NA)
  
  ALSvsHC_SigPro_Exist_Platform[ProS, "logFC_NONdel"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "NONdel"] == "YES", signif(ALSvsHC_logFC_NONdel_F[ProS],4), NA)
  ALSvsHC_SigPro_Exist_Platform[ProS, "padj_NONdel"] <- ifelse(ALSvsHC_SigPro_Exist_Platform[ProS, "NONdel"] == "YES", signif(result_tbl_adj_NONdel[ProS,"ALSvsHC"],4), NA)
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
  has_BEADdel <- ProS %in% colnames(BEADdel_ProF)
  has_NONdel  <- ProS %in% colnames(NONdel_ProF)
  
  # Olink vs BEADdel
  if (has_Olink && has_BEADdel) {
    print(
      plot_scatter(
        Olink_ProF, BEADdel_ProF,
        protein = ProS,
        label1 = "Olink",
        label2 = "BEADdel"
      )
    )
  }
  
  # Olink vs NONdel
  if (has_Olink && has_NONdel) {
    print(
      plot_scatter(
        Olink_ProF, NONdel_ProF,
        protein = ProS,
        label1 = "Olink",
        label2 = "NONdel"
      )
    )
  }
  
  # BEADdel vs NONdel
  if (has_BEADdel && has_NONdel) {
    print(
      plot_scatter(
        BEADdel_ProF, NONdel_ProF,
        protein = ProS,
        label1 = "BEADdel",
        label2 = "NONdel"
      )
    )
  }
  
}

dev.off()

platformList <- c("Olink", "BEADdel", "NONdel")

########################
########################
# Vocano Plot
######################## 
########################
Olink_eff <- read_df(paste0(outPath, "Olink_result_tbl_EffectSize.csv"))
BEADdel_eff <- read_df(paste0(outPath, "BEADdel_result_tbl_EffectSize.csv"))
NONdel_eff <- read_df(paste0(outPath, "NONdel_result_tbl_EffectSize.csv"))

effect_list <- list(
  Olink = Olink_eff,
  BEADdel = BEADdel_eff,
  NONdel = NONdel_eff
)

Olink_adjp <- read_df(paste0(outPath, "Olink_result_tbl_adjp.csv"))
BEADdel_adjp <- read_df(paste0(outPath, "BEADdel_result_tbl_adjp.csv"))
NONdel_adjp <- read_df(paste0(outPath, "NONdel_result_tbl_adjp.csv"))

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
### traits for testing
varAll_List <- c("ALSvsHC", "ALSvsDC", "DCvsHC", "BULBARvSpinal", "C9vsNonC9", "BODY_MASS_INDEX", "ALSFRSR_FINE_MOTOR", "ALSFRSR_GROSS_MOTOR", "ALSFRSR_RATE", "ECAS_SCORE")

### Ranking matrix calculation
Ranks_Olink_AllVar <- Calculate_Ranking_Metric(paste0(outPath, "Olink_result_tbl_p.csv"), 
                                               paste0(outPath, "Olink_result_tbl_EffectSize.csv"))

Ranks_BEADdel_AllVar <- Calculate_Ranking_Metric(paste0(outPath, "BEADdel_result_tbl_p.csv"), 
                                                 paste0(outPath, "BEADdel_result_tbl_EffectSize.csv"))

Ranks_NONdel_AllVar <- Calculate_Ranking_Metric(paste0(outPath, "NONdel_result_tbl_p.csv"), 
                                                paste0(outPath, "NONdel_result_tbl_EffectSize.csv"))

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

pdf(paste0(outPath, "DifferentialExpression_PathwayEnrichment.pdf"))

fileCt <- 1
for (pathBase in pathBaseList){
  for(ClinicVar in varAll_List){
    patternHere <- paste0(pathBase, "_", ClinicVar)
    makeBubble(intermediatePath, patternHere, 5, 5, Enrich_Rdat_files[fileCt])
    fileCt <- fileCt + 1
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
DCvsHC_logFC_BEADdel_F <- extract_logFC(BEADdel_merged, BEADdel_ProF_matchClinic, "DISEASE CONTROL", "HEALTHY CONTROL")
DCvsHC_logFC_NONdel_F <- extract_logFC(NONdel_merged, NONdel_ProF_matchClinic, "DISEASE CONTROL", "HEALTHY CONTROL")

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

ALSvsHC_common_DE <- intersect(ALSvsHC_limma_DE, ALSvsHC_Olink_DE)

library(VennDiagram)

venn.plot <- draw.pairwise.venn(
  area1 = length(ALSvsHC_limma_DE),
  area2 = length(ALSvsHC_Olink_DE),
  cross.area = length(ALSvsHC_common_DE),
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

DCvsHC_limma_DE <- rownames(DEGs)

