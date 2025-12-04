### This script generates mapping files for proteins and samples among MS and Olink data sets.
library(tidyverse)
library(magrittr)
library(dplyr)
library("xlsx")

inputPath <- "/home/rstudio/YD/MS/input/"

# Now I check the protein names among the three datasets
DEL_ProMap <- read_excel(paste0(inputPath, "ProteinMapping/BEADdepletions_ProteinMapping.xlsx"))
NON_ProMap <- read_excel(paste0(inputPath, "ProteinMapping/NONDEPLETED_ProteinMapping.xlsx"))
combinedP <- read.csv("/home/rstudio/YD/input/P230192_CSP1-5_CINO_and_C2I2N2O2_EXTENDED_NPX_2024-01-18.csv", sep = ";")
combinedP2 <- combinedP %>% filter(QC_Warning=="PASS" & Assay_Warning=="PASS" & Sample_Type=="SAMPLE" & (!grepl("control", Assay))) 
Olink_ProMap <- combinedP2[,c("UniProt", "Assay")] %>%  distinct() %>% mutate(ModifiedAssay = gsub("_", ".", Assay))
colnames(Olink_ProMap)[which(colnames(Olink_ProMap) == "Assay")] <- "OriginalAssay"

# write.table(Olink_ProMap,file = paste0(inputPath, "ProteinMapping/Olink_ProteinMapping.csv"))
write.xlsx(Olink_ProMap, row.names = FALSE, file = paste0(inputPath,"/ProteinMapping/Olink_ProteinMapping.xlsx"),
           sheetName="ProteinMapping", append=TRUE)

### clinical variables
ALS_ClinicVar1_F <- read.csv(paste0(inputPath, "20240827_minus_asymptomatic_at_risk_pseudo_anonymised_v2.csv"), sep=",")
ALS_ClinicVar2_F <- read.csv(paste0(inputPath, "20241014_asymptomatic_at_risk_binned_blinded.csv"), sep=",")
ALS_ClinicVarF <- rbind(ALS_ClinicVar1_F[intersect(colnames(ALS_ClinicVar1_F), colnames(ALS_ClinicVar2_F))], ALS_ClinicVar2_F[intersect(colnames(ALS_ClinicVar1_F), colnames(ALS_ClinicVar2_F))])
Olink_Frame <- read.csv("/home/rstudio/YD/input/imcm_team_als_csf_olink_cleaned_2024-08-29.csv", row.names = 1)
Olink_Frame %<>% select(-one_of(c("PlateID", "median_npx")))
DEL_Sample_Raw <-  read.table("/home/rstudio/YD/MS/intermediate/BEADdel/SAMPLE/ProDat_RmContam.csv")
NON_Sample_Raw <- read.table("/home/rstudio/YD/MS/intermediate/NONdel/SAMPLE/ProDat_RmContam.csv")

### Extract the raw sample IDs
Olink_Samp <- rownames(Olink_Frame)
DEL_Samp <- rownames(DEL_Sample_Raw)
NON_Samp <- rownames(NON_Sample_Raw)
ClinicVar <- ALS_ClinicVarF$CSF_OLINK_MANIFEST
ClinicVar <- ClinicVar[which(!ClinicVar=="")]

### Keep all the Olink samples in the final frame
keyID1 <- sub("C9_", "C9C_", Olink_Samp)
keyID2 <- sub("BM1_", "BX1_", keyID1)

Olink <- Olink_Samp
MS_NonDepleted <- unlist(sapply(keyID2, function(x){
  anyMap <- which(NON_Samp==x)
  if(length(anyMap) == 0){"NA"
  }else{NON_Samp[anyMap]}}))
MS_Depleted <- unlist(sapply(keyID2, function(x){
  anyMap <- which(DEL_Samp==x)
  if(length(anyMap) == 0){"NA"
  }else{DEL_Samp[anyMap]}}))
Clinical_Meta <- unlist(sapply(keyID1, function(x){
  if(x %in% c("BM2_0108_1", "BM2_0010_1", "BM2_0161_1")){
    x <- sub("M","X",x)}
  anyMap <- which(ClinicVar==x)
  if(length(anyMap) == 0){"NA"
  }else{ClinicVar[anyMap]}}))

Recommend_Consistent_Name <- keyID2
SampleMappingF1 <- cbind(Recommend_Consistent_Name, Olink, MS_NonDepleted, MS_Depleted, Clinical_Meta) %>% as.data.frame()

### Now add in MS samples which are not in Olink
MS_Sample_All <- unique(c(NON_Samp, DEL_Samp))
keyID3 <- MS_Sample_All[which(!(MS_Sample_All %in% SampleMappingF1$Recommend_Consistent_Name))]
Olink2 <- unlist(sapply(keyID3, function(x){
  anyMap <- which(Olink_Samp==x)
  if(length(anyMap) == 0){"NA"
  }else{Olink_Samp[anyMap]}}))
MS_NonDepleted2 <- unlist(sapply(keyID3, function(x){
  anyMap <- which(NON_Samp==x)
  if(length(anyMap) == 0){"NA"
  }else{NON_Samp[anyMap]}}))
MS_Depleted2 <- unlist(sapply(keyID3, function(x){
  anyMap <- which(DEL_Samp==x)
  if(length(anyMap) == 0){"NA"
  }else{DEL_Samp[anyMap]}}))
Clinical_Meta2 <- unlist(sapply(keyID3, function(x){
  if(x %in% c("BM2_0108_1", "BM2_0010_1", "BM2_0161_1")){
    x <- sub("M","X",x)}
  anyMap <- which(ClinicVar==x)
  if(length(anyMap) == 0){"NA"
  }else{ClinicVar[anyMap]}}))

Recommend_Consistent_Name2 <- keyID3
SampleMappingF2 <- cbind(Recommend_Consistent_Name2, Olink2, MS_NonDepleted2, MS_Depleted2, Clinical_Meta2) %>% as.data.frame()
colnames(SampleMappingF2) <- colnames(SampleMappingF1)

SampleMappingF <- rbind(SampleMappingF1,SampleMappingF2)

### Generate the sample mapping table
write.xlsx(SampleMappingF, row.names = FALSE, file = paste0(inputPath,"SampleMapping.xlsx"),
           sheetName="SampleMapping", append=TRUE)

# pay attention, we did not include 122 clinical samples which are not in Olink NOR MS in our mapping file.
# ClinicVar[which(!(ClinicVar %in% SampleMappingF1$Recommend_Consistent_Name))]

# additional exploration about the diff items in clinical meta
ALS_ID_F <-  ALS_ClinicVarF %>% filter(CSF_OLINK_MANIFEST != "" & SERUM_OLINK_MANIFEST != "" & CSF_OLINK_MANIFEST != SERUM_OLINK_MANIFEST) %>% select(CSF_OLINK_MANIFEST, SERUM_OLINK_MANIFEST)
whichDiff <- setdiff(ALS_ID_F$CSF_OLINK_MANIFEST, ALS_ID_F$SERUM_OLINK_MANIFEST)
DiffID <- which(ALS_ID_F$CSF_OLINK_MANIFEST != ALS_ID_F$SERUM_OLINK_MANIFEST)
DiffF <- ALS_ID_F[DiffID,]

Olink_csf_Frame <- read.csv("/home/rstudio/YD/input/imcm_team_als_csf_olink_cleaned_2024-08-29.csv", row.names = 1)
Olink_csf_Frame %<>% select(-one_of(c("PlateID", "median_npx")))
Olink_serum_Frame <- read.csv("/home/rstudio/YD/input/imcm_team_als_serum_olink_cleaned_2024-09-03.csv", row.names = 1)
Olink_serum_Frame %<>% select(-one_of(c("PlateID", "median_npx")))

which(rownames(Olink_csf_Frame) %in% DiffF$SERUM_OLINK_MANIFEST & rownames(Olink_csf_Frame) %in% DiffF$CSF_OLINK_MANIFEST)
which(rownames(Olink_serum_Frame) %in% DiffF$CSF_OLINK_MANIFEST & rownames(Olink_serum_Frame) %in% DiffF$CSF_OLINK_MANIFEST)
DiffF[which(DiffF$SERUM_OLINK_MANIFEST %in% rownames(Olink_csf_Frame) & DiffF$CSF_OLINK_MANIFEST %in% rownames(Olink_csf_Frame)),]
DiffF[which(DiffF$CSF_OLINK_MANIFEST %in% rownames(Olink_serum_Frame) & DiffF$SERUM_OLINK_MANIFEST %in% rownames(Olink_serum_Frame)),]
