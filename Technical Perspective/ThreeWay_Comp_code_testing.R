# This script is my testing for the code Comp_MS_Olink_Interactive.Rmd, Comp_Dataset_Call.R

# self written function interactiveRHO
testPro <- "APOB"
cor(Olink_Frame[intersect(rownames(Olink_Frame),rownames(NON_Sample_Keep)),testPro], NON_Sample_Keep[intersect(rownames(Olink_Frame),rownames(NON_Sample_Keep)),testPro])
cor(Olink_Frame[intersect(rownames(Olink_Frame),rownames(NON_Sample_Raw)),testPro], NON_Sample_Raw[intersect(rownames(Olink_Frame),rownames(NON_Sample_Raw)),testPro], use = "complete.obs")
cor(DEL_Sample_Keep[intersect(rownames(DEL_Sample_Keep),rownames(NON_Sample_Keep)),testPro], NON_Sample_Keep[intersect(rownames(DEL_Sample_Keep),rownames(NON_Sample_Keep)),testPro])
cor(DEL_Sample_Raw[intersect(rownames(DEL_Sample_Raw),rownames(NON_Sample_Raw)),testPro], NON_Sample_Raw[intersect(rownames(DEL_Sample_Raw),rownames(NON_Sample_Raw)),testPro], use = "complete.obs")

inputPath <- "/home/rstudio/YD/MS/input/"
# Now I check the protein names among the three datasets
DEL_ProMap <- read_excel(paste0(inputPath, "ProteinMapping/BEADdepletions_ProteinMapping.xlsx"))
NON_ProMap <- read_excel(paste0(inputPath, "ProteinMapping/NONDEPLETED_ProteinMapping.xlsx"))
combinedP <- read.csv("/home/rstudio/YD/input/P230192_CSP1-5_CINO_and_C2I2N2O2_EXTENDED_NPX_2024-01-18.csv", sep = ";")
combinedP2 <- combinedP %>% filter(QC_Warning=="PASS" & Assay_Warning=="PASS" & Sample_Type=="SAMPLE" & (!grepl("control", Assay))) 
Olink_ProMap <- combinedP2[,c("UniProt", "Assay")] %>%  distinct() %>% mutate(Assay1 = gsub("-", ".", Assay)) %>% mutate(Assay2 = gsub("_", ".", Assay))

# write.table(Olink_ProMap,file = paste0(inputPath, "ProteinMapping/Olink_ProteinMapping.csv"))

### "_" & "-", and previously MS ";", carefully map again
Olink_Frame <- read.csv("/home/rstudio/YD/input/imcm_team_als_csf_olink_cleaned_2024-08-29.csv", row.names = 1)
Olink_Frame %<>% select(-one_of(c("PlateID", "median_npx")))
Olink_ProName_Entrez <- colnames(Olink_Frame)
Olink_ProName_UniProt <- sapply(Olink_ProName_Entrez, function(x) {Olink_ProMap[which(Olink_ProMap$Assay1==x)[1], "UniProt"]}) ### pay attention to the "-" "_" in the names
Olink_ProName_UniProt <- Olink_ProName_UniProt[which(!(is.na(Olink_ProName_UniProt)))]

DEL_Sample <-  read.table("/home/rstudio/YD/MS/output/BEADdel/SAMPLE/QCed_Frame.csv")[,-1]
DEL_ProName_UniPro <- colnames(DEL_Sample)
DEL_ProName_Entrez <- sapply(DEL_ProName_UniPro, function(x){unique(DEL_ProMap[which(DEL_ProMap$Protein.Group==x),"Protein.Names"])[1]})
DEL_ProName_Card <- sapply(DEL_ProName_UniPro, function(x){unique(DEL_ProMap[which(DEL_ProMap$Protein.Group==x),"Genes"])[1]})

NON_Sample <- read.table("/home/rstudio/YD/MS/output/NONdel/SAMPLE/QCed_Frame.csv")[,-1]
NON_ProName_UniPro <- colnames(NON_Sample)
NON_ProName_Entrez <- sapply(NON_ProName_UniPro, function(x){unique(NON_ProMap[which(NON_ProMap$Protein.Group==x),"Protein.Names"])[1]})
NON_ProName_Card <- sapply(NON_ProName_UniPro, function(x){unique(NON_ProMap[which(NON_ProMap$Protein.Group==x),"Genes"])[1]})

set_list <- list(
  "Olink" = Olink_ProName_UniProt,
  "MSbd" = DEL_ProName_UniPro,
  "MS" = NON_ProName_UniPro
)

elements <- unique(unlist(set_list))  # Extract all unique elements from sets
binary_matrix <- as.data.frame(
  sapply(set_list, function(set) elements %in% set)  # Check presence across sets
)
rownames(binary_matrix) <- elements  # Assign row names as elements
colnames(binary_matrix) <- names(set_list)  # Assign column names as set names

ComplexUpset::upset(
  binary_matrix,
  intersect = names(set_list),  # Use the names of the sets as intersection columns
  name = "Proteins in Sets"  # Label for the UpSet plot
) + ggtitle("Mapping by UniProt ID")


set_list <- list(
  "Olink" = unique(Olink_ProName_Entrez),
  "MSbd" = unique(c(DEL_ProName_Card)),
  "MS" = unique(c(NON_ProName_Card))
)

elements <- unique(unlist(set_list))  # Extract all unique elements from sets
binary_matrix <- as.data.frame(
  sapply(set_list, function(set) elements %in% set)  # Check presence across sets
)
rownames(binary_matrix) <- elements  # Assign row names as elements
colnames(binary_matrix) <- names(set_list)  # Assign column names as set names

ComplexUpset::upset(
  binary_matrix,
  intersect = names(set_list),  # Use the names of the sets as intersection columns
  name = "Proteins in Sets"  # Label for the UpSet plot
) + ggtitle("Mapping by Gene Cards")

