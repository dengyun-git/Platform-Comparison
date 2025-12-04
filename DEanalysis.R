library(tidyverse)
library(magrittr)
library(readxl)
library(fgsea)
library(data.table)
library(MASS)

inputPath <- "/home/rstudio/workspace/Data Collection/"
intermediatePath <- "/home/rstudio/workspace/Data Collection/intermediate/"
outPath <- "/home/rstudio/workspace/Data Collection/output/"

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

Clinic_MetaF %<>%
  mutate(
    # C9 vs Non-C9 classification
    C9vsNonC9 = ifelse(
      PATHOGENIC_VARIANT %in% c("KIF5A", "MAPT", "SBMA", "SOD1", "SPG7", "TBK1"),
      "NEGATIVE",
      PATHOGENIC_VARIANT
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

Olink_merged <- Olink_ProF %>%
  rownames_to_column("ID") %>%
  left_join(Clinic_MetaF %>% rownames_to_column("ID"), by = "ID") %>%
  column_to_rownames("ID") 

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
CorSampReport <- Clinic_MetaF[corSamp_BEADdel_NONdel,]

# Define variables
categorical_vars <- c("BULBARvSpinal", "C9vsNonC9", "GENDER", "COHORT")
numeric_vars     <- c("AGE_AT_SYMPTOM_ONSET", "BODY_MASS_INDEX", "ALSFRSR_RATE", "ECAS_SCORE")

all_vars <- c(categorical_vars, numeric_vars)

# Initialize result vector
CorSamp_vs_AllSamp <- setNames(rep(NA, length(all_vars)), all_vars)

# Compare categorical variables using Fisher's Exact Test
for(var in categorical_vars){
  testF <- fisher.test(table(CorSampReport[,var]), table(Clinic_MetaF[,var]))
  CorSamp_vs_AllSamp[var] <- testF$p.value
}

# Compare numeric variables using Wilcoxon rank-sum test (non-parametric)
for(var in numeric_vars){
  testKS <- ks.test(table(CorSampReport[,var]), table(Clinic_MetaF[,var]))
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
varList1 = c("BULBARvSpinal", "C9vsNonC9", "BODY_MASS_INDEX", "ALSFRSR_FINE_MOTOR", "ALSFRSR_GROSS_MOTOR", "ALSFRSR_RATE", "ECAS_SCORE")
varList2 = c("AGE_AT_SYMPTOM_ONSET","GENDER", "COHORT")
varList = c(varList1,varList2)

### check whether the categorical variables are set as factor
for(id in c(1,2,9,10)){
  print(class(Clinic_MetaF[,varList[id]]))
} 

result_tbl_adj_BEADdel <- get_DE_Pvalue_Table(MS_BEADdel_ProF, MS_BEADdel_merged, c("AGE_AT_SYMPTOM_ONSET", "GENDER"), "/home/rstudio/workspace/Data Collection/output/", "BEADdel")
result_tbl_adj_NONdel <- get_DE_Pvalue_Table(MS_NONdel_ProF, MS_NONdel_merged, c("AGE_AT_SYMPTOM_ONSET", "GENDER"), "/home/rstudio/workspace/Data Collection/output/", "NONdel")
result_tbl_adj_Olink <- get_DE_Pvalue_Table(Olink_ProF, Olink_merged, c("AGE_AT_SYMPTOM_ONSET", "GENDER"), "/home/rstudio/workspace/Data Collection/output/", "Olink")

get_DE_Pvalue_Table <- function(ProExpF, Merged, covariates, outPath, whichPlatform){
  # --- Initialize results container --
  result_tbl_adj <- result_tbl <- matrix(NA, nrow=ncol(ProExpF), ncol=length(varList1))
  rownames(result_tbl_adj) <- rownames(result_tbl) <- colnames(ProExpF)
  colnames(result_tbl_adj) <- colnames(result_tbl) <- varList1
  
  # --- Loop over main predictors ---
  for(mainVar in varList1){
    
    cat("Running regression for predictor:", mainVar, "\n")
    
    # remove NA
    keepID <- which(!(is.na(Merged[, mainVar]) | Merged[, mainVar]=="" | is.null(Merged[, mainVar])))
    Merged_Here <- Merged[keepID, ]
    
    # Loop over proteins
    for(i in colnames(ProExpF)){
      ProS <- Merged_Here[,i]
      
      # Build formula: y ~ mainVar + covariates
      formula_str1 <- paste("ProS ~", paste(c(mainVar, covariates), collapse = " + "))
      formula_str2 <- paste("ProS ~", paste(covariates, collapse = " + "))
      
      fit1 <- lm(as.formula(formula_str1), data = Merged_Here)
      fit2 <- lm(as.formula(formula_str2), data = Merged_Here)
      ANOVAobj <- anova(fit1,fit2)
      
      # Extract coefficient for main predictor
      result_tbl[i,mainVar] <- ANOVAobj$`Pr(>F)`[2]
    }
    result_tbl_adj[,mainVar] <- p.adjust(result_tbl[,mainVar], method = "BH")
  }
  
  # Adjust p-values for multiple testing (optional)
  write.csv(result_tbl_adj, paste0(outPath, whichPlatform, "_result_tbl.csv"))
  write.csv(result_tbl_adj, paste0(outPath, whichPlatform, "_result_tbl_adj.csv"))
  
  return(result_tbl_adj)
}

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

# BEADdel vs NONdel
overlap_BEADdel_NONdel <- intersect(sig_BEADdel, sig_NONdel)

# BEADdel vs Olink
overlap_BEADdel_Olink <- intersect(sig_BEADdel, sig_Olink)

# NONdel vs Olink
overlap_NONdel_Olink <- intersect(sig_NONdel, sig_Olink)

overlap_all3 <- Reduce(intersect, list(sig_BEADdel, sig_NONdel, sig_Olink))

data.frame(
  Comparison = c("BEADdel vs NONdel", "BEADdel vs Olink", "NONdel vs Olink", "All three"),
  nProteins = c(length(overlap_BEADdel_NONdel),
                length(overlap_BEADdel_Olink),
                length(overlap_NONdel_Olink),
                length(overlap_all3))
)

write.csv(data.frame)

# for (confounder in confList){
#   
#   ClinicFrameHere <- ClinicFrame[!is.na(ClinicFrame[,confounder]),]
#   
#   if(confounder %in% confList1){
#     f1 <- as.formula(EndoLabel~1)
#     f2 <- as.formula(paste("EndoLabel~",confounder))
#   }else{
#     f1 <- as.formula(EndoLabel~cohort_name)
#     f2 <- as.formula(paste("EndoLabel~cohort_name+",confounder))
#   }
#   
#   ### tackle with the different cases of endotype numbers
#   if(length(unique(ClinicFrameHere$EndoLabel))==2){
#     fit1 <- glm(f1,data=ClinicFrameHere,family = binomial(link = "logit"))
#     fit2 <- glm(f2,data=ClinicFrameHere,family = binomial(link = "logit"))
#     ANOVAobj <- anova(fit1,fit2)
#     pContainer[confounder,1] = pchisq(ANOVAobj$Deviance[nrow(ANOVAobj)], df=ANOVAobj$Df[nrow(ANOVAobj)],lower.tail=FALSE)
#   }else{
#     fit1 <- nnet::multinom(f1,data=ClinicFrameHere)
#     fit2 <- nnet::multinom(f2,data=ClinicFrameHere)
#     ANOVAobj <- anova(fit1,fit2)
#     pContainer[confounder,1]= ANOVAobj$`Pr(Chi)`[nrow(ANOVAobj)]
#   } 
#   
#   if(confounder == "sf_iknee_proc_batch"){
#     fit3 <- lmerTest::lmer(as.numeric(EndoLabel) ~ get(confounder) + (1|cohort_name), data=ClinicFrameHere)
#     ANOVAobj <- anova(fit3)
#     pContainer[confounder,1]=ANOVAobj$`Pr(>F)`
#   }
#   
#   
#   for(clinicInd in c(2,5,8,11,12,13,14,3,4,9,10)){
#     clinic = clinics[clinicInd] ### define which clinic to explore this loop
#     
#     keepId <- which(!(is.na(CombinedFrameThis[,clinic]) | CombinedFrameThis[,clinic]=="NA")) 
#     CombinedFrameHere <- CombinedFrameThis[keepId,] ### remove NA entries before regression
#     
#     ### categorical clinical features
#     if(clinicInd %in% c(2,5,8,11,12,13,14)){
#       CombinedFrameHere[,clinic] <- as.factor(CombinedFrameHere[,clinic])
#       
#       ### bar chart displaying clinical feature distribution within endotypes
#       barDat = CombinedFrameHere[,c("EndoLabel",clinic)]
#       compT <- table(barDat)/rowSums(table(barDat))
#       addCol <- matrix(unlist(sapply(1:nrow(barDat),function(x){
#         endo = as.character(barDat[x,"EndoLabel"])
#         clin = as.character(barDat[x,clinic])
#         compT[endo,clin]})),ncol=1)
#       colnames(addCol) = "Composition"
#       barDat <- cbind(barDat,addCol)
#       
#       if(length(unique(EndoLabel))==2){
#         p1 <- ggplot(data = barDat) + geom_bar(aes(x=as.factor(get(clinic)),y=Composition,fill=as.factor(EndoLabel)),stat="identity",position=position_dodge()) +
#           ggtitle(paste0(clinics[clinicInd]," comparision for endotypes")) + labs(x="",y="Composition",fill="Endotype") + scale_fill_manual(name = "Endotype", labels = HighLowIntensityLabel,values = HighLowIntensityColor) +
#           theme(plot.title=element_text(size = 15,face="bold",hjust=0.5),legend.position="bottom",
#                 axis.text.x = element_text(size=12,face="bold",angle = 45,hjust = 1),axis.title.x =element_text(size = 12,face="bold"),
#                 legend.title = element_text(size = 13,face="bold"),legend.text = element_text(size = 12,face="bold"),
#                 axis.text.y =element_text(size=13,face="bold"),axis.title.y =element_text(size=12,face="bold"))
#       }else{ ### more than 2 endotypes
#         p1 <- ggplot(data = barDat) + geom_bar(aes(x=as.factor(get(clinic)),y=Composition,fill=as.factor(EndoLabel)),stat="identity",position=position_dodge()) +
#           ggtitle(paste0(clinics[clinicInd]," comparision for endotypes")) + labs(x="",y="Composition",fill="Endotype") + scale_fill_discrete(name = "Endotype") +
#           theme(plot.title=element_text(size = 15,face="bold",hjust=0.5),legend.position="bottom",
#                 axis.text.x = element_text(size=12,face="bold",angle = 45,hjust = 1),axis.title.x =element_text(size = 12,face="bold"),
#                 legend.title = element_text(size = 13,face="bold"),legend.text = element_text(size = 12,face="bold"),
#                 axis.text.y =element_text(size=13,face="bold"),axis.title.y =element_text(size=12,face="bold"))
#       }     
#       print(p1)
#       
#       ### an example plot of KL grade
#       # p1 <- ggplot(data = barDat) + geom_bar(aes(x=as.factor(get(clinic)),y=Composition,fill=as.factor(EndoLabel)),stat="identity",position=position_dodge()) +
#       #   ggtitle("Ordinal KL Grade Comparision for Endotypes") + labs(x="Ordinal KL grade", y="Composition",fill="endotype") + 
#       #   scale_fill_discrete(name = "Endotype", labels = c("High Intensity Endotype", "Low Intensity Endotype")) +
#       #   theme(plot.title=element_text(size = 15,face="bold",hjust=0.5),legend.position="bottom",
#       #         axis.text.x = element_text(size=15,face="bold"),axis.title.x =element_text(size = 12,face="bold"),
#       #         legend.title = element_text(size = 13,face="bold"),legend.text = element_text(size = 12,face="bold"),
#       #         axis.text.y =element_text(size=13,face="bold"),axis.title.y =element_text(size=12,face="bold"))
#       # print(p1)
#       
#       ### PCA plot, overlay with clinical features
#       PCAoverlayDat <- cbind(barDat,PCs[keepId,])
#       p1.1 <- ggplot(data=PCAoverlayDat) + geom_point(aes(x=PC1,y=PC2,color=as.factor(get(clinic)),shape=as.factor(EndoLabel))) +
#         labs(shape="endotype",color=clinic) +
#         theme(axis.text.x=element_blank(),axis.text.y=element_blank()) + ggtitle(paste0(clinics[clinicInd]," overlay on PCA"))
#       print(p1.1)
#       
#       ### umap, overlap with clinical features
#       df <- data.frame(x = umap[keepId,1],y =umap[keepId,2],barDat)
#       colnames(df)[c(1,2)] = c("umap.D1","umap.D2")
#       p1.2 <- ggplot(df, aes(x=umap.D1, y=umap.D2, color=as.factor(get(clinic)),shape=as.factor(EndoLabel))) + geom_point() + ggtitle(paste0(clinics[clinicInd]," overlay on umap")) + labs(color=clinics[clinicInd],shape="endotype") +
#         labs(shape="endotype",color=clinic) +
#         theme(plot.title=element_text(size=14,hjust=0.5),legend.text=element_text(size = 10), axis.title=element_text(size = 12))
#       print(p1.2)
#       
#     }else{### continuous clinical features
#       ### violin plot displaying clinical feature distribution within endotypes
#       violinDat = CombinedFrameHere[,c("EndoLabel",clinic)]
#       nn=length(unique(CombinedFrameHere$EndoLabel))
#       
#       if(length(unique(EndoLabel))==2){
#         p2 <- ggplot(data=violinDat,aes(x=as.factor(EndoLabel),y=get(clinic),fill=as.factor(EndoLabel))) + geom_violin(width=1) + geom_boxplot(width=0.1, color="grey", alpha=0.5) +
#           labs(title=(paste0(clinic," distribution for endotypes")),x="Endotype",fill="Endotype") + ylab("distribution") + scale_fill_manual(name = "Endotype", labels = HighLowIntensityLabel,values = HighLowIntensityColor) +
#           theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.text.x = element_blank(),legend.position="bottom") +
#           geom_vline(xintercept = 0.5:(nn+0.5),linetype="dotted")
#       }else{
#         p2 <- ggplot(data=violinDat,aes(x=as.factor(EndoLabel),y=get(clinic),fill=as.factor(EndoLabel))) + geom_violin(width=1) + geom_boxplot(width=0.1, color="grey", alpha=0.5) +
#           labs(title=(paste0(clinic," distribution for endotypes")),x="Endotype",fill="Endotype") + ylab("distribution") + scale_fill_discrete(name = "Endotype") +
#           theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.text.x = element_blank(),legend.position="bottom") +
#           geom_vline(xintercept = 0.5:(nn+0.5),linetype="dotted") 
#       }
#       
#       ### an example plot of womac pain
#       # p2 <- ggplot(data=violinDat,aes(x=as.factor(EndoLabel),y=get(clinic),fill=as.factor(EndoLabel))) + geom_violin(width=1) + geom_boxplot(width=0.1, color="grey", alpha=0.5) +
#       #   labs(title=(paste0("WOMAC pain score"," distribution for endotypes")),x="endotype",fill="endotype") + ylab("distribution") +
#       #   scale_fill_discrete(name = "Endotype", labels = c("High Intensity Endotype", "Low Intensity Endotype")) +
#       #   theme(plot.title=element_text(size = 15,face="bold",hjust=0.5),legend.position="bottom",
#       #         panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.text.x = element_blank(),
#       #         axis.title.x =element_text(size = 12,face="bold"),
#       #         legend.title = element_text(size = 13,face="bold"),legend.text = element_text(size = 12,face="bold"),
#       #         axis.text.y =element_text(size=13,face="bold"),axis.title.y =element_text(size=12,face="bold")) +
#       #   geom_vline(xintercept = 0.5:(nn+0.5),linetype="dotted")
#       print(p2)
#       
#       PCAoverlayDat <- cbind(violinDat,PCs[keepId,])
#       p2.1 <- ggplot(data=PCAoverlayDat) + geom_point(aes(x=PC1,y=PC2,color=get(clinic),shape=as.factor(EndoLabel))) +
#         theme(axis.text.x=element_blank(),axis.text.y=element_blank()) + ggtitle(paste0(clinic," overlay on PCA")) +
#         labs(color=clinic,shape="endotype")
#       print(p2.1)
#       
#       ### umap, overlap with clinical features
#       df <- data.frame(x = umap[keepId,1],y = umap[keepId,2],violinDat)
#       p2.2 <- ggplot(df, aes(x=x, y=y, color=get(clinic),shape=as.factor(EndoLabel))) + geom_point() + ggtitle(paste0(clinic," overlay on umap")) +
#         theme(plot.title=element_text(size=14,hjust=0.5),legend.text=element_text(size = 10), axis.title=element_text(size = 12),axis.text.x=element_blank(),axis.text.y=element_blank()) +
#         labs(color=clinic,shape="endotype",x="Dimension1",y="Dimension2")
#       print(p2.2)
#     }
#     
#     f1 <- as.formula(EndoLabel~cohort_name)
#     f2 <- as.formula(paste("EndoLabel~cohort_name+",clinic))
#     # f1 <- as.formula(EndoLabel~cohort_name+totalProHere)
#     # f2 <- as.formula(paste("EndoLabel~cohort_name+totalProHere+",clinic))
#     
#     ### tackle with different number of endotypes
#     if(length(unique(CombinedFrameHere$EndoLabel))==2){
#       CombinedFrameHere$EndoLabel <- as.factor(sapply(CombinedFrameHere$EndoLabel,function(x){if(x!=1){0}else{x}}))
#       fit1 <- glm(f1,family = binomial(link = "logit"),data=CombinedFrameHere)
#       fit2 <- glm(f2,family = binomial(link = "logit"),data=CombinedFrameHere)
#       ANOVAobj <- anova(fit1,fit2)
#       clinicP[clinic,1]= pchisq(ANOVAobj$Deviance[nrow(ANOVAobj)], df=ANOVAobj$Df[nrow(ANOVAobj)],lower.tail=FALSE)
#     }else{
#       CombinedFrameHere$EndoLabel <- as.factor(CombinedFrameHere$EndoLabel)
#       fit1 <- nnet::multinom(f1,data=CombinedFrameHere)
#       fit2 <- nnet::multinom(f2,data=CombinedFrameHere)
#       ANOVAobj <- anova(fit1,fit2)
#       clinicP[clinic,1]= ANOVAobj$`Pr(Chi)`[nrow(ANOVAobj)]
#     }
#   }
#   
#   clinicP2 <- cbind(clinicP,p.adjust(clinicP, method = "bonferroni"),p.adjust(clinicP, method = "BH"))
#   colnames(clinicP2)[c(2,3)] <- c("padj.bonferroni","padj.BH")
#   ### output the table: cilinic features vs p vlaues
#   write.table(clinicP2,paste0(outF1,"/ClinicPvalues.txt"),row.names = TRUE,col.names=TRUE,quote=FALSE,sep="\t")
# }  

# library(limma)
# library(dplyr)
# library(tibble)
# 
# # ---------------------------
# # 1. Filter metadata for groups of interest
# # ---------------------------
# meta <- Olink_merged %>% 
#   filter(GROUP %in% c("AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL") & CSF_LONGITUDINAL_NUMBER == "1") %>% 
#   mutate(
#     ID = CSF_OLINK_MANIFEST,
#     AGE_c = scale(as.numeric(AGE_AT_SYMPTOM_ONSET), center = TRUE, scale = FALSE), # centered age
#     GENDER = factor(GENDER),
#     GROUP = factor(GROUP)
#   )
# 
# # ---------------------------
# # 2. Check counts
# # ---------------------------
# table(meta$GROUP)
# table(meta$GENDER)
# 
# # ---------------------------
# # 3. Make sure expression matrix columns match metadata
# # ---------------------------
# expr <- as.matrix(Olink_ProF[meta$ID, ])
# 
# # ---------------------------
# # 4. Build design matrix (~0 + GROUP)
# # ---------------------------
# design <- model.matrix(~0 + GROUP + AGE_c + GENDER, data = meta)
# colnames(design) <- make.names(colnames(design))
# dim(design)
# 
# # ---------------------------
# # 5. Fit linear model
# # ---------------------------
# fit <- lmFit(expr, design)
# 
# # ---------------------------
# # 6. Define contrast ALS vs Healthy Control
# # ---------------------------
# contrast.matrix <- makeContrasts(
#   ALS_vs_Control = GROUPAMYOTROPHIC.LATERAL.SCLEROSIS - GROUPHEALTHY.CONTROL,
#   levels = design
# )
# 
# # ---------------------------
# # 7. Apply contrast and empirical Bayes
# # ---------------------------
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# 
# # ---------------------------
# # 8. Extract top differentially expressed proteins
# # ---------------------------
# results <- topTable(fit2, coef = "ALS_vs_Control", number = Inf)
# head(results)
# 
# # ---------------------------
# # 9. Filter significant DE proteins (optional)
# # ---------------------------
# DE <- subset(results, adj.P.Val < 0.05 & abs(logFC) > 1)
# nrow(DE)