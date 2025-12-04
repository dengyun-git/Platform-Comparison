library(tidyverse)
library(magrittr)
library(readxl)

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

### Assign the right class of the meta variables
Clinic_MetaF %<>% dplyr::mutate(GROUP = as.factor(GROUP), GENDER = as.factor(GENDER), WEAKNESS_SITE = as.factor(WEAKNESS_SITE))

### Merge the clinical meta and proteomic data

MS_BEADdel_merged <- MS_BEADdel_ProF %>%
  rownames_to_column("ID") %>%
  left_join(Clinic_MetaF %>% rownames_to_column("ID"), by = "ID") %>%
  column_to_rownames("ID") %>% mutate(GENDER = as.factor(GENDER))

MS_NONdel_merged <- MS_NONdel_ProF %>%
  rownames_to_column("ID") %>%
  left_join(Clinic_MetaF %>% rownames_to_column("ID"), by = "ID") %>%
  column_to_rownames("ID") %>% mutate(GENDER = as.factor(GENDER))

Olink_merged <- Olink_ProF %>%
  rownames_to_column("ID") %>%
  left_join(Clinic_MetaF %>% rownames_to_column("ID"), by = "ID") %>%
  column_to_rownames("ID") %>% mutate(GENDER = as.factor(GENDER))

########################
########################
# DE Between MS & OLINK
########################
########################
# Direct use paired t test

individualCor_BEADdel_Olink_Pro <- getIndividualPairTtest(MS_BEADdel_ProF, Olink_ProF, "Protein")

individualCor_BEADdel_Olink_Pro <- getIndividualRelation(MS_BEADdel_ProF, Olink_ProF, "Protein", "spearman")
individualCor_NONdel_Olink_Pro <- getIndividualRelation(MS_NONdel_ProF, Olink_ProF, "Protein", "spearman")
individualCor_BEADdel_NONdel_Pro <- getIndividualRelation(MS_BEADdel_ProF, MS_NONdel_ProF, "Protein", "spearman")

length(which(abs(individualCor_BEADdel_Olink_Pro)>0.5))
length(which(abs(individualCor_NONdel_Olink_Pro)>0.5))
length(which(abs(individualCor_BEADdel_NONdel_Pro)>0.5))

hist(individualCor_BEADdel_Olink, breaks=100)
hist(individualCor_NONdel_Olink, breaks=100)
hist(individualCor_BEADdel_NONdel, breaks=100)

length(intersect(colnames(MS_BEADdel_ProF), colnames(Olink_ProF)))
length(Reduce(intersect, list(colnames(MS_BEADdel_ProF), colnames(MS_NONdel_ProF), colnames(Olink_ProF))))

individualAdjP_BEADdel_Olink_Pro <- getIndividualPairTtest(MS_BEADdel_ProF, Olink_ProF, "Protein")

individualCor_BEADdel_Olink_Samp <- getIndividualRelation(MS_BEADdel_ProF, Olink_ProF, "Sample", "spearman")
individualCor_NONdel_Olink_Samp <- getIndividualRelation(MS_NONdel_ProF, Olink_ProF, "Sample", "spearman")
individualCor_BEADdel_NONdel_Samp <- getIndividualRelation(MS_BEADdel_ProF, MS_NONdel_ProF, "Sample", "spearman")

hist(individualCor_BEADdel_Olink_Samp, breaks=100)
hist(individualCor_NONdel_Olink_Samp, breaks=100)
hist(individualCor_BEADdel_NONdel_Samp, breaks=100)

individualAdjP_BEADdel_Olink_Samp <- getIndividualPairTtest(MS_BEADdel_ProF, Olink_ProF, "Sample")
commonSampID <- names(individualCor_BEADdel_Olink_Samp)[which(individualCor_BEADdel_Olink_Samp>0.05)]
View(Clinic_MetaF[commonSampID,])

# Compare with global, whether there is something important for the sample characteristics

# Direct compare the correlation coefficient

MS_BEADdel_DSvsHC <- testbinaryOutcome(outcome = "GROUP", covariates = c("AGE_AT_SYMPTOM_ONSET","GENDER"), proteins_matrix = MS_BEADdel_ProF, jointdata_frame_all = MS_BEADdel_merged,  subset=MS_BEADdel_merged$LONGITUDINAL_ENCOUNTER_NUMBER == "1" & MS_BEADdel_merged$GROUP %in% c("AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL"))

MS_NONdel_DSvsHC <- testbinaryOutcome(outcome = "GROUP", covariates = c("AGE_AT_SYMPTOM_ONSET","GENDER"), proteins_matrix = MS_NONdel_ProF, jointdata_frame_all = MS_NONdel_merged,  subset=MS_NONdel_merged$LONGITUDINAL_ENCOUNTER_NUMBER == "1" & MS_NONdel_merged$GROUP %in% c("AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL"))

Olink_DSvsHC <- testbinaryOutcome(outcome = "GROUP", covariates = c("AGE_AT_SYMPTOM_ONSET","GENDER"), proteins_matrix = Olink_ProF, jointdata_frame_all = Olink_merged,  subset= Olink_merged$LONGITUDINAL_ENCOUNTER_NUMBER == "1" & Olink_merged$GROUP %in% c("AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL"))

#######################
#######################
# DE Between Platforms
#######################
#######################

## 1) logistic regression models (binary outcome)
testbinaryOutcome <- function(outcome,covariates, proteins_matrix, jointdata_frame_all,subset){
  #function for testing log protein expression against outcome, with covariates included, in a logistic model
  #protein names to test are stored as the first column of proteins_list
  #protein abundance, outcome and covariates must be stored in jointdata_frame
  #output logistic model parameters (estimates, standards, t-stats and p-values)
  #stored in a named data frame of size Nprotein x 4 (i.e. proteins x {statistics}) 
  
  #initialize the named array for output
  model_stats <- array(NA,dim=c(ncol(proteins_matrix),7)) %>% as.data.frame() 
  dimnames(model_stats) <- list(colnames(proteins_matrix),c("estimate","std.error","z.value","p.value", "p_adj", "sig", "sig_2"))
  
  #Change the binary category to 0 and 1
  df <- jointdata_frame_all[subset, ]
  orig <- as.character(df[[outcome]])
  orig <- trimws(orig)
  df[,outcome] <- factor(ifelse(orig == unique(orig)[1], "1", "0"))
  jointdata_frame <- df
  
  #cycle through protein names
  for (i in 1:ncol(proteins_matrix)){
    
    # set the protein name
    protein_to_test <- names(proteins_matrix)[i]
    
    #construct the formula for the model
    fm_to_test1 <- formula(paste0(outcome, " ~ ", protein_to_test,"+",paste(covariates,collapse="+")))
    fm_to_test2 <- formula(paste0(outcome, " ~ ", paste(covariates,collapse="+")))
    
    #fit the model
    model1 <- glm(fm_to_test, data = jointdata_frame, family="binomial")
    model2 <- glm(fm_to_test, data = jointdata_frame, family="binomial")
    
    #extract out model statistics and store in a array
    model_summary <- summary(model)
    model_stats[i,1:4] <- model_summary$coef[protein_to_test,]
    print(paste("Iteration:", i))
    
  }
  
  ##apply a multiple testing correction
  model_stats[,5]   <- p.adjust(model_stats[,4],method="BH")
  
  ## at 0.05 significance level
  model_stats[,6]   <- ifelse(model_stats[,4]<= 0.05, 'yes', 'no')
  
  ## at FDR < 0.05 significance level
  model_stats[,7]  <- ifelse(model_stats[,5]<= 0.05, 'yes', 'no')
  
  model_stats %<>% arrange(p_adj) 
  
  #return array
  return(model_stats)
  
}


# “GROUP” Disease group i.e. ALS, PLS, HC, DC, AAR. Diagnostic biomarker theme. In particular, ALS versus healthy and disease control.
# 
# “PATHOGENIC_VARIANT” Gene mutations causing ALS. In particular, C9orf72 ALS versus negative ALS.
# 
# "WEAKNESS_SITE” Site of first weakness. In particular, bulbar versus non-bulbar.
# 
# "ALSFRSR_BULBAR"                  
# “ALSFRSR_FINE_MOTOR"
# "ALSFRSR_GROSS_MOTOR"
# "ALSFRSR_RESPIRATORY"             
# “ALSFRSR_SCORE"
# "ALSFRSR_RATE” Clinical rating scale of ALS disease severity. 

# "ECAS_ALS_SUBSCORE"
# "ECAS_NONALS_SUBSCORE" 
# “ECAS_SCORE” Clinical rating scale of cognitive impairment.
# 
# Disease group: ALS = 918, HC = 144, DC = 106, PLS = 68, FD = 4 (perhaps exclude FD due to the low sample size). – please note that many will have return visits and classification should either exclude these or take into account the repeated measures. I agree to the removal of FD given the low number (FTD?)
# PATHOGENIC_VARIANT: C9ORF72 = 94, NEGATIVE = 124, and the rest of the categories almost have fewer than 10 samples, so I suggest excluding them from the classification analysis. – David, what do you think about this? We could group together non-C9ORF72 present values as “other” i.e. non-C9. There is a question as to whether we can include those where we do not have a known genetic status (i.e. absent value) in the “other” group as well, but I don’ think that really works.
# WEAKNESS_SITE: BULBAR = 220, LOWER LIMB = 403, UPPER LIMB = 339, while other categories have only 5-30 samples, so I recommend not including them. – there is probably a secondary or even tertiary level to this i.e. bulbar and spinal (spinal including upper limb, lower limb, truncal, respiratory and generalised) and bulbar vs non-bulbar that we could use to classify. Again, David thoughts?

library(limma)
library(dplyr)
library(tibble)

# ---------------------------
# 1. Filter metadata for groups of interest
# ---------------------------
meta <- Olink_merged %>% 
  filter(GROUP %in% c("AMYOTROPHIC LATERAL SCLEROSIS", "HEALTHY CONTROL") & CSF_LONGITUDINAL_NUMBER == "1") %>% 
  mutate(
    ID = CSF_OLINK_MANIFEST,
    AGE_c = scale(as.numeric(AGE_AT_SYMPTOM_ONSET), center = TRUE, scale = FALSE), # centered age
    GENDER = factor(GENDER),
    GROUP = factor(GROUP)
  )

# ---------------------------
# 2. Check counts
# ---------------------------
table(meta$GROUP)
table(meta$GENDER)

# ---------------------------
# 3. Make sure expression matrix columns match metadata
# ---------------------------
expr <- as.matrix(Olink_ProF[meta$ID, ])

# ---------------------------
# 4. Build design matrix (~0 + GROUP)
# ---------------------------
design <- model.matrix(~0 + GROUP + AGE_c + GENDER, data = meta)
colnames(design) <- make.names(colnames(design))
dim(design)

# ---------------------------
# 5. Fit linear model
# ---------------------------
fit <- lmFit(expr, design)

# ---------------------------
# 6. Define contrast ALS vs Healthy Control
# ---------------------------
contrast.matrix <- makeContrasts(
  ALS_vs_Control = GROUPAMYOTROPHIC.LATERAL.SCLEROSIS - GROUPHEALTHY.CONTROL,
  levels = design
)

# ---------------------------
# 7. Apply contrast and empirical Bayes
# ---------------------------
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# ---------------------------
# 8. Extract top differentially expressed proteins
# ---------------------------
results <- topTable(fit2, coef = "ALS_vs_Control", number = Inf)
head(results)

# ---------------------------
# 9. Filter significant DE proteins (optional)
# ---------------------------
DE <- subset(results, adj.P.Val < 0.05 & abs(logFC) > 1)
nrow(DE)

# Disease group: ALS = 918, HC = 144, DC = 106, PLS = 68, FD = 4 (perhaps exclude FD due to the low sample size). – please note that many will have return visits and classification should either exclude these or take into account the repeated measures. I agree to the removal of FD given the low number (FTD?)
# PATHOGENIC_VARIANT: C9ORF72 = 94, NEGATIVE = 124, and the rest of the categories almost have fewer than 10 samples, so I suggest excluding them from the classification analysis. – David, what do you think about this? We could group together non-C9ORF72 present values as “other” i.e. non-C9. There is a question as to whether we can include those where we do not have a known genetic status (i.e. absent value) in the “other” group as well, but I don’ think that really works.
# WEAKNESS_SITE: BULBAR = 220, LOWER LIMB = 403, UPPER LIMB = 339, while other categories have only 5-30 samples, so I recommend not including them. – there is probably a secondary or even tertiary level to this i.e. bulbar and spinal (spinal including upper limb, lower limb, truncal, respiratory and generalised) and bulbar vs non-bulbar that we could use to classify. Again, David thoughts?

WEAKNESS_SITE: bulbarVSspinal
PATHOGENIC_VARIANT: C9VSnonC9

# ALSFRSR_rate – log transform – the key variable!
# ALSFRSR_total – 48 minus score then Box-Cox transformation. Raw values >48 are impossible and should be removed
# ECAS ALS score – 100 minus score then Box-Cox (I suspect – distribution will need to be inspected). Raw values >100 are impossible and should be removed
# ECAS total score – 134 minus score then Box-Cox (I suspect – distribution will need to be inspected). Raw values >134 are impossible and should be removed

confList1 = c("GROUP", "WEAKNESS_SITE", "PATHOGENIC_VARIANT", "BODY_MASS_INDEX", "ALSFRSR_BULBAR", "ALSFRSR_FINE_MOTOR", "ALSFRSR_GROSS_MOTOR", "ALSFRSR_RESPIRATORY", "ALSFRSR_RATE", "ECAS_ALS_SUBSCORE", "ECAS_SCORE")
confList2 = c("AGE_AT_SYMPTOM_ONSET","GENDER", "COHORT")
confList = c(confList1,confList2)

for(id in c(1,2,3,4,5,6,7,10,11)){
  ClinicFrame[,confList[id]] <- as.factor(ClinicFrame[,confList[id]])
} ### change some variables as factor

### p container -- independent variable confounder vs dependent variable endotype
pContainer = matrix(NA,nrow=length(confList),ncol=1)
rownames(pContainer) = c(confList)
colnames(pContainer) = "pvalue"

for (confounder in confList){
  
  ClinicFrameHere <- ClinicFrame[!is.na(ClinicFrame[,confounder]),]
  
  if(confounder %in% confList1){
    f1 <- as.formula(EndoLabel~1)
    f2 <- as.formula(paste("EndoLabel~",confounder))
  }else{
    f1 <- as.formula(EndoLabel~cohort_name)
    f2 <- as.formula(paste("EndoLabel~cohort_name+",confounder))
  }
  
  ### tackle with the different cases of endotype numbers
  if(length(unique(ClinicFrameHere$EndoLabel))==2){
    fit1 <- glm(f1,data=ClinicFrameHere,family = binomial(link = "logit"))
    fit2 <- glm(f2,data=ClinicFrameHere,family = binomial(link = "logit"))
    ANOVAobj <- anova(fit1,fit2)
    pContainer[confounder,1] = pchisq(ANOVAobj$Deviance[nrow(ANOVAobj)], df=ANOVAobj$Df[nrow(ANOVAobj)],lower.tail=FALSE)
  }else{
    fit1 <- nnet::multinom(f1,data=ClinicFrameHere)
    fit2 <- nnet::multinom(f2,data=ClinicFrameHere)
    ANOVAobj <- anova(fit1,fit2)
    pContainer[confounder,1]= ANOVAobj$`Pr(Chi)`[nrow(ANOVAobj)]
  } 
  
  if(confounder == "sf_iknee_proc_batch"){
    fit3 <- lmerTest::lmer(as.numeric(EndoLabel) ~ get(confounder) + (1|cohort_name), data=ClinicFrameHere)
    ANOVAobj <- anova(fit3)
    pContainer[confounder,1]=ANOVAobj$`Pr(>F)`
  }
  
  
  
  for(clinicInd in c(2,5,8,11,12,13,14,3,4,9,10)){
    clinic = clinics[clinicInd] ### define which clinic to explore this loop
    
    keepId <- which(!(is.na(CombinedFrameThis[,clinic]) | CombinedFrameThis[,clinic]=="NA")) 
    CombinedFrameHere <- CombinedFrameThis[keepId,] ### remove NA entries before regression
    
    ### categorical clinical features
    if(clinicInd %in% c(2,5,8,11,12,13,14)){
      CombinedFrameHere[,clinic] <- as.factor(CombinedFrameHere[,clinic])
      
      ### bar chart displaying clinical feature distribution within endotypes
      barDat = CombinedFrameHere[,c("EndoLabel",clinic)]
      compT <- table(barDat)/rowSums(table(barDat))
      addCol <- matrix(unlist(sapply(1:nrow(barDat),function(x){
        endo = as.character(barDat[x,"EndoLabel"])
        clin = as.character(barDat[x,clinic])
        compT[endo,clin]})),ncol=1)
      colnames(addCol) = "Composition"
      barDat <- cbind(barDat,addCol)
      
      if(length(unique(EndoLabel))==2){
        p1 <- ggplot(data = barDat) + geom_bar(aes(x=as.factor(get(clinic)),y=Composition,fill=as.factor(EndoLabel)),stat="identity",position=position_dodge()) +
          ggtitle(paste0(clinics[clinicInd]," comparision for endotypes")) + labs(x="",y="Composition",fill="Endotype") + scale_fill_manual(name = "Endotype", labels = HighLowIntensityLabel,values = HighLowIntensityColor) +
          theme(plot.title=element_text(size = 15,face="bold",hjust=0.5),legend.position="bottom",
                axis.text.x = element_text(size=12,face="bold",angle = 45,hjust = 1),axis.title.x =element_text(size = 12,face="bold"),
                legend.title = element_text(size = 13,face="bold"),legend.text = element_text(size = 12,face="bold"),
                axis.text.y =element_text(size=13,face="bold"),axis.title.y =element_text(size=12,face="bold"))
      }else{ ### more than 2 endotypes
        p1 <- ggplot(data = barDat) + geom_bar(aes(x=as.factor(get(clinic)),y=Composition,fill=as.factor(EndoLabel)),stat="identity",position=position_dodge()) +
          ggtitle(paste0(clinics[clinicInd]," comparision for endotypes")) + labs(x="",y="Composition",fill="Endotype") + scale_fill_discrete(name = "Endotype") +
          theme(plot.title=element_text(size = 15,face="bold",hjust=0.5),legend.position="bottom",
                axis.text.x = element_text(size=12,face="bold",angle = 45,hjust = 1),axis.title.x =element_text(size = 12,face="bold"),
                legend.title = element_text(size = 13,face="bold"),legend.text = element_text(size = 12,face="bold"),
                axis.text.y =element_text(size=13,face="bold"),axis.title.y =element_text(size=12,face="bold"))
      }     
      print(p1)
      
      ### an example plot of KL grade
      # p1 <- ggplot(data = barDat) + geom_bar(aes(x=as.factor(get(clinic)),y=Composition,fill=as.factor(EndoLabel)),stat="identity",position=position_dodge()) +
      #   ggtitle("Ordinal KL Grade Comparision for Endotypes") + labs(x="Ordinal KL grade", y="Composition",fill="endotype") + 
      #   scale_fill_discrete(name = "Endotype", labels = c("High Intensity Endotype", "Low Intensity Endotype")) +
      #   theme(plot.title=element_text(size = 15,face="bold",hjust=0.5),legend.position="bottom",
      #         axis.text.x = element_text(size=15,face="bold"),axis.title.x =element_text(size = 12,face="bold"),
      #         legend.title = element_text(size = 13,face="bold"),legend.text = element_text(size = 12,face="bold"),
      #         axis.text.y =element_text(size=13,face="bold"),axis.title.y =element_text(size=12,face="bold"))
      # print(p1)
      
      ### PCA plot, overlay with clinical features
      PCAoverlayDat <- cbind(barDat,PCs[keepId,])
      p1.1 <- ggplot(data=PCAoverlayDat) + geom_point(aes(x=PC1,y=PC2,color=as.factor(get(clinic)),shape=as.factor(EndoLabel))) +
        labs(shape="endotype",color=clinic) +
        theme(axis.text.x=element_blank(),axis.text.y=element_blank()) + ggtitle(paste0(clinics[clinicInd]," overlay on PCA"))
      print(p1.1)
      
      ### umap, overlap with clinical features
      df <- data.frame(x = umap[keepId,1],y =umap[keepId,2],barDat)
      colnames(df)[c(1,2)] = c("umap.D1","umap.D2")
      p1.2 <- ggplot(df, aes(x=umap.D1, y=umap.D2, color=as.factor(get(clinic)),shape=as.factor(EndoLabel))) + geom_point() + ggtitle(paste0(clinics[clinicInd]," overlay on umap")) + labs(color=clinics[clinicInd],shape="endotype") +
        labs(shape="endotype",color=clinic) +
        theme(plot.title=element_text(size=14,hjust=0.5),legend.text=element_text(size = 10), axis.title=element_text(size = 12))
      print(p1.2)
      
    }else{### continuous clinical features
      ### violin plot displaying clinical feature distribution within endotypes
      violinDat = CombinedFrameHere[,c("EndoLabel",clinic)]
      nn=length(unique(CombinedFrameHere$EndoLabel))
      
      if(length(unique(EndoLabel))==2){
        p2 <- ggplot(data=violinDat,aes(x=as.factor(EndoLabel),y=get(clinic),fill=as.factor(EndoLabel))) + geom_violin(width=1) + geom_boxplot(width=0.1, color="grey", alpha=0.5) +
          labs(title=(paste0(clinic," distribution for endotypes")),x="Endotype",fill="Endotype") + ylab("distribution") + scale_fill_manual(name = "Endotype", labels = HighLowIntensityLabel,values = HighLowIntensityColor) +
          theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.text.x = element_blank(),legend.position="bottom") +
          geom_vline(xintercept = 0.5:(nn+0.5),linetype="dotted")
      }else{
        p2 <- ggplot(data=violinDat,aes(x=as.factor(EndoLabel),y=get(clinic),fill=as.factor(EndoLabel))) + geom_violin(width=1) + geom_boxplot(width=0.1, color="grey", alpha=0.5) +
          labs(title=(paste0(clinic," distribution for endotypes")),x="Endotype",fill="Endotype") + ylab("distribution") + scale_fill_discrete(name = "Endotype") +
          theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.text.x = element_blank(),legend.position="bottom") +
          geom_vline(xintercept = 0.5:(nn+0.5),linetype="dotted") 
      }
      
      ### an example plot of womac pain
      # p2 <- ggplot(data=violinDat,aes(x=as.factor(EndoLabel),y=get(clinic),fill=as.factor(EndoLabel))) + geom_violin(width=1) + geom_boxplot(width=0.1, color="grey", alpha=0.5) +
      #   labs(title=(paste0("WOMAC pain score"," distribution for endotypes")),x="endotype",fill="endotype") + ylab("distribution") +
      #   scale_fill_discrete(name = "Endotype", labels = c("High Intensity Endotype", "Low Intensity Endotype")) +
      #   theme(plot.title=element_text(size = 15,face="bold",hjust=0.5),legend.position="bottom",
      #         panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.text.x = element_blank(),
      #         axis.title.x =element_text(size = 12,face="bold"),
      #         legend.title = element_text(size = 13,face="bold"),legend.text = element_text(size = 12,face="bold"),
      #         axis.text.y =element_text(size=13,face="bold"),axis.title.y =element_text(size=12,face="bold")) +
      #   geom_vline(xintercept = 0.5:(nn+0.5),linetype="dotted")
      print(p2)
      
      PCAoverlayDat <- cbind(violinDat,PCs[keepId,])
      p2.1 <- ggplot(data=PCAoverlayDat) + geom_point(aes(x=PC1,y=PC2,color=get(clinic),shape=as.factor(EndoLabel))) +
        theme(axis.text.x=element_blank(),axis.text.y=element_blank()) + ggtitle(paste0(clinic," overlay on PCA")) +
        labs(color=clinic,shape="endotype")
      print(p2.1)
      
      ### umap, overlap with clinical features
      df <- data.frame(x = umap[keepId,1],y = umap[keepId,2],violinDat)
      p2.2 <- ggplot(df, aes(x=x, y=y, color=get(clinic),shape=as.factor(EndoLabel))) + geom_point() + ggtitle(paste0(clinic," overlay on umap")) +
        theme(plot.title=element_text(size=14,hjust=0.5),legend.text=element_text(size = 10), axis.title=element_text(size = 12),axis.text.x=element_blank(),axis.text.y=element_blank()) +
        labs(color=clinic,shape="endotype",x="Dimension1",y="Dimension2")
      print(p2.2)
    }
    
    f1 <- as.formula(EndoLabel~cohort_name)
    f2 <- as.formula(paste("EndoLabel~cohort_name+",clinic))
    # f1 <- as.formula(EndoLabel~cohort_name+totalProHere)
    # f2 <- as.formula(paste("EndoLabel~cohort_name+totalProHere+",clinic))
    
    ### tackle with different number of endotypes
    if(length(unique(CombinedFrameHere$EndoLabel))==2){
      CombinedFrameHere$EndoLabel <- as.factor(sapply(CombinedFrameHere$EndoLabel,function(x){if(x!=1){0}else{x}}))
      fit1 <- glm(f1,family = binomial(link = "logit"),data=CombinedFrameHere)
      fit2 <- glm(f2,family = binomial(link = "logit"),data=CombinedFrameHere)
      ANOVAobj <- anova(fit1,fit2)
      clinicP[clinic,1]= pchisq(ANOVAobj$Deviance[nrow(ANOVAobj)], df=ANOVAobj$Df[nrow(ANOVAobj)],lower.tail=FALSE)
    }else{
      CombinedFrameHere$EndoLabel <- as.factor(CombinedFrameHere$EndoLabel)
      fit1 <- nnet::multinom(f1,data=CombinedFrameHere)
      fit2 <- nnet::multinom(f2,data=CombinedFrameHere)
      ANOVAobj <- anova(fit1,fit2)
      clinicP[clinic,1]= ANOVAobj$`Pr(Chi)`[nrow(ANOVAobj)]
    }
  }
  
  clinicP2 <- cbind(clinicP,p.adjust(clinicP, method = "bonferroni"),p.adjust(clinicP, method = "BH"))
  colnames(clinicP2)[c(2,3)] <- c("padj.bonferroni","padj.BH")
  ### output the table: cilinic features vs p vlaues
  write.table(clinicP2,paste0(outF1,"/ClinicPvalues.txt"),row.names = TRUE,col.names=TRUE,quote=FALSE,sep="\t")
  