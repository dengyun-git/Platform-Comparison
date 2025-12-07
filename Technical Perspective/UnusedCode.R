# Unused code but hope to record
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