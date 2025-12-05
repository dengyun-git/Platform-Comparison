# _____________________________________________________________________________________________________________________________________
# Function to generate DE results table using linear model based on base line samples
get_DE_Pvalue_Table <- function(ProExpF, Merged, covariateList, outPath, whichPlatform, varList1){
  # --- Initialize results container --
  result_tbl_EffectSize <- result_tbl_adj <- result_tbl_p <- matrix(NA, nrow=ncol(ProExpF), ncol=length(varList1))
  rownames(result_tbl_EffectSize) <- rownames(result_tbl_adj) <- rownames(result_tbl_p) <- colnames(ProExpF)
  colnames(result_tbl_EffectSize) <- colnames(result_tbl_adj) <- colnames(result_tbl_p) <- varList1
  
  # --- Loop over main predictors ---
  for(mainVar in varList1){
    
    if(mainVar %in% c("ALSvsHC", "DCvsHC", "ALSvsDC", "C9vsNonC9", "BODY_MASS_INDEX")){
      covariates <- covariateList[c(2,3)] ### HEALTH CONTROL only have records in AGE_AT_SAMPLING
    }else{covariates <- covariateList[c(1,3)]
    }
    
    cat("Running regression for predictor:", mainVar, "\n")
    
    # remove missing entries
    keepID <- vector(mode="list", length(covariates)+1)
    cnt <- 1
    for(j in c(mainVar, covariates)){
      keepID[[cnt]] <- which(!(is.na(Merged[, j]) | Merged[, j] == "NA" | Merged[, j]=="" | is.null(Merged[, j])))
      cnt <- cnt + 1  
    }
    
    Merged_Here <- Merged[Reduce(intersect, keepID), ]
    
    if(length(unique(Merged_Here[,mainVar])) > 1){
      # Loop over proteins
      for(i in colnames(ProExpF)){
        ProS <- Merged_Here[,i]
        
        formula_str1 <- paste("ProS ~", paste(covariates, collapse = " + "))
        formula_str2 <- paste("ProS ~", paste(c(mainVar, covariates), collapse = " + "))
        
        ### linear model with covariates as fixed effect
        fit1 <- lm(as.formula(formula_str1), data = Merged_Here) # formula_str1 is simpler than formula_str2
        fit2 <- lm(as.formula(formula_str2), data = Merged_Here)
        ANOVAobj <- anova(fit1,fit2)
        
        # Extract coefficient for main predictor
        result_tbl_p[i,mainVar] <- ANOVAobj$`Pr(>F)`[2]
        result_tbl_EffectSize[i,mainVar] <- ANOVAobj$F[2]
      }
      result_tbl_adj[, mainVar] <- p.adjust(result_tbl_p[, mainVar], method = "BH")   # Adjust p-values for multiple testing
    }else{
      result_tbl_adj[, mainVar] <- result_tbl_p[, mainVar] <- result_tbl_EffectSize[, mainVar] <- rep(NA, ncol(ProExpF)) # Filling NA due to mainVar less than 2 levels
    }
  }
  
  write.csv(result_tbl_p, paste0(outPath, whichPlatform, "_result_tbl_p.csv"))
  write.csv(result_tbl_EffectSize, paste0(outPath, whichPlatform, "_result_tbl_EffectSize.csv"))
  write.csv(result_tbl_adj, paste0(outPath, whichPlatform, "_result_tbl_adjp.csv"))
  
  return(result_tbl_adj)
}


# _____________________________________________________________________________________________________________________________________
# function to extrac logFC for proteins
extract_logFC <- function(mergedF, ProF_matchClinic, varA, varB) {
  
  # Use base R to get sample names
  SampGroupA <- mergedF$CSF_OLINK_MANIFEST[mergedF$GROUP == varA]
  SampGroupB <- mergedF$CSF_OLINK_MANIFEST[mergedF$GROUP == varB]
  
  # Compute log2 fold change (mean difference)
  logFC <- sapply(colnames(ProF_matchClinic), function(protein) {
    mean(ProF_matchClinic[SampGroupA, protein], na.rm = TRUE) -
      mean(ProF_matchClinic[SampGroupB, protein], na.rm = TRUE)
  })
  
  return(logFC)
}


# _____________________________________________________________________________________________________________________________________
# Function to generate DE results table using linear model based on all the samples

get_DE_Pvalue_Table_AllSamp <- function(ProExpF, Merged, covariateList, outPath, whichPlatform, varList1){
  # --- Initialize results container --
  result_tbl_EffectSize <- result_tbl_adj <- result_tbl_p <- matrix(NA, nrow=ncol(ProExpF), ncol=length(varList1))
  rownames(result_tbl_EffectSize) <- rownames(result_tbl_adj) <- rownames(result_tbl_p) <- colnames(ProExpF)
  colnames(result_tbl_EffectSize) <- colnames(result_tbl_adj) <- colnames(result_tbl_p) <- varList1
  
  # --- Loop over main predictors ---
  for(mainVar in varList1){
    
    if(mainVar %in% c("ALSvsHC", "ALSvsDC", "C9vsNonC9", "BODY_MASS_INDEX")){
      covariates <- covariateList[c(2,3)] ### HEALTH CONTROL only have records in AGE_AT_SAMPLING
    }else{covariates <- covariateList[c(1,3)]
    }
    
    cat("Running regression for predictor:", mainVar, "\n")
    
    # remove missing entries
    keepID <- vector(mode="list", length(covariates)+1)
    cnt <- 1
    for(j in c(mainVar, covariates)){
      keepID[[cnt]] <- which(!(is.na(Merged[, j]) | Merged[, j] == "NA" | Merged[, j]=="" | is.null(Merged[, j])))
      cnt <- cnt + 1  
    }
    
    Merged_Here <- Merged[Reduce(intersect, keepID), ]
    
    if(length(unique(Merged_Here[,mainVar])) > 1){
      # Loop over proteins
      for(i in colnames(ProExpF)){
        ProS <- Merged_Here[,i]
        
        formula_str <- formula(paste0(i, " ~ ", mainVar,"+",paste(covariates,collapse="+"), paste0(" + (1 | rownames(Merged_Here))")))
        
        #fit the model
        fit <- glmer(formula_str, data = Merged_Here, subset = subset, family="lm", 
                     control=glmerControl(optimizer="nloptwrap",optCtrl=list(algorithm = "NLOPT_LN_COBYLA", xtol_rel=1e-6,xtol_abs=1e-10)))
        
        model_summary <- summary(fit)
        result_tbl_p[i,1:4] <- model_summary$coef[ProS,]
        model_stats[i,8] <- ifelse(is.null(unlist(model_summary$optinfo$conv$lme4)),0,1) ## record if a given protein causes convergence issues
      }
      result_tbl_adj[, mainVar] <- p.adjust(result_tbl_p[, mainVar], method = "BH")   # Adjust p-values for multiple testing
    }else{
      result_tbl_adj[, mainVar] <- result_tbl_p[, mainVar] <- result_tbl_EffectSize[, mainVar] <- rep(NA, ncol(ProExpF)) # Filling NA due to mainVar less than 2 levels
    }
  }
  
  write.csv(result_tbl_p, paste0(outPath, whichPlatform, "_result_tbl_p.csv"))
  write.csv(result_tbl_EffectSize, paste0(outPath, whichPlatform, "_result_tbl_EffectSize.csv"))
  write.csv(result_tbl_adj, paste0(outPath, whichPlatform, "_result_tbl_adjp.csv"))
  
  return(result_tbl_adj)
}

# Function to apply Box-Cox transformation
boxcoxTranform <- function(VariableHere, maxScore) {
  # Adjust the variable by subtracting from maxScore to make all values positive
  adjustedVariable <<- maxScore - VariableHere
  zeroID <- which(adjustedVariable <= 0)
  
  # Handle zero values by adding a small positive constant
  if (length(zeroID) != 0) {
    adjustedVariable[zeroID] <<- runif(length(zeroID), min = 0.001, max = 0.01) * sort(adjustedVariable)[2]
  }
  
  # Perform Box-Cox transformation
  b <- boxcox(lm(adjustedVariable ~ 1))
  
  # Determine the best lambda
  lambda <- b$x[which.max(b$y)]
  new_VariableHere <- (adjustedVariable ^ lambda - 1) / lambda
  
  return(new_VariableHere)
}