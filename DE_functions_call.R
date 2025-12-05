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
# function to apply Box-Cox transformation
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

# _____________________________________________________________________________________________________________________________________
# Function to generate DE results table using GLM
get_DE_Pvalue_Table_GLM <- function(ProExpF, Merged, covariateList, outPath, whichPlatform, regMethod, varList1, use_random_YN){
  
  get_effect <- function(fit, protein_var) {
    # Safe version: returns NA if protein_var not in model coefficients
    coef_value <- NA
    
    if("merMod" %in% class(fit)) {        # lmer, glmer
      coef_names <- names(fixef(fit))
      if(protein_var %in% coef_names) {
        coef_value <- fixef(fit)[protein_var]
      }
      
    } else if("clmm" %in% class(fit)) {   # ordinal mixed model
      coef_names <- names(fixef(fit))
      if(protein_var %in% coef_names) {
        coef_value <- fixef(fit)[protein_var]
      }
      
    } else {                              # lm, glm, polr, multinom
      coef_names <- names(coef(fit))
      if(protein_var %in% coef_names) {
        coef_value <- coef(fit)[protein_var]
      }
    }
    
    # Return the coefficient (or NA if not found)
    return(coef_value)
  }
  
  # ---------- Initialize result tables ----------
  result_tbl_EffectSize <- result_tbl_adj <- result_tbl_p <-
    matrix(NA, nrow=ncol(ProExpF), ncol=length(varList1))
  
  rownames(result_tbl_EffectSize) <- rownames(result_tbl_adj) <- rownames(result_tbl_p) <- colnames(ProExpF)
  colnames(result_tbl_EffectSize) <- colnames(result_tbl_adj) <- colnames(result_tbl_p) <- varList1
  
  # ===============================================================
  # ---------- Loop over main phenotype variables ----------
  # ===============================================================
  for(mainVar in varList1){
    
    if(mainVar %in% c("ALSvsHC","DCvsHC","ALSvsDC","C9vsNonC9","BODY_MASS_INDEX")){
      covariates <- covariateList[c(2,3)]
    } else {
      covariates <- covariateList[c(1,3)]
    }
    
    cat("Running regression for predictor:", mainVar, "\n")
    
    # ---------- Remove samples with NA ----------
    keepID <- lapply(c(mainVar, covariates), function(v){
      which(!is.na(Merged[, v]) & Merged[, v] != "" & Merged[, v] != "NA")
    })
    Merged_Here <- Merged[Reduce(intersect, keepID), ]
    
    # =====================================================
    # ---------- LOOP over PROTEINS ----------
    # =====================================================
    for(i in seq_along(colnames(ProExpF))){
      protein_var <- colnames(ProExpF)[i]
      
      covar_str <- paste(covariates, collapse="+")
      
      if(use_random_YN){
        full_formula <- as.formula(
          paste0(mainVar, " ~ `", protein_var, "` + ", covar_str, " + (1|SubjectID_Random)")
        )
        null_formula <- as.formula(
          paste0(mainVar, " ~ ", covar_str, " + (1|SubjectID_Random)")
        )
      } else {
        full_formula <- as.formula(
          paste0(mainVar, " ~ `", protein_var, "` + ", covar_str)
        )
        null_formula <- as.formula(
          paste0(mainVar, " ~ ", covar_str)
        )
      }
      
      # =============== LINEAR MIXED MODEL ===============
      if(regMethod == "linear"){
        if(use_random_YN){
          fit <- lmer(full_formula, data=Merged_Here, REML=FALSE)
          null_fit <- lmer(null_formula, data=Merged_Here, REML=FALSE)
        } else {
          fit <- lm(full_formula, data=Merged_Here)
          null_fit <- lm(null_formula, data=Merged_Here)
        }
        
        # LRT
        LRT <- anova(null_fit, fit)
        result_tbl_p[i, mainVar] <- LRT$`Pr(>Chisq)`[2]
        result_tbl_EffectSize[i, mainVar] <- get_effect(fit, protein_var)
      }
      
      # =============== LOGISTIC MIXED MODEL ===============
      else if(regMethod == "logistic"){
        
        if(use_random_YN){
          fit <- glmer(full_formula, data=Merged_Here, family=binomial,
                       control=glmerControl(optimizer="bobyqa"))
          null_fit <- glmer(null_formula, data=Merged_Here, family=binomial,
                            control=glmerControl(optimizer="bobyqa"))
        } else {
          fit <- glm(full_formula, data=Merged_Here, family=binomial)
          null_fit <- glm(null_formula, data=Merged_Here, family=binomial)
        }
        
        LRT <- anova(null_fit, fit, test="LRT")
        result_tbl_p[i, mainVar] <- LRT$`Pr(>Chi)`[2]
        result_tbl_EffectSize[i, mainVar] <- get_effect(fit, protein_var)
      }
      
      # =============== MULTINOMIAL — NO RANDOM EFFECTS AVAILABLE ===============
      else if(regMethod == "multinomial"){
        
        fit <- multinom(full_formula, data=Merged_Here, trace=FALSE)
        null_fit <- multinom(null_formula, data=Merged_Here, trace=FALSE)
        
        LRT <- 2 * (logLik(fit) - logLik(null_fit))
        df <- attr(logLik(fit),"df") - attr(logLik(null_fit),"df")
        pval <- pchisq(LRT, df=df, lower.tail=FALSE)
        
        result_tbl_p[i, mainVar] <- pval
        
        # average effect across categories
        result_tbl_EffectSize[i, mainVar] <- get_effect(fit, protein_var)
      }
      
      # =============== ORDINAL MIXED MODEL ===============
      else if(regMethod == "ordinal"){
        
        if(use_random_YN){
          fit <- clmm(full_formula, data=Merged_Here, Hess=TRUE)
          null_fit <- clmm(null_formula, data=Merged_Here, Hess=TRUE)
        } else {
          fit <- polr(full_formula, data=Merged_Here, Hess=TRUE)
          null_fit <- polr(null_formula, data=Merged_Here, Hess=TRUE)
        }
        
        LRT <- anova(null_fit, fit)
        result_tbl_p[i, mainVar] <- LRT$`Pr(>Chisq)`[2]
        result_tbl_EffectSize[i, mainVar] <- get_effect(fit, protein_var)
      }
    } # end protein loop
    
    # -------------- BH correction for this variable --------------
    result_tbl_adj[, mainVar] <- p.adjust(result_tbl_p[, mainVar], method="BH")
    
  } # end var loop
  
  # save to files
  write.csv(result_tbl_p, paste0(outPath, whichPlatform, "_", regMethod, "_RandomEffect", use_random_YN, "_result_tbl_p.csv"))
  write.csv(result_tbl_EffectSize, paste0(outPath, whichPlatform, "_", regMethod, "_RandomEffect", use_random_YN, "_result_tbl_EffectSize.csv"))
  write.csv(result_tbl_adj, paste0(outPath, whichPlatform, "_", regMethod, "_RandomEffect", use_random_YN, "_result_tbl_adjp.csv"))
  
  return(result_tbl_adj)
}
