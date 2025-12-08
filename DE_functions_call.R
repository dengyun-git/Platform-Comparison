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
    pos_vals <- adjustedVariable[adjustedVariable > 0]
    if (length(pos_vals) >= 2) {
      adjustedVariable[zeroID] <<- runif(length(zeroID), min = 0.001, max = 0.01) * sort(pos_vals)[2]
    } else if (length(pos_vals) == 1) {
      adjustedVariable[zeroID] <<- runif(length(zeroID), min = 0.001, max = 0.01) * pos_vals
    } else {
      adjustedVariable[zeroID] <<- 0.001 * maxScore
    }
  }
  
  # Ensure strictly positive before Box-Cox
  adjustedVariable[adjustedVariable <= 0] <<- 0.001 * maxScore
  
  # Perform Box-Cox transformation
  b <- boxcox(lm(adjustedVariable ~ 1))
  
  # Determine the best lambda
  lambda <- b$x[which.max(b$y)]
  new_VariableHere <- (adjustedVariable ^ lambda - 1) / lambda
  
  remove(adjustedVariable)
  
  return(new_VariableHere)
}

# _____________________________________________________________________________________________________________________________________
# Function to generate DE results table using GLM
get_DE_Pvalue_Table_GLM <- function(ProExpF, Merged, covariateList, outPath, whichPlatform, regMethod, varList1, use_random_YN){
  
  # ---------- Helper function to extract effect ----------
  get_effect <- function(fit, protein_var) {
    coef_value <- NA
    
    if ("merMod" %in% class(fit)) {        # lmer, glmer
      coef_names <- names(fixef(fit))
      if (protein_var %in% coef_names) coef_value <- fixef(fit)[protein_var]
      
    } else if ("clmm" %in% class(fit)) {   # ordinal mixed model
      coef_names <- names(fixef(fit))
      if (protein_var %in% coef_names) coef_value <- fixef(fit)[protein_var]
      
    } else {                               # lm, glm, polr, multinom
      coef_names <- names(coef(fit))
      if (protein_var %in% coef_names) coef_value <- coef(fit)[protein_var]
    }
    
    return(coef_value)
  }
  
  # ---------- Initialize result tables ----------
  n_prot <- ncol(ProExpF)
  n_vars <- length(varList1)
  
  result_tbl_EffectSize <- result_tbl_adj <- result_tbl_p <- 
    matrix(NA, nrow = n_prot, ncol = n_vars)
  
  rownames(result_tbl_EffectSize) <- rownames(result_tbl_adj) <- rownames(result_tbl_p) <- colnames(ProExpF)
  colnames(result_tbl_EffectSize) <- colnames(result_tbl_adj) <- colnames(result_tbl_p) <- varList1
  
  # ---------- Loop over main phenotype variables ----------
  for (mainVar in varList1) {
    
    # Determine covariates (customize logic if needed)
    if (mainVar %in% c("ALSvsHC", "DCvsHC", "ALSvsDC", "C9vsNonC9", "BODY_MASS_INDEX")) {
      covariates <- covariateList[c(2,3)]
    } else {
      covariates <- covariateList[c(1,3)]
    }
    
    cat("Running regression for predictor:", mainVar, "\n")
    
    # ---------- Subset samples with valid values ----------
    keepID <- lapply(c(mainVar, covariates), function(v){
      which(!is.na(Merged[, v]) & Merged[, v] != "" & Merged[, v] != "NA")
    })
    keepID_final <- Reduce(intersect, keepID)
    if(length(keepID_final) < 3){  # too few samples to fit
      warning(paste("Too few samples for", mainVar, "- skipping"))
      next
    }
    Merged_Here <- Merged[keepID_final, ]
    
    # ---------- Loop over proteins ----------
    for (i in seq_len(ncol(ProExpF))) {
      protein_var <- colnames(ProExpF)[i]
      covar_str <- paste(covariates, collapse = "+")
      
      # ---------- Build formulas ----------
      if (use_random_YN) {
        full_formula <- as.formula(paste0(mainVar, " ~ `", protein_var, "` + ", covar_str, " + (1|SubjectID_Random)"))
        null_formula <- as.formula(paste0(mainVar, " ~ ", covar_str, " + (1|SubjectID_Random)"))
      } else {
        full_formula <- as.formula(paste0(mainVar, " ~ `", protein_var, "` + ", covar_str))
        null_formula <- as.formula(paste0(mainVar, " ~ ", covar_str))
      }
      
      # ---------- Regression ----------
      tryCatch({
        if (regMethod == "linear") {
          fit <- if (use_random_YN) lmer(full_formula, data = Merged_Here, REML = FALSE) else lm(full_formula, data = Merged_Here)
          null_fit <- if (use_random_YN) lmer(null_formula, data = Merged_Here, REML = FALSE) else lm(null_formula, data = Merged_Here)
          
          LRT <- anova(null_fit, fit)
          result_tbl_p[i, mainVar] <- LRT$`Pr(>F)`[2]
          result_tbl_EffectSize[i, mainVar] <- get_effect(fit, protein_var)
          
        } else if (regMethod == "logistic") {
          fit <- if (use_random_YN) {
            glmer(full_formula, data = Merged_Here, family = binomial, control = glmerControl(optimizer = "bobyqa"))
          } else {
            glm(full_formula, data = Merged_Here, family = binomial)
          }
          null_fit <- if (use_random_YN) {
            glmer(null_formula, data = Merged_Here, family = binomial, control = glmerControl(optimizer = "bobyqa"))
          } else {
            glm(null_formula, data = Merged_Here, family = binomial)
          }
          
          LRT <- tryCatch(anova(null_fit, fit, test = "LRT"), error = function(e) NA)
          result_tbl_p[i, mainVar] <-LRT$`Pr(>Chi)`[2] 
          result_tbl_EffectSize[i, mainVar] <- get_effect(fit, protein_var)
          
        } else if (regMethod == "multinomial") {
          fit <- multinom(full_formula, data = Merged_Here, trace = FALSE)
          null_fit <- multinom(null_formula, data = Merged_Here, trace = FALSE)
          LRT <- 2 * (logLik(fit) - logLik(null_fit))
          df <- attr(logLik(fit), "df") - attr(logLik(null_fit), "df")
          result_tbl_p[i, mainVar] <- pchisq(LRT, df = df, lower.tail = FALSE)
          result_tbl_EffectSize[i, mainVar] <- get_effect(fit, protein_var)
          
        } else if (regMethod == "ordinal") {
          fit <- if (use_random_YN) clmm(full_formula, data = Merged_Here, Hess = TRUE) else polr(full_formula, data = Merged_Here, Hess = TRUE)
          null_fit <- if (use_random_YN) clmm(null_formula, data = Merged_Here, Hess = TRUE) else polr(null_formula, data = Merged_Here, Hess = TRUE)
          
          LRT <- anova(null_fit, fit)
          result_tbl_p[i, mainVar] <- ifelse(length(LRT$`Pr(>Chisq)`) >= 2, LRT$`Pr(>Chisq)`[2], NA)
          result_tbl_EffectSize[i, mainVar] <- get_effect(fit, protein_var)
        }
      }, error = function(e) {
        warning(paste("Error in", mainVar, "protein", protein_var, ":", e$message))
      })
      
    } # end protein loop
    
    # ---------- BH correction ----------
    result_tbl_adj[, mainVar] <- p.adjust(result_tbl_p[, mainVar], method = "BH")
    
  } # end variable loop
  
  # ---------- Save results ----------
  write.csv(result_tbl_p, paste0(outPath, whichPlatform, "_", regMethod, "_RandomEffect", use_random_YN, "_result_tbl_p.csv"))
  write.csv(result_tbl_EffectSize, paste0(outPath, whichPlatform, "_", regMethod, "_RandomEffect", use_random_YN, "_result_tbl_EffectSize.csv"))
  write.csv(result_tbl_adj, paste0(outPath, whichPlatform, "_", regMethod, "_RandomEffect", use_random_YN, "_result_tbl_adjp.csv"))
  
  return(result_tbl_adj)
}


#---------------------------------------------------------------------
# Function to view distribution or composition of clinical variables
ViewClinincalVar <- function(cat_vars, num_vars, Clinic_MetaF_Here){
  # --- 1. Categorical variable plots ---
  cat_plots <- lapply(cat_vars, function(var){
    df <- Clinic_MetaF %>%
      filter(!is.na(.data[[var]]), .data[[var]] != "")
    
    ggplot(df, aes_string(x = var)) +
      geom_bar(fill = "#1f78b4", alpha = 0.7) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(size = 8, color = "black", angle = 20, hjust = 1),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        plot.title = element_blank()  # remove individual plot title
      ) +
      labs(x = var, y = "Count")
  })
  
  names(cat_plots) <- cat_vars
  
  # --- 2. Numerical variable plots ---
  num_plots <- lapply(num_vars, function(var){
    df <- Clinic_MetaF %>%
      filter(!is.na(.data[[var]]))
    
    ggplot(df, aes(x = .data[[var]])) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "#33a02c", alpha = 0.6) +
      geom_density(color = "red", linewidth = 1) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        plot.title = element_blank()  # remove individual plot title
      ) +
      labs(x = var, y = "Density")
  })
  
  names(num_plots) <- num_vars
  
  # --- 3. Combine all plots ---
  all_plots <- c(cat_plots, num_plots)
  
  # Combine using patchwork (example: 3 columns)
  combined_plot <- wrap_plots(all_plots, ncol = 4)
  
  return(combined_plot)
}

#---------------------------------------------------------------------
# Function to calculate ranking metrics for enrichment testing
read_df <- function(file) {
  df <- read.csv(file, check.names = FALSE)
  rownames(df) <- make.unique(df[,1])   # for duplicate protein names for MS
  df <- df[ , -1, drop = FALSE]        
  return(df)
}

Calculate_Ranking_Metric <- function(file1, file2){
  
  ### Load P-value files
  P_DataFrame <- read_df(file1)
  EF_DataFrame  <- read_df(file2)[rownames(P_DataFrame),] ### always pay attention to row names matching 
  
  ### Calculate ranking matrix per clinical variable: -log(p) * sign(effect size)
  Rank_Matrix_AllVar <- (-log(P_DataFrame)) * EF_DataFrame
  
  ### unit test here
  # all(rownames(P_DataFrame) == rownames(EF_DataFrame))
  # Rank_Matrix_AllVar[10,2] == (-log(P_DataFrame[10,2])) * EF_DataFrame[10, 2]
  
  return(Rank_Matrix_AllVar)
}

#_____________________________________________________________________________________________________________________________________
# Function to make bubble plots for pathway enrichment tests
makeBubble <- function(pathHere, patternHere, sizeHere1,sizeHere2, saveMessage){
  
  ### pay attention to the fileList order and resourceLabel order
  fileList <- list.files(path=pathHere, paste0(patternHere, ".*\\.txt$"))
  resourceLabel <- str_replace(fileList, ".*_(.*)\\.txt$", "\\1")
  
  ### read the data frame from each file.
  thisC = 1
  comSFDat=c()
  pThresh = 0.05
  for(fileSg in fileList){
    comSF <- read.csv(paste0(pathHere,fileSg),sep="\t",header=TRUE)[,c("pathway","pval","padj","log2err","ES", "NES", "size", "leadingEdge")] 
    
    comSF <- comSF %>%
      filter(padj < pThresh) %>%
      mutate(YNSignificant = if_else(padj < 0.05, "Y", "N")) %>%
      arrange(padj)
    
    
    addCol <- matrix(rep(resourceLabel[thisC],nrow(comSF)),ncol=1)
    
    comSF <- comSF %>% mutate(resource = addCol) %>% filter(!is.na(pval))
    
    comSFDat <- rbind(comSFDat,comSF)
    
    thisC = thisC+1
    remove(comSF)
  }
  
  comSFDat %<>% mutate(pathway = tolower(gsub("_", " ", gsub("^[^_]*_", "", pathway))))
  
  p.comp <- ggplot(data=comSFDat) + geom_point(aes(x= factor(resource,levels=resourceLabel),y=pathway,size=-log(pval),color = NES,shape=factor(YNSignificant,levels=c("Y","N")))) + xlab("") + ylab("") + 
    labs(color="NES",shape="padj<0.05",size="-log(pval)") + ggtitle(paste0("Enrichment test based on\n",patternHere, " database")) +
    scale_shape_manual(values=c(19,1),labels=c("Y","N")) +
    theme(plot.title=element_text(size = 8 ,face="bold",hjust=0.5),
          axis.text.x = element_text(angle = 20, vjust = 1, hjust=1,size= sizeHere1,face="bold"),axis.text.y = element_text(size = sizeHere2,face="bold"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))  
  
  
  save(comSFDat,file=paste0(pathHere,saveMessage))
  
  print(p.comp)
  return()
}

#_____________________________________________________________________________________________________________________________________
# Function to generate interactive enrichment testing plot
generate_interactive_pathway_plot <- function(file_path, titleMessage) {
  load(file = file_path)   # loads comSFDat
  
  # If no enriched pathways, return a blank plot
  if (!exists("comSFDat") || nrow(comSFDat) == 0) {
    return(
      plot_ly() %>%
        layout(
          title = titleMessage,
          xaxis = list(visible = FALSE),
          yaxis = list(visible = FALSE),
          annotations = list(
            text = "No significantly enriched pathways",
            x = 0.5, y = 0.5,
            showarrow = FALSE,
            font = list(size = 14)
          )
        )
    )
  }
  
  # Prepare hover text (collapse leadingEdge)
  comSFDat$hover_text <- sapply(comSFDat$leadingEdge, function(x) {
    if (is.null(x)) return("")
    paste(x, collapse = ", ")
  })
  
  # Build interactive plot
  ply <- plot_ly(
    comSFDat,
    type = 'scatter',
    mode = 'markers',
    x = ~factor(resource),
    y = ~pathway,
    marker = list(
      size = ~-log(pval),
      sizeref = 0.3,
      sizemode = 'area',
      color = ~NES,
      colorscale = 'RdYlGn',
      reversescale = FALSE,
      colorbar = list(title = 'NES'),
      sizebar = list(title = "-log(pvalue)")
    ),
    text = ~hover_text,
    hovertemplate = paste(
      "<b>Pathway:</b> %{y}<br>",
      "<b>Platform:</b> %{x}<br>",
      "<b>Proteins Driving Pathway:</b> %{text}<extra></extra>"
    )
  ) %>%
    layout(
      title = titleMessage,
      xaxis = list(title = "Platform", tickangle = -45,
                   tickfont = list(size = 10, family = "Arial Black")),
      yaxis = list(title = "", tickfont = list(size = 10, family = "Arial Black")),
      margin = list(b = 100)
    )
  
  return(ply)
}

#_____________________________________________________________________________________________________________________________________
# Function to make volcano plot
make_volcano <- function(df_eff, df_p, trait, platform,
                         p_threshold = 0.05,
                         effect_threshold = 0,
                         top_n = 10) {
  
  df <- data.frame(
    Protein = rownames(df_eff),
    Effect  = df_eff[[trait]],
    Pvalue  = df_p[[trait]]
  ) %>%
    mutate(
      logP = -log10(Pvalue),
      Significance = case_when(
        Pvalue < p_threshold & Effect > 0 ~ "Up",
        Pvalue < p_threshold & Effect < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Subset significant proteins (Up or Down)
  df_sig <- df %>% filter(Significance %in% c("Up", "Down"))
  
  # Select top N for labeling within Up and Down separately
  df_labels <- df_sig %>%
    group_by(Significance) %>%
    arrange(Pvalue) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  p <- ggplot(df, aes(x = Effect, y = logP)) +
    geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Up" = "firebrick", "Down" = "darkgreen", "NS" = "grey50")) +
    geom_vline(xintercept = effect_threshold, linetype = "dashed") +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
    geom_text_repel(
      data = df_labels,
      aes(label = Protein),
      size = 3,
      max.overlaps = Inf
    ) +
    theme_minimal(base_size = 14) +
    xlab(expression("Effect Size ("*beta*")")) +
    ylab("-log10(Adjusted p-value)") +
    ggtitle(paste0(trait, " (", platform, ")")) +
    theme(
      axis.title.x = element_text(size = 10, face = "bold", color = "black"),
      axis.title.y = element_text(size = 10, face = "bold", color = "black"),
      plot.title   = element_text(size = 12, face = "bold", color = "black", hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold", color = "black"),
      legend.text  = element_text(size = 9,  face = "bold", color = "black")
    )
  
  return(p)
}

make_volcano_panel <- function(trait, effect_list, pval_list) {
  
  p1 <- make_volcano(effect_list$Olink,  pval_list$Olink,  trait, "Olink")
  p2 <- make_volcano(effect_list$BEADdel, pval_list$BEADdel, trait, "BEADdel")
  p3 <- make_volcano(effect_list$NONdel,  pval_list$NONdel, trait, "NONdel")
  
  # Combine horizontally
  panel <- p1 + p2 + p3 + plot_layout(ncol = 3)
  
  return(panel)
}

make_volcano_plotly <- function(df_eff, df_p, trait, platform, p_threshold = 0.05, effect_threshold = 0) {
  
  df <- data.frame(
    Protein = rownames(df_eff),
    Effect  = df_eff[[trait]],
    Pvalue  = df_p[[trait]]
  ) %>%
    mutate(
      logP = -log10(Pvalue),
      Significance = case_when(
        Pvalue < p_threshold & Effect > 0 ~ "Up",
        Pvalue < p_threshold & Effect < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  plot_ly(
    df,
    x = ~Effect,
    y = ~logP,
    text = ~paste0("Protein: ", Protein,
                   "<br>Effect: ", round(Effect, 3),
                   "<br>Adjusted p-value: ", signif(Pvalue, 3),
                   "<br>Significance: ", Significance),
    type = "scatter",
    mode = "markers",
    color = ~Significance,
    colors = c("Up" = "firebrick", "Down" = "darkgreen", "NS" = "grey50"),
    marker = list(size = 8, opacity = 0.7)
  ) %>%
    layout(
      title = list(
        text = paste0(trait, " (", platform, ")"),
        font = list(size = 20, color = "black", family = "Arial", bold = TRUE),
        x = 0.5,
        y = 0.96  # lower the title slightly
      ),
      xaxis = list(
        title = list(text = "Effect Size (β)", font = list(size = 16, color = "black", family = "Arial", bold = TRUE)),
        tickfont = list(size = 14, color = "black", family = "Arial", bold = TRUE)
      ),
      yaxis = list(
        title = list(text = "-log10(Adjusted p-value)", font = list(size = 16, color = "black", family = "Arial", bold = TRUE)),
        tickfont = list(size = 14, color = "black", family = "Arial", bold = TRUE)
      ),
      legend = list(
        title = list(text = "Significance", font = list(size = 16, color = "black", family = "Arial", bold = TRUE)),
        font = list(size = 14, color = "black", family = "Arial", bold = TRUE)
      ),
      showlegend = TRUE
    )
}
