##---------------------------------------------------------------------------------------------------------##
## CREATE FUNCTIONS:
##---------------------------------------------------------------------------------------------------------##

## 1) linear regression models (continuous outcome) 
testContinuousOutcome <- function(outcome,covariates,proteins_list, proteins_matrix, jointdata_frame,subset){
  #function for testing log protein expression against outcome, with covariates included, in a linear model
  #protein names to test are stored as the first column of proteins_list
  #protein abundance, outcome and covariates must be stored in jointdata_frame
  #output linear model parameters (estimates, standards, t-stats and p-values)
  #stored in a named data frame of size Nprotein x 4 (i.e. proteins x {statistics}) 
  
  #initialize the named array for output
  model_stats <- array(NA,dim=c(nrow(proteins_list),7))
  dimnames(model_stats) <- list(proteins_list[,1],c("estimate","std.error","t.value","p.value", "p_adj", "sig", "sig_2"))
  
  #cycle through protein names
  for (i in 1:ncol(proteins_matrix)){
    
    # set the protein name and specify log transformation
    protein_to_test <- names(proteins_matrix)[i]
    
    #construct the formula for the model
    fm_to_test <- formula(paste0(outcome, " ~ ", protein_to_test,"+",paste(covariates,collapse="+")))
    
    #fit the model
    model <- lm(fm_to_test, data = jointdata_frame, subset = subset)
    
    #extract out model statistics and store in a array
    model_summary <- summary(model)
    model_stats[i,1:4] <- model_summary$coef[protein_to_test,]
    print(paste("Iteration:", i))
    
  }
  
  model_stats <- data.frame(model_stats)
  
  ##apply a multiple testing correction (Benjamini-Hochberg correction)
  model_stats[,5]   <- p.adjust(model_stats[,4],method="BH")
  
  ## at 0.05 significance level
  model_stats[,6]   <- ifelse(model_stats[,4]<= 0.05, 'yes', 'no')
  
  ## at FDR < 0.05 significance level
  model_stats[,7]  <- ifelse(model_stats[,5]<= 0.05, 'yes', 'no')
  
  #return array
  return(model_stats)
}               ## main fucntion for running linear regression models

## LINEAR MODEL: INCLUDING COHORT AS A RANDOM EFFECT
testContinuousOutcome_lme4 <- function(outcome,covariates,proteins_list, proteins_matrix, jointdata_frame,subset){
  #function for testing log protein expression against outcome, with covariates included, in a linear model
  #protein names to test are stored as the first column of proteins_list
  #protein abundance, outcome and covariates must be stored in jointdata_frame
  #output linear model parameters (estimates, standards, t-stats and p-values)
  #stored in a named data frame of size Nprotein x 4 (i.e. proteins x {statistics}) 
  #initialize the named array for output
  model_stats <- array(NA,dim=c(nrow(proteins_list),9))
  dimnames(model_stats) <- list(proteins_list[,1],c("estimate","std.error","df", "t.value","p.value", "p_adj", "sig", "sig_2","conv"))
  #cycle through protein names
  for (i in 1:ncol(proteins_matrix)){
    # set the protein name and specify log transformation
    protein_to_test <- names(proteins_matrix)[i]

    #construct the formula for the model
    fm_to_test <- formula(paste0(outcome, " ~ ", protein_to_test,"+",paste(covariates,collapse="+"), paste0(" + (1 | cohort_name)")))
    #fit the model
    model <- lmer(fm_to_test, data = jointdata_frame, subset = subset, REML=TRUE,
                  control=lmerControl(optimizer="nloptwrap",optCtrl=list(algorithm = "NLOPT_LN_COBYLA", xtol_rel=1e-6,xtol_abs=1e-10)))
    #extract out model statistics and store in a array
    model_summary <- summary(model)
    model_stats[i,1:5] <- model_summary$coef[protein_to_test,]
    model_stats[i,9] <- ifelse(is.null(unlist(model_summary$optinfo$conv$lme4)),0,1) ## record if a given protein causes convergence issues
    print(paste("Iteration:", i))
    
    
  }
  model_stats <- data.frame(model_stats)
  ##apply a multiple testing correction (Benjamini-Hochberg correction)
  model_stats[,6]   <- p.adjust(model_stats[,5],method="BH")
  ## at 0.05 significance level
  model_stats[,7]   <- ifelse(model_stats[,5]<= 0.05, 'yes', 'no')
  ## at FDR < 0.05 significance level
  model_stats[,8]  <- ifelse(model_stats[,6]<= 0.05, 'yes', 'no')
  
  #return array
  return(model_stats)
}     ## including cohort as random effect

## 2) logistic regression models (binary outcome)
testbinaryOutcome <- function(outcome,covariates, proteins_matrix, jointdata_frame,subset){
  #function for testing log protein expression against outcome, with covariates included, in a logistic model
  #protein names to test are stored as the first column of proteins_list
  #protein abundance, outcome and covariates must be stored in jointdata_frame
  #output logistic model parameters (estimates, standards, t-stats and p-values)
  #stored in a named data frame of size Nprotein x 4 (i.e. proteins x {statistics}) 
  
  #initialize the named array for output
  model_stats <- array(NA,dim=c(nrow(proteins_list),7))
  dimnames(model_stats) <- list(colnames(proteins_matrix),c("estimate","std.error","z.value","p.value", "p_adj", "sig", "sig_2"))
  
  #cycle through protein names
  for (i in 1:ncol(proteins_matrix)){
    
    # set the protein name
    protein_to_test <- names(proteins_matrix)[i]

    #construct the formula for the model
    fm_to_test <- formula(paste0(outcome, " ~ ", protein_to_test,"+",paste(covariates,collapse="+")))
    
    #fit the model
    model_2 <- glm(fm_to_test, data = jointdata_frame, subset = subset, family="binomial")
    
    #extract out model statistics and store in a array
    model_summary <- summary(model_2)
    model_stats[i,1:4] <- model_summary$coef[protein_to_test,]
    print(paste("Iteration:", i))
    
  }
  
  model_stats <- data.frame(model_stats)
  
  ##apply a multiple testing correction
  model_stats[,5]   <- p.adjust(model_stats[,4],method="BH")
  
  ## at 0.05 significance level
  model_stats[,6]   <- ifelse(model_stats[,4]<= 0.05, 'yes', 'no')
  
  ## at FDR < 0.05 significance level
  model_stats[,7]  <- ifelse(model_stats[,5]<= 0.05, 'yes', 'no')
  
  #return array
  return(model_stats)
  
  
}

## LOGISTIC MODEL: INCLUDING COHORT AS A RANDOM EFFECT
testbinaryOutcome_lme4 <- function(outcome,covariates,proteins_list, proteins_matrix, jointdata_frame,subset){
  #function for testing log protein expression against outcome, with covariates included, in a logistic model
  #protein names to test are stored as the first column of proteins_list
  #protein abundance, outcome and covariates must be stored in jointdata_frame
  #output logistic model parameters (estimates, standards, t-stats and p-values)
  #stored in a named data frame of size Nprotein x 4 (i.e. proteins x {statistics}) 
  
  #initialize the named array for output
  model_stats <- array(NA,dim=c(nrow(proteins_list),8))
  dimnames(model_stats) <- list(proteins_list[,1],c("estimate","std.error","z.value","p.value", "p_adj", "sig", "sig_2", "conv"))
  
  #cycle through protein names
  for (i in 1:ncol(proteins_matrix)){
    
    # set the protein name and specify log transformation
    protein_to_test <- names(proteins_matrix)[i]

    #construct the formula for the model
    fm_to_test <- formula(paste0(outcome, " ~ ", protein_to_test,"+",paste(covariates,collapse="+"), paste0(" + (1 | cohort_name)")))
    
    #fit the model
    model_2 <- glmer(fm_to_test, data = jointdata_frame, subset = subset, family="binomial", 
                     control=glmerControl(optimizer="nloptwrap",optCtrl=list(algorithm = "NLOPT_LN_COBYLA", xtol_rel=1e-6,xtol_abs=1e-10)))
    
    #extract out model statistics and store in a array
    model_summary <- summary(model_2)
    model_stats[i,1:4] <- model_summary$coef[protein_to_test,]
    model_stats[i,8] <- ifelse(is.null(unlist(model_summary$optinfo$conv$lme4)),0,1) ## record if a given protein causes convergence issues
    print(paste("Iteration:", i))
    
  }
  
  model_stats <- data.frame(model_stats)
  
  ##apply a multiple testing correction
  model_stats[,5]   <- p.adjust(model_stats[,4],method="BH")
  
  ## at 0.05 significance level
  model_stats[,6]   <- ifelse(model_stats[,4]<= 0.05, 'yes', 'no')
  
  ## at FDR < 0.05 significance level
  model_stats[,7]  <- ifelse(model_stats[,5]<= 0.05, 'yes', 'no')
  
  #return array
  return(model_stats)
  
  
}

## 3) ordinal regression models (ordered, categorical outcome)
testordinalOutcome <- function(outcome,covariates,proteins_list,proteins_matrix, jointdata_frame,subset,verbose=FALSE){
  #function for testing log protein expression against outcome, with covariates included, in a ordinal model
  #protein names to test are stored as the first column of proteins_list
  #protein abundance, outcome and covariates must be stored in jointdata_frame
  #output logistic model parameters (estimates, standards, t-stats and p-values) 
  #stored in a named data frame of size Nprotein x 4 (i.e. proteins x {statistics}) 
  
  #initialize the named array for output
  model_stats <- array(NA,dim=c(nrow(proteins_list),7))
  dimnames(model_stats) <- list(proteins_list[,1],c("estimate","std.error","t.value","p.value", "p_adj", "sig", "sig_2"))
  #cycle through protein names
  
  #fit the null model
  fm_to_test_null <- formula(paste0(outcome, " ~ ",paste(covariates,collapse="+")))
  null_model <- MASS::polr(fm_to_test_null, data = jointdata_frame, subset = subset, Hess = TRUE)
  
  for (i in 1:ncol(proteins_matrix)){
    
    
    # set the protein name and specify log transformation
    protein_to_test <- names(proteins_matrix)[i]

    if (verbose) cat(paste0("Testing protein: ",i,", ",protein_to_test,"..."))
    
    #construct the formula for the model
    fm_to_test <- formula(paste0(outcome, " ~ ", protein_to_test,"+",paste(covariates,collapse="+")))
    
    #fit the model
    model_3 <- MASS::polr(fm_to_test, data = jointdata_frame, subset = subset, Hess = TRUE)
    
    #extract out model statistics and store in a array
    model_summary <- summary(model_3)
    model_stats[i,1:3] <- model_summary$coef[protein_log_to_test,]
    model_stats[i,4] <- anova(model_3,null_model,test="Chisq")$Pr[2]
    print(paste("Iteration:", i))
    
    if (verbose) cat(paste0("done\n"))
    
  }
  
  model_stats <- data.frame(model_stats)
  
  ##apply a multiple testing correction
  model_stats[,5]   <- p.adjust(model_stats[,4],method="BH")
  
  ## at 0.05 significance level
  model_stats[,6]   <- ifelse(model_stats[,4]<= 0.05, 'yes', 'no')
  
  ## at FDR < 0.05 significance level
  model_stats[,7]  <- ifelse(model_stats[,5]<= 0.05, 'yes', 'no')
  
  #return array
  return(model_stats)
}

## ORDINAL LOGISTIC MODEL: INCLUDING COHORT AS A RANDOM EFFECT
testordinalOutcome_lme4 <- function(outcome,covariates,proteins_list,proteins_matrix, jointdata_frame,subset,verbose=FALSE){
  #function for testing log protein expression against outcome, with covariates included, in a ordinal model
  #protein names to test are stored as the first column of proteins_list
  #protein abundance, outcome and covariates must be stored in jointdata_frame
  #output logistic model parameters (estimates, standards, t-stats and p-values) 
  #stored in a named data frame of size Nprotein x 4 (i.e. proteins x {statistics}) 
  
  #initialize the named array for output
  model_stats <- array(NA,dim=c(nrow(proteins_list),7))
  dimnames(model_stats) <- list(proteins_list[,1],c("estimate","std.error","t.value","p.value", "p_adj", "sig", "sig_2"))
  #cycle through protein names
  
  #fit the null model
  fm_to_test_null <- formula(paste0(outcome, " ~ ",paste(covariates,collapse="+"), paste0(" + (1 | cohort_name)")))
  null_model <- clmm(fm_to_test_null, data = jointdata_frame, subset = subset, Hess = TRUE)
  
  for (i in 1:ncol(proteins_matrix)){
    
    
    # set the protein name and specify log transformation
    protein_to_test <- names(proteins_matrix)[i]

    if (verbose) cat(paste0("Testing protein: ",i,", ",protein_to_test,"..."))
    
    #construct the formula for the model
    fm_to_test <- formula(paste0(outcome, " ~ ", protein_to_test,"+",paste(covariates,collapse="+"), paste0(" + (1 | cohort_name)")))
    
    #fit the model
    model_3 <- clmm(fm_to_test, data = jointdata_frame, subset = subset, Hess = TRUE)
    
    #extract out model statistics and store in a array
    model_summary <- summary(model_3)
    model_stats[i,1:4] <- model_summary$coef[protein_to_test,] ## clmm generates a p-value 
    ##model_stats[i,4] <- anova(model_3,null_model,test="Chisq")$Pr[2,4]
    print(paste("Iteration:", i))
    
    if (verbose) cat(paste0("done\n"))
    
  }
  
  model_stats <- data.frame(model_stats)
  
  ##apply a multiple testing correction
  model_stats[,5]   <- p.adjust(model_stats[,4],method="BH")
  
  ## at 0.05 significance level
  model_stats[,6]   <- ifelse(model_stats[,4]<= 0.05, 'yes', 'no')
  
  ## at FDR < 0.05 significance level
  model_stats[,7]  <- ifelse(model_stats[,5]<= 0.05, 'yes', 'no')
  
  #return array
  return(model_stats)
}


## 4) Function for generating Volcano Plots 
plot_volcano_6             <- function(input, title, subtitle){
  
  #Subset the tibble to proteins that reach statistical significance (p_adj <0.05) and are the top 20 most significant (and |logFC| > 0.3) to label with Target name in the plot
  label <- input %>%
    dplyr::filter(input$sig_2 == "yes") %>% ## select only significant proteins at p_adj (<0.05)
    dplyr::slice_min(p_adj, n = 20)  %>%    ## select top 20 most significant proteins
    dplyr::slice(1:20)                             ## if multiple proteins share the same p-value, the first 20 in list will be labelled
  
  #Count the number of significantly regulated up and down proteins 
  number_upregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate > 0) %>%
    nrow()  
  
  number_downregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate < 0) %>%
    nrow()   
  
  annotations <- data.frame(
    xpos = c(-Inf,Inf),
    ypos =  c( Inf,Inf),
    annotateText = c(paste0(number_downregulated, " downregulated"), paste0(number_upregulated, " upregulated")),
    hjustvar = c(-0.25,1) ,
    vjustvar = c(0.9,0.9)) #<- adjust
  
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
  input$diffexpressed <- "No"
  # if adjusted pvalue < 0.05 and estimate >0, set as "UP"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate > 0] <- "Upregulated"
  
  # if adjusted pvalue < 0.05 and estimate <0, set as "DOWN"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate < 0] <- "Downregulated"
  
  ggplot(input, aes(x = estimate, y = -log10(p_adj))) + 
    geom_point(aes(col = diffexpressed), size = 1.5) + ## set datapoint colour based on grouping 
    geom_label_repel(data = label,
                     aes(x = estimate, y = -log10(p_adj),label = Target), size = 3.5) + ##label to 20 most significant data points with SOMAmer ID 
    geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4, colour = "grey43") +  ## add text labels to top of graph (counts of upregulated/downregulated proteins)
    xlab(expression("Log Odds Ratio (per SD)")) +              ## x-axis label                            
    ylab(expression("-log"[10]*"(Adjusted P-values)")) +       ## y-axis label
    labs(title = title, caption = subtitle) +                  ## title and caption headers
    theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(0,0,20,0)),  ## size, position and font style of title
          plot.caption = element_text(size = 10,hjust = 0, face = "italic")) +                          ## size, position and font style of caption
    geom_vline(xintercept = c(-1, 1), col = "grey", linetype = 'dashed') +                              ## vertical lines corresponding to effects greater / less than 1
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +                           ## horizontal line corresponding to the line of significance (log10(p_adj))
    scale_color_manual("Significant (padj ≤ 0.05):", values = c("#bb0c00", "#00AFBB", "grey"),
                       limits=c("Upregulated", "Downregulated","No"),
                       labels=c("Upregulated", "Downregulated", "Not Regulated"))+
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) ## legend_title and size
  
} ## top 20 most significant proteins
plot_volcano_6_top10updown <- function(input, title, subtitle){
  
  #Subset the tibble to proteins that reach statistical significance (p_adj <0.05) and are the top 20 most significant to label with Target name in the plot
  ## on ggplot, label only the top 20 most significant proteins (based on adjusted p-value)
  ## top 10 upregulated and top 10 downregulated
  
  label <- input 
  label <- as.data.frame(label)
  
  label <- label[label$sig_2 == "yes",]                                                   ## pull out only significantly regulated proteins based on adjusted p-value
  label_order <-  label[order(label$p_adj, decreasing = FALSE),]                          ## order data based on p_adj 
  label_order <-  label_order[!duplicated(label_order$TargetFullName),]                   ## remove duplicates, selecting only the first occurrence of duplicated SOMAmers that appear in list
  label_up    <-  head(label_order[label_order$estimate >0,], 15)                         ## select the top 15 most strongly upregulated proteins
  label_down  <-  head(label_order[label_order$estimate <0,], 15)                         ## select the top 15 most strongly downregulated proteins
  label <- rbind(label_up, label_down)                                                    ## combine lists of 15x upregulated & 15x downregulated proteins
  
  ## if multiple proteins share the same p-value, the first 30 in list will be labelled
  
  #Count the number of significantly regulated up and down proteins 
  number_upregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate > 0) %>%
    nrow()  
  
  number_downregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate < 0) %>%
    nrow()   
  
  annotations <- data.frame(
    xpos = c(-Inf,Inf),
    ypos =  c(Inf,Inf),
    annotateText = c(paste0(number_downregulated, " downregulated"), paste0(number_upregulated, " upregulated")),
    hjustvar = c(-0.25,1) ,
    vjustvar = c(0.9,0.9)) #<- adjust
  
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated
  input$diffexpressed <- "No"
  # if adjusted pvalue < 0.05 and estimate >0, set as "UP"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate > 0] <- "Upregulated"
  
  # if adjusted pvalue < 0.05 and estimate <0, set as "DOWN"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate < 0] <- "Downregulated"
  
  ggplot(input, aes(x = estimate, y = -log10(p_adj))) + 
    geom_point(aes(col = diffexpressed), size = 1.5) +                                                 ## set datapoint colour based on grouping 
    geom_label_repel(data = label,
                     aes(x = estimate, y = -log10(p_adj),label = Target), size = 3.5) +                ##label to 20 most significant data points with SOMAmer ID 
    geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4, colour = "grey43") +  ## add text labels to top of graph (counts of upregulated/downregulated proteins)
    xlab(expression("Log Odds Ratio (per SD)")) +                                                       ## x-axis label                            
    ylab(expression("-log"[10]*"(Adjusted P-values)")) +                                                ## y-axis label
    labs(title = title, caption = subtitle) +                                                           ## title and caption headers
    theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(0,0,20,0)),  ## size, position and font style of title
          plot.caption = element_text(size = 10,hjust = 0, face = "italic")) +                          ## size, position and font style of caption
    geom_vline(xintercept = c(-1, 1), col = "grey", linetype = 'dashed') +                              ## vertical lines corresponding to effects greater / less than 1
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +                           ## horizontal line corresponding to the line of significance (log10(p_adj))
    scale_color_manual("Significant (padj ≤ 0.05):", values = c("#bb0c00", "#00AFBB", "grey"),
                       limits=c("Upregulated", "Downregulated","No"),
                       labels=c("Upregulated", "Downregulated", "Not Regulated"))+
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  
  ## legend_title and size
  
} ## top 30 most significant proteins, top 10x upregulated & top 10x downregulated (not labeling repeats)
plot_volcano_7  <- function(input, title, subtitle){
  
  label <- input 
  label <- as.data.frame(label)
  
  label <- label[label$sig_2 == "no",]                                                   ## pull out only signifcaintly regulated proteins based on adjusted p-value
  label_order <-  label[order(label$p_adj, decreasing = FALSE),]                          ## order data based on p_adj 
  label_order <-  label_order[!duplicated(label_order$TargetFullName),]                   ## remove duplicates, selecting only the first occurrence of SOMAmer
  label_up    <-  head(label_order[label_order$estimate >0,], 10)                         ## select the top 10 most strongly upregulated proteins
  label_down  <-  head(label_order[label_order$estimate <0,], 10)                         ## select the top 10 most strongly downregulated proteins
  label <- rbind(label_up, label_down)                                                    ## combine lists of 10x upregulated & 10x downregulated proteins
  
  ## if multiple proteins share the same p-value, the first 20 in list will be labelled
  
  #Count the number of significantly regulated up and down proteins 
  number_upregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate > 0) %>%
    nrow()  
  
  number_downregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate < 0) %>%
    nrow()   
  
  annotations <- data.frame(
    xpos = c(-Inf,Inf),
    ypos =  c( Inf,Inf),
    annotateText = c(paste0(number_downregulated, " downregulated"), paste0(number_upregulated, " upregulated")),
    hjustvar = c(-0.25,1) ,
    vjustvar = c(1,1)) #<- adjust
  
  ##ggplot(input, aes(x = estimate, y = -log10(p_adj))) + geom_point(aes(col = p_adj < 0.05)) + geom_text_repel(data = label,aes(label = Target)) + theme_pubr(legend = "right") + scale_color_manual(values = c("grey", "red")) +
  ##  geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4, colour = "grey43") + 
  ##  xlab(expression("Log Odds Ratio (per SD)")) +
  ##  ylab(expression("-log"[10]*"(Adjusted P-values)")) +
  ##  labs(title = title, caption = subtitle) + 
  ##  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold", margin = margin(0,0,20,0)), 
  ##        plot.caption = element_text(size = 10,hjust = 0.1, face = "bold")) +
  ##  geom_hline(yintercept=-log10(0.05), col="red", linetype = "longdash") 
  
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
  input$diffexpressed <- "No"
  # if adjusted pvalue < 0.05 and estimate >0, set as "UP"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate > 0] <- "Upregulated"
  
  # if adjusted pvalue < 0.05 and estimate <0, set as "DOWN"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate < 0] <- "Downregulated"
  
  ggplot(input, aes(x = estimate, y = -log10(p_adj))) + 
    geom_point(aes(col = diffexpressed), size = 1.5) + ## set datapoint colour based on grouping 
    geom_label_repel(data = label %>% dplyr::slice_min(p_adj, n = 20),aes(x = estimate, y = -log10(p_adj),label = Target), size = 3.5) + ##label to 20 most significant data points with SOMAmer ID 
    geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4, colour = "grey43") +  ## add text labels to top of graph (counts of upregulated/downregulated proteins)
    xlab(expression("Log Odds Ratio (per SD)")) +              ## x-axis label                            
    ylab(expression("-log"[10]*"(Adjusted P-values)")) +       ## y-axis label
    labs(title = title, caption = subtitle) +                  ## title and caption headers
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold", margin = margin(0,0,20,0)),  ## size, position and font style of title
          plot.caption = element_text(size = 12,hjust = 0, face = "italic")) +                          ## size, position and font style of caption
    geom_vline(xintercept = c(-1, 1), col = "grey", linetype = 'dashed') +                              ## vertical lines corresponding to effects greater / less than 1
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +                           ## horizontal line corresponding to the line of significance (log10(p_adj))
    scale_color_manual("Significant (p_adj <0.05):", values = c("#bb0c00", "#00AFBB", "grey"),
                       limits=c("Upregulated", "Downregulated","No"),
                       labels=c("Upregulated", "Downregulated", "Not Regulated"))+ 
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
}                                   ## top 10 proteins if no significant proteins identified
plot_volcano_10 <- function(input, title, subtitle, discovery, replication, combined){
  
  ## merge all three working dataframes together (Discovery, Replication and Combined)
  df <- full_join(replication,discovery, by = "Row.names", suffix = c("","_Disc")) %>% 
    full_join(combined, by = "Row.names", suffix = c("_Rep", "_Comb"))
  
  label <- df 
  label <- as.data.frame(label)
  
  ## create a binary indicator variable that indicates whether a protein is replicated (based on statistical significance and direction of effect across Discovery & Replication)
  label$replicated <- ""
  label$replicated[label$sig_2_Disc == "yes" & label$sig_2_Rep == "yes" & label$estimate_Disc >0 & label$estimate_Rep >0 | 
                     label$sig_2_Disc == "yes" & label$sig_2_Rep == "yes" & label$estimate_Disc <0 & label$estimate_Rep <0] <- "Replicated"
  ## assign colour to group variable
  label$replicated <- ifelse(label$replicated == "Replicated", "orange", "White") 
  
  ## on ggplot, label only the top 20 most significant proteins (based on adjusted p-value)
  
  label <- label[label$sig_2_Comb == "yes",]                                       ## select out only significant SOMAmers by adjusted p-value
  label_order <-  label[order(label$p_adj_Comb, decreasing = FALSE),]              ## order data based on p_adj 
  label_order <-  label_order[!duplicated(label_order$TargetFullName_Comb),]       ## remove duplicates, selecting the first occurrence based on ordered list
  label_up    <-  head(label_order[label_order$estimate_Comb >0,], 15)             ## select the top 15 most strongly regulated proteins
  label_down  <-  head(label_order[label_order$estimate_Comb <0,], 15)             ## select the top 15 most strongly down-regulated proteins
  label       <- rbind(label_up, label_down)                                       ## combine lists of top 30 most strongly upregulated and downregulated proteins
  
  #Count the number of significantly regulated up and down proteins 
  number_upregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate > 0) %>%
    nrow()  
  
  number_downregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate < 0) %>%
    nrow()   
  
  annotations <- data.frame(
    xpos = c(-Inf,Inf),
    ypos =  c( Inf,Inf),
    annotateText = c(paste0(number_downregulated, " downregulated"), paste0(number_upregulated, " upregulated")),
    hjustvar = c(-0.25,1) ,
    vjustvar = c(0.9,0.9)) #<- adjust
  
  # Add a column to the input data frame to colour datapoints as to whether theyare UP- or DOWN- regulated 
  input$diffexpressed <- "No"
  # if adjusted pvalue < 0.05 and estimate_Comb >0, set as "UP"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate > 0] <- "Upregulated"
  
  # if adjusted pvalue < 0.05 and estimate_Comb <0, set as "DOWN"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate < 0] <- "Downregulated"
  
  ggplot(input, aes(x = estimate, y = -log10(p_adj))) + 
    geom_point(aes(col = diffexpressed), size = 1.5) +                                                  ## set datapoint color based on grouping 
    geom_label_repel(data = label,                                                                      ## now using the label subset, label proteins regulated in combined dataset and colour text boxes if replicated across Discovery & Replication
                     aes(x = estimate_Comb, y = -log10(p_adj_Comb),label = Target_Comb), size = 3.5, fill = label$replicated) + 
    geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4, colour = "grey43") +  ## add text labels to top of graph (counts of upregulated/downregulated proteins)
    xlab(expression("Log Odds Ratio (per SD)")) +                                                       ## x-axis label                            
    ylab(expression("-log"[10]*"(Adjusted P-values)")) +                                                ## y-axis label
    labs(title = title, caption = subtitle) +                                                           ## title and caption headers
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold", margin = margin(0,0,20,0)),  ## size, position and font style of title
          plot.caption = element_text(size = 12,hjust = 0, face = "italic")) +                          ## size, position and font style of caption
    geom_vline(xintercept = c(-1, 1), col = "grey", linetype = 'dashed') +                              ## vertical lines corresponding to effects greater / less than 1
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +                           ## horizontal line corresponding to the line of significance (log10(p_adj_Comb))
    scale_color_manual("Significant (padj ≤ 0.05):", values = c("#bb0c00", "#00AFBB", "grey"),
                       limits=c("Upregulated", "Downregulated","No"),
                       labels=c("Upregulated", "Downregulated", "Not Regulated"))+
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14), 
          axis.text.y = element_text(size = 12), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) ## legend_title and size
  
} ## top 20 most significant proteins in combined dataset & replicated across Discovery & Replication, top 10 upregulated & top 10 downregulated
plot_volcano_11 <- function(input, title, subtitle, model1, model2){
  
  ## merge two working dataframes together 
  df <- full_join(model1,model2, by = "Row.names", suffix = c("_model1","_model2"))
  
  label <- df 
  label <- as.data.frame(label)
  
  ## create a binary indicator variable that indicates whether a protein is replicated (based on statistical significance and direction of effect across Discovery & Replication)
  label$replicated <- ""
  label$replicated[label$sig_2_model1 == "yes" & label$sig_2_model2 == "yes" & label$estimate_model1 >0 & label$estimate_model2 >0 | 
                     label$sig_2_model1 == "yes" & label$sig_2_model2 == "yes" & label$estimate_model1 <0 & label$estimate_model2 <0] <- "Replicated"
  ## assign colour to group variable
  label$replicated <- ifelse(label$replicated == "Replicated",  "#D55E00", "black") 
  
  ## on ggplot, label only the top 20 most significant proteins (based on adjusted p-value)
  
  label       <- label[label$sig_2_model1 == "yes",]                                 ## select out only significant SOMAmers by adjusted p-value
  label_order <-  label[order(label$p_adj_model1, decreasing = FALSE),]              ## order data based on p_adj 
  label_order <-  label_order[!duplicated(label_order$TargetFullName_model1),]       ## remove duplicates, selecting the first occurrence based on ordered list
  label_up    <-  head(label_order[label_order$estimate_model1 >0,], 10)             ## select the top 10 most strongly regulated proteins
  label_down  <-  head(label_order[label_order$estimate_model1 <0,], 10)             ## select the top 10 most strongly down-regulated proteins
  label       <-  rbind(label_up, label_down)                                        ## combine lists of top 20 most strongly upregulated and downregulated proteins
  
  #Count the number of significantly regulated up and down proteins 
  number_upregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate > 0) %>%
    nrow()  
  
  number_downregulated <- input %>%
    dplyr::filter(p_adj < 0.05 & estimate < 0) %>%
    nrow()   
  
  annotations <- data.frame(
    xpos = c(-Inf,Inf),
    ypos =  c( Inf,Inf),
    annotateText = c(paste0(number_downregulated, " downregulated"), paste0(number_upregulated, " upregulated")),
    hjustvar = c(-0.25,1) ,
    vjustvar = c(0.9,0.9)) #<- adjust
  
  # Add a column to the input data frame to colour datapoints as to whether theyare UP- or DOWN- regulated 
  input$diffexpressed <- "No"
  # if adjusted pvalue < 0.05 and estimate_model1 >0, set as "UP"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate > 0] <- "Upregulated"
  
  # if adjusted pvalue < 0.05 and estimate_model1 <0, set as "DOWN"
  input$diffexpressed[input$p_adj < 0.05 & input$estimate < 0] <- "Downregulated"
  
  ggplot(input, aes(x = estimate, y = -log10(p_adj))) + 
    geom_point(aes(col = diffexpressed), size = 1.5) +                                                  ## set datapoint color based on grouping 
    geom_label_repel(data = label,                                                                      ## now using the label subset, label proteins regulated in combined dataset and colour text boxes if replicated across Discovery & Replication
                     aes(x = estimate_model1, y = -log10(p_adj_model1),label = Target_model1), size = 3.5, color = label$replicated,
                     segment.colour = "black",
                     ##fontface = 'bold',
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.5, "lines")) + 
    geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4, colour = "grey43") +  ## add text labels to top of graph (counts of upregulated/downregulated proteins)
    xlab(expression("Log Odds Ratio (per SD)")) +                                                       ## x-axis label                            
    ylab(expression("-log"[10]*"(Adjusted P-values)")) +                                                ## y-axis label
    labs(title = title, caption = subtitle) +                                                           ## title and caption headers
    theme(plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(0,0,20,0)),  ## size, position and font style of title
          plot.caption = element_text(size = 10,hjust = 0, face = "italic")) +                          ## size, position and font style of caption
    geom_vline(xintercept = c(-1, 1), col = "grey", linetype = 'dashed') +                              ## vertical lines corresponding to effects greater / less than 1
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +                           ## horizontal line corresponding to the line of significance (log10(p_adj_model1))
    scale_color_manual("Significant (padj ≤ 0.05):", values = c("#bb0c00", "#00AFBB", "grey"),
                       limits=c("Upregulated", "Downregulated","No"),
                       labels=c("Upregulated", "Downregulated", "Not Regulated"))+
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16), 
          axis.text.y = element_text(size = 14), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) ## legend_title and size
  
}                   ## top 20 most significant proteins comparing two separate models, top 10 upregulated & top 10 downregulated
plot_volcano_final <- function(combined_noCohort, title, subtitle, discovery, replication, combined_cohort){
  
  ## merge all three working dataframes together (Discovery, Replication and Combined)
  df <- full_join(combined_noCohort, discovery, by = "Row.names", suffix = c("","_Disc")) %>%
    full_join(replication, by = "Row.names", suffix = c("_combNoCohort", "_Rep")) %>%
    full_join(combined_cohort, by = "Row.names", suffix = c("_Disc", "_combCohort")) 
  
  label <- df 
  label <- as.data.frame(label)
  
  ## create a binary indicator variable that indicates whether a protein is replicated (based on statistical significance and direction of effect across Discovery & Replication)
  label$replicated <- ""
  label$replicated[label$sig_2_Disc == "yes" & label$sig_2_Rep == "yes" & label$estimate_Disc >0 & label$estimate_Rep >0 | 
                     label$sig_2_Disc == "yes" & label$sig_2_Rep == "yes" & label$estimate_Disc <0 & label$estimate_Rep <0] <- "Replicated"
  ## assign colour to group variable
  label$replicated <- ifelse(label$replicated == "Replicated", "orange", "grey") 
  
  ##View(select(label, estimate_Disc, p_adj_Disc, sig_2_Disc,estimate_Rep, p_adj_Rep, sig_2_Rep, replicated))
  ##View(select(filter(label, sig_2_Disc == "yes" & sig_2_Rep == "yes" & estimate_Disc >0 & estimate_Rep >0 | 
  ##                          sig_2_Disc == "yes" & sig_2_Rep == "yes" & estimate_Disc <0 & estimate_Rep <0), estimate_Disc, p_adj_Disc, sig_2_Disc,estimate_Rep, p_adj_Rep, sig_2_Rep, replicated))
  
  ## create a binary indicator variable that indicates whether a protein is 'robust' ie. significant in no cohort adjusted, combined dataset and significant in cohort adjusted, combined dataset 
  label$cohort_vs_noCohort <- ""
  label$cohort_vs_noCohort[label$sig_2_combNoCohort == "yes" & label$sig_2 == "yes" & label$estimate_combNoCohort >0 & label$estimate >0 | 
                             label$sig_2_combNoCohort == "yes" & label$sig_2 == "yes" & label$estimate_combNoCohort <0 & label$estimate <0] <- "Robust"
  ## assign colour to group variable
  label$cohort_vs_noCohort <- ifelse(label$cohort_vs_noCohort == "Robust", "blue", "darkgreen") 
  
  ##View(select(label, estimate, p_adj, sig_2,estimate_combNoCohort, p_adj_combNoCohort, sig_2_combNoCohort, cohort_vs_noCohort))
  ##View(select(filter(label, label$sig_2_combNoCohort == "yes" & label$sig_2 == "yes" & label$estimate_combNoCohort >0 & label$estimate >0 | 
  ##                          label$sig_2_combNoCohort == "yes" & label$sig_2 == "yes" & label$estimate_combNoCohort <0 & label$estimate <0), estimate, p_adj, sig_2,estimate_combNoCohort, p_adj_combNoCohort, sig_2_combNoCohort, cohort_vs_noCohort))
  
  ## on ggplot, label only the top 30 most significant proteins (based on adjusted p-value) in non-cohort adjusted analysis
  
  label <- label[label$sig_2_combNoCohort == "yes",]                                       ## select out only significant SOMAmers by adjusted p-value
  label_order <-  label[order(label$p_adj_combNoCohort, decreasing = FALSE),]              ## order data based on p_adj 
  label_order <-  label_order[!duplicated(label_order$TargetFullName_combNoCohort),]       ## remove duplicates, selecting the first occurrence based on ordered list
  label_up    <-  head(label_order[label_order$estimate_combNoCohort >0,], 15)             ## select the top 15 most strongly regulated proteins
  label_down  <-  head(label_order[label_order$estimate_combNoCohort <0,], 15)             ## select the top 15 most strongly down-regulated proteins
  label       <- rbind(label_up, label_down)                                       ## combine lists of top 10 most strongly upregulated and downregulated proteins
  
  #Count the number of significantly regulated up and down proteins 
  number_upregulated <- combined_noCohort %>%
    dplyr::filter(p_adj < 0.05 & estimate > 0) %>%
    nrow()  
  
  number_downregulated <- combined_noCohort %>%
    dplyr::filter(p_adj < 0.05 & estimate < 0) %>%
    nrow()   
  
  annotations <- data.frame(
    xpos = c(-Inf,Inf),
    ypos =  c( Inf,Inf),
    annotateText = c(paste0(number_downregulated, " downregulated"), paste0(number_upregulated, " upregulated")),
    hjustvar = c(-0.25,1) ,
    vjustvar = c(0.9,0.9)) #<- adjust
  
  # Add a column to the input data frame to colour datapoints as to whether theyare UP- or DOWN- regulated 
  combined_noCohort$diffexpressed <- "No"
  # if adjusted pvalue < 0.05 and estimate >0, set as "UP"
  combined_noCohort$diffexpressed[combined_noCohort$p_adj < 0.05 & combined_noCohort$estimate > 0] <- "Upregulated"
  
  # if adjusted pvalue < 0.05 and estimate <0, set as "DOWN"
  combined_noCohort$diffexpressed[combined_noCohort$p_adj < 0.05 & combined_noCohort$estimate < 0] <- "Downregulated"
  
  ggplot(combined_noCohort, aes(x = estimate, y = -log10(p_adj))) + 
    geom_point(aes(col = diffexpressed), size = 1.5) +                                                  ## set datapoint color based on grouping 
    geom_label_repel(data = label,                                                                      ## now using the label subset, label proteins regulated in combined dataset and colour text boxes if replicated across Discovery & Replication
                     aes(x = estimate_combNoCohort, y = -log10(p_adj_combNoCohort),label = Target_combNoCohort), size = 3.5, fill = label$replicated, color = label$cohort_vs_noCohort) + 
    geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4, colour = "grey43") +  ## add text labels to top of graph (counts of upregulated/downregulated proteins)
    xlab(expression("Log Odds Ratio (per SD)")) +                                                       ## x-axis label                            
    ylab(expression("-log"[10]*"(Adjusted P-values)")) +                                                ## y-axis label
    labs(title = title, caption = subtitle) +                                                           ## title and caption headers
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold", margin = margin(0,0,20,0)),  ## size, position and font style of title
          plot.caption = element_text(size = 12,hjust = 0, face = "italic")) +                          ## size, position and font style of caption
    geom_vline(xintercept = c(-1, 1), col = "grey", linetype = 'dashed') +                              ## vertical lines corresponding to effects greater / less than 1
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +                           ## horizontal line corresponding to the line of significance (log10(p_adj_Comb))
    scale_color_manual("Significant (padj ≤ 0.05):", values = c("#bb0c00", "#00AFBB", "grey"),
                       limits=c("Upregulated", "Downregulated","No"),
                       labels=c("Upregulated", "Downregulated", "Not Regulated"))+
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14), 
          axis.text.y = element_text(size = 12), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) ## legend_title and size
  
} ## top 30 most significant proteins in combined dataset & replicated across Discovery & Replication, top 10 upregulated & top 10 downregulated
## and robust (significant in cohort adjusted model also) with top 10 upregulated & top 10 downregulated
plot_volcano_RR    <- function(combined_noCohort, title, subtitle, discovery, replication, combined_cohort){
  
  ## merge all four working dataframes together (Discovery, Replication and Combined (cohort adjusted and not cohort adjusted))
  df <- full_join(combined_noCohort, discovery, by = "Row.names", suffix = c("","_Disc")) %>%
    full_join(replication, by = "Row.names", suffix = c("_combNoCohort", "_Rep")) %>%
    full_join(combined_cohort, by = "Row.names", suffix = c("_Disc", "_combCohort")) 
  
  label <- df 
  label <- as.data.frame(label)
  
  ## create an indicator variable that indicates whether a protein is replicated & robust, doesn't replicate and is not robust, and then either replicates or is robust
  ## 'replicated' = significant in both discovery and replication based on adjusted p-value with estimates in the same direction
  ## 'robust'     = significant in both combined analysis (not additionally adjusted for cohort) and combined analysis (additionally adjusted for cohort) based on adjusted p-values with estimates in the same direction
  label$replicated <- ""
  label$replicated[(label$sig_2_Disc == "yes" & label$sig_2_Rep == "yes" & label$estimate_Disc >0 & label$estimate_Rep >0 | 
                      label$sig_2_Disc == "yes" & label$sig_2_Rep == "yes" & label$estimate_Disc <0 & label$estimate_Rep <0) &
                     (label$sig_2 == "yes" & label$sig_2_combNoCohort == "yes" & label$estimate >0 & label$estimate_combNoCohort >0 | 
                        label$sig_2 == "yes" & label$sig_2_combNoCohort == "yes" & label$estimate <0 & label$estimate_combNoCohort <0) & label$replicated == ""]<- "Replicated"
  
  label$replicated[(label$sig_2_Disc == "no" | label$sig_2_Rep == "no") & (label$sig_2 == "no" | label$sig_2_combNoCohort == "no") & label$replicated == ""] <- "NotReplicated&Robust"
  
  label$replicated[((label$sig_2_Disc == "yes" & label$sig_2_Rep == "yes" & label$estimate_Disc >0 & label$estimate_Rep >0 | 
                       label$sig_2_Disc == "yes" & label$sig_2_Rep == "yes" & label$estimate_Disc <0 & label$estimate_Rep <0) |
                      (label$sig_2 == "yes" & label$sig_2_combNoCohort == "yes" & label$estimate >0 & label$estimate_combNoCohort >0 | 
                         label$sig_2 == "yes" & label$sig_2_combNoCohort == "yes" & label$estimate <0 & label$estimate_combNoCohort <0)) & label$replicated== ""]<- "Either"
  
  ## assign colour to group variable
  label$replicated[label$replicated == "Replicated"] <- "orange" 
  label$replicated[label$replicated == "NotReplicated&Robust"] <- "White"
  label$replicated[label$replicated == "Either"] <- "lightgreen" 
  
  label <- label[label$sig_2_combNoCohort == "yes",]                                       ## select out only significant SOMAmers by adjusted p-value
  label_order <-  label[order(label$p_adj_combNoCohort, decreasing = FALSE),]              ## order data based on p_adj 
  label_order <-  label_order[!duplicated(label_order$TargetFullName_combNoCohort),]       ## remove duplicates, selecting the first occurrence based on ordered list
  label_up    <-  head(label_order[label_order$estimate_combNoCohort >0,], 15)             ## select the top 15 most strongly regulated proteins
  label_down  <-  head(label_order[label_order$estimate_combNoCohort <0,], 15)             ## select the top 15 most strongly down-regulated proteins
  label       <- rbind(label_up, label_down)                                       ## combine lists of top 10 most strongly upregulated and downregulated proteins
  
  #Count the number of significantly regulated up and down proteins 
  number_upregulated <- combined_noCohort %>%
    dplyr::filter(p_adj < 0.05 & estimate > 0) %>%
    nrow()  
  
  number_downregulated <- combined_noCohort %>%
    dplyr::filter(p_adj < 0.05 & estimate < 0) %>%
    nrow()   
  
  annotations <- data.frame(
    xpos = c(-Inf,Inf),
    ypos =  c( Inf,Inf),
    annotateText = c(paste0(number_downregulated, " downregulated"), paste0(number_upregulated, " upregulated")),
    hjustvar = c(-0.25,1) ,
    vjustvar = c(0.9,0.9)) #<- adjust
  
  # Add a column to the input data frame to colour datapoints as to whether theyare UP- or DOWN- regulated 
  combined_noCohort$diffexpressed <- "No"
  # if adjusted pvalue < 0.05 and estimate >0, set as "UP"
  combined_noCohort$diffexpressed[combined_noCohort$p_adj < 0.05 & combined_noCohort$estimate > 0] <- "Upregulated"
  
  # if adjusted pvalue < 0.05 and estimate <0, set as "DOWN"
  combined_noCohort$diffexpressed[combined_noCohort$p_adj < 0.05 & combined_noCohort$estimate < 0] <- "Downregulated"
  
  ggplot(combined_noCohort, aes(x = estimate, y = -log10(p_adj))) + 
    geom_point(aes(col = diffexpressed), size = 1.5) +                                                  ## set datapoint color based on grouping 
    geom_label_repel(data = label,                                                                      ## now using the label subset, label proteins regulated in combined dataset and colour text boxes if replicated across Discovery & Replication
                     aes(x = estimate_combNoCohort, y = -log10(p_adj_combNoCohort),label = Target_combNoCohort), size = 3.5, fill = label$replicated, color = "black") + 
    geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 4, colour = "grey43", check_overlap = TRUE) +  ## add text labels to top of graph (counts of upregulated/downregulated proteins)
    xlab(expression("Log Odds Ratio (per SD)")) +                                                       ## x-axis label                            
    ylab(expression("-log"[10]*"(Adjusted P-values)")) +                                                ## y-axis label
    labs(title = title, caption = subtitle) +                                                           ## title and caption headers
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold", margin = margin(0,0,20,0)),  ## size, position and font style of title
          plot.caption = element_text(size = 10,hjust = 0, face = "italic")) +                          ## size, position and font style of caption
    geom_vline(xintercept = c(-1, 1), col = "grey", linetype = 'dashed') +                              ## vertical lines corresponding to effects greater / less than 1
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +                           ## horizontal line corresponding to the line of significance (log10(p_adj_Comb))
    scale_color_manual("Significant (padj ≤ 0.05):", values = c("#bb0c00", "#00AFBB", "grey"),
                       limits=c("Upregulated", "Downregulated","No"),
                       labels=c("Upregulated", "Downregulated", "Not Regulated"))+
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14), 
          axis.text.y = element_text(size = 12), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) ## legend_title and size
  
} ## top 30 most significant proteins in combined dataset & replicated across Discovery & Replication, top 10 upregulated & top 10 downregulated


## 5) Function for creating a matrix of statistically significant proteins with coefficients, p-values and gene map names
significance_matrix <- function(model_name, named.df){
  df <- merge(model_name, protein_gene_map, by = 0, all = FALSE)    ## merge regression output matrix with gene map (holds gene ID and names)
  row.names(df) <-df$Row.names
  table_1 <- table(model_name$sig_2[model_name$sig_2 == "yes"])     ## tabulate number of statistically significant proteins (adjusted p-value <=0.05)
  table_2 <- table(model_name$sig_2[model_name$estimate >0 & model_name$sig_2 == "yes"]) ## significant and up-regulated proteins
  table_3 <- table(model_name$sig_2[model_name$estimate <0 & model_name$sig_2 == "yes"]) ## significant and down-regulated proteins
  ## number of significant proteins (at p_adj <0.05)
  print(paste("Number of significant proteins (p_adj < 0.05): ",table_1))  
  ## number of upregulated proteins
  print(paste("Number of upregulated proteins: ",table_2))  
  ## number of downregulated proteins
  print(paste("Number of downregulated proteins: ",table_3))  
  
  assign(named.df, df, envir=.GlobalEnv) ## move dataframe to global working environment
}

## 6) Generating tables with proportions // summary statistics for continuous variables
tblFun <- function(x){
  tbl <- table(x, useNA = "always")
  res <- cbind(tbl,round(prop.table(tbl)*100,2))
  colnames(res) <- c('Count','Percentage')
  res
} ## including missingness
tblFun_2 <- function(x){
  tbl <- table(x, useNA = "no")
  res <- cbind(tbl,round(prop.table(tbl)*100,2))
  colnames(res) <- c('Count','Percentage')
  res
} ## excluding missingness
summary_stats <- function(dataset, variable){
  summary_stats <- dataset %>%
    filter(sf_iknee_qc_group == "0" & baseline == "1" ) %>%
    ##group_by(cohort_name) %>%
    summarise(
      QC_group = sum(sf_iknee_qc_group == "0" & baseline == "1"),
      count_with_data = sum(!is.na({{variable}})),
      mean = mean({{variable}}, na.rm = TRUE),
      sd = sd({{variable}}, na.rm = TRUE),
      median = median({{variable}}, na.rm = TRUE),
      IQR = IQR({{variable}}, na.rm = TRUE),
      min = min({{variable}}, na.rm = TRUE),
      max = max({{variable}}, na.rm = TRUE),
    )
  
  print(summary_stats)
  dataframe_name <- deparse(substitute(dataset))
  variable_name  <- deparse(substitute(variable))
  summary <- paste0("Dataset summarised: ", dataframe_name, ", for continuous variable: ", variable_name, " in baseline OA samples only")
  print(summary)
}

## 7) Generating plots comparing agreement between two different models (e.g. non-IPS regressed and IPS regressed):
compare_TPRvsnonTPR <- function(model1, model2){
  
  ## plot estimates/standard errors/p-values
  cols <- rep("grey", nrow(model1))                                  ## proteins non-significant in both models (based on adjusted p-value =<0.05)
  cols[model1$sig_2 == "yes" & model2$sig_2 == "yes"] <- "purple"    ## proteins significant in both models 
  cols[model1$sig_2 == "no"  & model2$sig_2 == "yes"] <- "red"       ## proteins significant in non-total signal regressed data only 
  cols[model1$sig_2 == "yes" & model2$sig_2 == "no"]  <- "blue"      ## proteins significant in total signal regressed data only 
  
  plot(model1$estimate, model2$estimate, 
       main = "Plot of Log Odds (per SD) from Model1 vs Model2 Datasets", 
       xlab = "Log Odds per SD (Model1)", 
       ylab = "Log Odds per SD (Model2)", 
       col=cols, pch = 19, cex.main = 0.9)
  points(model1$estimate[cols != "grey"], model2$estimate[cols != "grey"], col=cols[cols != "grey"], pch = 19)
  abline(lm(model2$estimate ~ model1$estimate),lty=2)
  abline(h=0,v=0)
  
  sum1 <- sum(model1$sig_2 == "yes" & model2$sig_2 == "yes")
  sum2 <- sum(model1$sig_2 == "no"  & model2$sig_2 == "no")
  sum3 <- sum(model1$sig_2 == "no"  & model2$sig_2 == "yes")
  sum4 <- sum(model1$sig_2 == "yes" & model2$sig_2 == "no")
  agreement <- ((sum1 + sum2)/6558)*100 ##% agreement based on number of non-significant and significant proteins in both models 
  
  print(paste("Significant (p_adj < 0.05) in both models: ",sum1, paste("(purple)"))) 
  print(paste("Non-Significant (p_adj > 0.05) in both models: ",sum2, paste("(grey)"))) 
  print(paste("Significant (p_adj < 0.05) in Model2 only: ",sum3, paste("(red)"))) 
  print(paste("Significant (p_adj < 0.05) in Model1 only: ",sum4, paste("(blue)"))) 
  print(paste("Agreement (%): ", round(agreement, digits = 1)))
  
}    ## comparing log-odds
compare_TPRvsnonTPR2 <- function(model1, model2){
  
  ## plot estimates/standard errors/p-values
  cols <- rep("grey", nrow(model1))                                  ## proteins non-significant in both models (based on adjusted p-value =<0.05)
  cols[model1$sig_2 == "yes" & model2$sig_2 == "yes"] <- "purple"    ## proteins significant in both models 
  cols[model1$sig_2 == "no"  & model2$sig_2 == "yes"] <- "red"       ## proteins significant in non-total signal regressed data only 
  cols[model1$sig_2 == "yes" & model2$sig_2 == "no"]  <- "blue"      ## proteins significant in total signal regressed data only 
  
  plot(-log10(model1$p_adj), -log10(model2$p_adj), 
       main = "Plot of Ajusted P-values from Model1 vs Model2 Datasets", 
       xlab = "Ajusted p-values (Model1)", 
       ylab = "Adjusted p-values (Model2)", 
       col=cols, pch = 19, cex.main = 0.9)
  points(model1$p_adj[cols != "grey"], model2$p_adj[cols != "grey"], col=cols[cols != "grey"], pch = 19)
  abline(lm(model2$p_adj ~ model1$p_adj),lty=2)
  abline(h=0,v=0)
  
  sum1 <- sum(model1$sig_2 == "yes" & model2$sig_2 == "yes")
  sum2 <- sum(model1$sig_2 == "no"  & model2$sig_2 == "no")
  sum3 <- sum(model1$sig_2 == "no"  & model2$sig_2 == "yes")
  sum4 <- sum(model1$sig_2 == "yes" & model2$sig_2 == "no")
  agreement <- ((sum1 + sum2)/6250)*100 ##% agreement based on number of non-significant and significant proteins in both models 
  
  print(paste("Significant (p_adj < 0.05) in both models: ",sum1, paste("(purple)"))) 
  print(paste("Non-Significant (p_adj > 0.05) in both models: ",sum2, paste("(grey)"))) 
  print(paste("Significant (p_adj < 0.05) in Model2 only: ",sum3, paste("(red)"))) 
  print(paste("Significant (p_adj < 0.05) in Model1 only: ",sum4, paste("(blue)"))) 
  print(paste("Agreement (%): ", round(agreement, digits = 1)))
  
}   ## comparing adjusted p-values

## 8) Saving volcano plots as PDFs
saveplot <- function(filename, path, plot1, plot2, plot3, plot4, plot5, plot6){
  Plotlist <- list(plot1,plot2, plot3, plot4, plot5, plot6)
  ggsave(filename = filename,
         path = c(paste0(path)),
         device = NULL,
         units = c("in"),
         width = 5,
         height = 4,
         dpi = 300, 
         scale = 2,
         marrangeGrob(Plotlist, nrow = 1, ncol = 1)
  )
}

## 9) exporting and saving dataframes
saving_dfs <- function(directory_folder, search_term1, search_term2, search_term3){
  # set the path to the output directory
  dir_path <- paste0(base_path, "/Replication/")
  output_dir_base <- paste0(dir_path,directory_folder)
  # set the string character to search for in the file name
  search_string1 <-  paste(c(search_term1))
  search_string2 <-  paste(c(search_term2))
  search_string3 <-  paste(c(search_term3))
  # loop through each data frame and save to a different folder
  for (i in seq_along(dfs)) {
    # get the data frame name
    ## df_name <- deparse(substitute(dfs[[i]]))
    df_name = names(dfs)[i]
    # check if the data frame name contains the search string
    if (grepl(search_string1, df_name) & grepl(search_string2, df_name) & grepl(search_string3, df_name)) {
      # create the output directory if it doesn't exist
      output_dir <- file.path(output_dir_base)
      if (!dir.exists(output_dir)) {
        dir.create(output_dir)}
      # set the file name for the output file
      output_file <- file.path(output_dir, paste0(df_name, ".csv"))
      # write the data frame to a CSV file
      write.csv(dfs[[i]], file = output_file, row.names = TRUE)
      # confirm that the file was saved
      if (file.exists(output_file)) {
        cat(sprintf("File for %s saved successfully!\n", df_name))}
    } else {
      cat(sprintf("Skipping %s because it doesn't contain the search string!\n", df_name))}
  }
  
}
## need to change the path in this function to align with folder structure shown below
saving_dfs_dataframes <- function(directory_folder, search_term1, search_term2, search_term3){
  # set the path to the output directory
  dir_path <- paste0(base_path, "/")
  output_dir_base <- paste0(dir_path,directory_folder)
  # set the string character to search for in the file name
  search_string1 <-  paste(c(search_term1))
  search_string2 <-  paste(c(search_term2))
  search_string3 <-  paste(c(search_term3))
  # loop through each data frame and save to a different folder
  for (i in seq_along(dfs)) {
    # get the data frame name
    ## df_name <- deparse(substitute(dfs[[i]]))
    df_name = names(dfs)[i]
    # check if the data frame name contains the search string
    if (grepl(search_string1, df_name) & grepl(search_string2, df_name) & grepl(search_string3, df_name)) {
      # create the output directory if it doesn't exist
      output_dir <- file.path(output_dir_base)
      if (!dir.exists(output_dir)) {
        dir.create(output_dir)}
      # set the file name for the output file
      output_file <- file.path(output_dir, paste0(df_name, ".R"))
      # write the data frame to a CSV file
      ##write.csv(dfs[[i]], file = output_file, row.names = FALSE)
      save(dfs[[i]], file = output_file)
      # confirm that the file was saved
      if (file.exists(output_file)) {
        cat(sprintf("File for %s saved successfully!\n", df_name))}
    } else {
      cat(sprintf("Skipping %s because it doesn't contain the search string!\n", df_name))}
  }
  
} ## this functions saves all dataframes as .R files

## 10) generate VennDiagrams assessing overlap across two models
VennDiagram <- function(model1, model2, label_model1, label_model2){
  
  ## compare the agreement between models
  test <- merge(model1, model2, by = 0, all = TRUE, suffix=c(".Model1",".Model2"))
  row.names(test) <- test$Row.names
  test <- test[, !names(test) %in% c("Row.names")]
  ## to remove row.names if duplicated through the merge function
  
  ## number of SOMAmers that are statistically significant (p_adj < 0.05) in both models
  value1 <- length(which(test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes"))
  
  ## number of SOMAmers that are statistically significant (p_adj < 0.05) and upregulated in both models
  value2 <- length(which(test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 >0 & test$estimate.Model2 >0))
  
  ## number of SOMAmers that are statistically significant (p_adj < 0.05) and downregulated in both models
  value3 <- length(which(test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 <0 & test$estimate.Model2 <0))
  
  ## number of SOMAmers that are statistically significant (p_adj < 0.05) but regulated in different directions
  value4 <- length(which(test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 <0 & test$estimate.Model2 >0 |
                           test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 >0 & test$estimate.Model2 <0))
  
  print(paste("Significant (p_adj < 0.05) in both models: N=",value1)) 
  print(paste("Significant (p_adj < 0.05) and upregulated in both models: N=",value2)) 
  print(paste("Significant (p_adj < 0.05) and downregulated in both models: N=",value3)) 
  print(paste("Significant (p_adj < 0.05) in both models but with opposing affects: N=",value4))
  
  ## print SOMAmer labels if significant in both models but in opposing directions
  ##print(list(test$TargetFullName.Model2[test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 <0 & test$estimate.Model2 >0 |
  ##                                     test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 >0 & test$estimate.Model2 <0]))
  
  ##View(select(filter(test, sig_2.Model1 == "yes" & sig_2.Model2 == "yes"), estimate.Model1, p_adj.Model1, sig_2.Model1, TargetFullName.Model1, Target.Model1, estimate.Model2, p_adj.Model2, sig_2.Model2, TargetFullName.Model2, Target.Model2))
  
  ## VENN DIAGRAM 
  x <- list(
    A = rownames(model1[(model1$sig_2 == "yes" & model1$estimate >0|
                           model1$sig_2 == "yes" & model1$estimate <0),]), 
    B = rownames(model2[(model2$sig_2 == "yes" & model2$estimate >0 |
                           model2$sig_2 == "yes" & model2$estimate <0),]))
  
  names(x) <- c(label_model1 = label_model1,label_model2 = label_model2)
  
  test <- ggvenn(
    x, 
    fill_alpha = 0.5,
    stroke_color = "black",
    stroke_linetype = "solid",
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4 
  )
  
  test <- test + 
    labs(caption = paste0("Of N = ",value1, " SOMAmers significantly regulated in both models, \n N = ", value4, " were regulated in opposing directions.")) + 
    theme(plot.caption = element_text(hjust = 0))
  
  return(test)
  
}
Venn_top50percent_overlap <- function(model_1, model_2, label_model1, label_model2){
  
  ##model 1
  df_sorted_m1 <- model_1[order(model_1$p_adj, decreasing = FALSE),] ## order entire dataframe by adjusted p-values (smallest to largest)
  df_sorted_m1 <- df_sorted_m1[df_sorted_m1$sig_2 == "yes",]         ## subset to include only proteins signficantly associated with outcome (p_adj)
  df_subset_m1 <- df_sorted_m1[seq(nrow(df_sorted_m1)/2),]           ## subset the 50% of dataframe (rows that fall in the first 50% of dataframe)
  
  ## model 2
  df_sorted_m2 <- model_2[order(model_2$p_adj, decreasing = FALSE),] 
  df_sorted_m2 <- df_sorted_m2[df_sorted_m2$sig_2 == "yes",]
  df_subset_m2 <- df_sorted_m2[seq(nrow(df_sorted_m2)/2),]
  
  
  ## compare the agreement between models
  test <- merge(df_subset_m1, df_subset_m2, by = 0, all = TRUE, suffix=c(".Model1",".Model2"))
  row.names(test) <- test$Row.names
  test <- test[, !names(test) %in% c("Row.names")]
  ## to remove row.names if duplicated through the merge function
  
  ## number of SOMAmers that are statistically significant (p_adj < 0.05) in both models
  value1 <- length(which(test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes"))
  
  ## number of SOMAmers that are statistically significant (p_adj < 0.05) and upregulated in both models
  value2 <- length(which(test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 >0 & test$estimate.Model2 >0))
  
  ## number of SOMAmers that are statistically significant (p_adj < 0.05) and downregulated in both models
  value3 <- length(which(test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 <0 & test$estimate.Model2 <0))
  
  ## number of SOMAmers that are statistically significant (p_adj < 0.05) but regulated in different directions
  value4 <- length(which(test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 <0 & test$estimate.Model2 >0 |
                           test$sig_2.Model2 == "yes" & test$sig_2.Model1 == "yes" & test$estimate.Model1 >0 & test$estimate.Model2 <0))
  
  print(paste("Significant (p_adj < 0.05) in both models: N=",value1)) 
  print(paste("Significant (p_adj < 0.05) and upregulated in both models: N=",value2)) 
  print(paste("Significant (p_adj < 0.05) and downregulated in both models: N=",value3)) 
  print(paste("Significant (p_adj < 0.05) in both models but with opposing affects: N=",value4))
  
  ##View(select(filter(test, sig_2.Model1 == "yes" & sig_2.Model2 == "yes"), estimate.Model1, p_adj.Model1, sig_2.Model1, TargetFullName.Model1, Target.Model1, estimate.Model2, p_adj.Model2, sig_2.Model2, TargetFullName.Model2, Target.Model2))
  
  ## VENN DIAGRAM 
  x <- list(
    A = rownames(df_subset_m1[(df_subset_m1$sig_2 == "yes" & df_subset_m1$estimate >0|
                                 df_subset_m1$sig_2 == "yes" & df_subset_m1$estimate <0),]), 
    B = rownames(df_subset_m2[(df_subset_m2$sig_2 == "yes" & df_subset_m2$estimate >0 |
                                 df_subset_m2$sig_2 == "yes" & df_subset_m2$estimate <0),]))
  
  names(x) <- c(label_model1 = label_model1,label_model2 = label_model2)
  
  test <- ggvenn(
    x, 
    fill_alpha = 0.5,
    stroke_color = "black",
    stroke_linetype = "solid",
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4 
  )
  
  test <- test + 
    labs(caption = paste0("Of N = ",value1, " SOMAmers significantly regulated in both models, \n N = ", value4, " were regulated in opposing directions.")) + 
    theme(plot.caption = element_text(hjust = 0))
  
  return(test)
} ## overlap for the top 50% most significant proteins

## 11) function for counting number of datapoints on ggplot      
give.n <- function(x){
  return(c(y = median(x)*1.01, label = length(x))) 
  # experiment with the multiplier to find the perfect position
} 

## 12) function for merging 4 dataframes together and exporting as .csv
# Define the function
join_and_process_matrices <- function(Disc_matrix, matrix_rep, matrix_comb_no_cohort, matrix_comb_cohort, output_path) {
  # Perform the joins
  df <- full_join(Disc_matrix[, c("Row.names", "estimate", "p.value", "p_adj", "sig", "sig_2")], 
                  matrix_rep[, c("Row.names", "estimate", "p.value", "p_adj", "sig", "sig_2")], by = "Row.names", suffix = c("","_Rep")) %>%
    full_join(matrix_comb_no_cohort[, c("Row.names", "estimate", "p.value", "p_adj", "sig", "sig_2")], by = "Row.names", suffix = c("_Disc", "_combNoCohort")) %>%
    full_join(matrix_comb_cohort[, c("Row.names", "estimate", "p.value", "p_adj", "sig", "sig_2", "Target", "TargetFullName", "EntrezGeneID", "EntrezGeneSymbol")], by = "Row.names", suffix = c("_Rep", "_combCohort"))
  
  # Rename columns
  names(df)[names(df) == "estimate"] <- "estimate_combCohort"
  names(df)[names(df) == "p.value"]  <- "p.value_combCohort"
  names(df)[names(df) == "p_adj"]    <- "p_adj_combCohort"
  names(df)[names(df) == "sig"]      <- "sig_combCohort"
  names(df)[names(df) == "sig_2"]    <- "sig_2_combCohort"
  
  # Select and reorder columns
  df <- df[, c("estimate_Disc", "estimate_Rep", "estimate_combNoCohort", "estimate_combCohort", 
               "p.value_Disc", "p.value_Rep", "p.value_combNoCohort", "p.value_combCohort", 
               "p_adj_Disc", "p_adj_Rep", "p_adj_combNoCohort", "p_adj_combCohort", 
               "sig_Disc", "sig_Rep", "sig_combNoCohort", "sig_combCohort", 
               "sig_2_Disc", "sig_2_Rep", "sig_2_combNoCohort", "sig_2_combCohort", 
               "Target", "TargetFullName", "EntrezGeneID", "EntrezGeneSymbol", "Row.names")]
  
  # Write to CSV
  write.csv(df, output_path, row.names = TRUE)
}