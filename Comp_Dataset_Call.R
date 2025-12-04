# This script includes selfwritten functions to be called by Comp_Dep_Non.R, which make comparisons between expMsets based on features and samples.

#__________________________________________________________________________________________________________
# function to plot VennDiagram to visualise the Overlap of proteins and samples
makeMyVenn <- function(vector1,vector2,cateMessage){
  venn.plot <- draw.pairwise.venn(
    area1 = length(vector1),
    area2 = length(vector2),
    cross.area = length(intersect(vector1, vector2)),
    category = cateMessage,
    fill = c("darkgreen", "firebrick"),
    alpha = 0.5,
    lty = "solid",
    cex = 1.5,
    cat.cex = 1.5,
    cat.pos = c(-30, 30),
    cat.dist = 0.05
  )
  return()
}

#__________________________________________________________________________________________________________
# function to investigate the relationship between individual same one between two datasets
getIndividualRelation <- function(expM1, expM2, SampOrPro, whichMethod){
  
  commonSamp <- intersect(rownames(expM1),rownames(expM2))
  commonPro <- intersect(colnames(expM1),colnames(expM2))
  
  expM1 <- expM1[commonSamp,commonPro]
  expM2 <- expM2[commonSamp,commonPro]
  
  if(SampOrPro=="Sample"){
    individualCor <- sapply(commonSamp, function(x){
      vec1 <- as.numeric(expM1[x, ])
      vec2 <- as.numeric(expM2[x, ])
      
      # Remove NA values
      valid_idx <- complete.cases(vec1, vec2)
      
      # Check if there are at least two non-NA values for correlation
      if (sum(valid_idx) > 1){
        cor(vec1[valid_idx], vec2[valid_idx], use = "complete.obs", method = whichMethod)
      }else{NA}  # Return NA if there aren't enough valid values
    })
    names(individualCor) <- commonSamp
    
  }else{
    individualCor <- sapply(commonPro, function(x){
      vec1 <- as.numeric(expM1[, x])
      vec2 <- as.numeric(expM2[, x])
      
      # Remove NA values
      valid_idx <- complete.cases(vec1, vec2)
      
      # Check if there are at least two non-NA values for correlation
      if (sum(valid_idx) > 1) {
        cor(vec1[valid_idx], vec2[valid_idx], use = "complete.obs", method = whichMethod)
      }else{NA}
    })
    names(individualCor) <- commonPro
  }
  
  #hist(individualCor, breaks=100, main=paste0("Histogram of ", whichMethod, " Correlation for ",SampOrPro))
  return(individualCor)
}

getIndividualPairTtest <- function(expM1, expM2, SampOrPro){
  
  commonSamp <- intersect(rownames(expM1),rownames(expM2))
  commonPro <- intersect(colnames(expM1),colnames(expM2))
  
  expM1 <- expM1[commonSamp,commonPro]
  expM2 <- expM2[commonSamp,commonPro]
  
  if(SampOrPro=="Sample"){
    individualP <- sapply(commonSamp, function(x){
      vec1 <- as.numeric(expM1[x, ])
      vec2 <- as.numeric(expM2[x, ])
      
      # Remove NA values
      valid_idx <- complete.cases(vec1, vec2)
      
      # Check if there are at least two non-NA values for correlation
      if (sum(valid_idx) > 1){
        t.test(vec1, vec2, paired = TRUE)$p.value
      }else{NA}  # Return NA if there aren't enough valid values
    })
    names(individualP) <- commonSamp
    
  }else{
    individualP <- sapply(commonPro, function(x){
      vec1 <- as.numeric(expM1[, x])
      vec2 <- as.numeric(expM2[, x])
      
      # Remove NA values
      valid_idx <- complete.cases(vec1, vec2)
      
      # Check if there are at least two non-NA values for correlation
      if (sum(valid_idx) > 1) {
        t.test(vec1, vec2, paired = TRUE)$p.value
      }else{NA}
    })
    names(individualP) <- commonPro
  }
  
  individualAdjp <- p.adjust(individualP, method = "BH")
  
  #hist(individualCor, breaks=100, main=paste0("Histogram of ", whichMethod, " Correlation for ",SampOrPro))
  return(individualAdjp)
}

getIndividualJS <- function(expM1, expM2, SampOrPro, titleMessage) {
  # Find common samples and proteins (features)
  commonSamp <- intersect(rownames(expM1), rownames(expM2))
  commonPro <- intersect(colnames(expM1), colnames(expM2))
  
  # Subset the matrices to only include common elements
  expM1 <- expM1[commonSamp, commonPro, drop = FALSE]
  expM2 <- expM2[commonSamp, commonPro, drop = FALSE]
  
  # Function to compute JS divergence for a given row or column
  computeJS <- function(vec1, vec2) {
    # Convert to numeric and normalize
    vec1 <- as.numeric(vec1)
    vec2 <- as.numeric(vec2)
    
    if (sum(vec1) == 0 || sum(vec2) == 0) {
      return(NA)  # Return NA if either distribution is all zeros
    }
    
    P <- vec1 / sum(vec1)
    Q <- vec2 / sum(vec2)
    
    return(distance(rbind(P, Q), method = "jensen-shannon"))
  }
  
  # Compute JS divergence
  if (SampOrPro == "Sample") {
    JS_divergence <- sapply(commonSamp, function(x) computeJS(expM1[x, ], expM2[x, ]))
  } else {
    JS_divergence <- sapply(commonPro, function(x) computeJS(expM1[, x], expM2[, x]))
  }
  
  # Plot histogram
  hist(JS_divergence, breaks = 100, main = paste0(titleMessage, "\n",signif(length(which(is.na(JS_divergence)))/length(JS_divergence), digits=2), " Totally Different Distribution"), col = "skyblue", xlab = "JS Divergence")
  
  return(JS_divergence)
}

#__________________________________________________________________________________________________________
# function to compare pairwise KL divergence
getKLM <- function(expM,SampOrPro){
  if(SampOrPro=="Sample"){
    # Normalize the distributions to ensure they sum to 1
    expM <- expM / rowSums(expM)
    
    # Initialize an empty matrix to store pairwise KL divergences
    n <- nrow(expM)
    kl_matrix <- matrix(NA, nrow = n, ncol = n)
    
    # Compute pairwise KL Divergence
    for (i in 1:n) {
      for (j in 1:n) {
        # KL divergence from row i to row j
        kl_matrix[i, j] <- distance(as.matrix(expM[i, , drop = FALSE]), 
                                    as.matrix(expM[j, , drop = FALSE]), 
                                    method = "kullback-leibler")
      }
    }
  }else{
    expM <- expM / colSums(expM)
    n <- ncol(expM)
    kl_matrix <- matrix(NA, nrow = n, ncol = n)
    for (i in 1:n) {
      for (j in 1:n) {
        # KL divergence from row i to row j
        kl_matrix[i, j] <- distance(as.matrix(expM[, i, drop = FALSE]), 
                                    as.matrix(expM[, j, drop = FALSE]), 
                                    method = "kullback-leibler")
      }
    }
  }
  return(kl_matrix)
}

#__________________________________________________________________________________________________________
# function to compare pairwise JS divergence
getJSM <- function(expM,SampOrPro){
  if(SampOrPro=="Sample"){
    # Normalize the distributions to ensure they sum to 1
    expM <- expM/rowSums(expM)
    
    # Initialize an empty matrix to store pairwise JS divergences
    n <- nrow(expM)
    js_matrix <- matrix(NA, nrow = n, ncol = n)
    
    # Compute pairwise JS Divergence
    for (i in 1:n) {
      for (j in 1:i) {
        P_Q <- rbind(expM[i,], expM[j,])
        js_matrix[i,j] <- js_matrix[j,i] <- distance(P_Q, method = "jensen-shannon")
      }
    }
  }else{
    expM <- expM/colSums(expM)
    n <- ncol(expM)
    js_matrix <- matrix(NA, nrow = n, ncol = n)
    for (j in 1:n) {
      for (i in 1:j){
        P_Q <- rbind(expM[,i], expM[,j])
        js_matrix[i,j] <- js_matrix[j,i] <- distance(P_Q, method = "jensen-shannon")}
    }
  }
  
  diag(js_matrix) <- 0
  
  return(js_matrix)
}

#__________________________________________________________________________________________________________________________________________________________
# function to make interactive plot showing correlation between proteins or samples based on two data sets where different normalization stages are applied
interactiveRHO <- function(ProFrame1, ProFrame2, ProFrame3, ProFrame4, ColorVar, proOrSamp, whichCorMethod, titleMessage, xlabMessage, ylabMessage, legendMessage){
  CorVec1 <- getIndividualRelation(ProFrame1, ProFrame2, proOrSamp, whichCorMethod)
  CorVec2 <- getIndividualRelation(ProFrame3, ProFrame4, proOrSamp, whichCorMethod)
  
  common_name_F1_F2 <- intersect(names(CorVec1), names(CorVec2))
  cor_Frame <- data.frame(
    F1_Correlation = CorVec1[common_name_F1_F2],
    F2_Correlation = CorVec2[common_name_F1_F2],
    Observation = common_name_F1_F2,
    ColorRange = ColorVar[common_name_F1_F2]
  )
  
  if(class(ColorVar) != "numeric"){cor_Frame %<>% mutate(ColorRange = ifelse(is.na(ColorRange),"no mapping", ColorRange))}
  
  min_val <- min(cor_Frame$F1_Correlation, cor_Frame$F2_Correlation, na.rm = TRUE) 
  max_val <- max(cor_Frame$F1_Correlation, cor_Frame$F2_Correlation, na.rm = TRUE) 
  
  if(proOrSamp == "Sample"){
    titlePlus <- paste0("\n",length(Reduce(intersect, list(rownames(ProFrame1), rownames(ProFrame2), rownames(ProFrame3), rownames(ProFrame4))))," Common Samples Displayed") 
  }else{
    titlePlus <- paste0("\n", length(Reduce(intersect, list(colnames(ProFrame1), colnames(ProFrame2), colnames(ProFrame3), colnames(ProFrame4))))," Common Proteins Displayed")
  }
  xlabPlus <- paste0("\n(Based on ", length(intersect(rownames(ProFrame1), rownames(ProFrame2)))," Samples and ", length(intersect(colnames(ProFrame1), colnames(ProFrame2)))," Proteins)")
  ylabPlus <- paste0("\n(Based on ", length(intersect(rownames(ProFrame3), rownames(ProFrame4)))," Samples and ", length(intersect(colnames(ProFrame3), colnames(ProFrame4)))," Proteins)")
  
  # For different types of color variable, use different color scheme when generating interactive plot
  plotly1 <- if(class(ColorVar) == "character"){
    plot_ly(cor_Frame, 
            x = ~F1_Correlation, 
            y = ~F2_Correlation, 
            type = 'scatter', 
            mode = 'markers', 
            color = ~ColorRange, 
            colors = c("darkgreen", "firebrick", "royalblue", "gold", "pink", "cyan", "darkorchid", "orange", "darkblue", "azure")[1:length(unique(ColorVar))],  
            marker = list(size = 4, opacity = 0.9),
            text = ~paste(Observation, "(", signif(F1_Correlation, digits = 2), ",", signif(F2_Correlation, digits = 2), ")"), 
            hoverinfo = 'text')
  }else{
    plot_ly(cor_Frame, 
            x = ~F1_Correlation, 
            y = ~F2_Correlation, 
            type = 'scatter', 
            mode = 'markers', 
            color = ~ColorRange, 
            colors = c("#F2F2FF", "#C8A2C8", "#A070A0", "purple"), 
            marker = list(size = 6, line = list(color = "black", width = 0.3), colorbar = list(title = legendMessage)),
            text = ~paste(Observation, "(", signif(F1_Correlation, digits = 2), ",", signif(F2_Correlation, digits = 2), ")\n", legendMessage, ":", signif(ColorRange, digits = 2)), 
            hoverinfo = 'text')
  }
  plotly2 <- plotly1%>%
    layout(
      title = list(text = paste0(titleMessage, titlePlus), font = list(size = 12)),
      xaxis = list(title = list(text = paste0(xlabMessage, xlabPlus), font = list(size = 10))),
      yaxis = list(title = list(text = paste0(ylabMessage, ylabPlus), font = list(size = 10))),
      legend = list(title = list(text = legendMessage), font = list(size = 10), itemsizing='constant'),
      xaxis = list(range = c(min_val - 0.1, max_val + 0.1), scaleanchor = "y"), ### guarantee the rounding values are still in the plot area 
      yaxis = list(range = c(min_val - 0.1, max_val + 0.1), scaleanchor = "x"),
      shapes = list(list(
        type = "line", 
        x0 = min_val, 
        x1 = max_val, 
        xref = "x",
        y0 = min_val, 
        y1 = max_val,
        yref = "y",
        line = list(color = "black",dash = "dash")))
    )
  
  return(plotly2)
}

#__________________________________________________________________________________________________________________________________________________________
# function to make heatmap based on adjacency matrix
getRowOrder <- function(proDFhere, whichMethod, titleMessage){
  adjM <-cor(proDFhere, method = whichMethod)
  
  plotObj <- pheatmap(
    mat = abs(adjM),
    cluster_rows = TRUE, 
    cluster_cols = TRUE, 
    show_rownames = FALSE,  
    show_colnames = FALSE,
    color = colorRampPalette(c("ivory", "firebrick"))(1000), 
    main = titleMessage,
    silent = TRUE
  )
  
  row_order <- plotObj$tree_row$order  # Order of rows
  col_order <- plotObj$tree_col$order  # Order of columns
  ordered_row_names <- rownames(adjM)[row_order]
  
  return(ordered_row_names)
}

PlotAdjM <- function(proDFhere, whichMethod, titleMessage){
  adjM <-cor(proDFhere, method = whichMethod)
  
  plotObj <- pheatmap(
    mat = abs(adjM),
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    show_rownames = FALSE,  
    show_colnames = FALSE,
    color = colorRampPalette(c("ivory", "firebrick"))(1000), 
    main = titleMessage,
    silent = TRUE
  )
  
  return(plotObj)
}

#__________________________________________________________________________________________________________
# function to make Scatter plot with marginal histogram
scatterHist <- function(ProFrame1, ProFrame2, ProFrame3, ProFrame4, proOrSamp, whichCorMethod, titleMessage, xMessage, yMessage){
  
  CorVec1 <- getIndividualRelation(ProFrame1, ProFrame2, proOrSamp, whichCorMethod)
  CorVec2 <- getIndividualRelation(ProFrame3, ProFrame4, proOrSamp, whichCorMethod)
  
  common_name_F1_F2 <- intersect(names(CorVec1), names(CorVec2))
  
  cor_Frame <- cbind(CorVec1[common_name_F1_F2], CorVec2[common_name_F1_F2]) %>% as.data.frame()
  colnames(cor_Frame) <- c("x","y")
  
  scatter_plot <- ggplot(cor_Frame, aes(x = x, y = y)) + geom_abline(slope=1, intercept=0, linetype = "dashed") + 
    xlim(-1, 1) + ylim(-1, 1) +
    geom_point(alpha = 0.4, color = "firebrick", size=0.5) +
    theme_minimal() +
    labs(title = titleMessage, x=xMessage, y=yMessage) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12)
    ) + coord_fixed(ratio = 1)
  
  scatter_with_hist <- ggExtra::ggMarginal(
    scatter_plot,
    type = "histogram",
    fill = "lightblue",
    bins = 30
  )
  
  return(scatter_with_hist)
}

#__________________________________________________________________________________________________________
# function to visualize sample overlap on UMAP based on two protein expression profile frame
UmapOverlap <- function(ProFrame1, ProFrame2){
  commonSamp <- intersect(rownames(ProFrame1), rownames(ProFrame2))
  commonPro <- intersect(colnames(ProFrame1), colnames(ProFrame2))
  
  umap_1 <- umap(ProFrame1[commonSamp,commonPro])
  umap_2 <- umap(ProFrame2[commonSamp,commonPro])
  
  # Convert the UMAP results to data frames for easier manipulation
  umap_1_df <- data.frame(umap_1$layout)
  umap_2_df <- data.frame(umap_2$layout)
  
  # Add observation names (assuming row names or some unique identifiers)
  umap_1_df %<>% mutate("Observation" = commonSamp)
  umap_2_df %<>% mutate("Observation" = commonSamp)
  colnames(umap_1_df)[c(1,2)] <- colnames(umap_2_df)[c(1,2)]  <- c("x","y")
  
  i=10
  j=100
  Samp1 <- commonSamp[i]
  Samp2 <- commonSamp[j]
  
  # Plot using plotly
  fig_UMAP_overlay <- plot_ly() %>%
    add_trace(
      x = umap_1_df$x, y = umap_1_df$y, type = 'scatter', mode = 'markers',
      marker = list(color = 'cyan', size = 4),
      text = umap_1_df$Observation, hoverinfo = 'text',
      name = 'UMAP based on DEL'
    ) %>%
    add_trace(
      x = umap_2_df$x, y = umap_2_df$y, type = 'scatter', mode = 'markers',
      marker = list(color = 'pink', size = 4),
      text = umap_2_df$Observation, hoverinfo = 'text',
      name = 'UMAP based on NON'
    ) %>%
    add_annotations(
      x = umap_1_df[Samp1,1], y = umap_1_df[Samp1,2],
      ax = umap_1_df[Samp2,1], ay = umap_1_df[Samp2,2],
      xref = "x", yref = "y", axref = "x", ayref = "y",
      text = "two samples in DEL", showarrow = TRUE,
      arrowcolor = "darkgreen"
    ) %>%
    add_annotations(
      x = umap_2_df[Samp1,1], y = umap_2_df[Samp1,2],
      ax = umap_2_df[Samp2,1], ay = umap_2_df[Samp2,2],
      xref = "x", yref = "y", axref = "x", ayref = "y",
      text = "two samples in NON", showarrow = TRUE,
      arrowcolor = "royalblue"
    ) %>%
    layout(
      title = 'Overlayed UMAP Plots',
      xaxis = list(title = 'Dimension 1'),
      yaxis = list(title = 'Dimension 2'),
      showlegend = TRUE,
      legend = list(
        orientation = "h",   # Horizontal legend
        x = 0.5,             # Centered horizontally
        y = -0.2,            # Positioned below the plot
        xanchor = "center",
        yanchor = "top"
      )
    )
  
  return(fig_UMAP_overlay)
}

UmapOverlap2 <- function(ProFrame1, ProFrame2, catlog1Plus, catlog2Plus){
  set.seed(123)
  
  commonSamp <- intersect(rownames(ProFrame1), rownames(ProFrame2))
  commonPro <- intersect(colnames(ProFrame1), colnames(ProFrame2))
  
  umap_1 <- umap(ProFrame1[commonSamp,commonPro])
  umap_2 <- umap(ProFrame2[commonSamp,commonPro])
  
  # Convert the UMAP results to data frames for easier manipulation
  umap_1_df <- data.frame(umap_1$layout)
  umap_2_df <- data.frame(umap_2$layout)
  
  # Add observation names (assuming row names or some unique identifiers)
  umap_1_df %<>% mutate("Observation" = commonSamp)
  umap_2_df %<>% mutate("Observation" = commonSamp)
  colnames(umap_1_df)[c(1,2)] <- colnames(umap_2_df)[c(1,2)]  <- c("x","y")
  
  # Plot using plotly
  fig_UMAP_overlay <- plot_ly() %>%
    add_trace(
      x = umap_1_df$x, y = umap_1_df$y, type = 'scatter', mode = 'markers',
      marker = list(color = 'cyan', size = 4),
      text = umap_1_df$Observation, hoverinfo = 'text',
      name = paste0('UMAP based on ', catlog1Plus)
    ) %>%
    add_trace(
      x = umap_2_df$x, y = umap_2_df$y, type = 'scatter', mode = 'markers',
      marker = list(color = 'pink', size = 4),
      text = umap_2_df$Observation, hoverinfo = 'text',
      name = paste0('UMAP based on ', catlog2Plus)
    ) %>%
    layout(
      title = list(text = paste0("Sample Visualisation on UMAP\n", length(commonSamp), " Samples ", length(commonPro), " Proteins"), font = list(size = 12)),
      xaxis = list(title = 'Dimension 1'),
      yaxis = list(title = 'Dimension 2'),
      showlegend = TRUE,
      legend = list(
        orientation = "h",   # Horizontal legend
        x = 0.5,             # Centered horizontally
        y = -0.2,            # Positioned below the plot
        xanchor = "center",
        yanchor = "top"
      )
    )
  
  return(fig_UMAP_overlay)
}

#__________________________________________________________________________________________________________
# function to plot histogram comparisons between actual and random generated data
plotCompHist <- function(correlation_truth,correlation2_random,whichMethod){
  df <- data.frame(
    value = c(correlation_truth, correlation2_random),
    group = rep(c("True", "Random"), each = length(correlation_truth))
  )
  
  # Plot overlapping histograms
  pObj <- ggplot(df, aes(x = value, fill = group)) +
    geom_histogram(alpha = 0.5, bins = 100, position = "identity") + 
    scale_fill_manual(values = c("darkgreen", "firebrick")) + 
    theme_minimal() +
    labs(title = paste0("Comparison of Correlatoins derived from True Data Set and Random Data Set\nbased on ", whichMethod),
         x = "Correlation Coefficient",
         y = "Count",
         fill = "Group") +
    theme(legend.position = "top")
  
  return(pObj)
}

