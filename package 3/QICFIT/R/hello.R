# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(Hmisc)
library(ggplot2)
library(reshape2)
library(gridExtra)

#' Spearman calculation using rcorr function for a single sample
#'
#'
#' @param query.df A data frame containing one sample in a single column with ENSEMBL IDs as row names.  ENSEMBL IDs should not have decimals.
#' @param ref.df A data frame containing the averaged GTEx reference transcriptomes with ENSEMBL IDs.  IDs should not have decimals.
#' @keywords
#' @export spearman.results A data frame of one column listing the spearman correlation coefficients between the queried sample and each reference tissue
#' @examples
#' spearman.corr.for.single.sample()

spearman.corr.for.single.sample <- function(query.df, ref.df) {   # Calculate the spearman correlation between the query.df and the references
  spearman.results <- data.frame(rep(0,ncol(ref.df)))                             # Generate empty results table
  row.names(spearman.results) <- names(ref.df)                                # Name empty results table
  names(spearman.results) <- names(query.df)

  for(ref.tissue.num in 1:ncol(ref.df)) {                                     # Loop through each reference
    ref.tissue.data <- ref.df[ref.tissue.num]                                 # Call the current tissue from the references
    tissue <- names(ref.tissue.data)
    ### Generate a data frame containing the two transcriptomes being compared
    ref.tissue.data <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE]         # Filter out missing values from tissue
    genes.present.in.ref <- row.names(ref.tissue.data)                                   # Declare the genes present in the reference
    genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(query.df))    # Declare genes in reference missing from query.df
    rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))                       # Generate a zero data frame the size of the missing rows
    row.names(rows.to.add) <- genes.missing.in.query                                       # Name rows after missing rows
    names(rows.to.add) <- names(query.df)                                           # Name column the same as the query
    query.df2 <- rbind(query.df,rows.to.add)                                # Use rbind to make a full data frame containing exactly the genes in the reference
    query.df2 <- query.df2[genes.present.in.ref, , drop = FALSE]           # Reorder query.df to match reference
    spearman.input <- cbind(query.df2, ref.tissue.data)                          # Bind query.df and reference
    result <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]]                     # Perform spearman calculation using rcorr
    result <- round(result[2], 5)                        # Round result
    spearman.results[tissue,] <- result                        # Add to results table
  }
  return(spearman.results)
}


#' Spearman calculation performed manually for a single sample
#'
#'
#' @param query.df A data frame containing one sample in a single column with ENSEMBL IDs as row names.  ENSEMBL IDs should not have decimals.
#' @param ref.df A data frame containing the averaged GTEx reference transcriptomes with ENSEMBL IDs.  IDs should not have decimals.
#' @param weighted If TRUE, delta squared values are weighted by rank.  Defaults to FALSE.
#' @keywords
#' @export spearman.results A data frame of one column listing the spearman correlation coefficients between the queried sample and each reference tissue
#' @examples
#' spearman.manual.for.single.sample()

cf <- function(counts) {return(counts*(counts^2-1)/12)}

spearman.manual.for.single.sample <- function(query.df, ref.df, weighted = FALSE) {   # Calculate the spearman correlation between the query.df and the references
  spearman.results <- data.frame(rep(0,ncol(reference.input)))                             # Generate empty results table
  row.names(spearman.results) <- names(ref.df)                                # Name empty results table
  names(spearman.results) <- names(query.df)

  for(ref.tissue.num in 1:ncol(ref.df)) {                                     # Loop through each reference
    ref.tissue.data <- ref.df[ref.tissue.num]                                 # Call the current tissue from the references
    tissue <- names(ref.tissue.data)
    ### Generate a data frame containing the two transcriptomes being compared
    ref.tissue.data <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE]         # Filter out missing values from tissue
    genes.present.in.ref <- row.names(ref.tissue.data)                                   # Declare the genes present in the reference
    genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(query.df))    # Declare genes in reference missing from query.df
    rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))                       # Generate a zero data frame the size of the missing rows
    row.names(rows.to.add) <- genes.missing.in.query                                       # Name rows after missing rows
    names(rows.to.add) <- names(query.df)                                           # Name column the same as the query
    query.df2 <- rbind(query.df,rows.to.add)                                # Use rbind to make a full data frame containing exactly the genes in the reference
    query.df2 <- query.df2[genes.present.in.ref, , drop = FALSE]           # Reorder query.df to match reference
    spearman.input <- cbind(query.df2, ref.tissue.data)                          # Bind query.df and reference

    ### Assign rank and delta rank
    spearman.input$sample_rank <- rank(-spearman.input[1])
    spearman.input$ref_rank <- rank(-spearman.input[2])
    spearman.input$d <- spearman.input$sample_rank - spearman.input$ref_rank
    spearman.input$d2 <- (spearman.input$d)^2

    if(weighted == TRUE) {
      ### Calculate weighted delta rank
      spearman.input$d2 <- spearman.input$d2 * ((-2*spearman.input$sample_rank/max(spearman.input$sample_rank))+2)
    }

    ### Define correction factor
    counts <- append(table(spearman.input$sample_rank), table(spearman.input$ref_rank))
    counts <- append(table(spearman.input[1]), table(spearman.input[2]))
    counts.counts <- table(counts)
    counts.df <- data.frame(counts.counts)
    counts.df$cf <- cf(as.numeric(counts.df[,1]))
    counts.df$cf_full <- counts.df[,2] * counts.df[,3]
    cf.sum <- sum(counts.df$cf_full)

    ### Calculate the correlation
    sum.di.sq <- sum(spearman.input$d2)
    sum.di.sq.corrected <- sum.di.sq + cf.sum
    result <- 1 - (6*sum.di.sq.corrected)/(nrow(spearman.input)*(nrow(spearman.input)^2-1))
    result <- round(result, 5)                        # Round result
    spearman.results[tissue,] <- result                        # Add to results table
  }
  return(spearman.results)
}


#' Spearman calculation for multiple samples
#'
#'
#' @param query.dataframe A data frame containing samples in a each column with ENSEMBL IDs as row names.  ENSEMBL IDs should not have decimals.
#' @param reference.dataframe A data frame containing the averaged GTEx reference transcriptomes with ENSEMBL IDs.  IDs should not have decimals.
#' @param method Options are 'corr', 'manual', and 'weighted'.  Defaults to 'corr'.
#' @param order.by.score If TRUE, outputs a table reordered by average score.  If FALSE, output table is in same order as reference data frame columns.  Defaults to TRUE
#' @keywords
#' @export spearman.results A data frame listing the spearman correlation coefficients between the queried samples and each reference tissue, with each sample in a column and each reference in a row
#' @examples
#' spearman.calc()

spearman.calc <- function(query.df, ref.df, method = 'corr', order.by.score = TRUE) {
  ### Initialize empty results table
  spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
  names(spearman.results) <- names(query.df)[1]                 # Name empty results table
  row.names(spearman.results) <- names(ref.df)        # Name empty results table

  ### Calculate Spearman correlations
  if(method == 'corr'){
    for(sample.number in 1:ncol(query.df)) {                      # Loop through all sample columns in query.df
      spearman.results[sample.number] <- spearman.corr.for.single.sample(query.df[sample.number],ref.df)
    }
  }
  if(method == 'manual'){
    for(sample.number in 1:ncol(query.df)) {                      # Loop through all sample columns in query.df
      spearman.results[sample.number] <- spearman.manual.for.single.sample(query.df[sample.number],ref.df)
    }
  }
  if(method == 'weighted'){
    for(sample.number in 1:ncol(query.df)) {                      # Loop through all sample columns in query.df
      spearman.results[sample.number] <- spearman.manual.for.single.sample(query.df[sample.number],ref.df, weighted = TRUE)
    }
  }

  spearman.results <- round(spearman.results*100,1)

  ### Reorder results
  if(order.by.score == TRUE) {
    row.means <- apply(spearman.results,1,mean)                  # Calculate the average score for a tissue across all 8 samples
    spearman.results <- spearman.results[order(row.means, decreasing = TRUE), , drop = FALSE] # Order the results from highest average tissue to lowest
  }
  return(spearman.results)
}


#' Plot spearman results
#'
#'
#' @param query.dataframe A data frame containing samples in a each column with ENSEMBL IDs as row names.  ENSEMBL IDs should not have decimals.
#' @param reference.dataframe A data frame containing the averaged GTEx reference transcriptomes with ENSEMBL IDs.  IDs should not have decimals.
#' @param method Options are 'corr', 'manual', and 'weighted'.  Defaults to 'corr'.
#' @param order.by.score If TRUE, outputs a table reordered by average score.  If FALSE, output table is in same order as reference data frame columns.  Defaults to TRUE
#' @keywords
#' @export spearman.results A data frame listing the spearman correlation coefficients between the queried samples and each reference tissue, with each sample in a column and each reference in a row
#' @examples
#' spearman.calc()

plot.spearman <- function(spearman.results) {
  plot.list <- list()
  for (col.num in 1:ncol(spearman.results)) {
    spearman.for.plot <- spearman.results[1:20,,drop = FALSE][col.num]
    spearman.for.plot <- spearman.for.plot[order(spearman.for.plot[1], decreasing = TRUE),, drop = FALSE]
    spear.m <- melt(as.matrix(spearman.for.plot))
    spear.m$Var1 <- factor(spear.m$Var1, levels = rev(as.character(row.names(spearman.for.plot))))
    spear.m$label <- paste0(spear.m$value, ' - ', spear.m$Var1)
    ## Plot heatmap
    plot <- ggplot(data = spear.m, aes(x = Var2, y = Var1, fill = value)) +
      geom_tile() +
      geom_text(aes(label = label), hjust = 0, nudge_x = -0.4, size = 3) +
      scale_fill_gradient(low = 'white', high = 'red') +
      #ggtitle(as.character(dataset.metadata[selection.number,][,1])) +
      theme(axis.text.x = element_text(angle = -90, hjust = 1),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            legend.position="none")
    plot.list[col.num] <- list(plot)
  }
  ### Convert plot list to list of graphic objects (grobs)
  grob.list <- lapply(plot.list, ggplotGrob)
  ### Plot grob.list
  return(grid.arrange(grobs = grob.list, ncol = ncol(spearman.results), top = tissue.name[1]))
}

#' Spearman gene responsibility assessment
#'
#'
#' @param query.df A data frame containing one sample in a single column with ENSEMBL IDs as row names.  ENSEMBL IDs should not have decimals.
#' @param ref.df A data frame containing the averaged GTEx reference transcriptomes with ENSEMBL IDs.  IDs should not have decimals.
#' @param weighted If TRUE, delta squared values are weighted by rank.  Defaults to FALSE.
#' @keywords
#' @export spearman.results A data frame of one column listing the spearman correlation coefficients between the queried sample and each reference tissue
#' @examples
#' spearman.manual.for.single.sample()

spearman.gene.assessor.single.sample <- function(query.df, ref.df, weighted = FALSE) {   # Calculate the spearman correlation between the query.df and the references
  spearman.results <- data.frame(rep(0,ncol(reference.input)))                             # Generate empty results table
  row.names(spearman.results) <- names(ref.df)                                # Name empty results table
  names(spearman.results) <- names(query.df)

  for(ref.tissue.num in 1:ncol(ref.df)) {                                     # Loop through each reference
    ref.tissue.data <- ref.df[ref.tissue.num]                                 # Call the current tissue from the references
    tissue <- names(ref.tissue.data)
    ### Generate a data frame containing the two transcriptomes being compared
    ref.tissue.data <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE]         # Filter out missing values from tissue
    genes.present.in.ref <- row.names(ref.tissue.data)                                   # Declare the genes present in the reference
    genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(query.df))    # Declare genes in reference missing from query.df
    rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))                       # Generate a zero data frame the size of the missing rows
    row.names(rows.to.add) <- genes.missing.in.query                                       # Name rows after missing rows
    names(rows.to.add) <- names(query.df)                                           # Name column the same as the query
    query.df2 <- rbind(query.df,rows.to.add)                                # Use rbind to make a full data frame containing exactly the genes in the reference
    query.df2 <- query.df2[genes.present.in.ref, , drop = FALSE]           # Reorder query.df to match reference
    spearman.input <- cbind(query.df2, ref.tissue.data)                          # Bind query.df and reference

    ### Assign rank and delta rank
    spearman.input$sample_rank <- rank(-spearman.input[1])
    spearman.input$ref_rank <- rank(-spearman.input[2])
    spearman.input$d <- spearman.input$sample_rank - spearman.input$ref_rank
    spearman.input$d2 <- (spearman.input$d)^2

    if(weighted == TRUE) {
      ### Calculate weighted delta rank
      spearman.input$d2 <- spearman.input$d2 * ((-2*spearman.input$ref_rank/max(spearman.input$ref_rank))+2)
    }

    gene.report <- spearman.input[order(spearman.input$d2, decreasing = TRUE),]
    # Add to results table
  }
  return(spearman.results)
}

#' Spearman gene responsibility assessment, pt 2: Bioconductor conversion and filtering
#'
#'
#' @param query.df A data frame containing one sample in a single column with ENSEMBL IDs as row names.  ENSEMBL IDs should not have decimals.
#' @param ref.df A data frame containing the averaged GTEx reference transcriptomes with ENSEMBL IDs.  IDs should not have decimals.
#' @param weighted If TRUE, delta squared values are weighted by rank.  Defaults to FALSE.
#' @keywords
#' @export spearman.results A data frame of one column listing the spearman correlation coefficients between the queried sample and each reference tissue
#' @examples
#' spearman.manual.for.single.sample()



