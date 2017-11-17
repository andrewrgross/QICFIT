### Spearman comparison using filtered data -- Andrew R Gross -- 2016/12/16
### This script uses GTEx data as the input for testing various strategies of recognizing samples and displaying results.

########################################################################
### Header
library(gplots)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(grid)
#install.packages("gridExtra")
library("gridExtra")
library("cowplot")

########################################################################
### Functions
spearman.calc <- function(sample, reference.input) {           # Calculate the spearman correlation between the sample and the references
  spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
  row.names(spearman.results) <- names(reference.input)        # Name empty results table
  names(spearman.results) <- names(sample)
  for(ref.tissue.num in 1:ncol(reference.input)) {             # Loop through each reference 
    ref.tissue.data <- reference.input[ref.tissue.num]         # Call the current tissue from the references
    tissue <- names(ref.tissue.data)
    ref.tissue.data.2 <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE] # Filter out missing values from tissue
    genes.present.in.ref <- row.names(ref.tissue.data.2)         # Declare the genes present in the reference
    genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(sample)) # Declare genes in reference missing from sample
    rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))  # Generate a zero data frame the size of the missing rows 
    row.names(rows.to.add) <- genes.missing.in.query           # Name rows after missing rows
    names(rows.to.add) <- names(sample)                        # Name column the same as the query 
    sample.2 <- rbind(sample,rows.to.add)                        # Use rbind to make a full data frame containing exactly the genes in the reference
    sample.3 <- sample.2[genes.present.in.ref, , drop = FALSE]     # Reorder sample to match reference
    spearman.input <- cbind(sample.3, ref.tissue.data.2)           # Bind sample and reference
    result <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]] # Perform spearman calculation
    result <- round(result[2] * 100, 1)                        # Round result
    spearman.results[tissue,] <- result                        # Add to results table
  }
  return(spearman.results)
}
benchmark.calc <- function(spearman.results, target) {
  benchmark.results <- data.frame(matrix(nrow = ncol(spearman.results), ncol = 4)) # Generate empty table for results
  names(benchmark.results) <- c('Correct', 'Score', 'Discernment', 'Runner.up')
  
  for(sample.num in 1:ncol(spearman.results)) {                                # Loop through each sample in the results table
    sample <- spearman.results[sample.num]                                     # Call the sample
    target.match <- row.names(sample)[1] == target                             # Check each metric
    target.score <- sample[1,]
    runner.up.score <- sample[2,]
    discernment <- target.score - runner.up.score
    benchmark.row <- c(target.match,target.score,discernment,runner.up.score)  # Concatenate metrics
    benchmark.results[sample.num,] <- benchmark.row
  }
  benchmark.results[1] <- as.logical(benchmark.results[,1])
  return(benchmark.results)
}
spearman.calc.for.multiple.samples <- function(samples, reference.input) {
  ### Initialize empty results table
  spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
  row.names(spearman.results) <- names(reference.input)        # Name empty results table
  
  ### Calculate Spearman correlations
  for(sample.number in 1:ncol(samples)) {                      # Loop through all 8 samples
    spearman.results[sample.number] <- spearman.calc(samples[sample.number],reference.input)
  }
  ### Reorder results
  row.means <- apply(spearman.results,1,mean)                  # Calculate the average score for a tissue across all 8 samples
  spearman.results <- spearman.results[order(row.means, decreasing = TRUE), , drop = FALSE] # Order the results from highest average tissue to lowest
}
### Plot single boxplot from results table function
plot.tissue.match.boxplot <- function(spearman.results, title) {
  
  ### Reshape data for compatibility with geom_boxplot
  spearman.results[9] <- row.names(spearman.results)                    # Add the tissues as a new row
  names(spearman.results)[9] <- 'ref'                                   # Label the new row
  spearman.t <- t(spearman.results[-9])                                 # Transpose the data frame, minus new row
  spearman.melt <- melt(spearman.t)                                     # Melt transposed data frame
  names(spearman.melt) <- c("sample","ref","value")                     # Rename columns of melted data frame
  spearman.melt$ref <- factor(spearman.melt$ref, levels = rev(row.names(spearman.results)), ordered = TRUE) # Assign order
  
  ### Generate data frame for color and label data
  spearman.for.color <- spearman.results                                # Duplicate spearman results
  spearman.for.color[1:8] <- apply(spearman.for.color[1:8], 1, mean)    # Reassign all values to mean of group
  names(spearman.for.color)[2] = 'color.val'
  spearman.t <- t(spearman.for.color[-9])                               # Transpose again
  spearman.melt.color <- melt(spearman.t)                               # Melt again
  spearman.melt$color <- spearman.melt.color$value                      # Copy repeating values as 'color column to main df
  
  g <- ggplot(data = spearman.melt, aes(x = ref, y = -value, fill = color)) +
    geom_boxplot() +
    coord_flip() +
    scale_fill_gradientn(colors = c('red','orange','yellow','white')) +
    scale_y_continuous(limits = c(-100,-50), position = 'right') +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = 'none', panel.background = element_rect(fill = 'grey97'),
          plot.title = element_text(size = 14, face = 'bold', margin = margin(5,0,5,0), hjust = 0.5)) +
    labs(title = row.names(spearman.results)[1],
         x = '',
         y = title) +
    geom_text(data = spearman.for.color, aes(x = ref, y = -color.val+8, fill = color.val, label = ref), 
              hjust = 0)
  
  return(g)
}
### Plot multiplot with all five filter levels
multiplot.spearman.results.list <- function(spearman.results.list) {
  boxplot.list <- list()
  filter.levels <- c('Full', 'Loose', 'Neutral', 'Tight', 'Supertight')
  level.num = 1
  for (spearman.results in spearman.results.list) {
    filter.level = filter.levels[level.num]
    level.num <- level.num + 1
    boxplot.list[[length(boxplot.list)+1]] <- plot.tissue.match.boxplot(spearman.results, filter.level)
  }
  print('All figures complete')
  multiboxplot <- plot_grid(boxplot.list[[1]], boxplot.list[[2]], boxplot.list[[3]], boxplot.list[[4]], boxplot.list[[5]], labels=c("F", "L", "N", "T", "ST"), ncol = 5, nrow = 1)
  #multiboxplot <- multiplot(boxplot.list[[1]], boxplot.list[[2]], boxplot.list[[3]], boxplot.list[[4]], boxplot.list[[5]], cols = 5)
  return(multiboxplot)
}
########################################################################
### Data input
### Import sample key -- This table includes the ID, Tissue type, Tissue group, and unique name for all 8555 GTEx samples
sample.key <- read.csv('Z:/Data/Andrew/reference_data/gtex/Sample.key.sorted.csv',header=TRUE)

### Open tables containing all tissue for each level -- These tables each contain the averages of each tissue, filtered for expression and sd
gtex.low.sd.loose <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.loose.csv', row.names = 1)
gtex.low.sd.neutral <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.neutral.csv', row.names = 1)
gtex.low.sd.tight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.tight.csv', row.names = 1)
gtex.low.sd.supertight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.supertight.csv', row.names = 1)

### Import full references  -- The average expression of each tissue without filtering
references <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.full.csv',header=TRUE,row.names=1)

### Import known samples  -- A table containing transcriptomes of 8 samples from each tissue in the GTEx collection
reference.transcriptomes <- read.csv('Z:/Data/Andrew/reference_data/gtex/individual-ref-transcriptomes-for-testing.csv', row.names = 1) # ~10 seconds

########################################################################
### Format

### Remove Description column
ref.full <- references[2:ncol(references)]
ref.loose <- gtex.low.sd.loose[2:ncol(gtex.low.sd.loose)]
ref.neutral <- gtex.low.sd.neutral[2:ncol(gtex.low.sd.neutral)]
ref.tight <- gtex.low.sd.tight[2:ncol(gtex.low.sd.tight)]
ref.supertight <- gtex.low.sd.supertight[2:ncol(gtex.low.sd.supertight)]

### Define a dataframe of just adult hypothalamus
aHT <- TPMdata[c(7,8,9,10,12)]

########################################################################
### Generate lists to select tissues and samples

#samples.list <- list(TPMdata[c(7,8,9,10,12)], TPMdata[5:6])
#names(samples.list) <- c('iMN_87iCTR', 'iMN_201iCTR')

### Generate a list of known trankscriptomes from each 8-sample collection
samples.list <- list()
known.samples.names <- c()

for(tissue.num in 1:50) {
  pos.start <- ((tissue.num-1)*8)+2
  pos.end <- pos.start
  referenceColumns <- reference.transcriptomes[pos.start:(pos.start+7)]
  samples.list[[tissue.num]] <- referenceColumns
  known.samples.names <- c(known.samples.names,names(referenceColumns)[1])
}
names(samples.list) <- known.samples.names

### Generate a list of filtered reference sets
references.list <- list(ref.full, ref.loose, ref.neutral, ref.tight, ref.supertight)
names(references.list) <- c('Ref.full', 'Ref.loose' ,'Ref.neutral', 'Ref. tight', 'Ref.supertight')

### Generate a list of tissues that can be specified as targets
targets <- names(ref.full)[-c(24,25,31)]


########################################################################
########################################################################
### Plot samples compared against references

########################################################################
### Group of samples against one reference set

### Specify a set of samples
names(samples.list)
selection = 10
samples <- samples.list[[selection]]
names(samples)[1]
(target <- targets[selection])
print(str(samples))

### Select a reference set (based on filter level)
names(references.list)
ref.set.num = 1
reference.input <- references.list[[ref.set.num]]              # Specify current reference list and its name
(reference.set.name <- names(references.list)[ref.set.num])    # Declare the name of the rererence filter set

### Define title
(title <- paste(target, 'vs', reference.set.name))

### Calculate spearman corr. for all samples in group against a chosen reference
spearman.results <- spearman.calc.for.multiple.samples(samples[1], reference.input)

### Generate a single plot for a single reference set
g <- plot.tissue.match.boxplot(spearman.results, title)
  
g

########################################################################
### Group of samples against all reference sets

### Specify a set of samples
names(samples.list)
selection = 10
samples <- samples.list[[selection]]
names(samples)[1]
(target <- targets[selection])
print(str(samples))

### Initialize empty results lists
spearman.results.list <- list()

### Loop through reference sets
for(ref.set.num in 1:length(references.list)) {                # Loop through all five reference filter levels
  
  reference.input <- references.list[[ref.set.num]]            # Specify current reference list and its name
  (reference.set.name <- names(references.list)[ref.set.num])  # Declare the name of the rererence filter set
  
  ### Initialize empty results table
  spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
  row.names(spearman.results) <- names(reference.input)        # Name empty results table
  names(spearman.results) <- names(samples)
  
  ### Calculate spearman corr. for all samples in group against a chosen reference
  spearman.results <- spearman.calc.for.multiple.samples(samples, reference.input)
  
  ### Add spearman results to list
  spearman.results.list[[ref.set.num]] <- spearman.results     # Add current Spearman results table of 8 samples compared to all tissues to a list
  names(spearman.results.list)[ref.set.num] <- reference.set.name # Name the new entry on the list with the tissue being assessed
}
### Generate multiplot

multiplot.boxplot <- multiplot.spearman.results.list(spearman.results.list)



### AUTOMATED LOOP #####################################################
########################################################################
### Generate a multiplot of each tissue

setwd("z:/Data/Andrew/QICFIT/Spearman results testing/")
benchmark.summary.list <- list()

for(selection in 40:length(samples.list)) {
  
  ########################################################################
  ### Spearman, run through loops
  
  ### Initialize empty results lists
  spearman.results.list <- list()
  benchmark.results.list <- list()
  
  ### Select samples --- Generate a list of samples to loop through
  samples <- samples.list[[selection]]
  names(samples)[1]
  (target <- targets[selection])
  
  ### Loop through reference sets
  for(ref.set.num in 1:length(references.list)) {
    
    ### Specify current reference list and its name
    reference.input <- references.list[[ref.set.num]]
    reference.set.name <- names(references.list)[ref.set.num]
    print(reference.set.name)  
    
    ### Initialize empty results table
    spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
    row.names(spearman.results) <- names(reference.input)        # Name empty results table
    names(spearman.results) <- names(samples)
    
    ### Calculate Spearman correlations
    for(sample.number in 1:ncol(samples)) {
      spearman.results[sample.number] <- spearman.calc(samples[sample.number],reference.input)
    }
    ### Reorder results
    row.means <- apply(spearman.results,1,mean)                  # Calculate the average score for a tissue across all 8 samples
    spearman.results <- spearman.results[order(row.means, decreasing = TRUE), , drop = FALSE] # Order the results from highest average tissue to lowest
    
    ### Add spearman results to list
    spearman.results.list[[ref.set.num]] <- spearman.results     # Add current Spearman results table of 8 samples compared to all tissues to a list
    names(spearman.results.list)[ref.set.num] <- reference.set.name # Name the new entry on the list with the tissue being assessed
    
    ### Benchmark results, add to list
    benchmark.results <- benchmark.calc(spearman.results, target) # Run the results through the benchmark function to 
    benchmark.results.list[[ref.set.num]] <- benchmark.results    # Add the benchmark results (four metrics) to the list
    names(benchmark.results.list)[ref.set.num] <- reference.set.name # Name the new entry on the list
  }
  ########################################################################
  ### Generate multiplot
  
  multiplot.boxplot <- multiplot.spearman.results.list(spearman.results.list)
  
  png(filename=paste0(target,".png"), 
      type="cairo",
      units="in", 
      width=22, 
      height=12, 
      pointsize=12, 
      res=120)
  print(multiplot.boxplot)
  dev.off()
}
names(benchmark.summary.list) <- names(samples.list)


















########################################################################
########################################################################
########################################################################
### Scratchwork

test <- plot_grid(g, g, labels=c("A", "B"), ncol = 2, nrow = 1)

########################################################################
### Summarize benchmark results from list

### Initialize summary data.frame, name
benchmark.summary <- data.frame(matrix(nrow = 1, ncol = 5))
names(benchmark.summary) <- c('Correct', 'Score.med', 'Score.sd', 'Discern.min', 'Reliability')

for(ref.set.num in 1:length(benchmark.results.list)) {         # Loop through all five filter levels
  benchmark.results <- benchmark.results.list[[ref.set.num]]   # Call benchmark results from each filter level
  reference.set.name <- names(benchmark.results.list)[ref.set.num]
  print(reference.set.name)
  
  ### Calculate each of the five metrics for each of the five filter levels
  correct <- length(which(benchmark.results$Correct == 1))/nrow(benchmark.results) * 100  # How many of the 8 matched the target?
  score.med <- median(benchmark.results$Score)                 # What was the median score for the top hit of each
  score.sd <- round(sd(benchmark.results$Score),3)             # What was the standard deviation across the scores
  discern.min <- min(benchmark.results$Discernment)            # How far was the closest second place score from the first?
  reliability <- length(which(benchmark.results$Score >= max(benchmark.results$Runner.up)))/nrow(benchmark.results) * 100 # Percent of top scores above the highest 2nd place score
  
  benchmark.summary.row <- c(correct, score.med, score.sd, discern.min, reliability) # Define the new row
  benchmark.summary[ref.set.num,] <- benchmark.summary.row     # Add the new row
  row.names(benchmark.summary)[ref.set.num] <- reference.set.name # Name new row
}
print(benchmark.summary)