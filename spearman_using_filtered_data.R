### Spearman comparison using filtered data -- Andrew R Gross -- 2016/12/16
### 

########################################################################
### Header
library(gplots)
library(Hmisc)

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

########################################################################
### Data input
### Import sample key
sample.key <- read.csv('Z:/Data/Andrew/reference_data/gtex/Sample.key.sorted.csv',header=TRUE)

### Open tables containing all tissue for each level
gtex.low.sd.loose <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.loose.csv', row.names = 1)
gtex.low.sd.neutral <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.neutral.csv', row.names = 1)
gtex.low.sd.tight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.tight.csv', row.names = 1)
gtex.low.sd.supertight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.supertight.csv', row.names = 1)

### Import full references
references <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.full.csv',header=TRUE,row.names=1)

### Import known samples
reference.transcriptomes <- read.csv('Z:/Data/Andrew/reference_data/gtex/individual-ref-transcriptomes-for-testing.csv', row.names = 1) # ~10 seconds

########################################################################
### Format

ref.full <- references[2:ncol(references)]
ref.loose <- gtex.low.sd.loose[2:ncol(gtex.low.sd.loose)]
ref.neutral <- gtex.low.sd.neutral[2:ncol(gtex.low.sd.neutral)]
ref.tight <- gtex.low.sd.tight[2:ncol(gtex.low.sd.tight)]
ref.supertight <- gtex.low.sd.supertight[2:ncol(gtex.low.sd.supertight)]

### Define a dataframe of just adult hypothalamus
aHT <- TPMdata[c(7,8,9,10,12)]

########################################################################
### Generate lists

#samples.list <- list(TPMdata[c(7,8,9,10,12)], TPMdata[5:6])
#names(samples.list) <- c('iMN_87iCTR', 'iMN_201iCTR')

### Group reference transcriptomes into a list
samples.list <- list()
known.samples.names <- c()
for(tissue.num in 1:40) {
  pos.start <- ((tissue.num-1)*8)+2
  pos.end <- pos.start
  referenceColumns <- reference.transcriptomes[pos.start:(pos.start+7)]
  samples.list[[tissue.num]] <- referenceColumns
  known.samples.names <- c(known.samples.names,names(referenceColumns)[1])
}
names(samples.list) <- known.samples.names

#samples.list <- list(reference.transcriptomes[122:124])

### Generate reference list
references.list <- list(ref.full, ref.loose, ref.neutral, ref.tight, ref.supertight)
names(references.list) <- c('Ref.full', 'Ref.loose' ,'Ref.neutral', 'Ref. tight', 'Ref.supertight')

### Define list of targets
targets <- names(ref.full)[-c(24,25,31)]

########################################################################
### Spearman, run through loops

### Initialize empty results lists
spearman.results.list <- list()
benchmark.results.list <- list()

### Select samples --- Generate a list of samples to loop through
selection <- 32
samples <- samples.list[[selection]]
names(samples)[1]
(target <- targets[selection])

### Loop through reference sets
for(ref.set.num in 1:length(references.list)) {
    # ref.set.num <- 2
  
  ### Specify current reference list and its name
  reference.input <- references.list[[ref.set.num]]
  
  ###reference.input <- references.list[[1]]
  reference.set.name <- names(references.list)[ref.set.num]
  print(reference.set.name)  

  ### Initialize empty results table
  spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
  row.names(spearman.results) <- names(reference.input)        # Name empty results table
  names(spearman.results) <- names(sample)
  
  ### Calculate Spearman correlations
  for(sample.number in 1:ncol(samples)) {
    print(sample.number)
    spearman.results[sample.number] <- spearman.calc(samples[sample.number],reference.input)
  }
  ### Reorder results
  row.means <- apply(spearman.results,1,mean)
  spearman.results <- spearman.results[order(row.means, decreasing = TRUE), , drop = FALSE]
  
  ### Add spearman results to list
  spearman.results.list[[ref.set.num]] <- spearman.results
  names(spearman.results.list)[ref.set.num] <- reference.set.name
  
  ### Benchmark results
  benchmark.results <- benchmark.calc(spearman.results, target)
  
  ### Add benchmark results to list
  benchmark.results.list[[ref.set.num]] <- benchmark.results
  names(benchmark.results.list)[ref.set.num] <- reference.set.name
}

### Check output
#str(spearman.results.list)
#str(benchmark.results.list)

########################################################################
### Summarize benchmark results

### Summarize all benchmark results from list
benchmark.summary <- data.frame(matrix(nrow = 1, ncol = 5))
names(benchmark.summary) <- c('Correct', 'Score.med', 'Score.sd', 'Discern.min', 'Reliability')

for(ref.set.num in 1:length(benchmark.results.list)) {
  benchmark.results <- benchmark.results.list[[ref.set.num]]   # Call the current benchmark results table
  reference.set.name <- names(benchmark.results.list)[ref.set.num]
  print(reference.set.name)
  
  correct <- length(which(benchmark.results$Correct == 1))/nrow(benchmark.results) * 100  # Define metrics
  score.med <- median(benchmark.results$Score)
  score.sd <- round(sd(benchmark.results$Score),3)
  discern.min <- min(benchmark.results$Discernment)
  reliability <- length(which(benchmark.results$Score >= max(benchmark.results$Runner.up)))/nrow(benchmark.results) * 100
  
  benchmark.summary.row <- c(correct, score.med, score.sd, discern.min, reliability)
  benchmark.summary[ref.set.num,] <- benchmark.summary.row
  row.names(benchmark.summary)[ref.set.num] <- reference.set.name
}
print(benchmark.summary)


########################################################################
### Generate heatmap

spearman.results <- spearman.results.list[[1]]

heatmap.2(as.matrix(spearman.results),
          #main = title, # heat map title
          cellnote = spearman.results,
          notecol = "gray40",
          density.info="none",  # turns off density plot inside color legend
          notecex=0.7,
          trace="none",         # turns off trace lines inside the heat map
          distfun=dist,
          lhei = c(1,5),
          margins=c(10,12),
          breaks = 30,
          Rowv="FALSE",
          dendrogram="none"     # only draw a row dendrogram
)











########################################################################
########################################################################
########################################################################
### Scratchwork

########################################################################
### Perform comparison for given sample

### Select sample (replace with loop)
sample <- reference.transcriptomes[122]
names(sample)
reference.input <- ref.loose
spearman.results <- data.frame(rep(0,53))
row.names(spearman.results) <- names(references[2:54])

### Perform Spearman comparison one-to-one for each column of reference input
for(ref.tissue.num in 1:ncol(reference.input)) {
  
  sample.temp <- sample
  
  #ref.tissue.num <- 16
  
  ### Generate a data frame with only genes expressed in that tissue
  ref.tissue.data <- reference.input[ref.tissue.num]
  tissue <- names(ref.tissue.data)
  print(tissue)
  ref.tissue.data <- ref.tissue.data[which(ref.tissue.data[,1] >0), , drop = FALSE]
  
  ### Declare the genes expressed in the reference tissue
  genes.present.in.ref <- row.names(ref.tissue.data)   # Declare a list of genes from reference

  ### Declare genes in reference missing from query
  genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(sample.temp)) 

  ### Add missing rows to sample
  rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))  # Generate a zero data frame of correct size 
  row.names(rows.to.add) <- genes.missing.in.query  # Name rows after missing rows
  names(rows.to.add) <- names(sample.temp)  # Name column the same as the query 
  sample.temp <- rbind(sample.temp,rows.to.add)  # Use rbind to make a full data frame containing exactly the genes in the reference
  
  ### Combine sample and reference
  sample.temp <- sample.temp[genes.present.in.ref, , drop = FALSE]
  spearman.input <- cbind(sample.temp, ref.tissue.data)
  spearman.result <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]]
  result.2 <- round(spearman.result[2] * 100, 1)
  
  #spearman.result <- round(spearman.result * 100, 0)
  print(result.2)
  spearman.results[tissue,] <- result.2
}

(spearman.results <- spearman.results[order(spearman.results, decreasing = TRUE), , drop = FALSE])




### Formatting for benchmarking of full pairwise calculation
spearman.results.ref.only <- as.matrix(spearman.results.trimmed[6:nrow(spearman.results.trimmed),])
spearman.results.samp.only <- sort(as.matrix(spearman.results.trimmed[1:5,]), decreasing = TRUE)
spearman.results.samp.only <- spearman.results.samp.only[6:length(spearman.results.samp.only)]

### Min spearman value; lower is better
min.spearman <- min(spearman.results.ref.only)
### Max value: higher is better
max.ref <- max(spearman.results.ref.only)
### Lowest score among top hit
min.val.for.match <- min(spearman.results.ref.only[1,])
### Range for top hit
range.match <- max.ref - min.val.for.match

### Second highest reference value: lower is better

### Minimum within samples; undefined
min.sample <- min(spearman.results.samp.only)
### Maximum within samples (discounting selves)
max.sample <- max(spearman.results.samp.only)
### Variance within samples
var.sample <- max.sample - min.sample



references.loose.pruned <- references.loose[!apply(references.loose, 1, function(x){all(is.na(x))}),]
references.neutral.pruned <- references.neutral[!apply(references.neutral, 1, function(x){all(is.na(x))}),]
references.tight.pruned <- references.tight[!apply(references.tight, 1, function(x){all(is.na(x))}),]
references.supertight.pruned <- references.supertight[!apply(references.loose, 1, function(x){all(is.na(x))}),]