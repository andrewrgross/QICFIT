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
    ref.tissue.data <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE] # Filter out missing values from tissue
    genes.present.in.ref <- row.names(ref.tissue.data)         # Declare the genes present in the reference
    genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(sample)) # Declare genes in reference missing from sample
    rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))  # Generate a zero data frame the size of the missing rows 
    row.names(rows.to.add) <- genes.missing.in.query           # Name rows after missing rows
    names(rows.to.add) <- names(sample)                        # Name column the same as the query 
    sample <- rbind(sample,rows.to.add)                        # Use rbind to make a full data frame containing exactly the genes in the reference
    sample <- sample[genes.present.in.ref, , drop = FALSE]     # Reorder sample to match reference
    spearman.input <- cbind(sample, ref.tissue.data)           # Bind sample and reference
    result <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]] # Perform spearman calculation
    result <- round(result[2] * 100, 1)                        # Round result
    spearman.results[tissue,] <- result                        # Add to results table
  }
  return(spearman.results)
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

### Import references
references <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.full.csv',header=TRUE,row.names=1)

### Query Data
### Normalized
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)
metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)


########################################################################
### Format

references.full <- references[2:ncol(references)]
references.loose <- gtex.low.sd.loose[2:ncol(gtex.low.sd.loose)]
references.neutral <- gtex.low.sd.neutral[2:ncol(gtex.low.sd.neutral)]
references.tight <- gtex.low.sd.tight[2:ncol(gtex.low.sd.tight)]
references.supertight <- gtex.low.sd.supertight[2:ncol(gtex.low.sd.supertight)]

### Define a dataframe of just adult hypothalamus
aHT <- TPMdata[c(7,8,9,10,12)]

### Combine samples with references
#shared.rows <- intersect(row.names(references.supertight),row.names(aHT))
#spearman.input <- cbind(aHT[shared.rows,],references.supertight[shared.rows,])

########################################################################
### Spearman, full pairwise

#spearman.results <- cor(spearman.input, method = 'spearman', use = 'complete.obs')
spearman.results <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]]
spearman.results <- round(spearman.results*100, 0)

spearman.results.trimmed <- spearman.results[,1:5]
order.spearman.mean <- apply(spearman.results.trimmed, 1, mean)
(spearman.results.trimmed <- spearman.results.trimmed[order(order.spearman.mean, decreasing = TRUE),])

########################################################################
### Spearman, Loop through references

samples <- aHT
references <- references.full

reference.set.name <- 'Ref.full'

spearman.results <- data.frame(rep(0,ncol(references)))
row.names(spearman.results) <- names(references)
names(spearman.results) <- names(samples)[1]

for(sample.number in 1:ncol(samples)) {
  print(sample.number)
  spearman.results[sample.number] <- spearman.calc(samples[sample.number],references)
}

row.means <- apply(spearman.results,1,mean)
spearman.results <- spearman.results[order(row.means, decreasing = TRUE),]
spearman.results.trimmed <- as.matrix(spearman.results)


########################################################################
### Generate heatmap

heatmap.2(spearman.results.trimmed,
          main = title, # heat map title
          cellnote = spearman.results.trimmed,
          notecol = "gray40",
          density.info="none",  # turns off density plot inside color legend
          notecex=0.7,
          trace="none",         # turns off trace lines inside the heat map
          distfun=dist,
#          col = 
          breaks = 30,
          margins =c(4,12),     # widens margins around plot
          Rowv="FALSE",
          dendrogram="none"     # only draw a row dendrogram
)

########################################################################
### Benchmark results

### For each sample (make loop)
target <- 'Brain...Hypothalamus'

benchmark.results <- data.frame(matrix(nrow = ncol(spearman.results), ncol = 4))
names(benchmark.results) <- c('Correct', 'Score', 'Discernment', 'Runner.up')

for(sample.num in 1:ncol(spearman.results)) {
  sample <- spearman.results[sample.num]
  target.match <- row.names(sample)[1] == target
  target.score <- sample[1,]
  runner.up.score <- sample[2,]
  discernment <- target.score - runner.up.score
  benchmark.row <- c(target.match,target.score,discernment,runner.up.score)
  benchmark.results[sample.num,] <- benchmark.row
}
benchmark.results[1] <- as.logical(benchmark.results[,1])
benchmark.results




### Summarize benchmark results
benchmark.summary <- data.frame(matrix(nrow = 1, ncol = 5))
names(benchmark.summary) <- c('Correct', 'Score.med', 'Score.sd', 'Discern.min', 'Reliability')

correct <- min(benchmark.results$Correct) == 1
score.med <- median(benchmark.results$Score)
score.sd <- round(sd(benchmark.results$Score),3)
discern.min <- min(benchmark.results$Discernment)
reliability <- length(which(benchmark.results$Score >= max(benchmark.results$Runner.up)))/nrow(benchmark.results) * 100

benchmark.summary.row <- c(correct, score.med, score.sd, discern.min, reliability)
row <- 1
benchmark.summary[row,] <- benchmark.summary.row
benchmark.summary[1] <- as.logical(benchmark.summary[,1])

row.names(benchmark.summary)[row] <- reference.set.name
benchmark.summary








########################################################################
########################################################################
########################################################################
### Scratchwork

########################################################################
### Perform comparison for given sample

### Select sample (replace with loop)
sample <- aHT[1]
reference.input <- references.full
spearman.results <- data.frame(rep(0,53))
row.names(spearman.results) <- names(references[2:54])

for(ref.tissue.num in 1:ncol(reference.input)) {
  sample.temp <- sample
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

spearman.results <- spearman.results[order(spearman.results, decreasing = TRUE), , drop = FALSE]


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