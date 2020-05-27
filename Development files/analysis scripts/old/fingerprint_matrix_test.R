### Generate and compare fingerprint matricies -- Andrew R Gross -- 2016-12-29
### A set of functions to generate fingerprint matricies and compare them

########################################################################
### Header
library(ggplot2)

########################################################################
### Functions
convertIDs <- function(dataframe) {                      # Remove the decimal from Ensembl IDs
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]         # Remove the decimal point and following digits from a string
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs                       # Replace the rownames with new ones
  return(dataframe)
}
addMedSD <- function(dataframe) {
  median <- apply(dataframe,1,median)
  sd <- apply(dataframe,1,sd)
  return(data.frame(dataframe,median,sd))  
}
sortByMed <- function(dataframe) {
  order <- order(dataframe$median,decreasing=TRUE)
  return(dataframe[order,])
}
########################################################################
### Data input
### Import each of the five transcriptome summary tables

### Open tables containing all tissue for each level
gtex.low.sd.supertight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.supertight.csv', row.names = 1)

### Import sample key
sample.key <- read.csv('Z:/Data/Andrew/reference_data/gtex/Sample.key.sorted.csv',header=TRUE)

### Import transcriptome location table
transcriptome.index <- read.csv("Z:/Data/Andrew/reference_data/gtex/transcriptome_tables_by_tissue/transcriptome.index.csv")
print(transcriptome.index[1])
selection <- 2
selected.tissue <- as.character(transcriptome.index[,1][selection])
print(selected.tissue)

### Import transcriptome of specific tissue
trans.df <- read.delim(as.character(transcriptome.index$file.locations[selection]))  # Load selected dataset (may take >10s)

### Load the housekeeping ids of interest
housekeeping.ids <- as.vector(read.table('Z:/Data/Andrew/QICFIT/housekeeping.ids.txt')[,1])

### Query Data
### Normalized
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)
metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)

### Define a dataframe of just adult hypothalamus
aHT <- TPMdata[c(7,8,9,10,12)]

########################################################################
### Format

ref.supertight <- gtex.low.sd.supertight[2:ncol(gtex.low.sd.supertight)]

########################################################################
### Generate fingerprint matricies for all tissues

reference.fingerprint.list <- list()

for(tissue in names(ref.supertight)) {
  print(tissue)
  ### Call the 36 housekeeping gene expression values
  housekeeping.vals <- tissue.exprs[housekeeping.ids, , drop = FALSE]
  housekeeping.vals <- 10^housekeeping.vals
  ### Call the top 100 genes of the current reference ******* TEMP *********
  tissue.exprs <- ref.supertight[tissue]
  tissue.exprs.sorted <- tissue.exprs[order(tissue.exprs, decreasing = TRUE), , drop = FALSE]
  top.100.ids <- row.names(tissue.exprs.sorted)[1:100]
  ### Generate the fingerprint matrix
  top.100.vals <- tissue.exprs[top.100.ids, , drop = FALSE]
  numerator <- matrix(10^top.100.vals)
  denominator <- t(1/hk)
  fingerprint.matrix <- numerator %*% denominator
  ### Add to list
  reference.fingerprint.list[[length(reference.fingerprint.list)+1]] <- list(top.100.ids, fingerprint.matrix)
  names(reference.fingerprint.list)[length(reference.fingerprint.list)] <- tissue
}

str(reference.fingerprint.list)

########################################################################
### Generate fingerprint matrix for each tissue for query

### Select query
query.exprs <- aHT[1]

### Call the 36 housekeeping gene expression values
housekeeping.vals <- query.exprs[housekeeping.ids, , drop = FALSE]

### Select tissue to compare to ******* TEMP *********
names(reference.fingerprint.list)
(tissue <- names(reference.fingerprint.list)[16])

### Call the top 100 genes of the current reference 
top.100.ids <- reference.fingerprint.list[[tissue]][[1]]

### Generate the fingerprint matrix
top.100.vals <- query.exprs[top.100.ids, , drop = FALSE]

numerator <- as.matrix(top.100.vals)
denominator <- t(1/housekeeping.vals)
fingerprint.matrix <- numerator %*% denominator

### find ratio of fingerprint matrices
current.reference.fingerprint <- reference.fingerprint.list[[tissue]][[2]]
comparison.matrix <- fingerprint.matrix/current.reference.fingerprint
comparison.matrix2 <- 10^(abs(log10(comparison.matrix)))-1
#round(comparison.matrix[1:10,1:8],2)
#comparison.matrix2[1:10,1:8]

### Generate a histogram of values
fingerprint.vector <- as.vector(comparison.matrix2)
fingerprint.vector <- sort(fingerprint.vector)
plot(density(fingerprint.vector))

summary(fingerprint.vector)

### Save summary to stats df

query.results <- data.frame(matrix(nrow = 53, ncol = 6))
names(query.results) <- c('Min','First','Median','Mean','Third','Max')
row.names(query.results) <- names(ref.supertight)
query.results[tissue,] <- summary

########################################################################
########################################################################
### Loop through all tissues
### Generate results df
query.results <- data.frame(matrix(nrow = 53, ncol = 6))
names(query.results) <- c('Min','First','Median','Mean','Third','Max')
row.names(query.results) <- names(ref.supertight)

### Select query
query.exprs <- aHT[1]
### Call the 36 housekeeping gene expression values
housekeeping.vals <- query.exprs[housekeeping.ids, , drop = FALSE]

for(tissue in names(reference.fingerprint.list)) {
  ### Call the top 100 genes of the current reference 
  top.100.ids <- reference.fingerprint.list[[tissue]][[1]]
  ### Generate the fingerprint matrix
  top.100.vals <- query.exprs[top.100.ids, , drop = FALSE]
  numerator <- as.matrix(top.100.vals)
  denominator <- t(1/housekeeping.vals)
  fingerprint.matrix <- numerator %*% denominator
  ### find ratio of fingerprint matrices
  current.reference.fingerprint <- reference.fingerprint.list[[tissue]][[2]]
  comparison.matrix <- fingerprint.matrix/current.reference.fingerprint
  comparison.matrix2 <- 10^(abs(log10(comparison.matrix)))-1
  ### Calculate stats
  fingerprint.vector <- as.vector(comparison.matrix2)
  fingerprint.vector <- sort(fingerprint.vector)
  ### Add to results df
  #print(summary(fingerprint.vector))
  query.results[tissue,] <- summary(fingerprint.vector)
}

query.results[order(query.results$Median),]
query.results[order(query.results$Min),]












### Examine columns and rows
row.means <- apply(comparison.matrix2,1,mean)
col.means <- apply(comparison.matrix2,2,mean)
row.means<- row.means[!is.na(row.means)]
col.means<- col.means[!is.na(col.means)]
row.mean.sum <- sum(row.means)
col.mean.sum <- sum(col.means)
(total <- row.mean.sum + col.mean.sum)











#####################################
testa <- fingerprint.matrix[1:5,1:5]
testb <- current.reference.fingerprint[2:6,2:6]

fingerprint.matrix[1:4,1:4]
current.reference.fingerprint[1:4,1:4]

fingerprint.heat <- round(log(fingerprint.matrix[10:20,1:4]),2)
heatmap.2(fingerprint.heat,cellnote = fingerprint.heat)





### Remove NA rows
ref.loose.pruned <- ref.loose[!apply(ref.loose, 1, function(x){all(is.na(x))}),]
ref.neutral.pruned <- ref.neutral[!apply(ref.neutral, 1, function(x){all(is.na(x))}),]
ref.tight.pruned <- ref.tight[!apply(ref.tight, 1, function(x){all(is.na(x))}),]
ref.supertight.pruned <- ref.supertight[!apply(ref.supertight, 1, function(x){all(is.na(x))}),]

ref.supertight.pruned.2 <- ref.supertight[!apply(ref.supertight, 1, function(x){any(is.na(x))}),]
temp <- 10^ref.supertight.pruned.2
