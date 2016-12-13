### Transcriptome Preprocessing -- Andrew R Gross -- 2016/12/02
### Make minor edits to formatting and begin identifying key genes

########################################################################
### Header
########################################################################

library(ggplot2)

########################################################################
### Functions
########################################################################
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
addGene <- function(dataframe) {
  genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) {
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}
########################################################################
### Data input
########################################################################

### Import transcriptome location table
transcriptome.index <- read.csv("Z:/Data/Andrew/reference_data/gtex/transcriptome_tables_by_tissue/transcriptome.index.csv")
print(transcriptome.index[1])
selection <- 3
print(transcriptome.index[,1][selection])

### Import transcriptome of specific tissue
#transcriptome <- read.delim('Z:/Data/Andrew/reference_data/gtex/transcriptome_tables_by_tissue/adrenalgland')
trans.df <- read.delim(as.character(transcriptome.index$file.locations[selection]))  # Load selected dataset (may take >10s)

### Import sample key
sample.key <- read.csv('Z:/Data/Andrew/reference_data/gtex/Sample.key.sorted.csv',header=TRUE)

########################################################################
### Formatting
########################################################################
### Format names

matched <- match(sample.key$Sample_ID,names(trans.df))           # Match the sample IDs in the sample key with IDs in the current file
matched.2 <- matched[!is.na(matched)]                            # Retain non-na matches

new.names <- as.vector(sample.key$unique.names[matched.2])       # Generate a list of new names
names(trans.df)[3:ncol(trans.df)] <- new.names                   # Assign the new names to the current samples

########################################################################
### Add med and SD

trans.med <- addMedSD(trans.df[3:ncol(trans.df)])                # Make new DF with median & sd
ncol(trans.med)

med.v.sd <- trans.med[146:147]                                   # Make new DF with just median & sd

med.v.sd$ratio <- med.v.sd$median/med.v.sd$sd                    # Add ratio column to med-sd df

########################################################################
### Log transform

med.v.sd <- log(med.v.sd+1)                                      # Log transform

########################################################################
### Plot log med. vs. log SD
########################################################################

plot(log(med.v.sd$median),log(med.v.sd$sd),pch=20,cex = 1)      # Plot crudely

ggplot(data = med.v.sd,aes(x = median, y = sd)) +
  geom_point() +
  geom_smooth(method = lm)
  
ggplot(data = med.v.sd,aes(x = median, y = sd)) +
  geom_point() +
  geom_abline(intercept = -1, slope = 1, col = 'blue') +
  geom_abline(intercept = -1.5, slope = 1.2, col = 'red') +
  geom_abline(intercept = -1.5, slope = 1.15, col = 'green')

########################################################################
### Highlight points beneath red line

med.v.sd$ratio <- med.v.sd$median/med.v.sd$sd
selected <- which((med.v.sd$median - 1.4)/med.v.sd$sd <= 1/1.15)
selected <- which((med.v.sd$median - 1.4)/med.v.sd$sd >= 1/1.15)

med.v.sd$selected <- FALSE
med.v.sd$selected[selected] <- TRUE

ggplot(data = med.v.sd,aes(x = median, y = sd)) +
  geom_point(aes(col = selected)) +
  geom_abline(intercept = -1.5, slope = 1.15, col = 'red')

### Filter for low SD genes

trans.df.filt <- trans.df[selected,]

### SCRATCH WORK
#trans.dist <- dist(med.v.sd$median)

plot(trans.dist)

hist(med.v.sd$median)
plot(hist(log2(med.v.sd$ratio),bi))

hist(log(med.v.sd$ratio,100),plot=TRUE,breaks=60)
