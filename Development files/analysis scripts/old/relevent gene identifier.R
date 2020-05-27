### Identify most relevant genes -- Andrew R Gross -- 2017-04-17
### I want to identify the genes with the greatest relevance to the low Spearman correlation, then screen out housekeeping genes

########################################################################
### Header
library(Hmisc)
library(biomaRt)
#ensembl = useMart(host="www.ensembl.org")
ensembl = useMart(host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL')
listDatasets(ensembl)
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
#attributes[grep('type', attributes[,1]),]
#filters[grep('ensembl', filters[,1]),]

########################################################################
### Functions
convertIDs <- function(dataframe) {
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs
  return(dataframe)
}

function(sample, reference.input) {           # Calculate the spearman correlation between the sample and the references
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

addGeneColumn <- function(dataframe) {
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

addType <- function(dataframe) {
  type <- getBM(attributes=c('ensembl_gene_id','gene_biotype'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  type <- type[match(row.names(dataframe),type[,1]),]
  Type <- c()
  for (rowNumber in 1:length(type[,1])) {
    newType <- type[rowNumber,][,2]
    Type <- c(Type, newType)
  }
  dataframe[length(dataframe)+1] <- as.factor(Type)
  names(dataframe)[ncol(dataframe)] <- "Type"
  return(dataframe)
}

addDesc <- function(dataframe) {
  desc <- getBM(attributes=c('ensembl_gene_id','description'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  desc <- desc[match(row.names(dataframe),desc[,1]),]
  Desc <- c()
  for (rowNumber in 1:length(desc[,1])) {
    newDesc <- desc[rowNumber,][,2]
    Desc <- c(Desc, newDesc)
  }
  dataframe[length(dataframe)+1] <- Desc
  names(dataframe)[ncol(dataframe)] <- "Description"
  return(dataframe)
}

########################################################################
### Import formatted data sets
### Fullwood 2015, GSE69360
df.fullwood <- read.table("Z:/Data/Andrew/reference_data/geo/GSE69360/GSE69360_RNAseq.counts.txt", sep = '\t', header = TRUE, row.names = 1)
### Housekeeping genes
housekeeping.genes <- as.character(read.csv('Z:/Data/Andrew/QICFIT/housekeeping.ids_3559.csv')[,1])
### 
df.yu <- read.csv('Z:/Data/Andrew/reference_data/geo/yu_GSE9440.csv', row.names = 1)
### 
df.han <- read.csv('Z:/Data/Andrew/reference_data/geo/han_GSE35108.csv', row.names = 1)
### Normalized
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

########################################################################
### Format

df.fullwood <- df.fullwood[6:24]
df.fullwood <- convertIDs(df.fullwood)

### Replace row and column names
sampleNames <- c("iHT_03iCTR","iHT_90iOBS","iHT_77iOBS","iHT_02iOBS","iMN_87iCTR","iMN_201iCTR","aHT_1662_S6","aHT_1838_S7","aHT_1843_S8","aHT_2266_S9","iHT_02iCTR_S13","aHT_2884_S11","21-Reference","iHT_87iCTR","iHT_201iCTR","iHT_25iCTR_S16","iHT_688iCTR_","iHT_80iCTR","iHT_74iOBS","iHT_03iOBS")
names(TPMdata) <- sampleNames
TPMdata <- convertIDs(TPMdata)
### Remove 02iOBS and both iMN
TPMdata[c("iHT_02iOBS","iMN_87iCTR","iMN_201iCTR","21-Reference")] <- NULL
sampleNames <- c("iHT_03iCTR","iHT_90iOBS","iHT_77iOBS","aHT_1662","aHT_1838","aHT_1843","aHT_2266","iHT_02iCTR","aHT_2884","iHT_87iCTR","iHT_201iCTR","iHT_25iCTR","iHT_688iCTR","iHT_80iCTR","iHT_74iOBS","iHT_03iOBS")

########################################################################
### Run test using Spearman
column.num <- 1
sample.data <- TPMdata[column.num]
sample.data <- df.fullwood[column.num]
sample.data <- df.yu[column.num]
sample.data <- df.han[column.num]
reference.df <- ref.full

#spearman.results <- spearman.calc(sample.data, reference.df)
sample <- sample.data
reference.input <- reference.df                              # Calculate the spearman correlation between the sample and the references
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
spearman.results <- spearman.results[order(spearman.results[,1], decreasing = TRUE),,drop = FALSE]
top.reference.hit <- row.names(spearman.results)[1]

### Generate a dataframe of genes and their rank difference between two lists
ref.tissue.data <- reference.input[top.reference.hit]         # Call the current tissue from the references
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

### Generate rank table
sample.rank <- rank(-spearman.input[1])
reference.rank <- rank(-spearman.input[2])
rank.diff <- (reference.rank - sample.rank)^2
rank.table <- cbind(spearman.input, sample.rank, reference.rank, rank.diff)
rank.table <- rank.table[order(rank.diff, decreasing = TRUE),]
rank.table$rank.diff.div <- rank.table$rank.diff/rank.table$reference.rank
rank.table <- rank.table[order(rank.table$rank.diff.div, decreasing = TRUE),]

### Add name, type, and description
gene.rank.results <- rank.table[c(1,4,6)][1:1000,]
gene.rank.results <- addGeneColumn(gene.rank.results)
gene.rank.results <- addType(gene.rank.results)
gene.rank.results <- gene.rank.results[!is.na(gene.rank.results$Type),]

### report results
summary(gene.rank.results,16)

### Filter results
protein.coding.genes <- gene.rank.results$Type == 'protein_coding'
protein.coding.results <- gene.rank.results[protein.coding.genes,]
protein.coding.results <- addDescription(protein.coding.results[-5])
new.desc <- c()
for(description in protein.coding.results$Description) {new.desc <- c(new.desc, strsplit(description, ' [', fixed = TRUE)[[1]][1])}
protein.coding.results$Description <- new.desc
head(protein.coding.results,30)


### Estimate Spearman value
n <- nrow(rank.table)
1 - (6*sum(rank.diff))/(n*(n^2-1))


sample.rank <- cbind(sample.data, rank(-sample.data, ties.method = 'average'))
closest.reference <- reference.df[16]
reference.rank <- cbind(closest.reference, rank(-closest.reference, ties.method = 'average'))


title <- names(spearman.results)[column.num]
#spearman.results <- spearman.results[1]
names(spearman.results) <- 'value'
spearman.results$ref <- row.names(spearman.results)

spearman.results <- spearman.results[order(spearman.results$value, decreasing = TRUE),]
spearman.results$ref <- factor(spearman.results$ref, levels = rev(spearman.results$ref))

### Plot
g <- ggplot(data = spearman.results, aes(x = ref, y = value, fill = value)) +
  geom_bar(stat = 'identity') +
  scale_fill_gradientn(colors = c('white','yellow','orange','red')) +
  coord_flip() +
  scale_y_continuous(position = 'right', limits = c(0,100)) +
  labs(title = row.names(spearman.results)[1],
       x = '',
       y = title)
g
