### Running QICFIT on E099 RNA-seq data   2018-11-30    Andrew R Gross

### Scratch work

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/Motor Neurons/E099 - PM - ALS v CTRL/')
counts.als <- read.csv('PM-5119--07--02--2018_COUNTS.csv', row.names = 1)
tpm.als <- read.csv('PM-5119--07--02--2018_TPM.csv', row.names = 1)
sample.names.als <- c('02iCTR','03iCTR','159iALS','172iCTR','372iALS_n1','372iALS_n2','372iALS_n3','395iCTR')

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/Motor Neurons/E0x - PP - D10 KO v CTRL/')
counts.ko <- read.csv('COUNTS.csv', row.names = 1)
sample.names.ko <- read.csv('Sample names.txt', header = FALSE)

### Reassign names
names(counts.als) <- sample.names.als
names(counts.ko) <- sample.names.ko[,1]



convert.ids <- function(dataframe, add.gene.name.column = TRUE) {
  ### This function will convert a row name consisting of a contactenated ensembl ID and gene to one or the other,
  ### based on the users instruction (2018-10-04)
  ensemblIDs <- c()                                           # Empty lists are initialized to receive IDs as they're created
  gene.names <- c()
  for (rowName in row.names(dataframe)) {                     # Loops through all rows in the data frame
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]                 # Splits the row name and declares the ensembl ID
    gene.name <- strsplit(rowName,"\\_")[[1]][2]                 # Splits the row name, declares the gene name
    ensemblIDs <- c(ensemblIDs, ensemblID)                       # Adds ensembl ID and gene name to appropriate lists
    gene.names <- c(gene.names, gene.name)
  }
  row.names(dataframe) <- ensemblIDs                          # assigns the new row names
  if(add.gene.name.column == TRUE) {
    dataframe$Gene <- gene.names
  }
  return(dataframe)                                           # Returns the data frame with new rows
}


counts.als <- convert.ids(counts.als)
counts.als <- counts.als[1:8]

counts.ko <- convert.ids(counts.ko)
counts.ko <- counts.ko[1:6]

samples.df <- counts.ko
