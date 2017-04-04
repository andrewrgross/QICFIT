### Transcriptome Preprocessing -- Andrew R Gross -- 2016/12/02
### Make minor edits to formatting and begin identifying key genes

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

########################################################################
### Data input
### Import sample key
sample.key <- read.csv('Z:/Data/Andrew/reference_data/gtex/Sample.key.sorted.csv',header=TRUE)

### Import transcriptome location table
transcriptome.index <- read.csv("Z:/Data/Andrew/reference_data/gtex/transcriptome_tables_by_tissue/transcriptome.index.csv")
print(transcriptome.index[1])

### Reference names, formatted
referenceNames <- c("Adipose, subcutaneous","Adipose, omentum","Adrenal gland","Aorta","Coronary artery","Tibial artery","Bladder","Amygdala","Anteriror cingulate nucleous","Caudate nucleous",
                    "Cerebellar hemisphere","Cerebellum","Cortex","Frontal cortex BA9","Hippocampus","Hypothalamus","Nucleus accumbens","Putamen",
                    "Spinal cord","Substantia nigra","Mammary","Lymphocyte","Fibroblast","Ectocervix","Endocervix","Colon, sigmoid","Colon, transverse","Gastroesophageal junction","Esophagus, mucosa",
                    "Esophagus, muscularis","Fallopian tube","Heart, Atrial","Heart, left ventricle","Kidney","Liver","Lung","Salvitory gland","Skeletal muscle","Tibial nerve","Ovary","Pancreas",
                    "Pituitary","Prostate","Skin, sun-hidden","Skin, sun-exposed","Small intestine","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole blood")

########################################################################
### Prep dataframe

selection <- 2
trans.df <- read.delim(as.character(transcriptome.index$file.locations[selection]))  # Load selected dataset (may take >10s)

### Rename rows by ensembl IDs w/o decimal
row.names(trans.df) <- trans.df$Name
trans.df <- convertIDs(trans.df)

testing.df <- trans.df[2]
row.names(testing.df) <- row.names(trans.df)

########################################################################
### Loop through all

for(selection in 1:nrow(transcriptome.index)) {
  selected.tissue <- as.character(transcriptome.index[,1][selection])
  currentReferenceName <- referenceNames[selection]
  #print(paste(selected.tissue,currentReferenceName))
  
  ### Import transcriptome of specific tissue
  trans.df <- read.delim(as.character(transcriptome.index$file.locations[selection]))  # Load selected dataset (may take >10s)
  print(as.character(transcriptome.index$file.locations[selection]))
  print(currentReferenceName)
  
  ########################################################################
  ### Formatting
  ### Validate references are from correct tissue
  #sample.key.filtered <- sample.key[sample.key$Tissue.Specific == selected.tissue,]    # Filter the sample key to include only selected tissue
  #matched <- match(names(trans.df),sample.key.filtered$Sample_ID)  # Match reference column names to filtered sample key
  #matched <- matched[!is.na(matched)]                              # Remove non-matches.  There should be just two.
  #nrow(sample.key.filtered) == length(matched)                     # Report whether the sample names match up with the sample IDs for the selected tissue
  
  ### Rename rows by ensembl IDs w/o decimal
  #row.names(trans.df) <- trans.df$Name
  #trans.df <- convertIDs(trans.df)
  
  ### Pare down to first 8 columns
  trans.df <- trans.df[3:10]
  
  ### Round to three decimals
  trans.df <- round(trans.df,3)
  
  ### Rename Colums
  tissueName <- substr(currentReferenceName,1,15)
  tissueName <- rep(tissueName,8)
  tissueNames <- make.unique(tissueName, sep = '-')
  
  names(trans.df) <- tissueNames

  ########################################################################
  ### Add pared down columns to master dataframe
  
  testing.df <- cbind(testing.df, trans.df)
}

########################################################################
### Remove empty rows

rowSums <- rowSums(testing.df[2:ncol(testing.df)])
rowsAbove2 <- (rowSums >2)

testing.df.filtered <- testing.df[rowsAbove2,]

########################################################################
### Output result


write.csv(testing.df.filtered,'Z:/Data/Andrew/reference_data/gtex/individual-ref-transcriptomes-for-testing.csv')

#transcriptome.index2 <- transcriptome.index[1][-c(24,25,31),]
#transcriptome.index2 <- c(transcriptome.index2,'Cervix')
### correct names
tissueNamesFull <- c('Gene')

for(tissueName in referenceNames) {
  tissueName <- substr(tissueName,1,15)
  currentTissueNames <- rep(tissueName,10)
  currentTissueNames <- make.unique(currentTissueNames, sep = '-')
  tissueNamesFull <- c(tissueNamesFull,currentTissueNames)
  
}

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
