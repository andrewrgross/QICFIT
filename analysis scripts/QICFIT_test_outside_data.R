### QICFIT test using outside data -- Andrew R Gross -- 2017-03-01
### This script is used to load and format data from other sources, then to test various recognition methods
### This script reads in expression tables; converts the IDs to ENSEMBL; Loops through each transcriptome and does a Spearman check against references; and outputs figures and tables reporting
### The similarity of each query to each of the GTEx references

########################################################################
### Header

library(Hmisc)
library(biomaRt)
#listMarts(host="www.ensembl.org")
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL')
#listDatasets(ensembl)
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
filters[grep('uniprot', filters[,1]),]
attributes[grep('pattern', attributes[,1])]

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

uniprot_gn.to.gene.name <- function(dataframe) {
  ### Generate initial conversion table
  conversion.df <- getBM(attributes=c('uniprot_gn','external_gene_name'), filters='uniprot_gn', values= row.names(dataframe), mart=ensembl)
  unique.ids <- unique(conversion.df$uniprot_gn)
  ### Identify the redundancy of each ID
  redundancy.df <- data.frame(unique.ids, 'no.of.results' = rep('',length(unique.ids)), stringsAsFactors = FALSE)
  for (row.num in 1:nrow(redundancy.df)) {
    id = redundancy.df[row.num,][,1]
    hits <- grep(id,conversion.df[,1])
    redundancy.df[row.num,][2] <- length(hits)
  }
  ### Generate a df containing single-result IDs and multiple-result IDs
  single.result.df <- data.frame('uniprot_gn' = c(), 'external_gene_name' = c(), stringsAsFactors = FALSE)
  multiple.result.df <- data.frame('uniprot_gn' = c(), 'external_gene_name' = c(), stringsAsFactors = FALSE)
  
  for (row.num in 1:nrow(redundancy.df)) {
    if (redundancy.df$no.of.results[row.num] == 1) {
      row.to.add <- conversion.df[grep(redundancy.df$unique.ids[row.num], conversion.df$uniprot_gn),]
      single.result.df <- rbind(single.result.df,row.to.add)
    }
    if (redundancy.df$no.of.results[row.num] > 1) {
      all.gene.names <- conversion.df$external_gene_name[grep(redundancy.df$unique.ids[row.num], conversion.df$uniprot_gn)]
      all.gene.names <- paste(all.gene.names, collapse = '-')
      row.to.add <- data.frame('uniprot_gn' = redundancy.df$unique.ids[row.num], 'external_gene_name' = all.gene.names)
      multiple.result.df <- rbind(multiple.result.df,row.to.add)    
    }
  }
  ### Generate a df containing missing IDs
  missing.ids <- setdiff(row.names(dataframe), conversion.df$uniprot_gn)
  missing.result.df <- data.frame('uniprot_gn' = missing.ids, 'external_gene_name' = missing.ids)
  ### Join dfs and reorder
  conversion.df <- rbind(single.result.df, multiple.result.df, missing.result.df)
  conversion.df <- conversion.df[match(row.names(dataframe),conversion.df$uniprot_gn),]
  ### Convert old IDs to new IDs
  converted.df <- dataframe
  row.names(converted.df) <- conversion.df$external_gene_name
  return(converted.df)
}

drop.decimals <- function(input.dataframe) {
  input.id <- input.dataframe[,1]
  temp.id <- strsplit(input.id, '.', fixed = TRUE)
  output.id <- c()
  for (id in temp.id) {output.id <- c(output.id, id[1])}
  row.names(input.dataframe) <- output.id
  input.dataframe <- input.dataframe[2:ncol(input.dataframe)]
  return(input.dataframe)
}
########################################################################
### Import formatted data sets, references

### Set working directory
setwd("Z:/Data/Andrew/reference_data/qicfit_ready/")

### Acquire list of available files
metadata.df <- read.csv("METADATA.csv", stringsAsFactors = FALSE)
(available.transcriptomes <- list.files("files with their original ids/"))

### SELECT INPUT FILE
selection.number <- 12
(input.file <- available.transcriptomes[selection.number])
#metadata.df[metadata.df$Author == "Uosaki",]

### Import from available file list
reference.data.df <- read.csv(paste0('files with their original ids/',input.file), stringsAsFactors = FALSE)

### Import full references  -- The average expression of each tissue without filtering
references <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.full.csv',header=TRUE,row.names=1)

########################################################################
### Convert IDs to ensembl

### Drop Decimal on ENSEMBL IDs
#reference.data.df.converted.full <- drop.decimals(reference.data.df)

### Identify correct input filter
names(reference.data.df)[1]
head(reference.data.df[1:4])
#filters[grep('u133', filters[,1]),]
#filters[grep('hugene', filters[,1]),]
length(filters[grep('gene', filters[,1]),])
filters[grep('illumina', filters[,1]),]
#attributes[grep('ucsc', attributes[,1]),]

### Declare the input, output IDs, and IDs to convert
input.id <- 'affy_hg_u133_plus_2'
#input.id <- 'affy_hugene_1_0_st_v1'
#input.id <- 'with_ucsc'
input.id <- 'wikigene_name'
#input.id <- 'entrezgene'
output.id <- 'ensembl_gene_id'
current.ids <- reference.data.df[,1][1:30]

### Generate initial conversion table
conversion.df <- getBM(attributes=c(input.id, output.id), filters=input.id, values= current.ids, mart=ensembl) # ~60 s to run
unique.ids <- unique(conversion.df[,1])

### Report stats on initial conversion table
head(unique.ids,30)
print(paste(round(length(unique.ids)/length(current.ids),2)*100, 'percent of IDs recognized'))

### Match data to convesion table and join
reference.data.df.reordered <- reference.data.df[match(conversion.df[,1],reference.data.df[,1]),]
reference.data.df.reordered.appended <- cbind(conversion.df,reference.data.df.reordered)

### Define list of ENsEMBL IDs and empty data frame to fill
#reference.data.df.reordered.appended <- reference.data.df.reordered.appended[1:10,]
ensembl.col <- reference.data.df.reordered.appended[,2]
unique.ensembl.ids <- unique(reference.data.df.reordered.appended[,2])
single.matches.df <- reference.data.df.reordered.appended[1,]

### Copy a unique row for each ENSEMBL ID to a new data frame ### ~5 minutes
for (id in unique.ensembl.ids) {
  matches <- grep(id, ensembl.col)
  temp.df <- reference.data.df.reordered.appended[matches,]
  new.row <- apply(temp.df, 2, max)
  single.matches.df <- rbind(single.matches.df, new.row)
}

### Rename rows based on ENSEMBL IDs; remove defunct rows/columns
single.matches.df <- single.matches.df[2:nrow(single.matches.df),]
row.names(single.matches.df) <- single.matches.df[,2]
reference.data.df.converted <- single.matches.df[4:ncol(single.matches.df)]

### Generate a list of missing ENSEMBL IDs
missing.ids <- setdiff(row.names(references),row.names(reference.data.df.converted))
absent.data.df <- data.frame(matrix(0,length(missing.ids), ncol(reference.data.df.converted)) )
names(absent.data.df) <- names(reference.data.df.converted)
row.names(absent.data.df) <- missing.ids

### Join empty rows to query data with converted ENSEMBL IDs
reference.data.df.converted.full <- rbind(reference.data.df.converted, absent.data.df)
reference.data.df.converted.full <- reference.data.df.converted.full[sort(row.names(reference.data.df.converted.full)),]


head(reference.data.df.converted.full,25)

########################################################################
### Write file with correct IDs

getwd()

(new.file.name <- paste0(substr(available.transcriptomes[selection.number], 0, nchar(available.transcriptomes[selection.number])-4), '_ID_CORR.csv'))

write.csv(reference.data.df.converted.full, new.file.name)

test <- read.csv(new.file.name, row.names = 1)





########################################################################
########################################################################

### Fullwood 2015, GSE69360
df.fullwood <- read.table("Z:/Data/Andrew/reference_data/geo/GSE69360/GSE69360_RNAseq.counts.txt", sep = '\t', header = TRUE, row.names = 1)

### Housekeeping genes
housekeeping.genes <- as.character(read.csv('Z:/Data/Andrew/QICFIT/housekeeping.ids_3559.csv')[,1])

df.fullwood

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


### Remove Housekeeping genes from references
test <- setdiff(row.names(ref.full),housekeeping.genes)
test2 <- ref.full[test,]

ref.full.pruned <- ref.full[setdiff(row.names(ref.full), housekeeping.genes),]
ref.loose.pruned <- ref.loose[setdiff(row.names(ref.loose), housekeeping.genes),]
ref.neutral.pruned <- ref.neutral[setdiff(row.names(ref.neutral), housekeeping.genes),]
ref.tight.pruned <- ref.tight[setdiff(row.names(ref.tight), housekeeping.genes),]
ref.supertight.pruned <- ref.supertight[setdiff(row.names(ref.supertight), housekeeping.genes),]

### Generate a list of filtered reference sets
references.list <- list(ref.full, ref.full.pruned, ref.loose, ref.loose.pruned, ref.neutral, ref.neutral.pruned, ref.tight, ref.tight.pruned, ref.supertight, ref.supertight.pruned)
names(references.list) <- c('Ref.full', 'ref.full.pruned', 'Ref.loose', 'ref.loose.pruned', 'Ref.neutral', 'ref.neutral.pruned', 'Ref. tight', 'ref.tight.pruned', 'Ref.supertight', 'ref.supertight.pruned')

########################################################################
### Run test using Spearman

column.num <- 4
sample.data <- TPMdata[column.num]
sample.data <- df.fullwood[column.num]
sample.data <- df.yu[column.num]
sample.data <- df.han[column.num]
reference.df <- ref.full
sample.data <- test[column.num]

spearman.results <- spearman.calc(sample.data, reference.df)
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


########################################################################
### Loop through all five reference cutoffs

column.num <- 11
### Initialize empty results lists
spearman.results.list <- list()
### Loop through reference sets
for(ref.set.num in 1:length(references.list)) {                # Loop through all five reference filter levels
  reference.input <- references.list[[ref.set.num]]            # Specify current reference list and its name
  (reference.set.name <- names(references.list)[ref.set.num])  # Declare the name of the rererence filter set
  spearman.results <- spearman.calc(df.fullwood[column.num], reference.input)
  title <- names(spearman.results)[column.num]
  names(spearman.results) <- 'value'
  spearman.results$ref <- row.names(spearman.results)
  ### Define order
  spearman.results <- spearman.results[order(spearman.results$value, decreasing = TRUE),]
  spearman.results$ref <- factor(spearman.results$ref, levels = rev(spearman.results$ref))
  ### Add spearman results to list
  spearman.results.list[[ref.set.num]] <- spearman.results     # Add current Spearman results table of 8 samples compared to all tissues to a list
  names(spearman.results.list)[ref.set.num] <- reference.set.name # Name the new entry on the list with the tissue being assessed
}
### Generate multiplot
barplot.list <- list()
filter.levels <- names(references.list)
level.num = 1
for (spearman.results in spearman.results.list) {
  filter.level = filter.levels[level.num]
  level.num <- level.num + 1
  new.plot  <- ggplot(data = spearman.results, aes(x = ref, y = value, fill = value)) +
    geom_bar(stat = 'identity') +
    scale_fill_gradientn(colors = c('white','yellow','orange','red')) +
    coord_flip() +
    scale_y_continuous(position = 'right') +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 7),axis.ticks.y = element_blank(),
          legend.position = 'none', panel.background = element_rect(fill = 'grey97'),
          plot.title = element_text(size = 14, face = 'bold', margin = margin(5,0,5,0), hjust = 0.5)) +
    labs(title = row.names(spearman.results)[1],
         x = '',
         y = filter.level) 
  barplot.list[[length(barplot.list)+1]] <- new.plot
}
print('All figures complete')
multibarplot <- plot_grid(barplot.list[[1]], barplot.list[[2]], barplot.list[[3]], barplot.list[[4]], barplot.list[[5]], labels=c("F", "L", "N", "T", "ST"), ncol = 5, nrow = 1)

multibarplot <- plot_grid(barplot.list[[1]], barplot.list[[2]], barplot.list[[3]], barplot.list[[4]], barplot.list[[5]], 
                          barplot.list[[6]], barplot.list[[7]], barplot.list[[8]], barplot.list[[9]], barplot.list[[10]],
                          ncol = 4, nrow = 2)
multibarplot <- plot_grid(barplot.list[[9]], barplot.list[[10]], ncol = 2, nrow = 1)


multibarplot
names(df.fullwood)[column.num]


### Breaking down Spearman correlation

sample <- df.fullwood[1]
reference.input <- ref.full

### Generate & format empty results df
spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
row.names(spearman.results) <- names(reference.input)        # Name empty results table
names(spearman.results) <- names(sample)

### Loop through each reference  tissue
for(ref.tissue.num in 1:ncol(reference.input)) {   
  #ref.tissue.num <- 1
  ref.tissue.data <- reference.input[ref.tissue.num]         # Call the current tissue from the references
  tissue <- names(ref.tissue.data)                           # Declare the name of the current tissue
  
  ### Add missing genes and assign value of zero
  ref.tissue.data.2 <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE] # Filter out missing values from tissue
  genes.present.in.ref <- row.names(ref.tissue.data.2)         # Declare the genes present in the reference
  genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(sample)) # Declare genes in reference missing from sample
  rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))  # Generate a zero data frame the size of the missing rows 
  row.names(rows.to.add) <- genes.missing.in.query           # Name rows after missing rows
  names(rows.to.add) <- names(sample)                        # Name column the same as the query 
  sample.2 <- rbind(sample,rows.to.add)                        # Use rbind to make a full data frame containing exactly the genes in the reference
  sample.3 <- sample.2[genes.present.in.ref, , drop = FALSE]     # Reorder sample to match reference
  spearman.input <- cbind(sample.3, ref.tissue.data.2)           # Bind sample and reference
  
  ### Perform Spearman calculation
  result <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]] # Perform spearman calculation
  result <- round(result[2] * 100, 1)                        # Round result
  spearman.results[tissue,] <- result                        # Add to results table
}

test2 <- cor(spearman.input, use="complete.obs", method = "spearman") 















spearman.calc <- function(sample, reference.input) {           # Calculate the spearman correlation between the sample and the references

  return(spearman.results)
}