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
    result <- round(result[2], 5)                        # Round result
    spearman.results[tissue,] <- result                        # Add to results table
  }
  return(spearman.results)
}
spearman.calc.for.multiple.samples <- function(samples, reference.input) {
  ### Initialize empty results table
  spearman.results <- data.frame(rep(0,ncol(reference.input))) # Generate empty results table
  names(spearman.results) <- names(samples)[1]                 # Name empty results table
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
  new.col <- ncol(spearman.results)+1 
  spearman.results[new.col] <- row.names(spearman.results)                    # Add the tissues as a new row
  names(spearman.results)[new.col] <- 'ref'                                   # Label the new row
  spearman.t <- t(spearman.results[-new.col])                                 # Transpose the data frame, minus new row
  spearman.melt <- melt(spearman.t)                                     # Melt transposed data frame
  names(spearman.melt) <- c("sample","ref","value")                     # Rename columns of melted data frame
  spearman.melt$ref <- factor(spearman.melt$ref, levels = rev(row.names(spearman.results)), ordered = TRUE) # Assign order
  
  ### Generate data frame for color and label data
  spearman.for.color <- spearman.results                                # Duplicate spearman results
  spearman.for.color[-new.col] <- apply(spearman.for.color[-new.col], 1, mean)    # Reassign all values to mean of group
  names(spearman.for.color)[2] = 'color.val'
  spearman.t <- t(spearman.for.color[-new.col])                               # Transpose again
  spearman.melt.color <- melt(spearman.t)                               # Melt again
  spearman.melt$color <- spearman.melt.color$value                      # Copy repeating values as 'color column to main df
  
  g <- ggplot(data = spearman.melt, aes(x = ref, y = value)) +
    geom_boxplot(fill = 'red') +
    coord_flip() +
    scale_y_continuous(limits = c(0.49,0.76), position = 'right', breaks = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75), labels = c('0.5','','0.6','','0.7',''), expand = c(0,0)) +
    theme(axis.ticks.y = element_blank(),
          legend.position = 'none', panel.background = element_rect(fill = 'grey99',size = 2, linetype = 'solid', color = 'black'),
          plot.title = element_text(size = 14, face = 'bold', margin = margin(5,0,5,0), hjust = 0.5)) +
    labs(title = title,
         x = 'Reference Tissues from GTEx',
         y = 'Spearman correlation') #+ geom_text(data = spearman.for.color, aes(x = ref, y = color.val-0.25, label = ref), hjust = 0)
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
### Open selected dataset

setwd('z:/Data/Andrew/reference_data/qicfit_ready/')                    # Sets directory

### Import dataset and sample metadata
samples.metadata <- read.csv('METADATA_samples.csv')                    # Import sample metadata
dataset.metadata <- read.csv('METADATA_datasets.csv')                   # Import dataset metadata
print(dataset.metadata[c(1,2,4,5,6)])                                   # Print main dataset stats

#data.frame(list.files())                                               # Lists available datasets
selection.number <- 5                                                   # Specifies file of interest
(sample.file <- as.character(dataset.metadata[,4][selection.number]))   # Defines file of interest
samples.df <- read.csv(sample.file, row.names = 1)                       # Imports file of interest as current sample.df

### Report dataset and sample stats
accession <- as.character(dataset.metadata$GEO[selection.number])       # Define the accession number of current dataset
sample.rows <- grep(accession, samples.metadata$Accession)              # Identify corresponding rows in sample metadata
print(dataset.metadata[selection.number,][c(1,2,4,5,6,7,8,9,10,12)])    # print dataset info
print(samples.metadata[sample.rows[1],][c(8,7,10,11,12)])               # print dataset info from samples
current.samples.metadata <- samples.metadata[grep(accession, samples.metadata$Accession),][c(1,2,3,4)]
print(current.samples.metadata)

########################################################################
### Filter out low or absent genes


### Try alternate reference sets
reference.input <- ref.full 
reference.input <- ref.loose 
reference.input <- ref.neutral 
reference.input <- ref.tight 
reference.input <- ref.supertight 

### Filter samples by max row value
gene.maxes <- apply(samples.df,1, max)
samples.df.filtered <- samples.df[gene.maxes >= 0,]

### Filter samples by reference
samples.df.filtered <- samples.df.filtered[row.names(reference.input),]




### Filter reference input
reference.input.filtered <- reference.input[row.names(samples.df.filtered),]
nrow(reference.input.filtered)

### Calc Spearman coeff for filtered data
spearman.results <- spearman.calc.for.multiple.samples(samples.df.filtered, reference.input.filtered)

### Report each sample's top 3
for (sample.selection in 1:12) {print(spearman.results[sample.selection][order(spearman.results[sample.selection], decreasing = TRUE),, drop = FALSE][1:3,, drop = FALSE])}

########################################################################
### Generate Spearman profile for each sample

### Calculate Spearman coefficents
spearman.results <- spearman.calc.for.multiple.samples(samples.df, reference.input)

### Format results for heatmap
(spearman.for.plot <- spearman.results[c(1,2,3,4,5,6,7,8,10,12,14,16,18,22,26,30,35,40,50),])
spear.m <- melt(as.matrix(spearman.for.plot))
spear.m$X1 <- factor(spear.m$X1, levels = rev(as.character(row.names(spearman.for.plot))))

### Plot heatmap
ggplot(data = spear.m, aes(x = X2, y = X1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 0)), size = 3) +
  scale_fill_gradient(low = 'white', high = 'red') +
  ggtitle(as.character(dataset.metadata[selection.number,][,1])) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 1))


### Subset data by tissue type
(available.tissues <- levels(factor(current.samples.metadata$Tissue)))
tissue.selection <- 3
(selected.samples.metadata <- current.samples.metadata[current.samples.metadata$Tissue == as.character(available.tissues[tissue.selection]),])
selected.samples <- as.character(selected.samples.metadata$Sample)
spearman.for.plot <- spearman.results[selected.samples]
spearman.for.plot <- spearman.for.plot[rev(order(apply(spearman.for.plot,1,mean))),]

### Format results for heatmap
(spearman.for.plot <- spearman.for.plot[c(1,2,3,4,5,6,7,8,10,12,14,16,18,22,26,30,35,40,50),])
spear.m <- melt(as.matrix(spearman.for.plot))
spear.m$X1 <- factor(spear.m$X1, levels = rev(as.character(row.names(spearman.for.plot))))


### View specific tissues
sample.selection <- 11
spearman.results[sample.selection][order(spearman.results[sample.selection], decreasing = TRUE),, drop = FALSE][1:4,, drop = FALSE]

### Report each sample's top 4
for (sample.selection in 1:ncol(spearman.results)) {print(spearman.results[sample.selection][order(spearman.results[sample.selection], decreasing = TRUE),, drop = FALSE][1:4,, drop = FALSE])}

#spearman.results <- spearman.calc(sample.df[1],reference.input)




