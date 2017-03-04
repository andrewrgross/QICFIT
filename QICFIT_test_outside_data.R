### QICFIT test using outside data -- Andrew R Gross -- 2017-03-01
### This script is used to load and format data from other sources, then to test various recognition methods


########################################################################
### Header



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

########################################################################
### Import

### Fullwood 2015, GSE69360
df.fullwood <- read.table("Z:/Data/Andrew/reference_data/geo/GSE69360/GSE69360_RNAseq.counts.txt", sep = '\t', header = TRUE, row.names = 1)

### Housekeeping genes
housekeeping.genes <- as.character(read.csv('Z:/Data/Andrew/QICFIT/housekeeping.ids_3559.csv')[,1])

########################################################################
### Format

df.fullwood <- df.fullwood[6:24]
df.fullwood <- convertIDs(df.fullwood)

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

column.num <- 5
reference.df <- ref.loose.pruned

spearman.results <- spearman.calc(df.fullwood[column.num], reference.df)
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