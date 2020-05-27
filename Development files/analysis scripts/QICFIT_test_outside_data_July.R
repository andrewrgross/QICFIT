### optimize analysis of fullwood -- Andrew R Gross -- 2017-06-05
### This script is used to compare Fullwood sample transcriptomes directly to GTEx transcriptomes in order to improve comparison methods.

########################################################################
### Header

library(Hmisc)
library(biomaRt)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(QICFIT, lib.loc = "C:/Users/grossar/Documents/GitHub/QICFIT/package 3")

########################################################################
### Functions

########################################################################
### Import Data

### Import dataset and sample metadata -- Dataset metadata contains filepath for all available samples
samples.metadata <- read.csv('Z:/Data - DEPRECIATED/Andrew/reference_data/qicfit_ready/METADATA_samples.csv')                    # Import sample metadata
dataset.metadata <- read.csv('Z:/Data - DEPRECIATED/Andrew/reference_data/qicfit_ready/METADATA_datasets.csv')                   # Import dataset metadata



### Import references  -- The average expression of each tissue without filtering
ref.df <- read.csv('Z:/Data - DEPRECIATED/Andrew/reference_data/gtex/sd.filtered.tables/gtex.full.csv',header=TRUE,row.names=1)[2:54]

### Import individual reference transcriptomes for testing (sets of 8 from GTEx)
gtex.transcriptomes.df <- read.csv('Z:/Data - DEPRECIATED/Andrew/reference_data/gtex/individual-ref-transcriptomes-for-testing.csv', row.names = 1)
individual.transcriptomes.df <- gtex.transcriptomes.df

########################################################################
### Select and import data from list of available datasets

### Print selection table
print(dataset.metadata[c(1,2,5,6,7,8,9,10)])

### Select data set
dataset.num <- 9              ###---+++---+++---### USER INPUT GOES HERE ###---+++---+++---###

### Load selected dataset
user.dataset <- read.csv(paste0('Z:/Data - DEPRECIATED/Andrew/reference_data/qicfit_ready/', dataset.metadata$file_corrected[dataset.num]), row.names = 1)
print(head(user.dataset,4))

########################################################################
### Formatting

### ???
samples.metadata <- samples.metadata[grep('Fullwood', samples.metadata$Author),]
(available.tissues <- levels(factor(samples.metadata$Tissue)))
print(dataset.metadata[c(1,2,4,5,6)])

### Select samples from within datasets

### Subset Fullwood data by tissue
fullwood.df <- read.csv('Z:/Data - DEPRECIATED/Andrew/reference_data/qicfit_ready/GSE69360-fullwood_ID_CORR.csv', row.names = 1)
fw.colon.df <- fullwood.df[grep('colon', samples.metadata$Tissue)]
fw.heart.df <- fullwood.df[grep('heart', samples.metadata$Tissue)]
fw.kidney.df <- fullwood.df[grep('kidney', samples.metadata$Tissue)]
fw.liver.df <- fullwood.df[grep('liver', samples.metadata$Tissue)]
fw.lung.df <- fullwood.df[grep('lung', samples.metadata$Tissue)]
fw.stomach.df <- fullwood.df[grep('stomach', samples.metadata$Tissue)]

### Subset GTEx samples by tissue
gtex.colon.s.df <- individual.transcriptomes.df[186:193]
gtex.colon.t.df <- individual.transcriptomes.df[194:201]
gtex.heart.a.df <- individual.transcriptomes.df[226:233]
gtex.heart.b.df <- individual.transcriptomes.df[234:241]
gtex.kidney.df <- individual.transcriptomes.df[242:249]
gtex.liver.df <- individual.transcriptomes.df[250:257]
gtex.lung.df <- individual.transcriptomes.df[258:265]
gtex.stomach.df <- individual.transcriptomes.df[354:361]

########################################################################
### Generate df lists to select from

query.list <- list(gtex.colon.s.df, gtex.colon.t.df, gtex.heart.a.df, gtex.heart.b.df, gtex.kidney.df, gtex.liver.df, gtex.lung.df, gtex.stomach.df, fw.colon.df, fw.heart.df, fw.kidney.df, fw.liver.df, fw.lung.df, fw.stomach.df)
names(query.list) <- c('gtex.colon.sig', 'gtex.colon.tran','gtex.heart a', 'gtex.heart b', 'gtex.kidney', 'gtex.liver', 'gtex.lung', 'gtex.stomach', 'fw.colon', 'fw.heart', 'fw.kidney', 'fw.liver', 'fw.lung', 'fw.stomach')

reference.input <- ref.df

########################################################################
### Generate Spearman profile for each sample

### Select transcriptome to analyze
data.frame(names(query.list))
tissue.selection <- 1
samples.df <- query.list[[tissue.selection]]; tissue.name <- names(query.list)[tissue.selection]
print(names(samples.df)[1])

### Calculate Spearman coefficents
#spearman.results <- spearman.calc.for.multiple.samples(samples.df, reference.input)
spearman.results <- spearman.calc(samples.df, reference.input, method = 'corr')
spearman.results <- spearman.calc(samples.df, reference.input, method = 'weighted')

### Plot function

plot <- plot.spearman(spearman.results)
########################################################################
### New plot type: heat map of individual samples
plot.list <- list()
for (column in 1:ncol(spearman.results)) {

spearman.for.plot <- spearman.results[1:20,][column]
spearman.for.plot <- spearman.for.plot[order(spearman.for.plot[1], decreasing = TRUE),, drop = FALSE]
spear.m <- melt(as.matrix(spearman.for.plot))
spear.m$Var1 <- factor(spear.m$Var1, levels = rev(as.character(row.names(spearman.for.plot))))
spear.m$label <- paste0(spear.m$value, '-', spear.m$Var1)

## Plot heatmap
plot <- ggplot(data = spear.m, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(low = 'white', high = 'red') +
  #ggtitle(as.character(dataset.metadata[selection.number,][,1])) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position="none")
plot.list[column] <- list(plot)
}

### Convert plot list to list of graphic objects (grobs)
grob.list <- lapply(plot.list, ggplotGrob)

### Plot grob.list
(plot.to.save <- grid.arrange(grobs = grob.list, ncol = ncol(spearman.results), top = tissue.name[1]))

########################################################################
### Calculate and plot using QICFIT package

spearman.results <- spearman.calc(samples.df, reference.df)



### Save plots
setwd("z:/Data/Andrew/QICFIT/spearman testing on Fullwood-June/")
png(filename=paste0(tissue.name,'_plot.png'),
    type="cairo",
    units="in",
    width=6,
    height=6,
    pointsize=12,
    res=120)
print(plot)
dev.off()

ggsave('test.png')

#grid.arrange(a, a, ncol=2, top = "Main title")

setwd("z:/Data/Andrew/QICFIT/spearman testing on Fullwood-June/")

for (tissue.selection in 1:length(query.list)) {
  samples.df <- query.list[[tissue.selection]]; tissue.name <- names(query.list)[tissue.selection]
  print(names(samples.df)[1])

  ### Calculate Spearman coefficents
  spearman.results <- spearman.calc(samples.df, reference.input, method = 'corr')
  spearman.results <- spearman.calc(samples.df, reference.input, method = 'weighted')

  ### Plot function
  plot <- plot.spearman(spearman.results)
  ggsave(file = paste0(tissue.name, '.png'), plot = plot)
}



########################################################################
########################################################################
########################################################################
### Rank-weighted Spearman calculation
### This script should perform a spearman calculation, but allow results to be weighted by their rank.  It should also fascilitate
### the process of identifying the genes responsible for diminished values

tissue.selection <- 2

### Select Fullwood tissue to analyze
samples.df <- fw.samples.list[[tissue.selection]]; tissue.name <- names(fw.samples.list)[tissue.selection]
#samples.df <- gtex.samples.list[[tissue.selection]]; tissue.name <- names(gtex.samples.list)[tissue.selection]
print(names(samples.df))


### Calculate Spearman coefficents
#spearman.results <- spearman.calc.for.multiple.samples(samples.df, reference.input)
#spearman.results <- spearman.calc.for.multiple.samples.manual(samples.df, reference.input)
spearman.results <- spearman.calc.for.multiple.samples.manual.weighted(samples.df, reference.input)

########################################################################
### New plot type: heat map of individual samples
plot.list <- list()
for (column in 1:ncol(spearman.results)) {
  spearman.for.plot <- spearman.results[1:20,][column]
  spearman.for.plot <- spearman.for.plot[order(spearman.for.plot[1], decreasing = TRUE),, drop = FALSE]
  spear.m <- melt(as.matrix(spearman.for.plot))
  spear.m$Var1 <- factor(spear.m$Var1, levels = rev(as.character(row.names(spearman.for.plot))))
  spear.m$label <- paste0(spear.m$value, '-', spear.m$Var1)
  ## Plot heatmap
  plot <- ggplot(data = spear.m, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient(low = 'white', high = 'red') +
    #ggtitle(as.character(dataset.metadata[selection.number,][,1])) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position="none")
  plot.list[column] <- list(plot)
}

### Convert plot list to list of graphic objects (grobs)
grob.list <- lapply(plot.list, ggplotGrob)
### Plot grob.list
(plot.to.save <- grid.arrange(grobs = grob.list, ncol = ncol(spearman.results), top = tissue.name[1]))







########################################################################
########################################################################

### Format results for heatmap
(spearman.for.plot <- spearman.results[c(1,2,3,4,5,6,7,8,10,12,14,16,18,22,26,30,35,40,50),])
spear.m <- melt(as.matrix(spearman.for.plot))
spear.m$Var1 <- factor(spear.m$Var1, levels = rev(as.character(row.names(spearman.for.plot))))

### Plot heatmap
ggplot(data = spear.m, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), size = 3) +
  scale_fill_gradient(low = 'white', high = 'red') +
  #ggtitle(as.character(dataset.metadata[selection.number,][,1])) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1))




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
