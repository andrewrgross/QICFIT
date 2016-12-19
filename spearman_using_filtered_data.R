### Spearman comparison using filtered data -- Andrew R Gross -- 2016/12/16
### 

########################################################################
### Header
library(ggplot2)

########################################################################
### Functions


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



aHT <- TPMdata[c(7,8,9,10,12)]

### Combine samples with references

shared.rows <- intersect(row.names(references.supertight),row.names(aHT))
spearman.input <- cbind(aHT[shared.rows,],references.supertight[shared.rows,])


title <- 'references.supertight.pruned'
#spearman.input <- references.full3
dim(spearman.input)
########################################################################
### Spearman

#spearman.results <- cor(spearman.input, method = 'spearman', use = 'complete.obs')
spearman.results <- rcorr(as.matrix(spearman.input), type = 'spearman')[[1]]
spearman.results <- round(spearman.results*100, 0)

spearman.results.trimmed <- spearman.results[,1:5]
order.spearman.mean <- apply(spearman.results.trimmed, 1, mean)
spearman.results.trimmed[order(order.spearman.mean, decreasing = TRUE),]


########################################################################
### Generate heatmap

heatmap.2(spearman.results,
          main = title, # heat map title
          cellnote = spearman.results,
          notecol = "gray40",
          density.info="none",  # turns off density plot inside color legend
          notecex=0.7,
          trace="none",         # turns off trace lines inside the heat map
          distfun=dist,
          margins =c(12,12),     # widens margins around plot
          dendrogram="none"     # only draw a row dendrogram
)






########################################################################
########################################################################
########################################################################
### Scratchwork

references.loose.pruned <- references.loose[!apply(references.loose, 1, function(x){all(is.na(x))}),]
references.neutral.pruned <- references.neutral[!apply(references.neutral, 1, function(x){all(is.na(x))}),]
references.tight.pruned <- references.tight[!apply(references.tight, 1, function(x){all(is.na(x))}),]
references.supertight.pruned <- references.supertight[!apply(references.loose, 1, function(x){all(is.na(x))}),]