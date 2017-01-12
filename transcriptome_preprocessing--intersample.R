### Transcriptome Preprocessing -- Andrew R Gross -- 2016/12/02
### Rank the consistency of gene expression levels between samples

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
gtex.full <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.full.csv',header=TRUE,row.names=1)
gtex.low.sd.loose <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.loose.csv', row.names = 1)
gtex.low.sd.neutral <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.neutral.csv', row.names = 1)
gtex.low.sd.tight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.tight.csv', row.names = 1)
gtex.low.sd.supertight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.supertight.csv', row.names = 1)

########################################################################
### Format

ref.full <- gtex.full[2:ncol(gtex.full)]
ref.loose <- gtex.low.sd.loose[2:ncol(gtex.low.sd.loose)]
ref.neutral <- gtex.low.sd.neutral[2:ncol(gtex.low.sd.neutral)]
ref.tight <- gtex.low.sd.tight[2:ncol(gtex.low.sd.tight)]
ref.supertight <- gtex.low.sd.supertight[2:ncol(gtex.low.sd.supertight)]

### Remove NA rows
ref.loose.pruned <- ref.loose[!apply(ref.loose, 1, function(x){all(is.na(x))}),]
ref.neutral.pruned <- ref.neutral[!apply(ref.neutral, 1, function(x){all(is.na(x))}),]
ref.tight.pruned <- ref.tight[!apply(ref.tight, 1, function(x){all(is.na(x))}),]
ref.supertight.pruned <- ref.supertight[!apply(ref.supertight, 1, function(x){all(is.na(x))}),]

ref.supertight.pruned.2 <- ref.supertight[!apply(ref.supertight, 1, function(x){any(is.na(x))}),]
temp <- 10^ref.supertight.pruned.2

########################################################################
### Calculate variance for each gene

trans.df <- ref.supertight.pruned

trans.df <- addMedSD(trans.df[3:ncol(trans.df)])                # Make new DF with median & sd


### Get list of housekeeping IDs
ids.of.interest <- row.names(ref.supertight.pruned.2)

genes.full <- gtex.low.sd.supertight[1]
genes.of.interest <- genes.full[ids.of.interest, , drop = FALSE]

row.names(temp) <- genes.of.interest[,1]

temp <- addMedSD(temp)
temp <- sortByMed(temp)

### Import transcriptome location table
transcriptome.index <- read.csv("Z:/Data/Andrew/reference_data/gtex/transcriptome_tables_by_tissue/transcriptome.index.csv")
print(transcriptome.index[1])
selection <- 2
selected.tissue <- as.character(transcriptome.index[,1][selection])
print(selected.tissue)

### Import transcriptome of specific tissue
trans.df <- read.delim(as.character(transcriptome.index$file.locations[selection]))  # Load selected dataset (may take >10s)

########################################################################
### Formatting
### Validate references are from correct tissue
sample.key.filtered <- sample.key[sample.key$Tissue.Specific == selected.tissue,]    # Filter the sample key to include only selected tissue
matched <- match(names(trans.df),sample.key.filtered$Sample_ID)  # Match reference column names to filtered sample key
matched <- matched[!is.na(matched)]                              # Remove non-matches.  There should be just two.
nrow(sample.key.filtered) == length(matched)                     # Report whether the sample names match up with the sample IDs for the selected tissue

### Rename rows by ensembl IDs w/o decimal
row.names(trans.df) <- trans.df$Name
trans.df <- convertIDs(trans.df)

########################################################################
### Add med and SD
trans.df <- addMedSD(trans.df[3:ncol(trans.df)])                # Make new DF with median & sd

########################################################################
### Plot, unfiltered
ggplot(data = log10(trans.df),aes(x = median, y = sd)) + geom_point(size = 0.5) + theme_light() +geom_smooth(method = lm)

########################################################################
### Filter by expression value
### Filter rows with any zeros
zero.mins <- which(apply(trans.df,1,min) == 0)                  # Generate list of rows with at least one occurence of no detection
trans.df <- trans.df[-zero.mins,]                               # Remove rows with one occurrence of no detection
summary(trans.df$median)                                        # Report the new range of values

### Log transform
trans.df.log <- log10(trans.df)                                 # Log transform

### Plot 2, after first filter
ggplot(data = trans.df.log,aes(x = median, y = sd)) + geom_point(size = 0.5) + theme_light() + geom_smooth(method = lm)

########################################################################
### Filter out lowest tenth percentile
tenth.perc <- quantile(trans.df.log$median, 0.1)                # Calculate the hightes expression within the lowest tenth percentile
trans.df.log <- trans.df.log[-which(trans.df.log$median <= tenth.perc),] # Remove all genes contained in lowest tenth percentile

### Prune all except med.v.sd
med.v.sd <- trans.df.log[(ncol(trans.df.log)-1):ncol(trans.df.log)] # Generate a data frame containing only the median and SD columns
summary(med.v.sd)

### Plot 3, after second filter
ggplot(data = med.v.sd, aes(x = median, y = sd)) +
  geom_point(size = 0.5) + theme_light() +
  geom_smooth(method = lm)

########################################################################
### Calculate slopes & intercepts for cutoff values
#med.v.sd.offset <- med.v.sd + 0.5                               # Generate a new data frame containing the median and SD with a 0.5 offset
offset <- min(med.v.sd$median)
med.v.sd.offset <- med.v.sd - offset                             # Generate a new data frame containing the median and SD with a 0.5 offset

### Calculate intercepts at 0.5
near.point.five <- med.v.sd.offset[med.v.sd.offset$median > 0.5 &med.v.sd.offset$median < 0.6,]
summary(near.point.five)

sd.75.at.point.five <- as.numeric(quantile(near.point.five$sd, 0.75)) # Green
sd.50.at.point.five <- as.numeric(quantile(near.point.five$sd, 0.50)) # Blue
sd.05.at.point.five <- as.numeric(quantile(near.point.five$sd, 0.05)) # Purple
sd.01.at.point.five <- as.numeric(quantile(near.point.five$sd, 0.01)) # Red

### Calculate intercepts at 2.5
near.2.5 <- med.v.sd.offset[med.v.sd.offset$median > 2.4 &med.v.sd.offset$median < 2.6,]
summary(near.2.5)

sd.99.at.2.5 <- as.numeric(quantile(near.2.5$sd, 0.99))               # Green
sd.95.at.2.5 <- as.numeric(quantile(near.2.5$sd, 0.95))               # Blue
sd.75.at.2.5 <- as.numeric(quantile(near.2.5$sd, 0.75))               # Purple
sd.50.at.2.5 <- as.numeric(quantile(near.2.5$sd, 0.50))               # Red

### Calculate slopes
slope.green <- (sd.99.at.2.5 - sd.75.at.point.five) / 2               # Green
slope.blue <- (sd.95.at.2.5 - sd.50.at.point.five) / 2                # Blue
slope.purple <- (sd.75.at.2.5 - sd.05.at.point.five) / 2              # Purple
slope.red <- (sd.50.at.2.5 - sd.01.at.point.five) / 2                 # Red

### Calculate intercepts at zero
intercept.green <- sd.75.at.point.five - slope.green * 0.5            # Green
intercept.blue <- sd.50.at.point.five - slope.blue * 0.5              # Blue
intercept.purple <- sd.05.at.point.five - slope.purple * 0.5          # Purple
intercept.red <- sd.01.at.point.five - slope.red * 0.5                # Red

### Plot 4, after adding cutoff lines
(final.plot <- ggplot(data = med.v.sd.offset,aes(x = median, y = sd)) +
  geom_point(size = 0.5) +
  geom_abline(intercept = intercept.green, slope = slope.green, col = 'green', size = 1) +
  geom_abline(intercept = intercept.blue, slope = slope.blue, col = 'blue', size = 1) +
  geom_abline(intercept = intercept.purple, slope = slope.purple, col = 'purple', size = 1) +
  geom_abline(intercept = intercept.red, slope = slope.red, col = 'red', size = 1) +
  annotate("text", label = nrow(med.v.sd.loose), x = 1, y = 3, size = 6, colour = "green4") +
  annotate("text", label = nrow(med.v.sd.neutral), x = 1.5, y = 2.6, size = 6, colour = "blue4") +
  annotate("text", label = nrow(med.v.sd.tight), x = 2.5, y = 0.9, size = 6, colour = "purple4") +
  annotate("text", label = nrow(med.v.sd.supertight), x = 3, y = 0.5, size = 6, colour = "red4") +
  labs(title = selected.tissue) +
  theme_light())

########################################################################
### Generate four filtered sets
med.v.sd.loose <- med.v.sd[(med.v.sd.offset$sd - intercept.green)/med.v.sd.offset$median <= slope.green,]
med.v.sd.neutral <- med.v.sd[(med.v.sd.offset$sd - intercept.blue)/med.v.sd.offset$median <= slope.blue,]
med.v.sd.tight <- med.v.sd[(med.v.sd.offset$sd - intercept.purple)/med.v.sd.offset$median <= slope.purple,]
med.v.sd.supertight <- med.v.sd[(med.v.sd.offset$sd - intercept.red)/med.v.sd.offset$median <= slope.red,]

########################################################################
### Add filtered median expression sets to each of four dataframes
### Open tables containing all tissue for each level
gtex.low.sd.loose <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.loose.csv', row.names = 1)
gtex.low.sd.neutral <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.neutral.csv', row.names = 1)
gtex.low.sd.tight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.tight.csv', row.names = 1)
gtex.low.sd.supertight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.supertight.csv', row.names = 1)

### Find the positions in the full tables of the genes remaining in the filtered tables
positions.loose <- match(row.names(med.v.sd.loose),row.names(gtex.low.sd.loose))
positions.neutral <- match(row.names(med.v.sd.neutral),row.names(gtex.low.sd.neutral))
positions.tight <- match(row.names(med.v.sd.tight),row.names(gtex.low.sd.tight))
positions.supertight <- match(row.names(med.v.sd.supertight),row.names(gtex.low.sd.supertight))

### Declare selected tissue
print(selected.tissue)

### Paste the new values into the full tables
gtex.low.sd.loose[selected.tissue][positions.loose,] <- med.v.sd.loose$median
gtex.low.sd.neutral[selected.tissue][positions.neutral,] <- med.v.sd.neutral$median
gtex.low.sd.tight[selected.tissue][positions.tight,] <- med.v.sd.tight$median
gtex.low.sd.supertight[selected.tissue][positions.supertight,] <- as.vector(med.v.sd.supertight$median)

### Save updated tables
write.csv(gtex.low.sd.loose, 'Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.loose.csv')
write.csv(gtex.low.sd.neutral, 'Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.neutral.csv')
write.csv(gtex.low.sd.tight, 'Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.tight.csv')
write.csv(gtex.low.sd.supertight, 'Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.supertight.csv')

### Update tracking file
tracking.file <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/tracking.file.csv',row.names = 1)
tracking.file <- rbind(tracking.file,c(selected.tissue, strftime(Sys.time(),"%a%b%d%H%M")))
write.csv(tracking.file,'Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/tracking.file.csv')

### Save plot
png(filename=paste0('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/plot.',selected.tissue,'.png'), 
    type="cairo",
    units="in", 
    width=8, 
    height=6, 
    pointsize=12, 
    res=100)
print(final.plot)
dev.off()



