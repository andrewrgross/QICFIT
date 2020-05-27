### Amalgamating Reference Tissues -- Andrew R Gross -- 2017-01-17
### Generate a dendrogram relating all the tissues, then generate amalgams of tissue groups to search.

########################################################################
### Header
########################################################################


########################################################################
### Import data
########################################################################

mtcars
mtcars.t <- t(mtcars)

### Import sample key
sample.key <- read.csv('Z:/Data/Andrew/reference_data/gtex/Sample.key.sorted.csv',header=TRUE)

### Open tables containing all tissue for each level
gtex.low.sd.loose <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.loose.csv', row.names = 1)
gtex.low.sd.neutral <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.neutral.csv', row.names = 1)
gtex.low.sd.tight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.tight.csv', row.names = 1)
gtex.low.sd.supertight <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.low.sd.supertight.csv', row.names = 1)

### Import references
references <- read.csv('Z:/Data/Andrew/reference_data/gtex/sd.filtered.tables/gtex.full.csv',header=TRUE,row.names=1)


########################################################################
### Format
########################################################################

### Reference names, formatted
referenceNames <- c("Adipose, subcutaneous","Adipose, omentum","Adrenal gland","Aorta","Coronary artery","Tibial artery","Bladder","Amygdala","Anteriror cingulate nucleous","Caudate nucleous",
                    "Cerebellar hemisphere","Cerebellum","Cortex","Frontal cortex BA9","Hippocampus","Hypothalamus","Nucleus accumbens","Putamen",
                    "Spinal cord","Substantia nigra","Mammary","Lymphocyte","Fibroblast","Ectocervix","Endocervix","Colon, sigmoid","Colon, transverse","Gastroesophageal junction","Esophagus, mucosa",
                    "Esophagus, muscularis","Fallopian tube","Heart, Atrial","Heart, left ventricle","Kidney","Liver","Lung","Salvitory gland","Skeletal muscle","Tibial nerve","Ovary","Pancreas",
                    "Pituitary","Prostate","Skin, sun-hidden","Skin, sun-exposed","Small intestine","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole blood")

ref.full <- references[2:ncol(references)]
ref.loose <- gtex.low.sd.loose[2:ncol(gtex.low.sd.loose)]
ref.neutral <- gtex.low.sd.neutral[2:ncol(gtex.low.sd.neutral)]
ref.tight <- gtex.low.sd.tight[2:ncol(gtex.low.sd.tight)]
ref.supertight <- gtex.low.sd.supertight[2:ncol(gtex.low.sd.supertight)]

names(ref.full) <- referenceNames
names(ref.supertight) <- referenceNames
names(ref.loose) <- referenceNames
names(ref.neutral) <- referenceNames
names(ref.tight) <- referenceNames

### Filter out empty rows

empty.rows <- rowSums(is.na(ref.loose))>=52
ref.loose <- ref.loose[!empty.rows,]
empty.rows <- rowSums(is.na(ref.neutral))>=52
ref.neutral <- ref.neutral[!empty.rows,]
empty.rows <- rowSums(is.na(ref.tight))>=52
ref.tight <- ref.tight[!empty.rows,]
empty.rows <- rowSums(is.na(ref.supertight))>=52
ref.supertight <- ref.supertight[!empty.rows,]


### Declare rows to keep from full

rows.loose <- row.names(ref.loose)
rows.neutral <- row.names(ref.neutral)
rows.tight <- row.names(ref.tight)
rows.supertight <- row.names(ref.supertight)

### Generate new filtered expression sets

ref.loose <- ref.full[rows.loose,]
ref.neutral <- ref.full[rows.neutral,]
ref.tight <- ref.full[rows.tight,]
ref.supertight <- ref.full[rows.supertight,]

### Filter full by average expression
row.sums <-rowSums(ref.full)
summary(row.sums)
rows.over100 <- row.sums > 800
ref.over100 <- ref.full[rows.over100,]
ref.over100 <- log(ref.over100+10)
nrow(ref.over100)

### Calculate and plot
input.data <- ref.supertight[8:20]
condition <- 'over 100, filtered'
ref.dist <- dist(t(as.matrix(input.data)))
ref.hclust <- hclust(ref.dist)
plot(ref.hclust, main = paste(condition))










### Define conditions
conditions <- c('complete', 'na = 0', 'remove empty rows', 'remove empty rows & na = 0', 'remove incomplete rows')

### Select base
base <- ref.supertight
baseName <- 'supertight'

### Define input.data, unedited
input.data <- base
(condition <- conditions[1])
nrow(input.data)

ref.dist <- dist(t(as.matrix(input.data)))
ref.hclust <- hclust(ref.dist)
plot(ref.hclust, main = paste(baseName, condition))


### Define input.data, na = 0
input.data <- base
input.data[is.na(input.data)] = 0
(condition <- conditions[2])
nrow(input.data)

ref.dist <- dist(t(as.matrix(input.data)))
ref.hclust <- hclust(ref.dist)
plot(ref.hclust, main = paste(baseName, condition))



### Define input.data, remove incomplete rows
incomplete.rows <- rowSums(is.na(base))>0
input.data <- base[!incomplete.rows,]
(condition <- conditions[5])
nrow(input.data)

ref.dist <- dist(t(as.matrix(input.data)))
ref.hclust <- hclust(ref.dist)
plot(ref.hclust, main = paste(baseName, condition))












main = 'tight, no na'

########################################################################
### Cluster
########################################################################

ref.dist <- dist(t(as.matrix(ref.st.ts.some.na)))
ref.dist <- dist(t(as.matrix(ref.st.ts)))

ref.hclust <- hclust(ref.dist)

########################################################################
### Plot
########################################################################

plot(ref.hclust,main = main)

########################################################################
### Generate amalgams
########################################################################

# Loop


########################################################################
### Save
########################################################################



########################################################################
### 
########################################################################





