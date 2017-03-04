### Manual Spearman calculation

### Declare data
sample.data <- sample
ref.data <- ref.full[1]

### Filter data
sample.data.filtered <- sample.data[sample.data$AA_Colon >1,,drop = FALSE]
ref.data.filtered <- ref.data[row.names(sample.data.filtered),, drop = FALSE]

### Join data
joined.data <- cbind(sample.data.filtered,ref.data.filtered)
joined.data$sample_rank <- rank(-joined.data$AA_Colon)
joined.data$ref_rank <- rank(-joined.data$Adipose...Subcutaneous)

### Define the difference in each row
joined.data$d <- joined.data$sample_rank - joined.data$ref_rank
joined.data$d2 <- (joined.data$d)^2

### Define correction factor
sample.ties <- nrow(joined.data)- length(unique(joined.data$sample_rank))
ref.ties <- nrow(joined.data)- length(unique(joined.data$ref_rank))
sample.cf <- (sample.ties^2 - sample.ties)/12
ref.cf <- (ref.ties^2 - ref.ties)/12

### Calculate the correlation
sum.di.sq <- sum(joined.data$d2)
sum.di.sq.corrected <- sum.di.sq + sample.cf + ref.cf
(p <- 1 - (6*sum.di.sq.corrected)/(nrow(joined.data)^3-nrow(joined.data)))

### Check work
cor(as.matrix(joined.data[1:2]),use="complete.obs", method = "spearman")

rcorr(as.matrix(joined.data[1:2]), type = 'spearman')[[1]] # Perform spearman calculation



sample.data <- sample.data[]


test <- head(sample.data.filtered)

head(ref.data,20)

### Weighted correlation

library('wCorr')


