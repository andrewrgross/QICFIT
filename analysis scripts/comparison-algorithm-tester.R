### Comparison algorithm tester -- Andrew R Gross -- 2017-01-18
### Open individual samples from each of the GTEx references for testing cell identity prediction methods

########################################################################
### Header
########################################################################

########################################################################
### Import data
########################################################################

### Import sample key
sample.key <- read.csv('Z:/Data/Andrew/reference_data/gtex/Sample.key.sorted.csv',header=TRUE)

### Import abreviated reference set
reference.transcriptomes <- read.csv('Z:/Data/Andrew/reference_data/gtex/individual-ref-transcriptomes-for-testing.csv', row.names = 1) # ~10 seconds



########################################################################
### Format
########################################################################


head(trans.df[1:5])
trans.df <- trans.df[3:12]
