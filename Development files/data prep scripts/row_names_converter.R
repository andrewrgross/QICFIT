### row names converter -- Andrew R Gross

### Rename rows in tissue once and for all

transcriptome.index <- read.csv("Z:/Data/Andrew/reference_data/gtex/transcriptome_tables_by_tissue/transcriptome.index.csv")
print(transcriptome.index[1])
selection <- 10
(selected.tissue <- as.character(transcriptome.index[,1][selection]))

for(selection in 10:nrow(transcriptome.index)) {
  selected.tissue <- as.character(transcriptome.index[,1][selection])
  print(paste('Loading',selected.tissue))
  trans.df <- read.delim(as.character(transcriptome.index$file.locations[selection]))  # Load selected dataset (may take >10s)
  print(paste(selected.tissue,'loaded at', strftime(Sys.time(),"%a%b%d%H%M")))
  ### Rename rows
  row.names(trans.df) <- trans.df$Name
  trans.df <- convertIDs(trans.df)
  trans.df <- trans.df[2:ncol(trans.df)]
  ### Save the renamed table
  write.csv(trans.df, paste0(as.character(transcriptome.index$file.locations[selection]),'_renamed.csv'))
}

