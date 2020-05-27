### QICFIT:

#spearman.corr.for.single.sample <- function(query.df, ref.df) {   # Calculate the spearman correlation between the query.df and the references


### Load data


TPMdata <- read.csv("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E196 - Pituitary/Input files/LP-6689--04--08--2019_TPM.csv", row.names = 1)

pituitary.profile <- read.csv("Z:/Data/Andrew/reference_data/gtex/transcriptome_tables_by_tissue/pituitary_renamed.csv",row.names = 1)

### Format

pituitary.profile <- convertIDs(pituitary.profile)[2:12]



TPMdata <- convertIDs(TPMdata)

g1 <- TPMdata[7:9]
g2 <- TPMdata[10:12]
g3 <- TPMdata[13:15]
g4 <- TPMdata[1:3]
g5 <- TPMdata[4:6]

names(g1) <- c('LP_1', 'LP_2', 'LP_3')
names(g2) <- c('LP_4', 'LP_5', 'LP_6')
names(g3) <- c('LP_7', 'LP_8', 'LP_9')
names(g4) <- c('LP_10', 'LP_11', 'LP_12')
names(g5) <- c('LP_13', 'LP_14', 'LP_15')

gtex.pituitary <- ''


query.df <- TPMdata
spearman.results.pit <- spearman.calc(pituitary.profile, ref.df)

query.df <- g3
spearman.results <- spearman.calc(query.df,ref.df = ref.df)

spearman.results.g1 <- spearman.calc(g1,ref.df = ref.df)

spearman.results.g2 <- spearman.calc(g2,ref.df = ref.df)

spearman.results.g3 <- spearman.calc(g3,ref.df = ref.df)

spearman.results.g4 <- spearman.calc(g4,ref.df = ref.df)

spearman.results.g5 <- spearman.calc(g5,ref.df = ref.df)

### Plot

plot.spearman(spearman.results = spearman.results.g3, title = 'g3')

plot.spearman(spearman.results = spearman.results.g4, title = 'g4')

plot.spearman(spearman.results = spearman.results.g5, title = 'g5')

plot.spearman(spearman.results = spearman.results.g1, title = 'g1')

plot.spearman(spearman.results = spearman.results.g2, title = 'g2')

plot.spearman(spearman.results.pit[1:6], title = "GTEX")

### Gene assessment
query.df = g3[1]
assessment.g3 <- spearman.gene.assessor.single.sample(query.df = query.df, ref.df=ref.df)



spearman.gene.assessor.single.sample <- function(query.df, ref.df, weighted = FALSE) {   # Calculate the spearman correlation between the query.df and the references
  spearman.results <- data.frame(rep(0,ncol(reference.input)))                             # Generate empty results table
  row.names(spearman.results) <- names(ref.df)                                # Name empty results table
  names(spearman.results) <- names(query.df)

  ref.tissue.num = 42

  for(ref.tissue.num in 1:ncol(ref.df)) {                                     # Loop through each reference
    ref.tissue.data <- ref.df[ref.tissue.num]                                 # Call the current tissue from the references
    print(tissue)
    tissue <- names(ref.tissue.data)
    ### Generate a data frame containing the two transcriptomes being compared
    ref.tissue.data <- ref.tissue.data[which(ref.tissue.data[,1]>0),,drop=FALSE]         # Filter out missing values from tissue
    genes.present.in.ref <- row.names(ref.tissue.data)                                   # Declare the genes present in the reference
    genes.missing.in.query <- setdiff(genes.present.in.ref, row.names(query.df))    # Declare genes in reference missing from query.df
    rows.to.add <- data.frame(rep(0,length(genes.missing.in.query)))                       # Generate a zero data frame the size of the missing rows
    row.names(rows.to.add) <- genes.missing.in.query                                       # Name rows after missing rows
    names(rows.to.add) <- names(query.df)                                           # Name column the same as the query
    query.df2 <- rbind(query.df,rows.to.add)                                # Use rbind to make a full data frame containing exactly the genes in the reference
    query.df2 <- query.df2[genes.present.in.ref, , drop = FALSE]           # Reorder query.df to match reference
    spearman.input <- cbind(query.df2, ref.tissue.data)                          # Bind query.df and reference

    ### Assign rank and delta rank
    spearman.input$sample_rank <- rank(-spearman.input[1])
    spearman.input$ref_rank <- rank(-spearman.input[2])
    spearman.input$d <- spearman.input$sample_rank - spearman.input$ref_rank
    spearman.input$d2 <- (spearman.input$d)^2

    if(weighted == TRUE) {
      ### Calculate weighted delta rank
      spearman.input$d2 <- spearman.input$d2 * ((-2*spearman.input$ref_rank/max(spearman.input$ref_rank))+2)
    }

    gene.report <- spearman.input[order(spearman.input$d2, decreasing = TRUE),]
    # Add to results table
  }
  return(spearman.results)
}
