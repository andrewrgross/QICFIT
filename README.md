# QICFIT

QICFIT stands for Quantitative Identity Classifier For Interrogating Transcriptomes.  It's an R package which allows a user to assess which tissue a sample resembles based on its transcriptome.

QICFIT works by performing a Spearman rank correlation between an RNA-seq transcriptome and a set of reference trasncriptomes from the Genotype-Tissue Expression (GTEx) project.  This means that QICFIT is scale-independent.  Since a Spearman correlation is based on rank, it can compare RNA-seq data between linear and log scales without issue.  The expression values themselve are irrelevant.  

Currently, QICFIT is an in-house tool used by the Sareen lab to assess the success of stem cell differentiation protocols in maturing pluripotenet stem cells into differentiated cells.  As such, I have no composed a user guide or tested QICFIT under a range of conditions and data types.  If you would like to use QICFIT or assist in its development, feel free to email me.

#
