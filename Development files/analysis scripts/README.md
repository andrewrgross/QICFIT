#QICFIT/analysis scripts

This folder contains scripts written in the process of developing QICFIT.

fingerprint_matrix_test

comparison-algorithm-tester

Manual spearman calculations
This was an experiment in reproducing spearman calculations which previously were completed by an R function in the package corr.  It was successful, although assigning a Spearman correlation when values tie is more challenging, so the dedicated function will continue to serve as my preferred means of calculating Spearman values.

qicfit_plotting_draft

QICFIT_test_outside_data
### QICFIT test using outside data -- Andrew R Gross -- 2017-03-01
This script is used to load and format data from other sources, then to test various recognition methods. This script reads in expression tables; converts the IDs to ENSEMBL; Loops through each transcriptome and does a Spearman check against references; and outputs figures and tables reportingThe similarity of each query to each of the GTEx references

qicfit_version_01

relevant gene identifier
### Identify most relevant genes -- Andrew R Gross -- 2017-04-17
I want to identify the genes with the greatest relevance to the low Spearman correlation, then screen out housekeeping genes.
This script is an evolution of the manual Spearman calculation method.  It applies components of a manual Spearman correlation calculation in order to present genes which were most responsible for the divergence of the sample and its closest match.

<<<<<<< HEAD:analysis scripts/README.md.txt

spearman_using_filtered_data
### Spearman comparison using filtered data -- Andrew R Gross -- 2016/12/16 
This script uses GTEx data as the input for testing various strategies of recognizing samples and displaying results.  It is an early attempt at calculating Spearman values which used filtered datasets, the use of which have since been indefinitely discontinued in favor of using full datasets.
=======
spearman_using_filtered_data
>>>>>>> d20730591d5e7f3fc62ed5bc9cd605dfdb6c6a48:analysis scripts/README.md
