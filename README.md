# Longitudinal Study of Stool-associated Microbial Taxa in Sibling Pairs with and without Autism

## Primary Workflow

Within the "/scripts" folder there a series of rmd files and associated html outputs. The "/results" folder houses most of the output from the scripts. Some of the files require that previous scripts are ran beforehand. The following is the order in which they were ran for the analysis:

### Filtering/Nomalization

The first in the series is the "filtering_normalization.rmd" which loads the phyloseq object generated from DADA2, performs a variety of filtering, and two normalizations (although DESEQ2 is the normalization chosen for final analysis)

### Metadata comparison between ASD and TD siblings

The "metadata_analysis.rmd" file contains the chi-square tests, wilcoxon ranked-sum tests, and mixed repeated measures anovas performed on the metadata. In addition, there are breakdown tables created for a variety of demographic information from the participant.s

### 16S Contrast analysis between ASD and TD siblings

The "16S_analysis.rmd" file contain the DESEQ2, Metagenomeseq, and ANCOM contrast analysis performed on the dataset. The 117 significant taxa for the longitudinal analysis were conglomerated into the "significant_taxa_timeseries_genus.tsv" file in "/results".

#### Visualization of significant taxa

The "visualize_taxa_abundance.rmd" file contains various plots of the taxa found signficant in this study. Many of the plots are quantified by relative abundance but note that this is for visualization purposes as the format/normalization of the abundances used depended on the method of contrast analysis. In addition, permanova used DESEQ2 normalized counts as shown in the next file within the series.

### Permanova and Analysis of potentially confounding variables

"confounding_variable_analysis.rmd" is the file in which we conducted our permanova analysis. In addition, variables significant in the permvanova were tested via Metagenomeseq. In this way, the 117 significant taxa were labeled with other variables they were associated with that had a signficant effect on the microbial community structure. 

### Random Forest Classifiers 

The "random_forest.rmd" holds the code used for our Random Forest machine learning analysis using the taxa found to be signficant in this study. All of the code and results from these classifiers can be found here. 

## Additional scripts

"normalization_choices.R" is additional comparison regarding the different normalizations.

"digital-phenotype" is an analysis of MARA scores and specific taxa. This is not discussed within the paper.

