# Longitudinal Study of Stool-associated Microbial Taxa in Sibling Pairs with and without Autism

## Primary Workflow

Within the "/scripts" folder there a series of rmd files and associated html outputs. The "/results" folder houses most of the output from the scripts. Some of the files require that previous scripts are ran beforehand. The following is the order in which they were ran for the analysis:


## 1. Filtering/Nomalization (filtering_normalization.Rmd)
**Heading in manuscript:** "Normalization and Taxa Filtration"

**Description:** Take raw 16s counts, normalize using DeSeq2, CSS, and total sum scaling. Remove taxa that are not present in more than 3% of samples. Remove taxa that vary significantly over time within individuals. Call clean_mapping to reformat survey answers into numerals. Please see normalization_choices.R to see analyses that were used to selecte DeSeq2 as the chosen normalization method.

**Output Files:** <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- data/Normalized/ps_* <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- data/Filtered/ps_* <br />

## 2. Metadata comparison between ASD and TD siblings (metadata_analysis.Rmd)
**Heading in manuscript:** "Demographics, diet, and lifestyle differences between cohorts"

**Description:** Contains the chi-square tests, wilcoxon ranked-sum tests, and mixed repeated measures anovas performed on the metadata. In addition, there are breakdown tables created for a variety of demographic information from the participants

**Output Files:** <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- results/metadata_analysis/*

## 3. 16S differential analysis between ASD and TD siblings (16s_analysis.Rmd)
**Heading in manuscript:** "11 ASVs are significantly associated with the ASD Phenotype, as determined by the union of at least two differential analysis methods"

**Description:** DESEQ2, Metagenomeseq, and ANCOM contrast analysis performed on the dataset. The 117 significant taxa for the longitudinal analysis were conglomerated into results/differential_abundance/significant_taxa_timeseries_genus.tsv. Differential abundance testing on ASVs agglomerated by genus can be found in results/differential_abundance/significant_taxa_all_genus.tsv. Results for differential abundance at each timepoint separately and per differential abundance test can be found in results/differential_abundance under the respective method name.

**Output Files:**  <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- results/differential_abundance/*

## 4. Visualization of significant taxa (visualize_taxa_abundance.Rmd)
**Heading in manuscript:** "11 ASVs are significantly associated with the ASD Phenotype, as determined by the union of at least two differential analysis methods"

**Description:** Plot differentially abundance taxa counts (detected in 16s_analysis.Rmd). Many of the plots are quantified by relative abundance but note that this is for visualization purposes as the format/normalization of the abundances used depended on the method of contrast analysis. In addition, permanova used DESEQ2 normalized counts as shown in the next file within the series.

**Output Files:** <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- None - manually saved

## 5. Permanova analysis (pcoa_diversity_plots.Rmd)
**Heading in manuscript:** "Demographics, diet, and lifestyle differences between cohorts"

**Description:** Permanova analysis to understand how beta-diversity (overall microbiome composition) relates to specific lifestyle variables.

**Output Files:**  <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- results/permanova/*

## 6. Analysis of potentially confounding variables (confounding_variable_analysis.Rmd)
**Heading in manuscript:** "11 ASVs are significantly associated with the ASD Phenotype, as determined by the union of at least two differential analysis methods"

**Description:** Variables significant in the permvanova were tested via Metagenomeseq to identify which ASVs were specifically related to these variables. In this way, the 117 significant taxa were labeled with other variables they were associated with that had a signficant effect on the microbial community structure. 

**Output Files:** <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- results/differential_abundance/significant_taxa_timeseries_confounders.csv <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- results/differential_abundance/significant_taxa_all_confounders.csv <br />

## 7.  Summarize effect size of feature sets (pcoa_classifier.Rmd)
**Heading in manuscript:** "Dietary/lifestyle, but not global microbiome compositional features, explain ASD phenotype"

**Description:** Perform a logistic regression on data with different input feature sets: 1. Basic age + sex model, 2. age + sex + lifestyle/dietary variables model, 3. age + sex + coordinates along principal components of microbiome model, 4. age + sex + lifestyle/dietary variables + coordinates along PC of microbiome model. Compares performance to models built using random noise vectors to replace actual features. Also produces a correlation plot between lifestyle features, and a correlation plot between coordinates on the PCs of the microbiome data and lifestyle variables. See scripts/gsea_biomarkers.ipynb for the GSEA analysis that produced Figure 3D

**Output Files:**  <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- results/summarize_effects/*  <br />

## 8. Anxiety analysis (anxiety_associations.Rmd)
**eading in manuscript:** "11 Taxa correlate with anxiety scores within and across individuals"

**Description:** Find associations between cASV counts and changes in anxiety within an individual across time using a spearman correlation test. 

**Output Files:**  <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- results/summarize_effects/*

## 9. Additional scripts
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - **Step1_SequencingSummary_MakePhyloseqObject_2020-03-23.html** provides a summary of the results of the 16S sequencing processing.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - **digital-phenotype :** specific ASVs associated with MARA scores (ASD severity)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - **build_tree.R :** Build phylogenetic tree used to calculate unifrac distances

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - **clean_mapping.R :** Turn survey responses into numerals

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - **rarefaction curve.R :** Analyze sequencing depth


# Figure Guide:

**Figure 1:** <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - Manually created infographic

**Figure 2:** <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - 16s_analysis.Rmd <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - visualize_taxa_abundance.Rmd <br />

**Figure 3:** <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - pcoa_classifier.Rmd <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - gsea_biomarkers.ipynb <br />

**Figure 4:** <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - anxiety_associations.Rmd <br />


**Supplementary File 1:** Methods description  <br />

**Supplementary File 2:** Lifestyle/dietary survey results per person  <br />

**Supplementary File 3:** 117 ASVs identified by differential abundance testing, method used to identify, taxonomic annotation for the ASV, and confounding lifestyle variables associated with those ASVs  <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - 16s_analysis.Rmd

**Supplementary File 4:** Centered log transform ASVs significantly associated with the ASD cohort. <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - visualize_taxa_abundance.Rmd with results from 16s_analysis.Rmd 

**Supplementary File 5:**  Differentially abundant genera in ASD <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - visualize_taxa_abundance.Rmd with results from 16s_analysis.Rmd


**Supplementary File 6:** Lifestyle and dietary features significantly associated with overall microbiome composition <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - pcoa_diversity_plots.Rmd


**Supplementary File 7:** PCoAs using bray curtis distance constrained by phenotype (ASD or TD) or lifestyle variables <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- pcoa_diversity_plots.Rmd


**Supplementary File 8:** PCoAs using weighted UniFrac distances constrained by phenotype (ASD or TD) or lifestyle variables <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - pcoa_diversity_plots.Rmd


**Supplementary File 9:** PCoAs using unweighted UniFrac distances constrained by phenotype (ASD or TD) or lifestyle variables <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - pcoa_diversity_plots.Rmd


**Supplementary File 10:** Diversity between ASD and TD, diversity correlations between ASD severity metric and age <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - pcoa_diversity_plots.Rmd 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - age_MARA_diversity_correlations.Rmd


**Supplementary File 11:** PCoA using bray curtis distance and constrained and colored by ASD associated behavioral trait <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - pcoa_diversity_plots.Rmd


**Supplementary File 12:** Full comparisions of all lifestyle variables between ASD and TD <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - metadata_analysis.Rmd


**Supplementary File 13:** PCoA constrained by 1. Bacteroides to Prevotella ratio 2. Firmicutes to Bacteroidetes ratio 3. Differences between ratios of phyla compared between ASD and TD <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - enterotypes.Rmd

**Supplementary File 14:** Intra-individual/ inter-individual bray curtis distances compared across normalization methods <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - normalization_choices.Rmd