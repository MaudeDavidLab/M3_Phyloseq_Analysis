library(dplyr)
library(ggplot2)
library(reshape2)
library(phyloseq)
library(rstatix)
library(ggpubr)
library(pbapply)
library(DESeq2)
library(DT)
library(metagenomeSeq)
setwd("C:/Users/ctata/Documents/Lab/M3_longitudinal/M3_Phyloseq_Analysis")
#Read in phyloseq objects
ps_deseq <- readRDS("Normalized/ps_DeSeq_pass_min_postDD_min0.03.rds")
ps_css <- readRDS("Normalized/ps_CSS_pass_min_postDD_min0.03.rds")
ps_no_norm <- readRDS("Normalized/ps_not_norm_comp_pass_min_postDD_min0.03.rds")

deseq_cut = .05
mtgseq_cut = .05

#Christine: Use freidman test instead of anova because the distributions of taxa are not normal across samples
#Friedman is the parameteric version of repeated measure anova
getASVsChangeSigOverTime <- function(ps){
  asv_table<-t(otu_table(ps))
  asv_table <- as.data.frame(asv_table)
  asv_table$timepoint <- ps@sam_data$Within.study.sampling.date
  asv_table$timepoint <- gsub(" ", "", asv_table$timepoint)
  asv_table$timepoint <- as.factor(asv_table$timepoint)
  asv_table$Host.Name <- ps@sam_data$Host.Name
  
  asv_long <- melt(asv_table, id = c("timepoint", "Host.Name"))
  change_over_timepoints_pvals <- pbsapply( unique(asv_long$variable), function(seq){
    asv_tmp <- asv_long[asv_long$variable == seq, ]
    asv_tmp <- select(asv_tmp, -c("variable"))
    res.fried <- asv_tmp %>% friedman_test(value ~ timepoint |Host.Name)
    return(res.fried$p)
  })
  asvs <- unique(asv_long$variable)[change_over_timepoints_pvals < .1]
  return(asvs)
}

#Set threshold for p value at .1 - If taxa are changing that much over 3 time points, their fluctuations are likely due to seasons or dietary
#changes rather than long standing gut-brain axis phenomenon
asv_deseq <- as.character(getASVsChangeSigOverTime(ps_deseq))
asv_css <- as.character(getASVsChangeSigOverTime(ps_css))
asv_intersect <- intersect(asv_deseq, asv_css)

#Which ASVs change significantly over time independent of normalization method? 58 of them.
intersect(asv_deseq, asv_css)


#remove these taxa from phyloseqs (since they are noise essentially)
keep = taxa_names(ps_deseq)[!(taxa_names(ps_deseq) %in% asv_deseq)]
ps_deseq_filt <- prune_taxa(keep, ps_deseq)

keep = taxa_names(ps_css)[!(taxa_names(ps_css) %in% asv_css)]
ps_css_filt <- prune_taxa(keep, ps_css)

keep = taxa_names(ps_no_norm)[!(taxa_names(ps_no_norm) %in% asv_intersect)]
ps_no_norm_filt <- prune_taxa(keep, ps_no_norm)


# DA taxa identification
#Identify bacterial and archaeal taxa (genera, species and strains) whose abundance is observed significantly more or less in the ASD

## DESeq

dir.create(paste0(output_data, 'DESeq/'))
###Run DESeq proper (not just the normalization but all of it)
runDESeq <- function(ps, dcut){
  diagdds = phyloseq_to_deseq2(ps, ~ phenotype) 
  #diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  res = results(diagdds, contrast = c("phenotype", "N", "A"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj < dcut,]
  if (dim(sig)[1] == 0) 
  {sigtab<- as.data.frame(1, row.names="nothing")
  colnames(sigtab) <- 'padj'}
  else 
  {
    sigtab <- data.frame(sig)
  }
  return(list(res, sigtab))
}


###Running analysis 
###split thedata based on the real 3 timepoints
P1<-prune_samples(sample_names(ps_no_norm_filt)[ps_no_norm_filt@sam_data$Within.study.sampling.date == "Timepoint 1"], ps_no_norm_filt)
P2<-prune_samples(sample_names(ps_no_norm_filt)[ps_no_norm_filt@sam_data$Within.study.sampling.date == "Timepoint 2"], ps_no_norm_filt)
P3<-prune_samples(sample_names(ps_no_norm_filt)[ps_no_norm_filt@sam_data$Within.study.sampling.date == "Timepoint 3"], ps_no_norm_filt)

#several significants 
deseq_res_P1 <- runDESeq(P1, deseq_cut) 
deseq_res_P2 <- runDESeq(P2, deseq_cut)
deseq_res_P3 <- runDESeq(P3, deseq_cut)

# print significant taxa
datatable(deseq_res_P1[[2]])
datatable(deseq_res_P2[[2]])
datatable(deseq_res_P3[[2]])
#"ASV_1669" present twice timepoint 1 and 3 


View(as.data.frame(ps_no_norm_filt@tax_table[intersect(rownames(deseq_res_P1[[2]]), rownames(deseq_res_P2[[2]])), ]))

# save
saveRDS(deseq_res_P1, file=paste0(output_data, "DESeq/deseq_res_P1_adjp", deseq_cut, ".Rda"))
saveRDS(deseq_res_P2, file=paste0(output_data, "DESeq/deseq_res_P2_adjp", deseq_cut, ".Rda"))
saveRDS(deseq_res_P3, file=paste0(output_data, "DESeq/deseq_res_P3_adjp", deseq_cut, ".Rda"))

#Working with time series
#according to the DeSeq vignette: design including the time factor, and then test using the likelihood ratio test as described
#the following section, where the time factor is removed in the reduced formula
runDESeq_time <- function(ps, dcut){
  diagdds = phyloseq_to_deseq2(ps, ~ phenotype + Within.study.sampling.date) 
  diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
  diagdds <- DESeq(diagdds,fitType="parametric", betaPrior = FALSE) 
  #resultsNames(diagdds): to determine the constrast
  res = results(diagdds, contrast = c("phenotype", "A", "N"))
  res$padj[is.na(res$padj)] = 1
  sig <- res[res$padj < dcut,]
  if (dim(sig)[1] == 0) 
  {sigtab<- as.data.frame(1, row.names="nothing")
  colnames(sigtab) <- 'padj'}
  else 
  {
    sigtab <- data.frame(sig)
  }
  return(list(res, sigtab))
}

#and this time when factoring in the interaction for longitudinal study! 
deseq_res_time <- runDESeq_time(ps_no_norm_filt, deseq_cut)
View(as.data.frame(ps_no_norm_filt@tax_table[rownames(deseq_res_time[[2]]), ]))
saveRDS(deseq_res_time, file=paste0(output_data, "DESeq/Deseq_time_interaction_adjp", deseq_cut, ".Rda"))


## metagenomeSeq

dir.create(paste0(output_data, 'metagenomseq/'))
###Run ZIG model fitting and prediction
run_metagenom_seq<-function(ps,maxit, mcut, time = F){
  p_metag<-phyloseq_to_metagenomeSeq(ps)
  #filtering at least 4 samples 
  p_metag= cumNorm(p_metag, p=0.75)
  normFactor =normFactors(p_metag)
  normFactor =log2(normFactor/median(normFactor) + 1)
  if(time == F){
    mod = model.matrix(~phenotype + normFactor, data = pData(p_metag))
  }
  if(time == T){
    mod = model.matrix(~phenotype + Within.study.sampling.date +normFactor, data = pData(p_metag))
  }
  
  settings =zigControl(maxit =maxit, verbose =FALSE)
  fit =fitZig(obj = p_metag, mod = mod, useCSSoffset = FALSE, control = settings)
  #Note: changed fit$taxa to fit@taxa in light of error (probably from newer metagenomeseq ver.)
  res_fit<-MRtable(fit, number = length(fit$taxa))
  res_fit_nonfiltered <- copy(res_fit)
  res_fit<-res_fit[res_fit$adjPvalues<mcut,]
  #finally remove the ones that are not with enough samples
  Min_effec_samp<-calculateEffectiveSamples(fit)
  Min_effec_samp<-Min_effec_samp[ names(Min_effec_samp)  %in% rownames(res_fit)] #####there is a bug here 
  #manually removing the ones with "NA"
  res_fit<-res_fit[grep("NA",rownames(res_fit), inv=T),]
  res_fit$Min_sample<-Min_effec_samp
  res_fit<-res_fit[res_fit$`+samples in group 0` >= Min_effec_samp & res_fit$`+samples in group 1` >= Min_effec_samp,]
  return(list(res_fit_nonfiltered, res_fit))
}


zig_res_P1 <- run_metagenom_seq(P1,30, mtgseq_cut) # 30The maximum number of iterations for the expectation-maximization algorithm
zig_res_P2 <- run_metagenom_seq(P2,30, mtgseq_cut)
zig_res_P3 <- run_metagenom_seq(P3,30, mtgseq_cut)

zig_res_all<- run_metagenom_seq(ps_no_norm_filt,30, mtgseq_cut, time = T)

# print significant taxa
datatable(zig_res_P1[[2]])
datatable(zig_res_P2[[2]])
datatable(zig_res_P3[[2]])
datatable(zig_res_all[[2]])

zig_res_P1_filtered <- data.frame(cbind(zig_res_P1[[2]], tax_table(P1)[rownames(zig_res_P1[[2]]),]))
zig_res_P1_filtered$enriched <- ifelse(zig_res_P1_filtered$phenotypeN < 0, "Aut", "Control")

zig_res_P2_filtered <- data.frame(cbind(zig_res_P2[[2]], tax_table(P2)[rownames(zig_res_P2[[2]]), ]))
zig_res_P2_filtered$enriched <- ifelse(zig_res_P2_filtered$phenotypeN < 0, "Aut", "Control")

zig_res_P3_filtered <- data.frame(cbind(zig_res_P3[[2]], tax_table(P3)[rownames(zig_res_P3[[2]]), ]))
zig_res_P3_filtered$enriched <- ifelse(zig_res_P3_filtered$phenotypeN < 0, "Aut", "Control")

zig_res_all_filtered <- data.frame(cbind(zig_res_all[[2]], tax_table(ps_no_norm_filt)[rownames(zig_res_all[[2]]), ]))
zig_res_all_filtered$enriched <- ifelse(zig_res_all_filtered$phenotypeN < 0, "Aut", "Control")


saveRDS(zig_res_P1, file=paste0(output_data, "metagenomseq/zig_res_P1_adjp", mtgseq_cut, ".rds"))
saveRDS(zig_res_P2, file=paste0(output_data, "metagenomseq/zig_res_P2_adjp", mtgseq_cut, ".rds"))
saveRDS(zig_res_P3, file=paste0(output_data, "metagenomseq/zig_res_P3_adjp", mtgseq_cut, ".rds"))

#do we have any in ESV in common?
#Reduce(intersect, list(rownames(zig_res_P1_filtered),rownames(zig_res_P2_filtered),rownames(zig_res_P3_filtered)))
#rownames(zig_res_P1[[2]])[which(rownames(zig_res_P1[[2]]) %in% c(rownames(zig_res_P2[[2]]), rownames(zig_res_P3[[2]])))]
#rownames(zig_res_P2[[2]])[which(rownames(zig_res_P2[[2]]) %in% rownames(zig_res_P3[[2]]))]

metag_res<-cbind(rownames(zig_res_all[[2]]), as(tax_table(ps_no_norm_filt)[rownames(zig_res_all[[2]]), ], "matrix"))
metag_res_esv <-metag_res
rownames(metag_res) <- NULL
metag_res


#New ANCOM 2.1

runAncom <- function(ps){
  metada_ps<-sample_data(ps)
  metada_ps<-as.data.frame(metada_ps)
  metada_ps$HostName <- as.character(metada_ps$Host.Name)
  metada_ps$phenotype <- as.character(metada_ps$phenotype)
  
  metada_ps<-tibble(metada_ps$HostName, metada_ps$phenotype, rows = rownames(metada_ps))
  colnames(metada_ps) <- c("HostName", "phenotype", "Biospecimen.Barcode")
  res_table_filt<-otu_table(ps)
  res_table_filt<-as.data.frame(res_table_filt)
  
  prepro<-feature_table_pre_process(feature_table = res_table_filt, meta_data = metada_ps, sample_var = "Biospecimen.Barcode", group_var = "phenotype", out_cut = 0.05, zero_cut = 0.90, lib_cut = 1000, neg_lb = FALSE)
  
  feature_table = prepro$feature_table # Preprocessed feature table
  meta_data = prepro$meta_data # Preprocessed metadata
  struc_zero = prepro$structure_zeros # Structural zero info
  
  main_var = "phenotype"; p_adj_method = "BH"; alpha = 0.05
  adj_formula = NULL; rand_formula = NULL
  ancom_res <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                     alpha, adj_formula, rand_formula)#, control)
  return(ancom_res)
}

ancom.all.filt.0.05 <- runAncom(ps_no_norm_filt)
dir.create(path = "ANCOM")
saveRDS(ancom.all.filt.0.05, file=paste0("ANCOM/ancom_res", mtgseq_cut, ".rds"))
#no significant taxa from ANCOM



#Trying by each timepoint
#Timepoint 3
ps3<-subset_samples(ps_no_norm_filt, Within.study.sampling.date == "Timepoint 3" )
ancom.all.filt.0.05.3 <- runAncom(ps3)
#dir.create(path = "ANCOM")
saveRDS(ancom.all.filt.0.05.3, file=paste0(output_data, "ANCOM/ancom_res3_", mtgseq_cut, ".rds"))

#Format the Detected Table w/ Taxa 
sum(ancom.all.filt.0.05.3$out$detected_0.6)
# No sig


#Timepoint 2

ps2<-subset_samples(filtered_ps003, Within.study.sampling.date == "Timepoint 2" )
metada_ps<-sample_data(ps2)
metada_ps<-as.data.frame(metada_ps)
metada_ps$HostName <- as.character(metada_ps$Host.Name)
metada_ps$phenotype <- as.character(metada_ps$phenotype)

metada_ps<-tibble(metada_ps$HostName, metada_ps$phenotype, rows = rownames(metada_ps))
colnames(metada_ps) <- c("HostName", "phenotype", "Biospecimen.Barcode")
res_table_filt<-otu_table(ps2)
res_table_filt<-as.data.frame(res_table_filt)

prepro<-feature_table_pre_process(feature_table = res_table_filt, meta_data = metada_ps, sample_var = "Biospecimen.Barcode", group_var = "phenotype", out_cut = 0.05, zero_cut = 0.90, lib_cut = 1000, neg_lb = FALSE)

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = "phenotype"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
ancom.all.filt.0.05.2 <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                               alpha, adj_formula, rand_formula)#, control)
#dir.create(path = "ANCOM")
saveRDS(ancom.all.filt.0.05.2, file=paste0(output_data, "ANCOM/ancom_res2_", mtgseq_cut, ".rds"))

#Format the Detected Table w/ Taxa 
ancom.all.filt.0.05.2$detected
# No sig



#Timepoint 1

ps1<-subset_samples(filtered_ps003, Within.study.sampling.date == "Timepoint 1" )
metada_ps<-sample_data(ps1)
metada_ps<-as.data.frame(metada_ps)
metada_ps$HostName <- as.character(metada_ps$Host.Name)
metada_ps$phenotype <- as.character(metada_ps$phenotype)

metada_ps<-tibble(metada_ps$HostName, metada_ps$phenotype, rows = rownames(metada_ps))
colnames(metada_ps) <- c("HostName", "phenotype", "Biospecimen.Barcode")
res_table_filt<-otu_table(ps1)
res_table_filt<-as.data.frame(res_table_filt)

prepro<-feature_table_pre_process(feature_table = res_table_filt, meta_data = metada_ps, sample_var = "Biospecimen.Barcode", group_var = "phenotype", out_cut = 0.05, zero_cut = 0.90, lib_cut = 1000, neg_lb = FALSE)

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

main_var = "phenotype"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
ancom.all.filt.0.05.1 <- ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
                               alpha, adj_formula, rand_formula)#, control)
#dir.create(path = "ANCOM")
saveRDS(ancom.all.filt.0.05.1, file=paste0(output_data, "ANCOM/ancom_res1_", mtgseq_cut, ".rds"))

#Format the Detected Table w/ Taxa 
ancom.all.filt.0.05.1$detected
# No sig



```

## Results table summarizing significant taxa

```{r resultstable}
#Deseq results

#Assemble Deseq results into one table
DesP1<-tibble(rownames(deseq_res_P1[[2]]), rep("DESeq2_P1", length(rownames(deseq_res_P1[[2]]))))
colnames(DesP1) <- c("ASV", "Method+Data")

DesP2<-tibble(rownames(deseq_res_P2[[2]]), rep("DESeq2_P2", length(rownames(deseq_res_P2[[2]]))))
colnames(DesP2) <- c("ASV", "Method+Data")

DesP3<-tibble(rownames(deseq_res_P3[[2]]), rep("DESeq2_P3", length(rownames(deseq_res_P3[[2]]))))
colnames(DesP3) <- c("ASV", "Method+Data")

DesAll<-tibble(des_res[,1], rep("DESeq2_All", length(des_res[,1])))
colnames(DesAll) <- c("ASV", "Method+Data")

deseq_allsig<-rbind(DesP1,DesP2, DesP3, DesAll)

#Find ones found multiple times (either in indiviudual timepoints or all together) and add them back in with correct labeling

deseqcount<-plyr::ddply(deseq_allsig, "ASV", transform, count = length(ASV))
duplic<-tibble(deseqcount$ASV[which(deseqcount$count == 2)], deseqcount$Method.Data[which(deseqcount$count == 2)])
a<-tibble(duplic$`deseqcount$ASV[which(deseqcount$count == 2)]`[1], "DESeq2_P3+All")
colnames(a) <- c("ASV", "Method+Data")

b<-tibble(duplic$`deseqcount$ASV[which(deseqcount$count == 2)]`[3], "DESeq2_P1+P3")
colnames(b) <- c("ASV", "Method+Data")

to_add_backin<-rbind(a,b)

deseq_allsig<-deseq_allsig[-which(deseq_allsig$ASV == deseqcount$ASV[which(deseqcount$count == 2)][1]),]
deseq_allsig<-deseq_allsig[-which(deseq_allsig$ASV == deseqcount$ASV[which(deseqcount$count == 2)][3]),]

deseq_allsig<-rbind(deseq_allsig, to_add_backin)
deseq_allsig<-cbind(deseq_allsig, as(tax_table(filtered_ps003)[deseq_allsig$ASV, ], "matrix"))
saveRDS(deseq_allsig, "DESeq/full_res_table_deseq")

deseq_allsig.print <- deseq_allsig
deseq_allsig.print$ASV <- NULL
rownames(deseq_allsig.print) <- NULL
deseq_allsig.print


#MTG results 

MtgP1<-tibble(rownames(zig_res_P1[[2]]), rep("Mtg_P1", length(rownames(zig_res_P1[[2]]))))
colnames(MtgP1) <- c("ASV", "Method+Data")

MtgP2<-tibble(rownames(zig_res_P2[[2]]), rep("Mtg_P2", length(rownames(zig_res_P2[[2]]))))
colnames(MtgP2) <- c("ASV", "Method+Data")

MtgP3<-tibble(rownames(zig_res_P3[[2]]), rep("Mtg_P3", length(rownames(zig_res_P3[[2]]))))
colnames(MtgP3) <- c("ASV", "Method+Data")

MtgAll<-tibble(metag_res[,1], rep("Mtg_All", length(metag_res[,1])))
colnames(MtgAll) <- c("ASV", "Method+Data")

Mtg_allsig<-rbind(MtgP1,MtgP2, MtgP3, MtgAll)

Mtgcount<-plyr::ddply(Mtg_allsig, "ASV", transform, count = length(ASV))
duplic<-tibble(Mtgcount$ASV[which(Mtgcount$count >= 2)], Mtgcount$Method.Data[which(Mtgcount$count >= 2)])


replace <-paste(duplic$`Mtgcount$Method.Data[which(Mtgcount$count >= 2)]`[which(duplic$`Mtgcount$ASV[which(Mtgcount$count >= 2)]` ==  unique(duplic$`Mtgcount$ASV[which(Mtgcount$count >= 2)]`)[1])], collapse = "")
tmp<-tibble(unique(duplic$`Mtgcount$ASV[which(Mtgcount$count >= 2)]`)[1], replace)
colnames(tmp) <- c("ASV", "Method+Data")
replacetab <-tmp

for (i in 2:length(unique(duplic$`Mtgcount$ASV[which(Mtgcount$count >= 2)]`))) {
  replace <-paste(duplic$`Mtgcount$Method.Data[which(Mtgcount$count >= 2)]`[which(duplic$`Mtgcount$ASV[which(Mtgcount$count >= 2)]` ==  unique(duplic$`Mtgcount$ASV[which(Mtgcount$count >= 2)]`)[i])], collapse = "")
  tmp<-tibble(unique(duplic$`Mtgcount$ASV[which(Mtgcount$count >= 2)]`)[i], replace)
  colnames(tmp) <- c("ASV", "Method+Data")
  replacetab<-rbind(replacetab, tmp)
}

Mtg_allsig<-Mtg_allsig[-which(Mtg_allsig$ASV %in% replacetab$ASV),]
Mtg_allsig<-rbind(Mtg_allsig, replacetab)

Mtg_allsig<-cbind(Mtg_allsig, as(tax_table(filtered_ps003)[Mtg_allsig$ASV, ], "matrix"))
saveRDS(Mtg_allsig, "metagenomseq//full_res_table_mtg")

Mtg_allsig.print <- Mtg_allsig
Mtg_allsig.print$ASV <- NULL
rownames(Mtg_allsig.print) <- NULL
Mtg_allsig.print

#None for ANCOM, so will concatenate mtg and des into one table

fullsigtab<-rbind(Mtg_allsig.print, deseq_allsig.print)
fullsigtab_esv<-rbind(Mtg_allsig, deseq_allsig)
fullsigtab
