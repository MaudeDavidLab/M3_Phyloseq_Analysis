library(caret)
library(randomForest)
library(tibble)
library(ROCR)
library(dplyr)
library(reshape2)
library(phyloseq)
library(glmnet)
library(stringr)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggplot2)
library(pROC)
library(e1071)


#Write fasta file using only the 11 biomarkers from the longitudinal study
sig_taxa <- read.delim("~/Lab/M3_Phyloseq_Analysis/results/significant_taxa_timeseries.tsv")
biomarkers <- sig_taxa$asv[grepl(",", sig_taxa$method)][1:10] #11th is the t__186843 that we can't find a corresponding ASV to

headers <- paste0(">ASV", seq(1, 10), "\n")
seqs <- biomarkers
fasta <- paste0(headers, seqs, collapse = "\n")
write.table(fasta, "~/Lab/M3_Phyloseq_Analysis/data/11_biomarker.fasta")



#Write fasta of all sequences in pilot
ps <- readRDS("~/Lab/Microbiome_16S_mbio-master/data/ps_noDuplicates.RDS")
seqs <- taxa_names(ps)
headers <- paste0(">ASV", seq(1, length(seqs)), "\n")
fasta <- paste0(headers, seqs, collapse = "\n")
write.table(fasta, "~/Lab/M3_Phyloseq_Analysis/data/pilot_seqs.fasta")


#Keep only biomarkers in phyloseq object
blast_hits <- read.delim("~/Lab/M3_Phyloseq_Analysis/data/pilot/pilot_11biomarkers_100iden.tsv", header=FALSE)
asvids <- paste0("ASV", seq(1, ntaxa(ps)))
pred_sequences <- taxa_names(ps)[asvids %in% blast_hits$V1]


#Normalize taxa counts
deSeqNorm <- function(ps){
  library(DESeq2)
  ps_dds <- phyloseq_to_deseq2(ps, ~ Treatment )
  ps_dds <- estimateSizeFactors(ps_dds, type = "poscounts")
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  abund <- abund + abs(min(abund)) #don't allow deseq to return negative counts
  ps_deSeq <- phyloseq(otu_table(abund, taxa_are_rows = T), sample_data(ps), tax_table(ps), phy_tree(ps))
  return(ps_deSeq)
}


filterTaxaByPrevolence <- function(ps, percentSamplesPresentIn){
  prevalenceThreshold <- percentSamplesPresentIn * nsamples(ps)
  toKeep <- apply(data.frame(otu_table(ps)), 1, function(taxa) return(sum(taxa > 0) > prevalenceThreshold))
  ps_filt <- prune_taxa(toKeep, ps)
  return(ps_filt)
}

prevFiltThresh = .03
ps_filt <- filterTaxaByPrevolence(ps, prevFiltThresh)
ps <- deSeqNorm(ps_filt)




#Make fields numeric
keep_variables <- c("Multivitamin", "Probiotic", "VIT.B", "Vit.D", "OtherSupp", "bowel.freq", "bowel.qual", "C.section",
                    "WholeGrain", "Fermented.veg", "dairy", "Fruit", "Veg", "Meat.eggs", "OliveOil", "Milk..Cheese",
                    "Frozen.D", "Red.Meat", "Highfat.red.meat", "ChickenTurkey", "Seafood", "Sugar", "SweetBev", "Dog", "Cat", "Pair", "Treatment")

nas <- apply(ps@sam_data, 1, function(x) return(sum(is.na(x))))
keep <- nas < 20
ps <- prune_samples(keep, ps)

df_asv <- t(ps@otu_table)
df_metadata <- ps@sam_data[ , keep_variables]
df_metadata$Multivitamin <- ifelse(df_metadata$Multivitamin == "A", 0, 1)
df_metadata$OtherSupp <- ifelse(df_metadata$OtherSupp == "A", 0, 1)
df_metadata$bowel.qual <- ifelse(df_metadata$bowel.qual == "A", 0, 1)
df_metadata$C.section <- ifelse(df_metadata$C.section == "A", 0, 1)
df_metadata$dairy <- ifelse(df_metadata$dairy == "Dairy_int", 0, 1)
df_metadata$Dog <- ifelse(df_metadata$Dog == "A", 0, 1)
df_metadata$Cat <- ifelse(df_metadata$Cat == "A", 0, 1)
df_metadata$Treatment <- as.factor(df_metadata$Treatment)






rand_forest <- function(pred_sequences, training_ids, ps, df_metadata, include_taxa = T, include_meta = F){
  
  phen <- df_metadata$Treatment
  fam_id <- df_metadata$Pair
  data <- NULL
  if(include_taxa){
    ps <- prune_taxa(pred_sequences, ps )
    data <- t(otu_table(ps))
    
    #taxa names for plot labels
    species_names <- as.character(ps@tax_table[, 7])
    genus_names <- as.character(ps@tax_table[, 6])
    family_names <- as.character(ps@tax_table[, 5])
    taxa <- species_names
    #taxa[taxa == "s__unclassified"] = genus_names[taxa == "s__unclassified"]
    #taxa[taxa == "g__unclassified"] = family_names[taxa == "g__unclassified"]
    colnames(data) <- taxa
    
  }
  if(include_meta){
    Lifestyle_Variables <- data.frame(df_metadata)
    Lifestyle_Variables <- Lifestyle_Variables %>% select(-c(Treatment))

    if(sum(is.na(ps@sam_data))){
      Lifestyle_Variables <- rfImpute(Pair~., data = Lifestyle_Variables)
    }
    
    
    if(is.null(data)){
      data <- Lifestyle_Variables
    }else{
      data <- data.frame(cbind(data, Lifestyle_Variables))
    }
  }
  data <- data.frame(phen, data)
  data$phen <- ordered(data$phen, levels = c("Control", "Aut"))
  
  validate <- data[-training_ids,]
  training <- data[training_ids,]
  
  
  #to pick best parameters
  control <- trainControl(method='repeatedcv', 
                          number=3, 
                          repeats=3)
  
  tree_depths <- round(seq(2, ncol(data), by = 2))
  tunegrid <- expand.grid(.mtry= tree_depths) #mtry is the depth of each decision tree
  
  rf <- train(phen ~., 
              data = training, 
              method='rf', 
              metric='Accuracy', 
              tuneGrid=tunegrid, 
              trControl=control)
  
  
  mtry_best = as.numeric(rf$bestTune)
  print(paste0("Tree depth: ", mtry_best))
  
  pred_votes <- c()
  outputlist <- list()
  for(i in 1:2){
    set.seed(i)
    AR.classify <- randomForest(phen~., data = training, ntree = 128, mtry = mtry_best, importance = TRUE)
    rf <- AR.classify
    OOB.votes <- predict(rf, validate[,-1], type="prob")
    
    pred_votes <-  OOB.votes[,2]
    pred.obj <- prediction(pred_votes, data[names(pred_votes), ]$phen)
    
    tp <- (pred_votes > 0.5) & data[names(pred_votes), ]$phen == "Aut"
    
    #roc (tpr / fpr) perforamnce 
    roc_perf <- performance(pred.obj,"tpr", "fpr") #Calculate the AUC value
    auc_roc <- performance(pred.obj, "auc")@y.values[[1]]
    print(auc_roc)
    
    
    perf_obj <-list()
    perf_obj[[1]] <- auc_roc
    perf_obj[[2]] <- roc_perf
    perf_obj[[3]] <- rf
    #labels_ordered<-ordered(data[names(pred_votes), ]$phen, levels = c("Control", "Aut"))
    #perf_obj[[4]] <- ci.auc(labels_ordered, pred_votes)
    outputlist[[i]] <- perf_obj
  }
  return(outputlist)
}



getPerformance <- function(ps, df_metadata, include_taxa, include_meta, pred_sequences = NA, num_folds = 7, seed = 0){
  values <- NULL
  perf_objs <- list()
  seed=3
  
  num_folds <- 10
  set.seed(seed)
  folds = num_folds
  fam_id <- sample_data(ps)$Pair
  folds_by_family <- groupKFold(fam_id, folds)
  
  
  for(i in 1:folds){
    print(i)
    perf_obj <- rand_forest(pred_sequences = pred_sequences, training_ids = folds_by_family[[i]], ps = ps, df_metadata = df_metadata, include_taxa = include_taxa, include_meta = include_meta)
    i = i + 1
    perf_objs[[i]] <- perf_obj
    tmp_values <- unlist(lapply(perf_obj, function(x) return(x[[1]])))
    print(paste("Mean values across forests with same parameters: ", mean(tmp_values)))
    print(paste("SD values across forests with same parameters: ", sd(tmp_values)))
    values <- c(values, mean(tmp_values))
  }
  return(list(values, perf_objs))
}


# Metadata only
tmp <- getPerformance(ps, df_metadata, include_taxa = F, include_meta = T, pred_sequences = NA)
meta_values <- tmp[[1]]
meta_perf_objs <- tmp[[2]]
metadf<-data.frame('Lifestyle_Variables_only' = meta_values) #,  'meta_null_dist_no_db_sig' = meta_null_nodoublesig)
metadf<-melt(metadf)

summary(meta_values)
saveRDS(metadf, file = "results/pilot_ml/mlmetadf.rds")
saveRDS(meta_perf_objs, file = "results/pilot_ml/meta_perf_objs.rds")


# Metadata + biomarkers
tmp <- getPerformance(ps, df_metadata, include_taxa = T, include_meta = T, pred_sequences = pred_sequences)
biomarker_values <- tmp[[1]]
biomarker_perf_objs <- tmp[[2]]
biodf <- data.frame('Lifestyle_Variables_and_11_Biomarkers' = biomarker_values)
biodf <- melt(biodf)

summary(biomarker_values)
saveRDS(biodf, file = "results/pilot_ml/biodf.rds")
saveRDS(biomarker_perf_objs, file = "results/pilot_ml/biomarker_perf_objs.rds")


print(meta_values)
print(biomarker_values)
wilcox.test(meta_values, biomarker_values, paired = T)
boxplot(meta_values, biomarker_values)

#maybe instead of biomarkers exactly, go for those of the same genus
#mandate exact 100 match




#Biomarker only
tmp <- getPerformance(ps, df_metadata, include_taxa = T, include_meta = F, pred_sequences = pred_sequences)
bio_only_values <- tmp[[1]]
bio_only_perf_objs <- tmp[[2]]

bio_only_df <- data.frame('Biomarkers_only' = bio_only_values )
bio_only_df <- melt(bio_only_df)


saveRDS(bio_only_df, file = "results/pilot_ml/bio_only_df.rds")
saveRDS(bio_only_perf_objs, file = "results/pilot_ml/bio_only_perf_objs.rds")

print(meta_values)
print(biomarker_values)
print(bio_only_values )
wilcox.test(meta_values, bio_only_values, paired = T)





## 11 random raxa
rand_values <- c()
for(i in 1:5){
  pred_sequences_rand <- sample(taxa_names(ps)[-which(taxa_names(ps) %in% pred_sequences)], length(pred_sequences))
  tmp <- getPerformance(ps, df_metadata, include_taxa = T, include_meta = F, pred_sequences = pred_sequences_rand)
  tmp_values <- tmp[[1]]
  rand_values <- c(rand_values, tmp_values)
}

taxa_rand_df <- data.frame('Random_taxa' = rand_values )
taxa_rand_df <- melt(taxa_rand_df)

saveRDS(taxa_rand_df, file = "../results/ml/taxa_rand_df.rds")









