library(DECIPHER)
library(phangorn)
library(raxml)
setwd("C:/Users/kitikomp/Documents/Lab/M3_Phyloseq_Analysis/scripts")


sequences <- rownames(ps@otu_table)
sequences <- gsub("[^ATCG]", "", sequences)
sequences <- unique(sequences)

alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

gtr <- raxml(phang.align, 
             m="GTRGAMMAIX", # model
             f="a", # best tree and bootstrap
             p=1234, # random number seed
             x=2345, # random seed for rapid bootstrapping
             N=2, # number of bootstrap replicates
             file="alignment", # name of output files )
)


dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR$tree, save_file_name)

save_file_name <- "../Tree/ps_no_norm_age_bf_file_complete_family_tree.rds"
ps_not_norm_comp <- readRDS(paste0("../data/ps_no_norm_age_bf_filt_complete_family.rds"))


tree <- readRDS(file ="../Tree/tree_root")
sequences %in% tree$tip.label
