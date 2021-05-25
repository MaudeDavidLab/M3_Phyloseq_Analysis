sig <- read.csv("~/Lab/M3_Phyloseq_Analysis/Significant_res_deseq_q0.05_mtgseq_q0.05.csv")
single_timepoints <- sig[is.na(sig$deseq_timeseries_adjp) & any(!is.na(sig$deseq_P1_adjp), !is.na(sig$deseq_P2_adjp), !is.na(sig$deseq_P3_adjp)), ]
single_timepoints <- single_timepoints[is.na(single_timepoints$mtgseq_timeseries_adjp), ]


single_timepoints <- data.frame(ps_deseq@tax_table)[single_timepoints$rn, ]

write.csv(single_timepoints, "../results/ASVs_single_timepoint.csv", quote = F)