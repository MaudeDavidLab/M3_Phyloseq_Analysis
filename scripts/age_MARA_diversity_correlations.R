library(vegan)
library(phyloseq)
library(DT)
library(ggplot2)
library(data.table)
library(pander)
library(picante)

output_results = "../results/"
ps_deseq <- readRDS(paste0(output_results, "Filtered/ps_deseq_friedfilt.rds"))
ps_css <- readRDS(paste0(output_results,"Filtered/ps_css_friedfilt.rds"))
ps_no_norm <- readRDS(paste0(output_results,"Filtered/ps_no_norm_friedfilt.rds"))
#ps_use <- readRDS(paste0("../data/ps_no_norm_age_bf_filt_complete_family.rds"))

ps_use <- readRDS("../data/ps_nonorm_nobreastfed_fullseq.rds")

#Plot alpha diversity comparison between groups
ER <- estimate_richness(ps_not_norm_comp, measures=c("Observed", "Chao1", "Shannon"))
ER <- cbind(ER, sample_data(ps_not_norm_comp)[row.names(ER), c("phenotype", "Family.group.ID", "Within.study.sampling.date")])
faiths <- pd(t(ps_use@otu_table), ps_use@phy_tree, include.root = T)
ER$Faith <- as.numeric(faiths$PD)

ER <- data.table(ER, keep.rownames = TRUE)
ER_long <- melt(ER, id.vars=c('rn', 'phenotype', "Family.group.ID", "Within.study.sampling.date"))

# plot
ggplot(data=ER_long[variable!='se.chao1'], aes(x=phenotype, y=value, fill=phenotype))+
  geom_boxplot(width=0.7, outlier.colour='white')+
  geom_jitter(size=1, position=position_jitter(width=0.1))+
  xlab('')+ylab('')+
  scale_fill_manual(values=sgColorPalette)+
  theme_minimal()+facet_wrap(~variable, scales='free')

# run t-test to check significance
ttest=NULL
for(alphad in c('Observed', 'Chao1', 'Shannon', 'Faith')){
  ttest_res=t.test(value ~ phenotype, data=ER_long[ER_long$variable==alphad], var.equal=TRUE)
  ttest=rbindlist(list(ttest, data.table(alpha_index=alphad, pvalue=ttest_res$p.value)))
}

pander(ttest)







# Correlate alpha diversity with age

ER <- estimate_richness(ps_use, measures=c("Observed", "Chao1", "Shannon"))
age <- ps_use@sam_data$Age..months.


df <- data.frame(shannon = ER$Shannon, chao = ER$Chao1, faithsPD = faiths, age = age/12, phenotype = ps_use@sam_data$phenotype, MARA = ps_use@sam_data$Mobile.Autism.Risk.Assessment.Score)
df <- df[df$phenotype == "A", ]



color = "#2E4A76"
text_size = 15
size = 5
my_theme =theme(text= element_text(size = text_size),
                      axis.text.x = element_text(size = text_size),
                      axis.text.y = element_text(size = text_size),
                      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                      strip.text = element_text(size = text_size))

test <- cor.test(df$MARA, df$chao, method = "spearman")
p_mara_chao <- ggplot(data = df, aes(x = MARA, y = chao), color = color) +
  my_theme +
  ylab("chao1") +
  ylim(0, 450)+
  geom_point() +
  geom_smooth(method = 'lm')+
  geom_text(label = paste0("R^2 = ", round(test$estimate, 3), "  p = ", round(test$p.value, 5)), x = -4.75, y = 440, color = color, size = size)
  
  
test <- cor.test(df$MARA, df$shannon, method = "spearman")          
p_mara_shannon <- ggplot(data = df, aes(x = MARA, y = shannon), color = color) +
  my_theme +
  ylim(2, 5)+
  geom_point() +
  geom_smooth(method = 'lm')+
  geom_text(label = paste0("R^2 = ", round(test$estimate, 3), "  p = ", round(test$p.value, 6)), x = -5.2, y = 5, color = color, size = size)

test <- cor.test(df$MARA, df$faithsPD.PD, method = "spearman")          
p_mara_faiths <- ggplot(data = df, aes(x = MARA, y = faithsPD.PD), color = color) +
  my_theme +
  ylim(9, 45)+
  ylab("Faith's PD")+
  geom_point() +
  geom_smooth(method = 'lm')+
  geom_text(label = paste0("R^2 = ", round(test$estimate, 3), "  p = ", round(test$p.value, 6)), x = -4.3, y = 43, color = color, size = size)

test <- cor.test(df$chao, df$age, method = "spearman") 
p_age_chao <- ggplot(data = df, aes(x = age, y = chao), color = color) + 
  my_theme +
  ylim(0, 450)+
  ylab("chao1") +
  geom_point() +
  geom_smooth(method = 'lm')+
  geom_text(label = paste0("R^2 = ", round(test$estimate, 3), "  p = ", round(test$p.value, 5)), x = 4.7, y = 450, color = color, size = size)

test <- cor.test(df$shannon, df$age, method = "spearman") 
p_age_shannon <- ggplot(data = df, aes(x = age, y = shannon), color = color) + 
  my_theme +
  ylim(2, 5)+
  geom_point() +
  geom_smooth(method = 'lm')+
  geom_text(label = paste0("R^2 = ", round(test$estimate, 3), "  p = ", round(test$p.value, 5)), x = 4.5, y = 4.9, color = color, size = size)

test <- cor.test(df$faithsPD.PD, df$age, method = "spearman") 
p_age_faiths <- ggplot(data = df, aes(x = age, y = faithsPD.PD), color = color) + 
  my_theme +
  ylim(9, 45) +
  ylab("Faith's PD") +
  geom_point() +
  geom_smooth(method = 'lm')+
  geom_text(label = paste0("R^2 = ", round(test$estimate, 3), "  p = ", round(test$p.value, 5)), x = 4.5, y = 43, color = color, size = size)


test <- cor.test(df$MARA, df$age, method = "spearman")          
p_age_mara <- ggplot(data = df, aes(x = age, y = MARA), color = color) + 
  my_theme +
  ylim(-10, 5)+
  geom_point() +
  geom_smooth(method = 'lm')+
  geom_text(label = paste0("R^2 = ", round(test$estimate, 3), "  p = ", round(test$p.value, 5)), x = 4.7, y = 5, color = color, size = size)

#test <- cor.test(df$chao, df$shannon, method = "spearman") 
#p6 <- ggplot(data = df, aes(x = chao, y = shannon), color = color) + 
#  my_theme +
#  ylim(2, 5)+
##  xlab("chao1") +
#  geom_point() +
#  geom_smooth(method = 'lm')+
#  geom_text(label = paste0("R^2 = ", round(test$estimate, 3), "  p = ", round(test$p.value, 6)), x = 165, y = 5, color = color, size = size)


library(cowplot)
p <- plot_grid(p_mara_chao, p_mara_shannon, p_mara_faiths, p_age_chao, p_age_shannon, p_age_faiths, p_age_mara, labels = c("A", "B", "C", "D", "E", "F", "G"))
ggsave("../results/diversity_correlations.png", device = "png", plot = p, width = 12, height = 7, unit = "in", dpi = 400)




