tmp <- readRDS("../data/ps_no_norm_age_bf_filt_complete_family.rds")

hist(colSums(tmp@otu_table), breaks = 15)
median(colSums(tmp@otu_table))
mean(colSums(tmp@otu_table))

color <- c("#01B5BB","#84CF04")
color_use <- color[ifelse(tmp@sam_data$phenotype == "A", 1, 2)]




text_size = 15
size = 5
my_theme =theme(text= element_text(size = text_size),
                axis.text.x = element_text(size = text_size),
                axis.text.y = element_text(size = text_size),
                axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                strip.text = element_text(size = text_size))

p <- rarecurve(t(tmp@otu_table), step = 10000, col = color_use, label = F)

ggsave("../results/rarefaction_curve.pdf", device = "pdf", plot = p, width = 5, height = 5, unit = "in", dpi = 400)