library("Seurat")
library("tidyverse")
library("cowplot")
library("patchwork")

pigeon_brain_obj <- readRDS('/data/input/Files/shanhongcheng/complete_rds/Pigeon.scgcn.rds')

Idents(pigeon_brain_obj) <- pigeon_brain_obj$`mosttype`

main_type_anno <- c( 
               'Astro' = 'Astrocytes', 'CR' = 'Neuron', 'CT_SUB' = 'Neuron', 'DG' = 'Neuroblast',
               'Endo' = 'Endothelial', 'Endothelial_mural' = 'Endothelial_Mural', 'Granule' = 'Granule cell', 'L5_ET' = 'Neuron',
               'L6_IT' = 'Neuron', 'Lamp5' = 'Neuron', 'Micro-PVM' = 'Micro', 'NP_SUB' = 'Neuron',
               'Oligo' = 'Oligodendrocytes', 'OPC' = 'OPCs', 'Pvalb' = 'Neuron', 'Sncg' = 'Neuron',
               'Sst' = 'Neuron', 'Vip' = 'Neuron', 'VLMC' = 'VLMCs'
)
pigeon_brain_obj <- RenameIdents(pigeon_brain_obj, main_type_anno)
pigeon_brain_obj$class_annotation <- Idents(pigeon_brain_obj)

# Define the regions to retain
regions_to_keep <- c("L_Arcopallium", "L_Hyperpallium", "L_Mesorpallium", "L_duanyu", "NCL")

# Extract these regions using subset()
pigeon_Tele_obj <- subset(pigeon_brain_obj, subset = Tissue %in% regions_to_keep)

# Check results
table(pigeon_Tele_obj$Tissue)


cell_to_keep <- c("Neuron", "Neuroblast")

# Extract these cell types using subset()
pigeon_Tele_neu <- subset(pigeon_Tele_obj, subset = class_annotation %in% cell_to_keep)

# Check results
table(pigeon_Tele_neu$class_annotation)

Idents(pigeon_Tele_neu) <- pigeon_Tele_neu$`Tissue`
pigeon_Tele_neu.markers <- FindAllMarkers(pigeon_Tele_neu, only.pos = TRUE)
markers <- pigeon_Tele_neu.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    dplyr::filter(p_val_adj < 0.05)
write.csv(markers, file = "/data/work/pigeon_brain/pigeon_Tele_neu_Region_DEG.csv")

morpho_genes <- c("NPY",
                  "NR2F2", "SATB2", "NRN1", "TUBB3", "MAP2",
                  "EPHA5", "SLIT2",
                  "ZBTB18", 
                  "DCX"
)

pigeon_Tele_neu <- AddModuleScore(pigeon_Tele_neu, features = list(morpho_genes),
                                  name = "morpho_genes_Score", ctrl = 100)
p1 <- ggplot()

# Get data for plotting
vln_data <- FetchData(pigeon_Tele_obj, vars = c("Tissue", "morpho_genes_Score1"))
vln_data$morpho_genes_Score1_normalized <- (vln_data$morpho_genes_Score1 - min(vln_data$morpho_genes_Score1)) / 
                                          (max(vln_data$morpho_genes_Score1) - min(vln_data$morpho_genes_Score1))

library(ggsignif)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)

vln_data <- vln_data %>%
  mutate(Tissue = factor(Tissue))

# Global ANOVA
anova_res <- vln_data %>%
  anova_test(morpho_genes_Score1_normalized ~ Tissue)

# Tukey HSD
posthoc <- vln_data %>%
  tukey_hsd(morpho_genes_Score1_normalized ~ Tissue)

# Old versions of rstatix: do not support x= argument — remove it
posthoc <- posthoc %>%
  add_y_position(fun = "max", step.increase = 0.12)

y_max <- max(vln_data$morpho_genes_Score1_normalized, na.rm = TRUE)


p1 <- ggplot(vln_data, aes(x = Tissue, y = morpho_genes_Score1_normalized, fill = Tissue)) +
  geom_violin(width = 1, trim = FALSE) +
  scale_fill_manual(values = c("L_Arcopallium" = "#FFA500", "L_Hyperpallium" = "#A020F0",
                               "L_Mesorpallium" = "#0000FF", "L_duanyu" = "#FF0000", "NCL" = "#FFC0CB")) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  
  # Global ANOVA
  stat_compare_means(method = "anova", label = "p.format",
                     label.y = y_max * 1.10) +
  
  # Pairwise comparisons (Tukey)
  stat_pvalue_manual(
    posthoc,
    label = "p.adj.signif",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    hide.ns = TRUE,
    inherit.aes = FALSE   # Important: do not inherit global aes(Tissue, y = ...)
)

print(p1)
ggsave("/data/work/pigeon_brain/morpho_genes_Score1_normalized.pdf", width = 6, height = 5)
