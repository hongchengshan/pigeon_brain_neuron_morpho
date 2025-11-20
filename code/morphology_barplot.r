library(dplyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(tidyverse)
library(rstatix)
library(scales)

# ---------- Configuration (fixed) ----------
input_file <- "/data/work/pigeon_brain/pca.xlsx"
grouping_candidates <- c("region")   
palette <- c(
  "Hp"    = "#FFA500",
  "HA"    = "#A020F0",
  "MesoP" = "#0000FF",
  "SubP"  = "#FF0000",
  "NCL"   = "#FFC0CB"
)
# -------------------------------------

# Read data & detect grouping column
df_all <- read_excel(input_file) %>% as.data.frame()
grouping_var <- NULL
for(g in grouping_candidates) if(g %in% colnames(df_all)) { grouping_var <- g; break }
if(is.null(grouping_var)) stop("Grouping column 'region' or 'Tissue' not found.")

cat("Using grouping column:", grouping_var, "\n")

# ====== Variable : area / length / Terminal_Pts / bifurcation / branch ======
varname <- "length"
df_raw <- df_all %>%
  select(all_of(grouping_var), value = all_of(varname)) %>%
  filter(!is.na(value)) %>%
  mutate(!!grouping_var := as.factor(.data[[grouping_var]]))

# Checks
if(nrow(df_raw) == 0) stop("length: no data")
if(length(unique(df_raw$value)) <= 1) stop("length: no variance")
df_raw <- df_raw %>%
  mutate(region = factor(region))

# Global ANOVA
anova_res <- df_raw %>%
  anova_test(value ~ region)

# Tukey HSD
posthoc <- df_raw %>%
  tukey_hsd(value ~ region)

# Old versions of rstatix: do not support x= argument — remove it
posthoc <- posthoc %>%
  add_y_position(fun = "max", step.increase = 0.12)

y_max <- max(df_raw$value, na.rm = TRUE)

# Colors
my_cols <- c("Hp" = "#E9BBAF", 
             "HA" = "#BDDEBA",
             "MesoP" = "#B7E0EA",
             "SubP" = "#E98F9E",
             "NCL" = "#FFC0CB")

# Compute means & SD (for barplot)
summary_df <- df_raw %>%
  group_by(region) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE))

# Get upper y-axis limit to avoid overlap with significance marks
y_max <- max(df_raw$value, na.rm = TRUE)

# ---------------- Plot bar chart ----------------
p1 <- ggplot() +
  # Bars
  geom_col(data = summary_df, aes(x = region, y = mean, fill = region), width = 0.7) +
  # Error bars
  geom_errorbar(data = summary_df, aes(x = region, ymin = mean - sd, ymax = mean + sd),
                width = 0.1, linewidth = 0.6) +
  # Jittered raw points
  geom_jitter(data = df_raw, aes(x = region, y = value),
              width = 0.15, height = 0, size = 1.5, alpha = 0.7, inherit.aes = FALSE) +
  scale_fill_manual(values = my_cols) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +

  # Global ANOVA
  stat_compare_means(data = df_raw, aes(x = region, y = value),
                     method = "anova", label = "p.format",
                     label.y = y_max * 1.02) +

  # Pairwise comparison significance (must use posthoc table from df_raw)
  stat_pvalue_manual(
    posthoc,
    label = "p.adj.signif",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    hide.ns = TRUE,
    inherit.aes = FALSE
  )

print(p1)

ggsave("/data/work/pigeon_brain/morpho_length_barplot.pdf", p1, width = 5, height = 8)
