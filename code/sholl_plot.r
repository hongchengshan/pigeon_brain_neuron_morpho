library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(scales)

input_file <- "/data/work/pigeon_brain/sholl/sholl_Hp.xlsx"

input_file2 <- "/data/work/pigeon_brain/sholl/sholl_NCL.xlsx"

input_file3 <- "/data/work/pigeon_brain/sholl/sholl_HA.xlsx"

input_file4 <- "/data/work/pigeon_brain/sholl/sholl-SubP.xlsx"

input_file5 <- "/data/work/pigeon_brain/sholl/sholl_MesoP.xlsx"

df_all <- read_excel(input_file) %>% as.data.frame()

df_all <- read_excel(input_file2) %>% as.data.frame()

df_all <- read_excel(input_file3) %>% as.data.frame()

df_all <- read_excel(input_file4) %>% as.data.frame()

df_all <- read_excel(input_file5) %>% as.data.frame()

library(tidyr)
df_long <- pivot_longer(
  df_all,
  cols = -Radius,
  names_to = "Column",
  values_to = "Value"
)

# 查看结果
print(head(df_long, 20))

df_long_with_region1 <- df_long %>%
  mutate(region = "Hp")

df_long_with_region2 <- df_long %>%
  mutate(region = "NCL")

df_long_with_region3 <- df_long %>%
  mutate(region = "HA")

df_long_with_region4 <- df_long %>%
  mutate(region = "SubP")

df_long_with_region5 <- df_long %>%
  mutate(region = "MesoP")

combined_df <- bind_rows(df_long_with_region, df_long_with_region2,df_long_with_region3,df_long_with_region4,df_long_with_region5)

# check
combined_df

total_neurons <- 64
summary_group <- combined_df %>%
  group_by(region, Radius) %>%
  summarise(
    total_neurons = n(),                        # total number of neurons in each group
    n_reaching = sum(Value > 0),                # number of neurons reaching this radius
    mean_including0 = mean(Value),              # mean including zeros
    sem_including0 = sd(Value) / sqrt(total_neurons),
    mean_excluding0 = ifelse(n_reaching>0, mean(Value[Value>0]), NA),
    sem_excluding0 = ifelse(n_reaching>1, sd(Value[Value>0]) / sqrt(n_reaching), NA)
  ) %>% ungroup()
palette <- c(
 "Hp" = "#E9BBAF", 
             "HA" = "#BDDEBA",
             "MesoP" = "#B7E0EA",
             "SubP" = "#E98F9E",
             "NCL" = "#FFC0CB"
)

p <- ggplot(summary_group, aes(x = Radius, y = mean_including0, color = region)) +
  geom_ribbon(aes(ymin = mean_including0 - sem_including0,
                  ymax = mean_including0 + sem_including0,
                  fill = region),
              alpha = 0.25, color = NA) +
  geom_line(size = 1) +
  labs(x = "Radius (µm)", y = "Mean crossings",
       title = "Sholl Analysis (Group-wise Mean ± SEM)") +
  theme_classic() +
scale_color_manual(values = palette) +
scale_fill_manual(values = palette)
print(p)
ggsave("/data/work/pigeon_brain/morpho_sholl.pdf", p, width = 8, height = 5)