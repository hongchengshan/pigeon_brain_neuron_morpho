# --------------- Required packages ---------------
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)  # Can be commented out if point labels are not needed

# --------------- Read and rename columns (according to your actual column names) ---------------
df <- read_excel("col_analysis.xlsx") %>%
  rename(
    group = `...1`,         # first column
    x = area,               # second column
    y = `gene module`       # third column (use backticks because of space)
  )

# --------------- Transformation (z-score) ---------------
df <- df %>%
  mutate(
    x_t = as.numeric(scale(x)[,1]),
    y_t = as.numeric(scale(y)[,1])
  )

# --------------- Fit overall linear model (based on transformed data) ---------------
lm_mod <- lm(y_t ~ x_t, data = df)
sm <- summary(lm_mod)

slope <- coef(lm_mod)["x_t"]
intercept <- coef(lm_mod)["(Intercept)"]
r2 <- sm$r.squared
p_slope <- coef(sm)["x_t", "Pr(>|t|)"]

# Format the text to display
label_text <- sprintf("y = %.3fx + %.3f\nR² = %.3f, p = %.3g",
                      slope, intercept, r2, p_slope)

# Compute annotation position (place it near the top-left inside the plot)
x_pos <- min(df$x_t, na.rm = TRUE) + 0.05 * diff(range(df$x_t, na.rm = TRUE))
y_pos <- max(df$y_t, na.rm = TRUE) - 0.05 * diff(range(df$y_t, na.rm = TRUE))

# --------------- Plot: points colored by group, with overall regression line (black) ---------------
p <- ggplot(df, aes(x = x_t, y = y_t)) +
  # Group points: different colors
  geom_point(aes(color = factor(group)), size = 3, alpha = 0.9) +
  # (Optional) Fit a separate dashed regression line for each group: uncomment below
  # geom_smooth(aes(group = factor(group), color = factor(group)), method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.8, inherit.aes = TRUE) +
  # Overall regression line (must use data=df to avoid per-group lines)
  geom_smooth(data = df, aes(x = x_t, y = y_t), method = "lm",
              se = TRUE, color = "black", fill = "grey80", linewidth = 1.1, inherit.aes = FALSE) +
  # Point labels (optional)
  # geom_text_repel(aes(label = group), size = 3.5) +
  # Annotate regression equation / R2 / p-value
  annotate("text", x = x_pos, y = y_pos, label = label_text, hjust = 0, vjust = 1, size = 4.5) +
  theme_classic(base_size = 14) +
  labs(
    x = "x (z-scored)",
    y = "y (z-scored)",
    color = "Group",
    title = "Scatter plot with overall regression line and group colors"
  ) +
  theme(legend.position = "right") +
  scale_color_manual(values = c("#A020F0","#FFA500","#0000FF","#FFC0CB","#FF0000"))

# Display and save
print(p)
ggsave("regression_overall_line.pdf", p, width = 7, height = 5)
