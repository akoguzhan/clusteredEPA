# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)

# Input data
df <- data.frame(
  Tobs = c(20, 50, 100, 200),
  known = c(0.045,
            0.033,
            0.027,
            0.026
  ),
  known_cond = c(0.059,
                 0.043,
                 0.036,
                 0.043
  ),
  split = c(0.151,
            0.154,
            0.179,
            0.381
  ),
  split_cond = c(0.134,
                 0.138,
                 0.194,
                 0.348
  ),
  selective = c(0.049,
                0.039,
                0.038,
                0.048
  ),
  selective_cond = c(0.049,
                     0.035,
                     0.037,
                     0.045
  )
)

# Prepare tidy data
df_uncond <- df %>%
  select(Tobs, known, split, selective) %>%
  pivot_longer(cols = -Tobs, names_to = "Test", values_to = "RejectionRate")

df_cond <- df %>%
  select(Tobs, known_cond, split_cond, selective_cond) %>%
  rename_with(~gsub("_cond", "", .x)) %>%
  pivot_longer(cols = -Tobs, names_to = "Test", values_to = "RejectionRate")

# Define color and shape palettes
color_palette <- c("known" = "#1b9e77", "split" = "#d95f02", "selective" = "#7570b3")
shape_palette <- c("known" = 16, "split" = 17, "selective" = 15)  # circle, triangle, square
label_names <- c("known" = "Predetermined", "split" = "Split Sample", "selective" = "Selective Inference")

# Plot for Unconditional
p_uncond <- ggplot(df_uncond, aes(x = Tobs, y = RejectionRate, color = Test, shape = Test)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette, labels = label_names) +
  scale_shape_manual(values = shape_palette, labels = label_names) +
  labs(x = "T", y = "Rejection Rate") +
  ylim(0, 0.4) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Plot for Conditional
p_cond <- ggplot(df_cond, aes(x = Tobs, y = RejectionRate, color = Test, shape = Test)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = color_palette, labels = label_names) +
  scale_shape_manual(values = shape_palette, labels = label_names) +
  labs(x = "T", y = "Rejection Rate") +
  ylim(0, 0.4) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Export each plot to PDF
ggsave(here("Tools", "clusteredEPA", "Results", "unconditional_CEPA_size_sbreaks.pdf"), plot = p_uncond, width = 5.5, height = 4)
ggsave(here("Tools", "clusteredEPA", "Results", "conditional_CEPA_size_sbreaks.pdf"), plot = p_cond, width = 5.5, height = 4)
