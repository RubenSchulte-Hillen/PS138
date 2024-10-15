

setwd("~/PS138")
load("~/PS138/ampvis_subsets.Rdata")
library(dplyr)
library(tidyr)
library(ggplot2)

# Select the ASVs you want to plot
selected_asvs <- c("asv24", "asv9", "asv13", "asv44", "asv1236")

# Extract abundance data
abundance_data <- ampvis_hel_nofobs$abund %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ASV") %>%
  filter(ASV %in% selected_asvs) %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "Abundance")

# Extract metadata
metadata <- ampvis_hel_nofobs$metadata %>%
  select(sample_id = `CTD #`, depth_category)

# Prepare the data
plot_data <- abundance_data %>%
  left_join(metadata, by = "sample_id") %>%
  mutate(depth_category = factor(depth_category, levels = c(
    "MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "chl-max", "25", "50", 
    "100", "200", "500", "1000", "1500", "2000", "3000", 
    "20m from bottom", "5m from bottom", "bottom"
  )))

# Create the plot
p <- ggplot(plot_data, aes(x = depth_category, y = Abundance, fill = ASV)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = 21, outlier.size = 1.5) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major.x = element_blank()) +
  labs(x = "Depth Category", 
       y = "Abundance (log scale)", 
       title = "Selected ASV Abundances Across Depth Categories",
       fill = "ASV") +
  scale_fill_brewer(palette = "Set1")

# Print the plot
print(p)

plot_data_mean <- plot_data %>%
  group_by(ASV, depth_category) %>%
  summarize(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# Create the line plot
p2 <- ggplot(plot_data_mean, aes(x = depth_category, y = mean_abundance, color = ASV, group = ASV)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major.x = element_blank()) +
  labs(x = "Depth Category", 
       y = "Mean Abundance (log scale)", 
       title = "Selected ASV Abundances Across Depth Categories",
       color = "ASV") +
  scale_color_brewer(palette = "Set1")

# Print the plot
print(p2)