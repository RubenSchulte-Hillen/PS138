library(dplyr)
library(tidyr)
load( "ampvis_subsets.Rdata")

tax_table <- ampvis.bac_nofobs$tax

# Function to calculate percentage of unclassified taxa
calc_unclassified_percentage <- function(column) {
  total <- length(column)
  unclassified <- sum(grepl("uc$", column, ignore.case = TRUE))
  return((unclassified / total) * 100)
}

# Calculate percentages for each taxonomic level
unclassified_percentages <- sapply(tax_table, calc_unclassified_percentage)

# Create a data frame with the results
result_df <- data.frame(
  Taxonomic_Level = names(unclassified_percentages),
  Unclassified_Percentage = round(unclassified_percentages, 1)  # Round to 1 decimal place
) %>%
  filter(Taxonomic_Level != "OTU")  # Remove the OTU row
# Define the correct order of taxonomic levels
tax_order <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Reorder the data frame
result_df <- result_df %>%
  mutate(Taxonomic_Level = factor(Taxonomic_Level, levels = tax_order)) %>%
  arrange(Taxonomic_Level)

print(result_df)
# Optionally, you can create a bar plot to visualize the results
library(ggplot2)

ggplot(result_df, aes(x = Taxonomic_Level, y = Unclassified_Percentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(Unclassified_Percentage, "%")), 
            vjust = -0.5, color = "black", size = 3.5) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(title = "Percentage of Unclassified ASVs at Each Taxonomic Level",
       x = "Taxonomic Level",
       y = "Percentage Unclassified") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("unclassified_taxa_percentages.pdf", width = 10, height = 6)


