# Load required libraries
library(dplyr)
library(tibble)
library(mixOmics)
library(cluster)
library(proxy)
library(reshape2)
library(gtools)
library(tidyr)
library(ShortRead)
library(ape)
library(ampvis2)

# Load the saved objects
load("ampvis_subsets.Rdata")

# Define a list of object names
ampvis.bac_CTD <- amp_filter_samples(ampvis.bac_nofobs,
                                     gear %in% c("CTD"),
                                     normalise = FALSE )

ampvis.bac_ICE <- amp_filter_samples(ampvis.bac_nofobs,
                                     gear %in% c("si_corer_9cm", "hand_pump"),
                                     normalise = FALSE )

transect_amp <- amp_filter_samples(transect_amp,
                                     gear %in% c("CTD"),
                                     normalise = FALSE )

station_amp <- amp_filter_samples(station_amp,
                                     gear %in% c("CTD"),
                                     normalise = FALSE )
#for different subsets change the following dependend on  salinity and temprature data is available
#ice samples didnt have temp and sal data during MSc
# line : 276 and 288 in or exclude salinity and temp into analyisis
#depending on how many clusters: choose : line 105 if k=(num of branches) or h=(dist of branches) 

object_names <- c("ampvis.bac_nofobs")
object_names <- c("transect_amp")#,"station_amp")#,"ampvis.bac_CTD")
object_names <- c("ampvis.bac_ICE")#,"ampvis.bac_CTD",)
object_names <- c("ampvis.bac_ICE")#,"ampvis.bac_CTD",)
object_names <- c("ampvis.bac_CTD")#,"ampvis.bac_CTD",)
obj_name <- "transect_amp"
obj_name <- "station_amp"
obj_name <- "ampvis.bac_CTD"
dev.off()
# Loop through each object
for (obj_name in object_names) {
  # Create a new folder to save the plots and data
  plot_folder <- paste0("sPLS_", obj_name, "_Final_final")
  dir.create(plot_folder, showWarnings = FALSE)
  
 
  
  # Extract metadata and abundance data
  ENV <- get(obj_name)$metadata
  ASV <- get(obj_name)$abund
  TAX <- get(obj_name)$tax
  ENV$depth_category[ENV$depth_category == "25"] <- "chl-max"
  # Match rows
  asv <- t(ASV) %>% as.data.frame()
  '  ENV$`depth [m]` <- replace(ENV$`depth [m]`, ENV$depth_category == "Ice-TS", -2)
  ENV$`depth [m]` <- replace(ENV$`depth [m]`, ENV$depth_category == "Ice-BS", -1)
  ENV$`depth [m]`[ENV$`depth [m]` == "chlmax"] <- "25"
  ENV$`depth [m]`[ENV$`depth [m]` == "bottom"] <- "4000"
  ENV$`depth [m]`[ENV$`depth [m]` == "4216-9"] <- "4217"
  ENV$`depth [m]`[ENV$`depth [m]` == "4233-4"] <- "4233"
  ENV$`depth [m]`[ENV$`depth [m]` == "3873.6"] <- "3874"
  ENV$`depth [m]`[ENV$`depth [m]` == "under-ice water"] <- "2"
  ENV$`depth [m]`[ENV$`depth [m]` == "meltpond"] <- "-3"
  ENV$`depth [m]`[ENV$`depth [m]` == "ice-edge"] <- "0"
  ENV$`depth [m]`[ENV$`depth [m]` == "3873.6"] <- "3874"
  unique(ENV$`depth [m]`)
  '
  meta <- ENV %>%
    remove_rownames %>%
    column_to_rownames("CTD #") %>%
    dplyr::select(c("lon", "lat", "NO3+NO2_c_mean",
                    "PO4_c_mean", "Si(OH)4_c_mean",
                    "NO2_c_mean", "NH4_c_mean",
                    "TDN_c_mean", "TDP_c_mean",
                    "Sal00" ,"Fix_Temp_C",))# "depth [m]"))
  
  asv <- asv[row.names(meta), ]
  ' meta$`depth [m]` <- as.numeric(meta$`depth [m]`)'
  
  # Calculate sPLS
  PLS <- spls(asv, meta, ncomp = 3, tol = 1e-06, max.iter = 100, near.zero.var = F, mode = "regression")
  
  # Find significant ASVs
  cim_res <- cim(PLS, comp = 1:3, cutoff = 0.4, dist.method = c("correlation", "correlation"), clust.method = c("complete", "complete"), mapping = "XY")
  
  # Subset significant ASVs
  mat <- cim_res$mat
  subset <- row.names(mat)
  asv <- asv[names(asv) %in% subset]
  
  # Final PLS with subset
  PLS <- spls(asv, meta, ncomp = 3, tol = 1e-06, max.iter = 100, near.zero.var = F, mode = "regression")
  cim_res <- cim(PLS, comp = 1:3, dist.method = c("correlation", "correlation"), clust.method = c("complete", "complete"), mapping = "XY")
  
  # Extract correlations & clusters
  corr <- data.frame(cim_res$mat) %>% mutate(across(everything(), round, 2))
  hc <- hclust(dist(cim_res$mat, method = "correlation"), "complete")
  cutree <- cutree(hc, k=8)
  length(unique(cutree))
  # Define colors
  col.cluster <- c("C1" = "plum4", "C2" = "aquamarine2",
                   "C3" = "palevioletred1", "C4" = "black",
                   "C5" = "gray55", "C6" = "yellow3",
                   "C7" = "dodgerblue4", "C8" = "greenyellow",
                   "C9" = "orangered4", "C10" = "orange1",
                   "C11" = "mediumorchid3", "C12" = "darkcyan",
                   "C13" = "red", "C14" = "green")
  
  # Create a color vector for ASVs based on cluster assignments
  asv_colors <- col.cluster[cutree]
  
  # Plot clusters
  tree <- as.phylo(hc)
  pdf(file = file.path(plot_folder, paste0("Tree_", obj_name, ".pdf")))
  plot(tree, tip.color = asv_colors, cex = 0.3, main = "Cluster Tree")
  dev.off()
  
  # Prepare data for DIABLO analysis
  data <- list(ASV = asv, ENV = meta)
  Y <- as.factor(ENV$depth_category)  # Assuming 'depth_category' is the class label
  
  # Set parameters for DIABLO analysis
  ncomp <- 3
  list.keepX <- list(ASV = c(25, 25, 25), ENV = c(9, 9, 9))  # Adjusted to match the number of columns in ENV
  design <- matrix(0.1,
                   ncol = length(data),
                   nrow = length(data),
                   dimnames = list(names(data),
                                   names(data)))
  diag(design) <- 0
  

  
  # Plot the sPLS scatter plot
  pdf(file = file.path(plot_folder, paste0("sPLS_scatter_", obj_name, ".pdf")))
  plotIndiv(PLS, comp = c(1, 2), ind.names = FALSE, ellipse = TRUE, title = "sPLS Scatter Plot")
  plotVar(PLS, comp = c(1, 2), var.names = c(FALSE, TRUE), title = obj_name)
  dev.off()
  
  
  
  pdf(file = file.path(plot_folder, paste0("sPLS_scatter_arrow_", obj_name, ".pdf")))
  plotIndiv(PLS, comp = c(1, 2), ind.names = FALSE, ellipse = TRUE, title = "sPLS Scatter Plot")
  plotVar(PLS, comp = c(1, 2),
          var.names = c(FALSE, FALSE),
          abline = FALSE,
          title = obj_name)
  dev.off()
  # Reformat clusters
  clusters <- cbind(corr, cutree)
  clusters <- data.frame(cbind(rownames(clusters),
                               clusters[, "cutree"])) %>% dplyr::rename(asv = X1, cluster = X2) %>% mutate(cluster = paste0("C", cluster))
  
  # Save cluster information as a text file
  write.table(clusters, file = file.path(plot_folder,
                                         paste0("clusters_", obj_name, ".txt")),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Plot and save the cim plot
  cim(PLS, dist.method = c("correlation", "correlation"),
      clust.method = c("complete", "complete"),
      margins = c(5, 10),
      row.names = F,
      title = "Correlation Integral Mapping",
      mapping = "XY",
      color = colorRampPalette(c("pink4", "pink2", "floralwhite", "lightskyblue2", "skyblue4"))(20),
      save = "pdf", name.save = file.path(plot_folder, paste0("sPLS_heat_", obj_name, ".pdf")))
  
  # Combine cluster 
  tax <- TAX[match(clusters$asv, rownames(TAX)), ] %>%
    add_column(cluster = clusters$cluster) %>% 
    mutate(cluster = case_when(cluster %in% c("C6") ~ "C6", TRUE ~ cluster)) %>% 
    rownames_to_column("asv")
  

  
  # Reassign to clustCol
  clustCol <- tibble(cluster = names(col.cluster)) %>%
    left_join(tax) %>% dplyr::slice(mixedorder(asv)) %>% drop_na()
  
  # Add new colors for combined clusters
  clustCol$col <- col.cluster[match(clustCol$cluster, names(col.cluster))]
  
  # ASVs per cluster
  clustNum <- clustCol %>% group_by(cluster) %>% tally()
  
  # Export cluster info for SI table
  clustCol %>% dplyr::slice(mixedorder(cluster)) %>%
    dplyr::select(-c("Kingdom", "Species", "col")) %>%
    write.table(file = file.path(plot_folder,
                                 paste0("ASV_PLS-clust_", obj_name, ".txt")),
                sep = "\t", row.names = F, col.names = T, quote = F)
  
  # Export correlations for SI table
  corr %>% rownames_to_column("asv") %>%
    dplyr::slice(mixedorder(asv)) %>%
    write.table(file = file.path(plot_folder,
                                 paste0("ASV_PLS-corr_", obj_name, ".txt")),
                sep = "\t", row.names = F, col.names = T, quote = F)
  
  # More unclassifieds in EGC/deep?
  clustComp <- clustCol %>% group_by(cluster) %>%
    summarize(count = sum(grepl("uc", Genus))) %>% 
    left_join(clustNum) %>% mutate(frac = count / n * 100)
  
  # Export for SI table
  clustComp %>% dplyr::slice(mixedorder(cluster)) %>%
    relocate(cluster, n, count, frac) %>%
    dplyr::rename(`#ASV` = n, `#unclassified` = count, `%unclassified` = frac) %>% 
    write.table(file = file.path(plot_folder, paste0("ASV_PLS-comp_", obj_name, ".txt")),
                sep = "\t", row.names = F, col.names = T, quote = F)
  
  # SIGNATURE GENERA
  ASV.hel <- as.data.frame(apply(ASV, 2, function(x) sqrt(x / sum(x))))
  clustAsv <- ASV.hel %>% rownames_to_column("asv")
  
  # Filter low-abundant ASVs
  clustTax <- ASV.hel %>% rownames_to_column("asv") %>%
    reshape2::melt(id.vars = c("asv")) %>% group_by(asv) %>%
    summarise(value = max(value)) %>% left_join(clustCol) %>%
    drop_na() %>% filter(value > 0.05) %>% left_join(clustAsv) %>% 
    dplyr::select(-c(Kingdom, Phylum, Order, Family, Species)) %>% 
    reshape2::melt() %>% 
    aggregate(value ~ cluster + Class + Genus + variable, data = ., FUN = sum) %>%
    group_by(Genus) %>% top_n(1, value) %>% distinct(value, .keep_all = T) %>%
    mutate(cluster = factor(cluster, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11")))
  
  # Sort genera by cluster
  clustTax$Genus <- factor(clustTax$Genus, levels = rev(unique(clustTax$Genus[order(clustTax$cluster)])))
  clustTax <- ASV.hel %>% 
    rownames_to_column("asv") %>%
    reshape2::melt(id.vars = c("asv")) %>% 
    group_by(asv) %>%
    summarise(value = max(value)) %>% 
    left_join(clustCol) %>% 
    drop_na() %>% 
    filter(value > 0.05) %>% 
    left_join(clustAsv) %>% 
    dplyr::select(-c(Kingdom, Phylum, Order, Family, Species)) %>% 
    reshape2::melt() %>% 
    group_by(cluster, Class, Genus) %>%
    summarise(value = sum(value), .groups = 'drop') %>%
    mutate(cluster = factor(cluster, levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11")))
  # Filter most uc taxa / low-abundant
  # Filter ambiguities (e.g. Magnetospiraceae-uc/Magnetospira)
  # Export size 
  ggplot(data = clustTax) +
    geom_point(aes(x = cluster, y = Genus, size = value, colour = cluster), stat = "identity", position = "identity") +
    scale_colour_manual(values = col.cluster) +
    scale_size_continuous(range = c(1, 5), breaks = c(0.3, 0.6, 0.9), name = "Maximum abundance") +
    guides(fill = T, color = "none") +
    scale_y_discrete(limits = rev) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.ticks = element_blank(), axis.text.y = element_text(size = 10, colour = "black"),
          legend.position = "top")
  ggsave(file.path(plot_folder, paste0("signature_genera_", obj_name, ".pdf")), width = 5, height = 10)
  
  # Percent contribution of signature-pops -- layer/region
  ASV.rel <- as.data.frame(apply(ASV, 2, function(x) x / sum(x) * 100))
  clustFrac <- ASV.rel %>% rownames_to_column("asv") %>% 
    reshape2::melt(id.vars = c("asv")) %>% left_join(clustCol) %>%
    mutate(type = case_when(is.na(cluster) ~ "unassigned", TRUE ~ "signature")) %>%
    group_by(variable, type) %>% summarize(sum = sum(value)) %>%
    left_join(ENV, by = c("variable" = "CTD #")) %>%
    group_by("lon", "lat", "NO3+NO2_c_mean",
             "PO4_c_mean", "Si(OH)4_c_mean",
             "NO2_c_mean", "NH4_c_mean",
             "TDN_c_mean", "TDP_c_mean",
             "Sal00" ,"Fix_Temp_C") %>% 
    summarize(mean = mean(sum)) %>% ungroup()
  
  # Percent contribution of signature-pops -- year/region
  clustYear <- ASV.rel %>% rownames_to_column("asv") %>%
    reshape2::melt(id.vars = c("asv")) %>% left_join(clustCol) %>%
    group_by(variable, cluster) %>% summarize(sum = sum(value)) %>%
    left_join(ENV, by = c("variable" = "CTD #")) %>% 
    group_by("lon", "lat", "NO3+NO2_c_mean",
             "PO4_c_mean", "Si(OH)4_c_mean",
             "NO2_c_mean", "NH4_c_mean",
             "TDN_c_mean", "TDP_c_mean",
          "Sal00" ,"Fix_Temp_C")%>% 
    summarize(mean = mean(sum)) %>% ungroup()
  
  # Export for SI table
  clustFrac %>% arrange('depth [m]') %>% 
    dplyr::rename(`Relative abundance of signature clusters (%)` = mean) %>% 
    write.table(file = file.path(plot_folder, paste0("signature_clusters_abundance.txt")),
                sep = "\t", row.names = F, col.names = T, quote = F)
  
  
  library(gplots)
  library(dplyr)
  n_clusters <- as.numeric(length(unique(cutree)))
  # Extract the matrix from CIM results
  mat <- cim_res$mat
  
  # Use the original clustering
  clusters <- cutree(hc, k = n_clusters)
  
  # Ensure col.cluster has enough colors
  if(length(col.cluster) < n_clusters) {
    col.cluster <- colorRampPalette(col.cluster)(n_clusters)
  }
  
  # Assign colors to clusters using the original col.cluster
  cluster_colors <- col.cluster[1:n_clusters]
  names(cluster_colors) <- sort(unique(clusters))
  
  # Create side colors based on the original clusters
  side_colors_ordered <- cluster_colors[as.character(clusters)]
  
  # Create the heatmap
  png(file = file.path(plot_folder, paste0("cim_plot_cluster_col_", obj_name, ".png")), width = 1500, height = 850, res = 100)
  
  heatmap.2(mat,
            Rowv = as.dendrogram(hc),  # Use the original dendrogram
            Colv = FALSE,
            dendrogram = "row",
            col = colorRampPalette(c("darkred", "pink2", "floralwhite", "lightskyblue2", "darkblue"))(20),
            trace = "none",
            margins = c(10, 10),
            main = "Correlation Integral Mapping",
            cexRow = 0.4,
            cexCol = 1.1,
            keysize = 1,
            key.title = "Correlation",
            key.xlab = "Correlation Value",
            density.info = "none",
            lhei = c(1.5, 4),  
            lwid = c(1.5, 4, 1),
            RowSideColors = side_colors_ordered,
            labRow = FALSE)
  
  # Add legend for clusters
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("right", 
         legend = names(cluster_colors), 
         fill = cluster_colors, 
         title = "Clusters", 
         cex = 1.2,  # Increased text size for legend
         inset = c(0.015, 0),  # Adjusted position
         xpd = TRUE,  # Allow plotting outside the plot region
         title.cex = 1.4)  # Increased text size for legend title
  
  dev.off()
  # Step 1: Create a data frame with correlations and cluster assignments
  cor_df <- as.data.frame(mat) %>%
    mutate(Cluster = clusters) %>%
    pivot_longer(cols = -Cluster, names_to = "EnvFactor", values_to = "Correlation")
  
  # Step 2: Calculate mean correlation for each cluster and environmental factor
  mean_cors <- cor_df %>%
    group_by(Cluster, EnvFactor) %>%
    summarize(MeanCorrelation = mean(Correlation, na.rm = TRUE), .groups = 'drop') %>%
    arrange(Cluster, desc(abs(MeanCorrelation)))
  
  # Step 3: Rank correlations within each cluster
  ranked_cors <- mean_cors %>%
    group_by(Cluster) %>%
    mutate(Rank = row_number()) %>%
    arrange(Cluster, Rank)
  
  # Step 4: Print the results
  for (i in sort(unique(clusters))) {
    cat("\nCluster", i, ":\n")
    cluster_data <- ranked_cors[ranked_cors$Cluster == i, c("EnvFactor", "MeanCorrelation", "Rank")]
    print(cluster_data, row.names = FALSE)
    cat("\n")
  }
  
  # Step 5: Optionally, you can save these results to a file
  write.csv(ranked_cors, file = file.path(plot_folder, paste0("cluster_env_correlations_", obj_name, ".csv")), row.names = FALSE)  ############################
  ####################################
  library(dplyr)
  library(tidyr)
  
  # Step 1: Prepare the data
  asv_data <- as.data.frame(t(ASV))  # Transpose ASV so that samples are rows and ASVs are columns
  cluster_data <- data.frame(ASV = names(clusters), Cluster = clusters)
  
  # Step 2: Calculate mean abundance for each ASV in each depth category
  asv_depth_abundance <- asv_data %>%
    rownames_to_column("CTD") %>%
    pivot_longer(cols = -CTD, names_to = "ASV", values_to = "Abundance") %>%
    left_join(ENV %>% dplyr::select(!!sym("CTD #"), depth_category), by = c("CTD" = "CTD #")) %>%
    inner_join(cluster_data, by = "ASV") %>%
    group_by(Cluster, depth_category) %>%
    summarize(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = 'drop')
  
  # Step 3: Find the depth category with the highest mean abundance for each cluster
  cluster_dominant_depth <- asv_depth_abundance %>%
    group_by(Cluster) %>%
    slice_max(order_by = MeanAbundance, n = 1) %>%
    ungroup()
  
  # Step 4: Calculate the total abundance for each cluster in each depth category
  cluster_depth_distribution <- asv_depth_abundance %>%
    group_by(Cluster) %>%
    mutate(TotalAbundance = sum(MeanAbundance),
           RelativeAbundance = MeanAbundance / TotalAbundance) %>%
    ungroup()
  
  # Step 5: Print the results
  for (i in sort(unique(clusters))) {
    cat("\nCluster", i, ":\n")
    cat("Dominant depth category:", cluster_dominant_depth$depth_category[cluster_dominant_depth$Cluster == i], "\n")
    cat("Depth distribution:\n")
    print(cluster_depth_distribution %>% 
            filter(Cluster == i) %>% 
            dplyr::select(depth_category, MeanAbundance, RelativeAbundance) %>% 
            arrange(desc(RelativeAbundance)),
          row.names = FALSE)
    cat("\n")
  }
  
  library(dplyr)
  library(tidyr)
  
  cleaned_data <- cluster_depth_distribution %>%
    mutate(
      # Round MeanAbundance to 2 decimal places
      MeanAbundance = round(as.numeric(MeanAbundance), 2),
      # Recalculate TotalAbundance
      TotalAbundance = sum(MeanAbundance),
      # Recalculate and round RelativeAbundance to 4 decimal places
      RelativeAbundance = round(MeanAbundance / TotalAbundance, 4)
    ) %>%
    arrange(desc(RelativeAbundance))
  
  # Print the cleaned and recalculated data
  print(cleaned_data, n = Inf)
  
  # Verify that RelativeAbundance sums to 1
  cat("\nSum of RelativeAbundance:", sum(cleaned_data$RelativeAbundance))
  
  # Save the results to a CSV file with proper formatting
  write.csv(cleaned_data, 
            file = file.path(plot_folder, paste0("cluster_depth_distribution_", obj_name, "_rounded.csv")), 
            row.names = FALSE)

  # Group by cluster, sort by value, and select top 20 for each cluster
  top_20_per_cluster <- clustTax %>%
    group_by(cluster) %>%
    arrange(desc(value)) %>%
    slice_head(n = 25) %>%
    ungroup()
  
  # Sort the entire dataframe by cluster and then by value
  top_20_per_cluster <- top_20_per_cluster %>%
    arrange(cluster, desc(value))
  # Round the 'value' column to 2 decimal places
  top_20_per_cluster$value <- round(top_20_per_cluster$value, 2)
  
  # Write the result to a CSV file
  write.table(top_20_per_cluster,
              file = file.path(plot_folder, paste0("top_20_taxa_per_cluster_", obj_name, "_rounded.txt")),
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)
  # If you want to see the result in the console
  print(top_20_per_cluster, n = Inf)
  ###################
 
  ############
  
  # First, ensure depth_category is a factor with levels in the correct order
  depth_order <- c("MP", "Ice-TS", "Ice-BS", "UIW", "2", "10", "chl-max", "50", "100", "200",
                   "500", "1000", "1500", "2000", "3000", "20m from bottom", "5m from bottom", "bottom")
  cluster_depth_distribution$depth_category <- factor(cluster_depth_distribution$depth_category, levels = rev(depth_order))
  # Convert cluster to factor
  cluster_depth_distribution$Cluster <- as.factor(cluster_depth_distribution$Cluster)
  
  # Create the plot
  p <- ggplot(cluster_depth_distribution, aes(x = depth_category, y = RelativeAbundance, group = Cluster)) +
    geom_line(aes(color = Cluster), linewidth = 1) +
    geom_point(aes(color = Cluster), size = 3) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    facet_wrap(~ Cluster, scales = "fixed", ncol = 4) +
    scale_color_manual(values = cluster_colors) +  # Use your predefined color palette
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "lightgrey"),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8),
      plot.margin = unit(c(1, 1, 1, 4), "lines")
    ) +
    labs(x = "Depth Category", 
         y = "Relative Abundance", 
         title = "Cluster Relative Abundances Across Depth Categories") +
    coord_flip()
  
  # Print the plot
  print(p)
  
  # Save the plot
  ggsave(file.path(plot_folder, paste0("cluster_relative_abundance_by_depth_", obj_name, ".pdf")), 
         p, width = 12, height = 8, dpi = 300)
  
}

  