#!/usr/bin/env Rscript
# Upload.R
#
# This script performs spatial network analysis on glioma immune cell data.
# It reads in processed data, computes a global layout for immune cell types,
# calculates pairwise distances between cells within each image, aggregates them
# per condition (defined by hypoxia and IDH status), constructs network graphs,
# and generates primary as well as supplementary figures.
#
# Dependencies: tidyverse, ggnetwork, igraph, ggraph, tidygraph, circlize,
#               ggforce, readxl, reshape2, pheatmap
#
# Usage: Set the working directory and adjust file paths as needed. Then run the script.
#


# --- 1. Load Required Libraries ---
# Load a set of R packages for various tasks including data manipulation, plotting,
# network analysis, and more. Each package adds specific capabilities:
library(tidyverse)    # Provides a suite of tools including ggplot2, dplyr, readr, etc.
library(ggnetwork)    # Helps with converting network objects for plotting with ggplot2.
library(igraph)       # Offers comprehensive network analysis and graph manipulation.
library(ggraph)       # Extends ggplot2 for visualizing graphs and networks.
library(tidygraph)    # Optional: offers a tidy approach to handling graphs.
library(circlize)     # Useful for circular layout calculations and chord diagram plotting.
library(ggforce)      # Provides additional visualizations such as convex hulls (for clustering).
library(readxl)       # Enables reading Excel files.
library(reshape2)     # Supplies functions for reshaping data (melting/wcasting matrices).
library(pheatmap)     # Used to create heatmaps for visualizing distance matrices.

# --- 2. Read and Merge Data ---
# Set the appropriate working directory before running this script.
# Read sample information from an Excel file which includes metadata about images.
Sample_Info <- read_excel("Data/SampleInfo.xlsx")

# Read processed cell data stored in an RDS file containing combined measurements.
combined_data <- readRDS("Data/Processed/Combined.RDS")

# Filter the combined data to include only those images representing either 'Astrocytoma_Diffuse'
# or 'GBM' tumors as identified in the sample information. Then, add a new column 'IDH_Status'
# based on the tumor type.
all_gbm_data <- combined_data %>%
  filter(Image %in% subset(Sample_Info, Tumor %in% c("Astrocytoma_Diffuse", "GBM"))$Image) %>%
  mutate(IDH_Status = case_when(
    Image %in% subset(Sample_Info, Tumor == "Astrocytoma_Diffuse")$Image ~ "IDHmut",
    Image %in% subset(Sample_Info, Tumor == "GBM")$Image ~ "IDHwt"))

# --- 3. Compute a Global Layout for Cell Types ---
# Determine the unique cell types from the data and create a corresponding nodes data frame.
global_cell_types <- unique(all_gbm_data$Classification)
global_nodes <- data.frame(name = global_cell_types,
                           label = global_cell_types,
                           stringsAsFactors = FALSE)

# For the global network layout, we start with an empty data frame for edges because here
# we only want to establish node positions based on the cell types.
empty_edges <- data.frame(from = character(0),
                          to = character(0),
                          stringsAsFactors = FALSE)

# Create an igraph object using the nodes and no edges. This allows us to compute a fixed layout.
global_graph <- graph_from_data_frame(d = empty_edges, vertices = global_nodes)

# Compute a layout for the network using the Fruchterman-Reingold algorithm.
global_layout_df <- create_layout(global_graph, layout = "fr")

# --- 4. Define Functions ---

# Function: compute_avg_distances()
# This function calculates the pairwise median and standard deviation of distances between cells 
# of different cell types. It accepts the data frame and several column names, and it limits 
# the number of computed pairs to avoid excessive calculations.
compute_avg_distances <- function(df,
                                  x_col = "Centroid X µm",
                                  y_col = "Centroid Y µm",
                                  class_col = "Classification",
                                  max_pairs = 1e6) {
  # Get unique cell types from the provided classification column.
  cell_types <- unique(df[[class_col]])
  # Initialize an empty list to store the results.
  results_list <- list()
  
  # Loop over each combination of cell type pairs.
  for (ct1 in cell_types) {
    for (ct2 in cell_types) {
      # Extract coordinate matrices for cell type ct1 and ct2.
      A <- df %>% filter(.data[[class_col]] == ct1) %>%
        dplyr::select(all_of(x_col), all_of(y_col)) %>% as.matrix()
      B <- df %>% filter(.data[[class_col]] == ct2) %>%
        dplyr::select(all_of(x_col), all_of(y_col)) %>% as.matrix()
      
      # If either cell type has no cells, record NA values.
      if (nrow(A) == 0 || nrow(B) == 0) {
        results_list[[length(results_list) + 1]] <- data.frame(
          cell_type_1 = ct1, cell_type_2 = ct2,
          mean_dist = NA, sd_dist = NA,
          stringsAsFactors = FALSE)
        next
      }
      
      # When comparing cells of the same type, consider only one half of the symmetric distance matrix.
      if (ct1 == ct2) {
        if(nrow(A) < 2){
          results_list[[length(results_list) + 1]] <- data.frame(
            cell_type_1 = ct1, cell_type_2 = ct2,
            mean_dist = NA, sd_dist = NA,
            stringsAsFactors = FALSE)
          next
        }
        total_pairs <- choose(nrow(A), 2)
        if (is.na(total_pairs)) total_pairs <- Inf
        # If number of cell pairs exceeds the maximum allowed, sample a fixed number of pairs.
        if (total_pairs > max_pairs) {
          sample_n <- max_pairs
          idx1 <- sample(1:nrow(A), sample_n, replace = TRUE)
          idx2 <- sample(1:nrow(A), sample_n, replace = TRUE)
          dists <- sqrt((A[idx1, 1] - A[idx2, 1])^2 + (A[idx1, 2] - A[idx2, 2])^2)
        } else {
          # Compute the full distance matrix and then extract the upper triangular values.
          dist_mat <- as.matrix(dist(A))
          dists <- dist_mat[upper.tri(dist_mat)]
        }
      } else {
        # For two different cell types, compute the pairwise distances between each cell in A and B.
        total_pairs <- nrow(A) * nrow(B)
        if (is.na(total_pairs)) total_pairs <- Inf
        if (total_pairs > max_pairs) {
          sample_n <- max_pairs
          idxA <- sample(1:nrow(A), sample_n, replace = TRUE)
          idxB <- sample(1:nrow(B), sample_n, replace = TRUE)
          dists <- sqrt((A[idxA, 1] - B[idxB, 1])^2 + (A[idxA, 2] - B[idxB, 2])^2)
        } else {
          # Use vectorized computation to calculate the squared Euclidean distances.
          A_sq <- rowSums(A^2)
          B_sq <- rowSums(B^2)
          dist_sq <- outer(A_sq, B_sq, "+") - 2 * (A %*% t(B))
          # Guard against negative values due to numerical errors.
          dist_sq[dist_sq < 0] <- 0
          dists <- sqrt(as.vector(dist_sq))
        }
      }
      
      # Append the computed median and standard deviation of distances to the results list.
      results_list[[length(results_list) + 1]] <-
        data.frame(cell_type_1 = ct1,
                   cell_type_2 = ct2,
                   mean_dist = median(dists, na.rm = TRUE),
                   sd_dist = sd(dists, na.rm = TRUE),
                   stringsAsFactors = FALSE)
    }
  }
  # Combine all the small data frames into one final results data frame.
  results <- bind_rows(results_list)
  return(results)
}


# Function: create_network_plot()
# This function creates a network plot using pairwise distance data. It:
# - Removes self-loop entries and missing values.
# - Optionally inverts distances so that shorter distances have higher weights.
# - Constructs an igraph object, assigns node positions (using a global layout if available),
#   and performs community detection.
# - Highlights top edges and visualizes clusters with convex hulls.
create_network_plot <- function(dist_results,
                                region_label,
                                idh_label,
                                invert_dist = FALSE,
                                global_layout = NULL,
                                annotation_text = NULL) {
  # Filter edges: remove any with NA values and self-joins (i.e., same cell type).
  edges <- dist_results %>%
    filter(!is.na(mean_dist)) %>%
    filter(cell_type_1 != cell_type_2)
  
  tiny_epsilon <- 1e-5  # Prevent division by zero when inverting.
  if (invert_dist) {
    # Invert distances: shorter physical distances result in higher connection weights.
    edges <- edges %>% mutate(weight = 1 / (mean_dist + tiny_epsilon))
    edge_label_title <- "1 / Mean Dist"
  } else {
    # Use the mean distance directly as the weight.
    edges <- edges %>% mutate(weight = mean_dist)
    edge_label_title <- "Mean Dist"
  }
  
  # Build an undirected graph from the edges data frame.
  g <- graph_from_data_frame(d = edges, directed = FALSE)
  
  # If a fixed global layout is provided, assign x and y positions to each node.
  if (!is.null(global_layout)) {
    node_positions <- global_layout %>% filter(name %in% V(g)$name)
    node_positions <- node_positions[match(V(g)$name, node_positions$name), ]
    V(g)$x <- node_positions$x
    V(g)$y <- node_positions$y
  }
  
  # Set up an abbreviation lookup for certain cell types; adjust this mapping as needed.
  abbr_map <- c("Macrophage" = "M1_mac",
                "Microglia"   = "M2_mic",
                "DendriticCell" = "DC")
  V(g)$abbrev <- ifelse(V(g)$name %in% names(abbr_map),
                        abbr_map[V(g)$name],
                        V(g)$name)
  
  # Perform community detection using the Louvain method to find clusters of highly connected nodes.
  clust <- cluster_louvain(g)
  V(g)$cluster <- as.factor(membership(clust))
  
  # Prepare the base plot. Use a fixed layout if provided; otherwise, use the Fruchterman-Reingold layout.
  if (!is.null(global_layout)) {
    p <- ggraph(g, layout = "manual", x = V(g)$x, y = V(g)$y)
  } else {
    p <- ggraph(g, layout = "fr")
  }
  
  # Add edges to the plot with properties reflecting connection weights.
  p <- p + geom_edge_link(aes(width = weight, color = weight),
                          show.legend = TRUE, alpha = 0.8)
  
  # Highlight the top three edges (those with the smallest mean distances) by labeling them.
  top_edges <- edges %>% arrange(mean_dist) %>% head(3)
  if(nrow(top_edges) > 0) {
    # Get current node positions for later annotation.
    vertex_pos <- data.frame(name = V(g)$name,
                             x = V(g)$x,
                             y = V(g)$y,
                             stringsAsFactors = FALSE)
    # Merge positions for both source (from) and target (to) nodes in the top edges.
    top_edges <- top_edges %>%
      rename(from = cell_type_1, to = cell_type_2) %>%
      left_join(vertex_pos, by = c("from" = "name")) %>%
      rename(x_from = x, y_from = y) %>%
      left_join(vertex_pos, by = c("to" = "name")) %>%
      rename(x_to = x, y_to = y) %>%
      mutate(mid_x = (x_from + x_to) / 2,
             mid_y = (y_from + y_to) / 2)
    
    # Add text labels at the midpoints of these top edges.
    p <- p + geom_text(data = top_edges,
                       aes(x = mid_x, y = mid_y, label = round(mean_dist, 1)),
                       size = 4, color = "black")
  }
  
  # Draw convex hulls around clusters that have at least three nodes to visually highlight them.
  cluster_df <- data.frame(name = V(g)$name,
                           x = V(g)$x,
                           y = V(g)$y,
                           cluster = V(g)$cluster)
  hull_cluster_df <- cluster_df %>%
    group_by(cluster) %>%
    filter(n() >= 3) %>%
    slice(chull(x, y))  # chull() computes the convex hull for points in each cluster.
  
  if(nrow(hull_cluster_df) > 0){
    p <- p + geom_polygon(data = hull_cluster_df,
                          aes(x = x, y = y, group = cluster),
                          fill = NA, color = "black", linetype = "dashed",
                          inherit.aes = FALSE)
  }
  
  # Add nodes (cells) to the plot as points with labels using the previously defined abbreviations.
  p <- p + geom_node_point(fill = "lightgray", size = 5, shape = 21, color = "black") +
    geom_node_text(aes(label = abbrev), repel = TRUE, size = 6, family = "sans") +
    scale_edge_width(trans = "reverse", range = c(1, 6), name = edge_label_title) +
    scale_edge_color_gradient(trans = "reverse", low = "lightblue", high = "red", name = edge_label_title) +
    labs(title = paste("Network of Cell Types by Average Distance\n",
                       "Region:", region_label, "| IDH:", idh_label))
  
  # Optionally add an annotation label in the top right corner of the plot.
  if (!is.null(annotation_text)) {
    x_annot <- max(V(g)$x, na.rm = TRUE)
    y_annot <- max(V(g)$y, na.rm = TRUE)
    p <- p + annotate("label", x = x_annot, y = y_annot,
                      label = annotation_text, size = 3,
                      fill = "white", color = "black")
  }
  
  # Finalize the plot with a graph theme and adjust the legend position.
  p <- p + theme_graph(base_size = 14, base_family = "sans") +
    theme(legend.position = "right")
  
  return(p)
}


# Function: compute_raw_distances()
# This function computes the detailed distribution of pairwise distances for a given pair of cell types.
# It optionally limits the number of cells sampled to make the computation more efficient.
compute_raw_distances <- function(df, celltype1, celltype2,
                                  x_col = "Centroid X µm",
                                  y_col = "Centroid Y µm",
                                  class_col = "Classification",
                                  sample_size = 100) {
  # Retrieve coordinates for the two specified cell types.
  A <- df %>% filter(.data[[class_col]] == celltype1) %>%
    dplyr::select(all_of(x_col), all_of(y_col)) %>% as.matrix()
  B <- df %>% filter(.data[[class_col]] == celltype2) %>%
    dplyr::select(all_of(x_col), all_of(y_col)) %>% as.matrix()
  
  # If there are more cells than the sample size, randomly sample a subset.
  if(nrow(A) > sample_size) {
    A <- A[sample(seq_len(nrow(A)), sample_size), , drop = FALSE]
  }
  if(nrow(B) > sample_size) {
    B <- B[sample(seq_len(nrow(B)), sample_size), , drop = FALSE]
  }
  
  if(celltype1 == celltype2){
    # For the same cell type, calculate distances using the upper triangle of the distance matrix.
    if(nrow(A) < 2) return(numeric(0))
    dist_mat <- as.matrix(dist(A))
    dist_vals <- dist_mat[upper.tri(dist_mat)]
  } else {
    # For different cell types, compute every pairwise distance between A and B.
    dist_vals <- sqrt((outer(A[,1], B[,1], "-")^2) +
                        (outer(A[,2], B[,2], "-")^2))
    dist_vals <- as.vector(dist_vals)
  }
  return(dist_vals)
}

# --- 5. Main Execution Block ---
# Prepare an empty list to store computed distance metrics for different conditions.
all_chord_results <- list()

# Identify unique values for regions and IDH statuses from the dataset.
unique_regions <- unique(all_gbm_data$Parent)
unique_idh     <- unique(all_gbm_data$IDH_Status)

# Loop through each combination of region and IDH status.
for (this_region in unique_regions) {
  for (this_idh in unique_idh) {
    # For each condition, get the subset of cells and then compute distances within each image
    df_sub <- all_gbm_data %>%
      filter(Parent == this_region, IDH_Status == this_idh)
    
    if(nrow(df_sub) == 0) next
    
    # Get the unique images within this condition:
    image_list <- unique(df_sub$Identifier)
    
    # Initialize a list to store the distance tables for each image:
    dist_results_list <- list()
    
    # Loop over images
    for (img in image_list) {
      df_img <- df_sub %>% filter(Identifier == img)
      # Only compute if there is enough data in this image.
      if(nrow(df_img) > 0) {
        # Compute pairwise distances within this image
        dist_img <- compute_avg_distances(df_img)
        # Record the Image identifier (for future reference)
        dist_img$Identifier <- img
        dist_results_list[[img]] <- dist_img
      }
    }
    
    # Combine the per-image results into a single data frame:
    dist_results_combined <- bind_rows(dist_results_list)
    
    # Now aggregate across images for each cell type pair.
    # Here we take the median of the per-image median distances as a robust central estimate.
    dist_results_aggregated <- dist_results_combined %>%
      group_by(cell_type_1, cell_type_2) %>%
      summarize(mean_dist = mean(mean_dist, na.rm = TRUE),
                sd_dist = mean(sd_dist, na.rm = TRUE)) %>%  # Alternatively, use median for sd if desired.
      ungroup()
    
    # Save the aggregated results under a unique key:
    key_name <- paste(this_region, this_idh, sep = "_")
    all_chord_results[[key_name]] <- dist_results_aggregated
    
    # Optionally, print a message:
    message("\n--- Aggregated Distance Results for ", key_name, " ---")
    print(dist_results_aggregated)
  }
}

# Save incase computationally heavy #
#saveRDS(all_chord_results,"/Data/Processed/all_chord_results.RDS")

# --- 6. Generate Primary Figure: Network Plots ---
## Incase loading from here is preferred
#all_chord_results <- readRDS("/Data/Processed/all_chord_results.RDS")

# Open a PDF device to save network plot figures.
pdf("plots/CellType_Distance_Networks.pdf", width = 10, height = 10, family = "sans")

# Loop over each condition to generate a network plot.
for (key_name in names(all_chord_results)) {
  # Extract region and IDH information from the key name.
  parts <- unlist(strsplit(key_name, "_"))
  if (length(parts) == 2) {
    region_label <- parts[1]
    idh_label <- parts[2]
  } else if (length(parts) >= 3) {
    region_label <- paste(parts[1:2], collapse = "_")
    idh_label <- parts[3]
  } else {
    region_label <- key_name
    idh_label <- "Unknown"
  }
  
  dist_df <- all_chord_results[[key_name]]
  
  # Create the network plot for the current condition.
  plot_network <- create_network_plot(dist_df,
                                      region_label = region_label,
                                      idh_label = idh_label,
                                      invert_dist = FALSE,
                                      global_layout = global_layout_df)
  print(plot_network)
}
dev.off()  # Close the PDF device.

# --- 7. Generate Supplementary Figures ---
# Open a second PDF device to save additional supplementary visualizations (heatmaps, violin plots etc.).
pdf("plots/Supplementary_Visualizations_AllConditions.pdf", width = 10, height = 8)

for (key_name in names(all_chord_results)) {
  # Parse the key name to extract region and IDH labels.
  parts <- unlist(strsplit(key_name, "_"))
  if (length(parts) == 2) {
    region_label <- parts[1]
    idh_label <- parts[2]
  } else if (length(parts) >= 3) {
    region_label <- paste(parts[1:2], collapse = "_")
    idh_label <- parts[3]
  } else {
    region_label <- key_name
    idh_label <- "Unknown"
  }
  
  # Subset the dataset for the current condition.
  df_condition <- all_gbm_data %>%
    filter(Parent == region_label, IDH_Status == idh_label)
  
  dist_results <- all_chord_results[[key_name]]
  
  # Check if the number of unique cell types is manageable for plotting.
  unique_cells <- unique(c(dist_results$cell_type_1, dist_results$cell_type_2))
  if(length(unique_cells) > 50) {
    warning("Too many unique cell types (", length(unique_cells), "). Consider subsetting the data.")
    next
  }
  
  # 7.1 Pairwise Distance Heatmap
  # Reshape the distance data into a wide format suitable for a heatmap.
  heatmap_data <- dist_results %>%
    dplyr::select(cell_type_1, cell_type_2, mean_dist) %>%
    pivot_wider(names_from = cell_type_2,
                values_from = mean_dist,
                values_fill = list(mean_dist = NA)) %>%
    as.data.frame()
  
  # Use the cell types as row names.
  heatmap_data <- column_to_rownames(heatmap_data, var = "cell_type_1")
  
  # Generate the heatmap with specified cell sizes and styling.
  pheatmap(heatmap_data, 
           border_color = "white",
           main = paste("Heatmap ", key_name),
           cellwidth = 15,
           cellheight = 15,
           scale = "none",
           display_numbers = FALSE,
           color = colorRampPalette(c("pink", "white", "skyblue"))(100))
  
  # 7.2 Distance Distribution Violin/Box Plot
  # Create a data frame to store distribution information for each cell pair.
  pairs_df <- unique(dist_results[, c("cell_type_1", "cell_type_2")])
  df_violin <- data.frame(distance = numeric(0), pair = character(0))
  
  # Compute the raw distance distribution for each cell-type pair.
  for(i in 1:nrow(pairs_df)) {
    ct1 <- as.character(pairs_df[i, "cell_type_1"])
    ct2 <- as.character(pairs_df[i, "cell_type_2"])
    raw_dist <- compute_raw_distances(df_condition, ct1, ct2)
    if(length(raw_dist) > 0) {
      temp_df <- data.frame(distance = raw_dist,
                            pair = rep(paste(ct1, "-", ct2), length(raw_dist)),
                            stringsAsFactors = FALSE)
      df_violin <- rbind(df_violin, temp_df)
    }
  }
  
  # Only plot if data exists.
  if(nrow(df_violin) > 0){
    p_violin <- ggplot(df_violin, aes(x = pair, y = distance, fill = pair)) +
      geom_violin(trim = FALSE, alpha = 0.7) +  # Plot the full distribution as a violin plot.
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Overlay a boxplot for summary stats.
      labs(title = paste("Distance Distribution for", key_name),
           x = "Cell Pair", y = "Distance (µm)") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, hjust = 1))
    print(p_violin)
  }
  
  gc()  # Trigger garbage collection to free up memory.
}
dev.off()  # Close the supplementary figures PDF device.
