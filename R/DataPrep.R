#!/usr/bin/env Rscript
# DataPrep.R
#
# This script reads and concatenates multiple CSV files from a specified folder,
# adds TMA labels to each file, merges the data into a single data frame, updates
# the hypoxia annotation (Parent column) based on specific conditions, and saves
# the final processed data as an RDS file.
#
# Dependencies: readr, dplyr
# Usage: Modify the folder_path variable below as needed and run the script.

# Load necessary libraries
library(readr)  # For reading CSV files
library(dplyr)  # For data manipulation

# Function to read and concatenate CSV files from a given folder path.
concatenate_csv_files <- function(folder_path) {
  # List all CSV files in the folder
  csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Initialize an empty list to store data frames
  df_list <- list()
  
  # Read each CSV file and store in the list
  for (file in csv_files) {
    df_list[[file]] <- read_csv(file)
  }
  
  # Assign TMA labels to each file (adjust these if needed)
  df_list[[1]]$TMA <- "Allglioma"
  df_list[[2]]$TMA <- "Astro 15"
  df_list[[3]]$TMA <- "Astro 16"
  df_list[[4]]$TMA <- "Astro 21"
  df_list[[5]]$TMA <- "Astro 22"
  df_list[[6]]$TMA <- "Astro 22"
  df_list[[7]]$TMA <- "Astro 23"
  df_list[[8]]$TMA <- "Astro 24"
  
  # Concatenate all data frames into a single data frame
  combined_df <- bind_rows(df_list)
  
  return(combined_df)
}

# Set the folder path containing the raw CSV files (update this path as needed)
folder_path <- "~/Data/Raw"

# Read and combine the CSV files
combined_data <- concatenate_csv_files(folder_path)

# Update the Parent column: standardize the hypoxia labels based on conditions
combined_data <- combined_data %>%
  mutate(Parent = case_when(
    Parent %in% c("Normoxia", "Annotation (Non hypoxic)") ~ "Normoxia",
    Parent %in% c("Mild_Hypoxia", "Annotation (Low hypoxia)") ~ "Mild_Hypoxia",
    Parent %in% c("Hypoxia", "Annotation (High hypoxia)") ~ "Hypoxia",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Parent))  # Remove rows with NA Parent

# Create a unique identifier combining TMA and Image
combined_data$Identifier <- paste0(combined_data$TMA, "_", combined_data$Image)

# Optional: inspect the unique Identifiers
print(unique(combined_data$Identifier))

# Save the processed combined data as an RDS file
saveRDS(combined_data, "~/Data/Processed/Combined.RDS")
