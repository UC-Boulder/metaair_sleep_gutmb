# Header ----
# Title: Sleep and gut microbiota in MetaAIR PCoA Plots
# Alderete_metaair_sleep.R code used to create analytical dataset that includes metadata and sleep variables
# Author: Tanya Alderete
# Date: Dec 2023

# Loading Packages ----
library(tidyverse); library(dplyr); library(lubridate); 
library(tibble); library(ggplot2); library(qwraps2);
library(magrittr); library(gridExtra); library(vegan);
library(ape); library(ggrepel)

# Create theme for graphs ---
mytheme <- theme(axis.line = element_line(size = 0.5, colour = "black"),
                 panel.background = element_rect(fill = "white"))

mytheme.2 <- theme(axis.line = element_line(size = 2, colour = "black"),
                   panel.background = element_rect(fill = "white"))

# Clear workspace ----
rm(list=ls())

# Set filepaths ----
meta <-"~/Git/metaair_sleep_gutmb/"
meta2 <-"/Volumes/ADORLab/Lab Projects/MetaAIR/DADA2/Split_Tables/"

# Importing Data ----
metaair <- read.csv(paste0(meta,"metaair_sleep_output.csv"))
species <- read.csv(paste0(meta2,"species_MMDADA2_summary.csv"))

# Extract lab ID ----
species$X <- sub("^[^\\.]*\\.([^\\.]*)\\..*$", "\\1", species$X)

length(species$X %in% metaair$id_lab)

# Limit species df to those with sleep data ----
trimmed_species <- species %>%
  semi_join(metaair, by = c("X" = "id_lab"))

# Limit metadata df to those with species data ----
trimmed_metaair <-  metaair %>%
  semi_join(species, by = c("id_lab" = "X"))

# Set the ID column ('X' in this case) as row names
row.names(trimmed_species) <- trimmed_species$X
trimmed_species <- trimmed_species %>% select(-X)
trimmed_species[] <- lapply(trimmed_species, as.numeric)

# Create distance Matrix ----
bray_curtis_matrix <- vegdist(trimmed_species, method = "bray")

# To view the matrix
print(bray_curtis_matrix)

# Create PCoA ----
pcoa_result <- pcoa(bray_curtis_matrix, correction = "none")

# Extracting PCoA scores ----
scores <- as.data.frame(pcoa_result$vectors)

# Ensure identifiers match ----
rownames(scores)
scores <- rownames_to_column(scores, var = "id_lab")

scores$id_lab <- as.character(scores$id_lab)
trimmed_metaair$id_lab <- as.character(trimmed_metaair$id_lab)

# Check for identifiers in group_info not in pcoa_scores
setdiff(trimmed_metaair$id_lab, scores$id_lab)

# Check for identifiers in pcoa_scores not in group_info
setdiff(scores$id_lab, trimmed_metaair$id_lab)

# Sort
scores <- scores[order(scores$id_lab), ]
trimmed_metaair <- trimmed_metaair[order(trimmed_metaair$id_lab), ]

# Merge the data frames ----
merged_data <- merge(scores, trimmed_metaair, by = "id_lab", all = FALSE)

# Check the first few rows of the merged data frame
head(merged_data)

# Make factor variables ----
merged_data$snoring <- as.factor(merged_data$snoring)
merged_data$SleepDebtShort <- as.factor(merged_data$SleepDebtShort)

# Basic PCoA plot Snoring ----
pcoa_plot <- ggplot(merged_data, aes(x = Axis.1, y = Axis.2, color=snoring)) +
  geom_point() +  # Plots the points
  xlab("PCoA 1") + ylab("PCoA 2") +
  theme_minimal() +
  geom_text_repel(aes(label = merged_data$id_lab)) # Adds labels to points

# Display the plot
print(pcoa_plot)

# Basic PCoA plot SleepDebtShort ----
pcoa_plot <- ggplot(merged_data, aes(x = Axis.1, y = Axis.2, color=SleepDebtShort)) +
  geom_point() +  # Plots the points
  xlab("PCoA 1") + ylab("PCoA 2") +
  theme_minimal() +
  geom_text_repel(aes(label = merged_data$id_lab)) # Adds labels to points

# Display the plot
print(pcoa_plot)