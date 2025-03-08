##Reading in all data files
Fung.abund <- read.csv("/cloud/project/QB2025_Guevara/FINAL/data/MAT_fungal_abundances.csv")
Fung.enzymes <- read.csv("/cloud/project/QB2025_Guevara/FINAL/data/MAT_enzymes.csv")
Fung.ectocounts <- read.csv("/cloud/project/QB2025_Guevara/FINAL/data/MAT_ectocounts.csv")
Fung.seq.enzymes <- read.csv("/cloud/project/QB2025_Guevara/FINAL/data/MAT_sequences_enzymes.csv")
soils <- read.csv("/cloud/project/QB2025_Guevara/FINAL/data/MAT_soils.csv")

#To slim down the data to only columns we will need to incidence matrix
library(dplyr)
slim <- Fung.abund %>%
  select(Plot, Treatment, Species, Relative.Abundance)
slim
#Changing data so that each species is a column, values are the relative.abundance, and rows are each plot/site
library(tidyr)
Fung.abund <- slim %>%
  pivot_wider(names_from = Species, values_from = Relative.Abundance, values_fill = 0)
# Convert nonzero values to 1 while keeping 0s as is
Fung.abund <- Fung.abund %>%
  mutate(across(-c(Plot, Treatment), ~ ifelse(. > 0, 1, 0)))
#Creating complete incidence matrix 
Fung.abund.incidence <- Fung.abund %>%
  arrange(Treatment)
##Counting fungal presence
fungal_presence_counts <- Fung.abund.incidence %>%
  mutate(Presence_Count = rowSums(across(-c(Plot, Treatment))))


##ANOVA for difference in # of species
plot.incid.diff <- aov(formula = Presence_Count ~ Treatment, data = fungal_presence_counts)
summary(plot.incid.diff)

# Create a bar plot showing the number of species present per plot
ggplot(fungal_presence_counts, aes(x = reorder(Plot, -Presence_Count), y = Presence_Count, fill = Treatment)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Fungal spp. Presence/Plot",
       x = "Plot",
       y = "Number of Species Present") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size, bold, and center align
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.title.x = element_text(size = 16, face = "bold"),  # Increase x-axis label size
    axis.title.y = element_text(size = 16, face = "bold")   # Increase y-axis label size
  ) +
  scale_fill_brewer(palette = "Set2")

ggplot(fungal_presence_counts, aes(x = Treatment, y = Presence_Count)) +
  geom_boxplot(aes(fill = Treatment), color = "black") +
  labs(title = "Spp. Richness by Treatment",
       x = "Treatment", 
       y = "Spp. Richness") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##Now, we will go ahead and create clustering tree
install.packages("ggdendro")
library(ggdendro)

library(tibble)

# Convert tibble to a standard data frame to avoid row name issues
Fung.abund.incidence <- as.data.frame(Fung.abund.incidence)

# Set row names as Plot names and remove 'Plot' and 'Treatment' columns
rownames(Fung.abund.incidence) <- Fung.abund.incidence$Plot
fungal_matrix <- as.matrix(Fung.abund.incidence[, -c(1, 2)])  # Keep only presence/absence data

# Compute Jaccard distance matrix
jaccard_dist <- vegdist(fungal_matrix, method = "jaccard")

# Perform hierarchical clustering
jaccard_clust <- hclust(jaccard_dist, method = "average")  # UPGMA clustering

# Convert clustering object to dendrogram format
dendro <- as.dendrogram(jaccard_clust)
dendro_data <- ggdendro::dendro_data(dendro)

# Ensure correct plot labels
dendro_data$labels$label <- rownames(fungal_matrix)[jaccard_clust$order]

ggplot(dendro_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(
    data = dendro_data$labels,
    aes(x = x, y = y, label = label),
    hjust = 1, angle = 90, size = 5, nudge_y = -0.2  # Adjust branch label position
  ) +
  theme_minimal(base_size = 16) +  # Increase overall text size
  labs(
    title = "Fungal Community Clustering (Jaccard Similarity)",
    x = "",
    y = "Distance (Jaccard)"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 20)),  # Add space below title
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(1, 1, 2, 1, "cm"),  # Adjust overall plot margins
    axis.title.x = element_text(margin = margin(t = 20, b = 20)),  # Add space around x-axis title
    axis.title.y = element_text(margin = margin(r = 20))  # Add space around y-axis title
  ) +
  coord_cartesian(clip = "off")  # Prevent text cutoff

# Count the number of Root.tip.IDs per Treatment
treatment_counts <- Fung.seq.enzymes %>%
  group_by(Treatment) %>%
  summarise(Count = n_distinct(Root.tip.ID))
##There were 165 root tips in control group, 161 root tips in fertilized group, 161 tips in warmed, and 159 in warmed.fertilized
# Count the number of Root.tip.IDs per Plot
plot_counts <- Fung.seq.enzymes %>%
  group_by(Plot) %>%
  summarise(Count = n_distinct(Root.tip.ID))

# View the result
print(plot_counts)

# Getting actual abundances per plot
Abundances <- Fung.seq.enzymes %>%
  group_by(Plot, OTU) %>%
  summarise(Abundance = n_distinct(Root.tip.ID)) %>%
  arrange(Plot, OTU) %>%
  filter(!is.na(OTU))  # Remove rows with NA in Abundance

Abundances_transposed <- Abundances %>%
  pivot_wider(
    names_from = OTU,           # Make OTUs the column names
    values_from = Abundance,    # Fill values with Abundance
    values_fill = list(Abundance = 0)  # Fill missing values with 0
  )

# Select unique Plot-Treatment pairs from Fung.seq.enzymes
treatment_data <- Fung.seq.enzymes %>%
  select(Plot, Treatment) %>%
  distinct()  # Ensure only one row per Plot

# Merge the transposed Abundances matrix with the treatment data
merged_data <- Abundances_transposed %>%
  left_join(treatment_data, by = "Plot") %>%
  select(Plot, Treatment, everything())  # Ensure correct column order

# View the result
print(merged_data)



long_data <- merged_data %>%
  pivot_longer(cols = -c(Plot, Treatment), names_to = "Species", values_to = "Count")

# Summarize total individuals per treatment
summary_data <- long_data %>%
  group_by(Treatment, Plot) %>%
  summarise(Total_Individuals = sum(Count, na.rm = TRUE), .groups = 'drop')

# Perform ANOVA
anova_result <- aov(Total_Individuals ~ Treatment, data = summary_data)

# Print summary of ANOVA
summary(anova_result)

# Perform Tukey's HSD test for pairwise comparisons (optional)
TukeyHSD(anova_result)

# Extract the first two columns (Plot and Treatment)
first_two_columns <- merged_data[, 1:2]

# Extract the remaining columns (OTUs)
otu_columns <- merged_data[, -c(1, 2)]

# Convert OTU column names to numeric and sort them
sorted_otu_columns <- otu_columns[, order(as.numeric(names(otu_columns)))]

# Combine the first two columns with the sorted OTU columns
merged_data_sorted <- cbind(first_two_columns, sorted_otu_columns)

# View the sorted dataframe
print(merged_data_sorted)



# Calculate the Bray-Curtis dissimilarity matrix
bray_curtis_matrix <- vegdist(merged_data, method = "bray")

# View the Bray-Curtis matrix
print(bray_curtis_matrix)

# Perform hierarchical clustering
cluster_result <- hclust(bray_curtis_matrix, method = "ward.D2")  # You can use "ward.D", "complete", etc.

# Plot the dendrogram
plot(cluster_result, main = "Cluster Dendrogram (Bray-Curtis)", xlab = "", sub = "", cex = 0.9, labels = rownames(merged_data))

# Optional: Customize the dendrogram with colors and clusters
library(dendextend)
clusters <- cutree(cluster_result, k = 4)  # Cut into 4 clusters
dend <- as.dendrogram(cluster_result)
dend <- color_branches(dend, k = 4)  # Color branches by clusters

# Add labels to the dendrogram
labels(dend) <- rownames(merged_data)[order.dendrogram(dend)]

# Plot the customized dendrogram
plot(dend, main = "Cluster Dendrogram (Bray-Curtis)", ylab = "Squared Bray-Curtis Distance", sub = "", cex = 0.9)


library(dplyr)
library(vegan)

# Remove non-OTU columns and ungroup to avoid issues
otu_data <- merged_data_sorted %>%
  ungroup() %>%  # Remove any grouping that may have been applied
  select(-Plot, -Treatment) %>%
  mutate(across(everything(), as.numeric))  # Ensure all columns are numeric

# Compute Shannon Diversity Index for each plot
merged_data_sorted$Shannon_Index <- diversity(as.matrix(otu_data), index = "shannon")

# View results
print(merged_data_sorted[, c("Plot", "Treatment", "Shannon_Index")])

# Ensure data is properly ordered by Treatment
merged_data_sorted <- merged_data_sorted %>%
  arrange(Treatment, Plot)  # Sorting by Treatment first, then Plot

# Convert Plot to a factor to maintain order in the plot
merged_data_sorted$Plot <- factor(merged_data_sorted$Plot, levels = merged_data_sorted$Plot)

# Create the histogram (bar plot)
ggplot(merged_data_sorted, aes(x = Plot, y = Shannon_Index, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black") +  # Bar plot with black borders
  theme_minimal() +
  labs(x = "Plot", y = "Shannon Diversity Index", title = "Shannon Diversity Index by Plot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for readability
  scale_fill_manual(values = c("Control" = "blue", "Fertilized" = "green", 
                               "Warmed" = "orange", "Warmed.Fertilized" = "red"))  # Custom colors

# Perform ANOVA
shannon_anova <- aov(Shannon_Index ~ Treatment, data = merged_data_sorted)

# View ANOVA summary
summary(shannon_anova)
TukeyHSD(shannon_anova)

# Create the boxplot with the means
ggplot(merged_data_sorted, aes(x = Treatment, y = Shannon_Index)) +
  geom_boxplot(aes(fill = Treatment), color = "black") +
  labs(title = "Boxplot of Shannon Index by Treatment",
       x = "Treatment", 
       y = "Shannon Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
