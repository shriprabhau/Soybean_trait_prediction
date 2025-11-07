# load packages
library(tidyverse)
library(ggplot2)

# read in data
pca <- read.table("full_prt_filtered_results.eigenvec", header = FALSE)
eigenval <- scan("full_prt_filtered_results.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
pca$ind <- sub("_.*", "", pca$ind)

# first convert to percentage variance explained
pve <- data.frame(pve = eigenval/sum(eigenval)*100)



# Load species metadata from a CSV file
metadata <- read.table("../../info.txt", sep="\t", header = TRUE)

# Merge PCA results with species metadata
pca_merged <- merge(pca, metadata, by = "ind")

head(pca_merged)
                    
#### If you want to add labels####
pca_plot <- ggplot(pca_merged, aes(x = PC1, y = PC2, color = Group, label = Country)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +  # Add individual labels
  labs(title = "PCA Plot with Species type",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()
print(pca_plot)

#### PCA plot for grouping under varieties ####
pca_plot <- ggplot(pca_merged, aes(x = PC1, y = PC2, color = Country, shape = Group)) +
  geom_point(size = 3) +
    # Add individual labels
  labs(title = "PCA Plot with varieties",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

print(pca_plot)

### PCA plot to combine both varieties and country ###
pca_plot <- ggplot(pca_merged, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  # Add individual labels
  labs(title = "PCA Plot with varieties",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()
print(pca_plot)

## make the legends bigger

pca_plot <- ggplot(pca_merged, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  labs(title = "PCA Plot with varieties",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  # Bigger legend title
    legend.text  = element_text(size = 14), # Bigger legend labels
    axis.title   = element_text(size = 13),  # Axis titles
    axis.text    = element_text(size = 12)
  )

print(pca_plot)

ggsave("PCA_cluster.png", 
       dpi = 300, 
       width = 8, 
       height = 6, 
       units = "in")
##################################
# PLOT PLINK PCA eigenvals file #
##################################

# Load the eigenvalues from a file
eigenval <- read.table("full_prt_filtered_results.eigenval", header = FALSE)

# Calculate the percentage of variance explained by each component
percent_variance <- (eigenval$V1 / sum(eigenval$V1)) * 100
cumulative_variance <- cumsum(percent_variance)

# Create a data frame with principal component numbers, eigenvalues, and percentages
eigen_data <- data.frame(
  PC = 1:nrow(eigenval),
  Eigenvalue = eigenval$V1,
  PercentVariance = round(percent_variance, 2) # rounding to 2 decimal places
)

# Generate the bar plot

ggplot(eigen_data, aes(x = factor(PC), y = Eigenvalue)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = paste0(PercentVariance, "%")), 
            vjust = -0.5, size = 3.0) + # Adding percentage labels on top of bars
  labs(title = "PCA Eigenvalues with Percentage Variance Explained",
       x = "Principal Component",
       y = "Eigenvalue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Plot only the first 25 values ##

# Filter the data to include only the first X-number of PCs
eigen_data_top10 <- eigen_data[1:25, ]


# Generate the bar plot

ggplot(eigen_data_top10, aes(x = factor(PC), y = Eigenvalue)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = paste0(PercentVariance, "%")), 
            vjust = -0.5, size = 3.0) + # Adding percentage labels on top of bars
  labs(title = "PCA Eigenvalues (Top 25 PCs) with Percentage Variance Explained",
       x = "Principal Component",
       y = "Eigenvalue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Have the labels one below one below so they are not mushed together


# Create alternating up/down positions for labels
eigen_data_top10$label_pos <- rep(c(-0.5, 1.5), length.out = nrow(eigen_data_top10))

ggplot(eigen_data_top10, aes(x = factor(PC), y = Eigenvalue)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = paste0(PercentVariance), vjust = label_pos),
            size = 4.0, color = "black") + 
  labs(title = "PCA Eigenvalues (Top 25 PCs) with Percentage Variance Explained",
       x = "Principal Component",
       y = "Eigenvalue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("PCA_Eigenvalues.png", 
       dpi = 300, 
       width = 8, 
       height = 6, 
       units = "in")
ggsave("PCA_Eigenvalues.pdf", 
       width = 8, 
       height = 6, 
       units = "in")
