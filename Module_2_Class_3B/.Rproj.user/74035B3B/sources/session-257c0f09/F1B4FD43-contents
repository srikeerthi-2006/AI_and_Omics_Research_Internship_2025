# =====================================================================
#               AI and Omics Research Internship (2025)
# =====================================================================
#             Module II: Introduction to Genomics Data Analysis
# ---------------------------------------------------------------------
#                     Microarray Data Analysis - GSE285348
# =====================================================================

#######################################################################
#### 0. Install and Load Required Packages ####
#######################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery","affy","dplyr"))

# Load Required Libraries
library(GEOquery)             # Download GEO datasets
library(affy)                 # Pre-processing of Affymetrix microarray data
library(dplyr)                # Data manipulation

#######################################################################
#### 1. Downloading dataset (GSE285348) ####
#######################################################################

# -------------------------------------
#### Download Series Matrix Files ####
# -------------------------------------

gse_data <- getGEO("GSE285348", GSEMatrix = TRUE)

# Extract data
expression_data <- exprs(gse_data[[1]])
feature_data <- fData(gse_data[[1]])
phenotype_data <- pData(gse_data[[1]])

# Check missing values
sum(is.na(phenotype_data$source_name_ch1)) 

# --------------------------------------
#### Download Raw Data (CEL files)
# --------------------------------------
# --------------------------------------------------------------------
# Dataset: GSE285348 uses Clarion_S_Human arrays
# since Bioconductor annotation package 'clariomshumancdf' is not available
#        for current Bioconductor version (3.21)
#Hence, I have used Used pre-normalized expression data from GEO repository,
#           which has undergone standard RMA normalization pipeline.
# Reference: GEO automatically processes all data through standardized
#           Affymetrix pipelines including RMA normalization.
# --------------------------------------------------------------------

#######################################################################
#### 2. Quality Control (QC) Using Pre-normalized Data ####
#######################################################################

# ---------------------------------------------------
#### QC Before Processing ####
# ---------------------------------------------------

# Use pre-normalized data from GEO (already processed with RMA)
processed_data <- as.data.frame(expression_data)

# Basic QC plots
boxplot(expression_data, 
        main = "GSE285348 - Expression Distribution", 
        col = "lightblue",
        las = 2)

hist(expression_data, 
     main = "GSE285348 - Density Plot",
     xlab = "Expression Values")

#metrics
print("=== QC METRICS ===")
print(paste("Number of samples:", ncol(expression_data)))
print(paste("Number of probes:", nrow(expression_data)))
print(paste("Data dimensions:", dim(expression_data)))

# Sample intensity statistics
sample_means <- colMeans(expression_data)
print("Sample mean intensities:")
print(sample_means)

#######################################################################
#### 3. Filter Low-Variance Transcripts ####
#######################################################################

# Calculate median intensity
row_median <- apply(expression_data, 1, median)

# Visualize distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "GSE285348 - Median Intensity Distribution",
     xlab = "Median Intensity")

# Set a threshold to remove low variance probes
threshold <- quantile(row_median, 0.2)  # Remove bottom 20% low intensity probes
abline(v = threshold, col = "red", lwd = 2) 
legend("topright", legend = paste("Threshold =", round(threshold, 2)), col = "red", lwd = 2)

# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- expression_data[indx, ] 

print("=== FILTERING RESULTS ===")
print(paste("Original number of probes:", nrow(expression_data)))
print(paste("Probes after filtering:", nrow(filtered_data)))
print(paste("Probes removed:", nrow(expression_data) - nrow(filtered_data)))

#filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

# Overwrite processed data with filtered dataset
processed_data <- as.data.frame(filtered_data)

#######################################################################
#### 4. Phenotype Data Preparation ####
#######################################################################

# Check phenotype structure
print("=== PHENOTYPE DATA ===")
print(head(phenotype_data$source_name_ch1))
print(unique(phenotype_data$source_name_ch1))

# Define experimental groups based on metadata
print(colnames(phenotype_data))

print(unique(phenotype_data$characteristics_ch1.2))

print(unique(phenotype_data$characteristics_ch1))
print(unique(phenotype_data$characteristics_ch1.1))

print(unique(phenotype_data$description))
print(unique(phenotype_data$title))

# Create groups based on PWS vs control
groups <- factor(phenotype_data$description,
                 levels = c("not PWS,not obese", "not PWS,obese", "PWS,not obese", "PWS,obese"),
                 labels = c("control_NW", "control_obese", "PWS_NW", "PWS_obese"))

print("Experimental groups:")
print(table(groups))


#######################################################################
#### 5. Save Processed Data ####
#######################################################################

# Save the filtered expression data
write.csv(processed_data, "GSE285348_filtered_expression.csv")

# Save phenotype data with groups
phenotype_data$experimental_group <- groups
write.csv(phenotype_data, "GSE285348_phenotype_data.csv")


print("=== PROCESSING COMPLETE ===")
print(paste("Final dataset:", nrow(processed_data), "probes x", ncol(processed_data), "samples"))
