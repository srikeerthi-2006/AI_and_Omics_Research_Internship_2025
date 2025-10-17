# =====================================================================
#               AI and Biotechnology / Bioinformatics
# =====================================================================

# ---------------------------------------------------------------------
#              AI and Omics Research Internship (2025)
# ---------------------------------------------------------------------
#             Module II: Introduction to Genomics Data Analysis
# ---------------------------------------------------------------------
#                     Microarray Data Analysis
# =====================================================================


library(GEOquery)
library(affy)
library(dplyr)


gse_data <- getGEO("GSE285348", GSEMatrix = TRUE)
platform <- annotation(gse_data[[1]])
print(paste("Platform annotation needed:", platform))


BiocManager::install("clariomdhumantranscriptcluster.db", ask = FALSE)
library(clariomdhumantranscriptcluster.db)





gse_data <- getGEO("GSE285348", GSEMatrix = TRUE)
expression_data <- exprs(gse_data[[1]])
phenotype_data <- pData(gse_data[[1]])


row_median <- apply(expression_data, 1, median)
threshold <- quantile(row_median, 0.2)
indx <- row_median > threshold 
filtered_data <- expression_data[indx, ]

# Create groups
groups <- factor(phenotype_data$description,
                 levels = c("not PWS,not obese", "not PWS,obese", "PWS,not obese", "PWS,obese"),
                 labels = c("control_NW", "control_obese", "PWS_NW", "PWS_obese"))

# Save the data in current project
processed_data <- as.data.frame(filtered_data)
write.csv(processed_data, "GSE285348_filtered_expression.csv")
write.csv(phenotype_data, "GSE285348_phenotype_data.csv")

print("=== DATA REGENERATED ===")
print(paste("Processed data dimensions:", dim(processed_data)))
print("Groups:")
print(table(groups))




probe_ids <- rownames(processed_data)


gene_symbols <- mapIds(
  clariomdhumantranscriptcluster.db,  # Your platform-specific database
  keys = probe_ids,                   # Input probe IDs
  keytype = "PROBEID",                # Probe ID key type
  column = "SYMBOL",                  # Desired annotation column
  multiVals = "first"                 # Return first match if multiple exist
)


gene_map_df <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = ".")

print("=== PROBE-GENE MAPPING ===")
print(paste("Total probes mapped:", nrow(gene_map_df)))
print(paste("Probes with gene symbols:", sum(!is.na(gene_map_df$SYMBOL))))
print(paste("Probes without gene symbols:", sum(is.na(gene_map_df$SYMBOL))))


duplicate_summary <- gene_map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

print("Multiple probes per gene:")
print(head(duplicate_summary))







processed_data_df <- processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)


processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))


expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)


averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)

print("=== AFTER HANDLING DUPLICATES ===")
print(paste("Original probes with gene symbols:", nrow(processed_data_df)))
print(paste("Unique genes after averaging:", nrow(averaged_data)))
print(paste("Reduction in features:", nrow(processed_data_df) - nrow(averaged_data)))


data_matrix <- as.matrix(averaged_data)

print("Final data matrix ready for differential expression analysis:")
print(dim(data_matrix))



phenotype_data$simple_group <- ifelse(grepl("PWS", phenotype_data$description), "PWS", "Control")
simple_groups <- factor(phenotype_data$simple_group, levels = c("Control", "PWS"))

print("=== SIMPLIFIED GROUPS ===")
print(table(simple_groups))


design <- model.matrix(~0 + simple_groups)
colnames(design) <- levels(simple_groups)

print("Design matrix:")
print(design[1:10,])  # Show first 10 samples





library(limma)


design_final <- model.matrix(~0 + groups)
colnames(design_final) <- c("control_NW", "control_obese", "PWS_NW", "PWS_obese")


fit_final <- lmFit(data_matrix, design_final)


contrast_final <- makeContrasts(
  PWS_vs_Control = (PWS_NW + PWS_obese)/2 - (control_NW + control_obese)/2,
  levels = design_final
)


fit_contrast_final <- contrasts.fit(fit_final, contrast_final)
fit_ebayes_final <- eBayes(fit_contrast_final)

#results
deg_results_final <- topTable(fit_ebayes_final, 
                              coef = "PWS_vs_Control",
                              number = Inf,
                              adjust.method = "BH")


deg_results_final$threshold <- ifelse(
  deg_results_final$adj.P.Val < 0.05 & deg_results_final$logFC > 1, "Upregulated",
  ifelse(deg_results_final$adj.P.Val < 0.05 & deg_results_final$logFC < -1, "Downregulated", "No")
)


upregulated <- sum(deg_results_final$threshold == "Upregulated")
downregulated <- sum(deg_results_final$threshold == "Downregulated")

print("=== FINAL RESULTS ===")
print(paste("Upregulated genes:", upregulated))
print(paste("Downregulated genes:", downregulated))
print(paste("Total DEGs:", upregulated + downregulated))
print("DONE! Assignment complete.")








install.packages("pheatmap")
library(pheatmap)
library(ggplot2)


if(!dir.exists("Results")) dir.create("Results")


write.csv(deg_results_final, file = "Results/DEGs_Results.csv")
write.csv(subset(deg_results_final, threshold == "Upregulated"), file = "Results/Upregulated_DEGs.csv")
write.csv(subset(deg_results_final, threshold == "Downregulated"), file = "Results/Downregulated_DEGs.csv")

# 2. Create Volcano plot
png("Results/volcano_plot.png", width = 2000, height = 1500, res = 300)
ggplot(deg_results_final, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: PWS vs Control",
       x = "log2 Fold Change", y = "-log10(Adjusted P-value)")
dev.off()

# 3. Create Heatmap of top 25 DEGs
deg_updown <- subset(deg_results_final, threshold %in% c("Upregulated", "Downregulated"))
top_25_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 25)
heatmap_data <- data_matrix[top_25_genes, ]

png("Results/heatmap_top25_DEGs.png", width = 2000, height = 1500, res = 300)
pheatmap(heatmap_data,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         fontsize_row = 8,
         main = "Top 25 Differentially Expressed Genes\nPWS vs Control")
dev.off()

print("=== ALL FILES SAVED ===")
print("✓ DEGs_Results.csv")
print("✓ Upregulated_DEGs.csv") 
print("✓ Downregulated_DEGs.csv")
print("✓ volcano_plot.png")
print("✓ heatmap_top25_DEGs.png")





#required packages
library(limma)
library(ggplot2)
library(pheatmap)


library(clariomdhumantranscriptcluster.db)


gene_symbols <- mapIds(
  clariomdhumantranscriptcluster.db,
  keys = top_25_genes,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)


rownames(heatmap_data) <- gene_symbols


png("heatmap_top25_DEGs_FIXED.png", width = 2000, height = 1500, res = 300)
pheatmap(heatmap_data, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         fontsize_row = 8,
         main = "Top 25 DEGs: PWS vs Control")
dev.off()

print("FIXED HEATMAP CREATED: 'heatmap_top25_DEGs_FIXED.png'")
