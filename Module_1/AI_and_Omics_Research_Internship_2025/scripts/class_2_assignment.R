#load data
raw_data1 <- read.csv("C:/Users/eggce/oneDrive/AI_Omics_Internship_2025/Module_1/AI_and_Omics_Research_Internship_2025/raw_data/DEGs_Data_1.csv")
raw_data2 <- read.csv("C:/Users/eggce/oneDrive/AI_Omics_Internship_2025/Module_1/AI_and_Omics_Research_Internship_2025/raw_data/DEGs_Data_2.csv")

# For raw_data1
na_count_1 <- sum(is.na(raw_data1$padj))
cat("Missing 'padj' values in raw_data1:", na_count_1, "\n")

# For raw_data2
na_count_2 <- sum(is.na(raw_data2$padj))
cat("Missing 'padj' values in raw_data2:", na_count_2, "\n")


# Define the function to classify genes
classify_gene <- function(logFC, padj) {
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# Add status column to raw_data1
raw_data1$status <- mapply(classify_gene, raw_data1$logFC, raw_data1$padj)

# Add status column to raw_data2
raw_data2$status <- mapply(classify_gene, raw_data2$logFC, raw_data2$padj)

cat("Gene classification complete for both datasets.\n")


# Summary for raw_data1
cat("Summary for DEGs_data_1.csv:\n")
print(table(raw_data1$status))

# Summary for raw_data2
cat("\nSummary for DEGs_data_2.csv:\n")
print(table(raw_data2$status))

# Save processed data to CSV files
write.csv(raw_data1, file.path(output_dir, "Processed_DEGs_data_1.csv"), row.names = FALSE)
write.csv(raw_data2, file.path(output_dir, "Processed_DEGs_data_2.csv"), row.names = FALSE)
cat("Processed CSV files saved to Results folder.\n")

# Save all required objects 
save(raw_data1, raw_data2, classify_gene, file = "yourname_Class_2_Assignment.RData")
cat("RData file 'yourname_Class_2_Assignment.RData' created.\n")

# Save the cleaned data as CSV

write.csv(raw_data1, "Cleaned_DEGs_Data.csv", row.names = FALSE)
getwd()
shell.exec("Cleaned_DEGs_Data.csv")
write.csv(raw_data2, "Cleaned_DEGs_Data.csv", row.names = FALSE)
getwd()
shell.exec("Cleaned_DEGs_Data.csv")