# 1. FIRST, RUN THIS TO CLEAR ANY ERRORS
rm(list = ls())  # Clears all objects from environment

# 2. LOAD DATA PROPERLY (choose one method)

# METHOD A: Using file.choose() - interactive selection
raw_data <- read.csv(file.choose())  # Navigate to CORRECT file (patient_info.csv)

# METHOD B: Using exact path (fix the double extension)
raw_data <- read.csv("C:/Users/eggce/OneDrive/Ai_omics_internship_2025/class_1c/raw_data/patient_info.csv")

# 3. VERIFY SUCCESSFUL LOAD
head(raw_data)  # Should show first 6 rows
str(raw_data)   # Check column structure
View(raw_data)  # Spreadsheet view - must use object name, not filename

# 1. Convert character columns to factors
clean_data <- raw_data
clean_data$gender <- as.factor(clean_data$gender)
clean_data$diagnosis <- as.factor(clean_data$diagnosis)
clean_data$smoker <- as.factor(clean_data$smoker)

# 2. Verify changes
str(clean_data)  # Check gender/diagnosis/smoker now show as 'Factor'
View(clean_data)  # Visual confirmation

# 3. Save cleaned data (optional)
write.csv(clean_data, "clean_patient_data.csv", row.names = FALSE)