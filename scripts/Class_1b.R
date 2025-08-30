#create sub folders
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

#load data
patient_data <- read.csv("")

#converting to factors
patient_data$gender<- factor(patient_data$gender)
patient_data$diagnosis <- factor(patient_data$diagnosis,
                                 levels = c("Normal","Cancer"))  #levels in order
patient_data$smoker <- factor(patient_data$smoker,
                              levels = c("No","Yes"))

#binary smoker variable
patient_data$smoker_binary <- ifelse(patient_data$smoker == "Yes",1,0)

str(patient_data)
View(patient_data)

#save cleaned data
write.csv(patient_data,"clean_data/patient_info_clean.csv",row.names = FALSE)
save.image("srikeerthi_Class_1b_Assignment.RData")




