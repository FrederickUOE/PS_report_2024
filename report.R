#libraries used
library(dplyr)
library(stringr)
library(readxl)
library(ggplot2)
library(broom)

#question: Do the biomarker levels at in- clusion for patients with high VAS (≥ 5) differ from those for patients with low VAS (< 5)?

#reading in the data
covariates <- read_excel("data/covariates.xlsx")
biomarkers <- read_excel("data/biomarkers.xlsx")

#First inspection of the datasets 
head(covariates)
head(biomarkers)

#summary and check for missing values
summary(covariates)
summary(biomarkers)

str(covariates)
str(biomarkers)

#get only the inclusion data from the biomarkers
week0_data <- biomarkers[grep("0weeks", biomarkers$Biomarker), ]

#add a new column for patientID to the new data frame and extract the id from the Biomarker column
week0_data <- week0_data %>%
  mutate(PatientID = as.numeric(str_extract(Biomarker, "^[^-]+")))

#combining the datasets to the new dataset called patient_data
patient_data <- week0_data %>%
  right_join(covariates, by = "PatientID")

#control new dataframe
head(patient_data)
summary(patient_data)
str(patient_data)

#inspect missing values
na_rows <- patient_data %>%
  filter(rowSums(is.na(.)) > 0)
na_patient_ids <- na_rows$PatientID
na_patient_ids

#removal of entry with patientID 40, since it is missing inclusion data 
patient_data <- patient_data %>%
  filter(PatientID != 40)

#checking for outliers for biomarkers
biomarker_names <- c("IL-8", "VEGF-A", "OPG", "TGF-beta-1","IL-6", "CXCL9", "CXCL1", "IL-18", "CSF-1")

for(biomarker in biomarker_names){
  p <- ggplot(patient_data, aes( x = factor(1), y = .data[[biomarker]])) +
    geom_boxplot() + ggtitle(paste("Boxplot for", biomarker)) +
    xlab("") +
    ylab(biomarker)
  print(p)
}

#checking outliers for VAS score at inclusion
ggplot(patient_data, aes( x = factor(1), y = .data[["VAS-at-inclusion"]])) +
  geom_boxplot() + ggtitle("Boxplot for VAS at inclusion") +
  xlab("") +
  ylab(("VAS at inclusion"))

#compute some basic stats for each biomarker and the VAS score 
# calculate mean, median, standard deviation, interquartile range, minimum, maximum
basic_stats_biomarkers <- patient_data %>%
  select(all_of(biomarker_names)) %>%
  summarise(across(everything(), list(
    mean = ~mean(., na.rm = TRUE),
    median = ~median(., na.rm = TRUE),
    sd = ~sd(., na.rm = TRUE),
    IQR = ~IQR(., na.rm = TRUE),
    min = ~min(., na.rm = TRUE),
    max = ~max(., na.rm = TRUE),
    n = ~sum(!is.na(.))
  )))

basic_stats_vas <- patient_data %>%
  select(all_of("VAS-at-inclusion")) %>%
  summarise(across(everything(), list(
    mean = ~mean(., na.rm = TRUE),
    median = ~median(., na.rm = TRUE),
    sd = ~sd(., na.rm = TRUE),
    IQR = ~IQR(., na.rm = TRUE),
    min = ~min(., na.rm = TRUE),
    max = ~max(., na.rm = TRUE),
    n = ~sum(!is.na(.))
  )))


#scatter plot to look further into the relationship of the VAS score and the biomarkers
#creating a plot for each biomarker against the VAS Score
for(biomarker in biomarker_names) {
  p <- ggplot(patient_data, aes(x = .data$'VAS-at-inclusion', y = .data[[biomarker]])) + 
    geom_point(alpha = 0.5) + 
    geom_smooth(method = "lm", color = "blue", se = FALSE) + 
    ggtitle(paste("Scatter Plot of", biomarker, "vs VAS Score at inclusion")) +
    xlab("VAS Score") + 
    ylab(biomarker)
  print(p)
}

#create histrograms

#biomarkers
for(biomarker in biomarker_names) {
  p <- ggplot(patient_data, aes(x = .data[[biomarker]])) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "blue", color = "black") +
    geom_density(alpha = .2, fill = "#FF6666") +
    ggtitle(paste("Histogram and Density Plot for", biomarker))
    
  print(p)
}

#VAS score at inclusion
ggplot(patient_data, aes(x = .data[["VAS-at-inclusion"]])) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "blue", color = "black") +
  geom_density(alpha = .2, fill = "#FF6666") +
  ggtitle("Histogram and Density Plot for VAS at inclusion" )

  
#Statistical test for normality (biomarkers)
for(biomarker in biomarker_names) {
  print(paste(biomarker, "Shapiro-Wilk Test:"))
  print(shapiro.test(biomarkers[[biomarker]]))
}

#categorize patients based on VAS Score
patient_data <- patient_data %>%
  mutate(VAS_group = if_else(`VAS-at-inclusion` >= 5, "high", "low"))


#Hypothesis Testing (Mann-Whitney U test) for all biomarkers
test_results <- list()
for(biomarker in biomarker_names){
  formula <- as.formula(paste(patient_data[biomarker], "~ VAS_group"))
  test_results[[biomarker]] <- wilcox.test(formula, data = patient_data)
}

test_results$`IL-8`
test_results$`VEGF-A`
test_results$OPG
test_results$`TGF-beta-1`
test_results$`IL-6`
test_results$CXCL9
test_results$CXCL1
test_results$`IL-18`
test_results$`CSF-1`

#Calculating the median by VAS Group
medians_by_group <- patient_data %>%
  group_by(VAS_group) %>%
  summarise(across(all_of(biomarker_names), ~ median(., na.rm = TRUE), .names = "median_{.col}"))

print(medians_by_group)


# Regression model 

#splitting the data

biomarker_names <- c("`IL-8`", "`VEGF-A`", "OPG", "`TGF-beta-1`","`IL-6`", "CXCL9", "CXCL1", "`IL-18`", "`CSF-1`")
set.seed(747)
training_indices <- sample(1:nrow(patient_data), size = 0.8 * nrow(patient_data))

training_data <- patient_data[training_indices, ]
test_data <- patient_data[-training_indices, ]

#fitting the model 
predictors <- c(biomarker_names, "Age", "`Sex (1=male, 2=female)`", "`Smoker (1=yes, 2=no)`")
model_formula_string <- paste("`Vas-12months` ~", paste(predictors, collapse = " + "))
model_formula <- as.formula(model_formula_string)
model <- lm(model_formula, data = training_data)

summary(model)

#present fitted parameters
parameter_table <- broom::tidy(model)
write.csv(parameter_table, "model_summary.csv", row.names = FALSE)

# out of sample evaluation 
test_predictions <- predict(model, newdata = test_data)

test_data$predicted_VAS_Score_12months <- test_predictions

test_data$prediction_error <- test_data$`Vas-12months` - test_data$predicted_VAS_Score_12months

write.csv(test_data, "test_data.csv", row.names = FALSE)

summary_stats_predictions <- summary(test_data$prediction_error)
print(summary_stats_predictions)

MAE <- mean(abs(test_data$prediction_error), na.rm = TRUE)

print(MAE)

