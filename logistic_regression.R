##############################################
# Multivariable Logistic Regression Analysis
# Author: Zhichao Liang
# Description:
# This script performs multivariable logistic regression
# to evaluate the association between intestinal diseases
# and Alzheimer's disease (AD), adjusting for key covariates.
##############################################

# Set CRAN mirror
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Load required packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("broom")) install.packages("broom")
library(tidyverse)
library(broom)

# Load the dataset
data <- read.csv("cross_section_data.csv", stringsAsFactors = FALSE)

# Recode and ensure correct variable types
data <- data %>%
  mutate(
    AD = as.factor(AD),
    Sex = as.factor(Sex),
    Education = as.factor(Education),
    Racitizen = as.factor(Racitizen),
    Cancer = as.factor(Cancer),
    Diabetes = as.factor(Diabetes),
    Hypertension = as.factor(Hypertension),
    Heart_attack = as.factor(Heart_attack),
    Cerebrovascular = as.factor(Cerebrovascular),
    Health = as.factor(Health),
    Weight = as.numeric(Weight)
  )

# Convert intestinal disease variables to factors
intestinal_vars <- c("BCs", "IBD", "IBS", "CST", "IM", "Intestinal_disease")
data[intestinal_vars] <- lapply(data[intestinal_vars], as.factor)

# Descriptive summary of AD prevalence across intestinal diseases
disease_summary <- data %>%
  group_by(BCs, IBD, IBS, CST, IM, Intestinal_disease) %>%
  summarise(
    Total = n(),
    AD_1_Count = sum(as.numeric(as.character(AD)) == 1, na.rm = TRUE),
    AD_1_Percentage = round(AD_1_Count / Total * 100, 2),
    .groups = "drop"
  )

# Print disease-wise AD summary
print(disease_summary)

# Fit logistic regression model
logit_model <- glm(
  AD ~ BCs + IBD + IBS + CST + IM +
    Age + Sex + Racitizen + Education +
    Cancer + Diabetes + Hypertension + Heart_attack +
    Cerebrovascular + Health + Weight,
  data = data,
  family = binomial(link = "logit")
)

# Summarize model
summary(logit_model)

# Extract odds ratios and confidence intervals
logit_results <- tidy(logit_model, conf.int = TRUE) %>%
  mutate(
    OR = exp(estimate),
    OR_lower = exp(conf.low),
    OR_upper = exp(conf.high)
  ) %>%
  select(term, OR, OR_lower, OR_upper, p.value)

# Print results
print(logit_results, n = Inf)

# Save results to CSV
write.csv(logit_results, "logistic_results.csv", row.names = FALSE)

cat("\n--- Logistic regression results saved to 'logistic_results.csv' ---\n")
