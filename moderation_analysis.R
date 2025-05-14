###############################################
# Moderation Analysis for Cross-Sectional Data
# Author: Zhichao Liang
# Description:
# This script performs moderation analysis to examine
# whether certain variables moderate the relationship 
# between intestinal diseases and AD.
###############################################

# Load required packages
if (!require("dplyr")) install.packages("dplyr")
if (!require("broom")) install.packages("broom")
library(dplyr)
library(broom)

# Load dataset
data <- read.csv("cross_section_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Ensure Age is numeric
data$Age <- as.numeric(data$Age)

# Create binary variable for any intestinal disease
data <- data %>%
  mutate(
    Any_Intestinal_Disease = ifelse(
      BCs == 1 | IBD == 1 | IBS == 1 | CST == 1 | IM == 1, 1, 0
    )
  )

# Convert selected variables to factor
data <- data %>%
  mutate(
    Cerebrovascular = factor(Cerebrovascular),
    Heart_attack = factor(Heart_attack),
    Hypertension = factor(Hypertension),
    Education = factor(Education)
  )

# Discretize Age into two groups
data <- data %>%
  mutate(
    Age_Group = case_when(
      Age < 60 ~ 1,
      Age >= 60 ~ 2,
      TRUE ~ NA_real_
    ),
    Age_Group = factor(Age_Group)
  )

# Define moderators to test
moderators <- c("Age_Group", "Education", "Heart_attack", "Hypertension")

# Initialize result list
results <- list()

# Loop over each moderator
for (moderator in moderators) {
  cat("\n--- Analyzing Moderator:", moderator, "---\n")
  
  # Other moderators (for adjustment)
  other_moderators <- setdiff(moderators, moderator)
  
  # Construct covariate formula
  covariates <- paste(c("Sex", "Racitizen", "Smoking", "Cancer", "Diabetes", "Cerebrovascular",
                        "Weight", "Health", other_moderators), collapse = " + ")
  
  # Full model formula
  model_formula <- as.formula(paste("AD ~ Any_Intestinal_Disease +", covariates))
  
  # Group by moderator and fit models
  group_results <- data %>%
    group_by(!!sym(moderator)) %>%
    group_modify(~ {
      subset_data <- .x
      
      # Calculate summary counts
      stats <- subset_data %>%
        summarise(
          AD_Intestinal_Disease_1 = sum(Any_Intestinal_Disease == 1 & AD == 1, na.rm = TRUE),
          Total_N_Intestinal_Disease_1 = sum(Any_Intestinal_Disease == 1, na.rm = TRUE),
          AD_Intestinal_Disease_0 = sum(Any_Intestinal_Disease == 0 & AD == 1, na.rm = TRUE),
          Total_N_Intestinal_Disease_0 = sum(Any_Intestinal_Disease == 0, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Fit logistic regression model and extract OR
      model_results <- tryCatch({
        model <- glm(model_formula, data = subset_data, family = binomial)
        tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
          filter(term == "Any_Intestinal_Disease") %>%
          select(estimate, conf.low, conf.high, p.value) %>%
          rename(OR = estimate, CI_Lower = conf.low, CI_Upper = conf.high, P_Value = p.value)
      }, error = function(e) {
        data.frame(OR = NA, CI_Lower = NA, CI_Upper = NA, P_Value = NA)
      })
      
      # Combine counts and model result
      cbind(stats, model_results)
    }) %>%
    ungroup() %>%
    mutate(Moderator = moderator)
  
  # Store per-moderator results
  results[[moderator]] <- group_results
}

# Combine all results into one table
final_results <- bind_rows(results, .id = "Moderator")

# Output to CSV
write.csv(final_results, "Moderation_Analysis_Results.csv", row.names = FALSE)

cat("\n--- Results saved to 'Moderation_Analysis_Results.csv' ---\n")
