########################################
# Descriptive Statistics for Cross-Sectional and Follow-Up Data
# Author: Zhichao Liang
# Description: This script loads cross-sectional and follow-up datasets,
# computes summary statistics for continuous and categorical variables,
# and writes the results to separate CSV files for publication.
########################################

# Set CRAN mirror
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Load required package
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

# --------------------- Define summary functions ---------------------

# Function to summarize continuous variables
summarize_continuous <- function(df, vars) {
  df %>%
    summarise(across(all_of(vars), list(
      median = ~median(., na.rm = TRUE),
      Q1 = ~quantile(., 0.25, na.rm = TRUE),
      Q3 = ~quantile(., 0.75, na.rm = TRUE)
    ), .names = "{.col}_{.fn}"))
}

# Function to summarize categorical variables
summarize_categorical <- function(df, exclude_vars) {
  cat_vars <- setdiff(names(df), exclude_vars)
  result_list <- lapply(cat_vars, function(var) {
    counts <- table(df[[var]])
    percentages <- round(100 * prop.table(counts), 1)
    result <- data.frame(
      Variable = var,
      Category = names(counts),
      Count = as.vector(counts),
      Percentage = as.vector(percentages)
    )
    return(result)
  })
  do.call(rbind, result_list)
}

# --------------------- Cross-Sectional Data ---------------------

cross_data <- read.csv("cross_section_data.csv", stringsAsFactors = FALSE)
cross_continuous_vars <- c("Age", "BMI", "Income", "Sleep_duration")
cross_exclude <- c("start_time", "end_time", cross_continuous_vars)

# Compute summaries
cross_summary_cont <- summarize_continuous(cross_data, cross_continuous_vars)
cross_summary_cat <- summarize_categorical(cross_data, cross_exclude)

# Save to CSV
write.csv(cross_summary_cont, "summary_continuous_cross_section.csv", row.names = FALSE)
write.csv(cross_summary_cat, "summary_categorical_cross_section.csv", row.names = FALSE)

# --------------------- Follow-Up Data ---------------------

followup_data <- read.csv("follow_up_data.csv", stringsAsFactors = FALSE)
followup_continuous_vars <- c("Age", "BMI", "Income", "Sleep_duration",
                              "physical_activity", "diet_vriation", "summer_outdoor", "winter_outdoor", "medication")
followup_exclude <- c("start_time", "end_time", followup_continuous_vars)

# Compute summaries
followup_summary_cont <- summarize_continuous(followup_data, followup_continuous_vars)
followup_summary_cat <- summarize_categorical(followup_data, followup_exclude)

# Save to CSV
write.csv(followup_summary_cont, "summary_continuous_followup.csv", row.names = FALSE)
write.csv(followup_summary_cat, "summary_categorical_followup.csv", row.names = FALSE)

# Done
cat("Summary statistics for both datasets exported successfully.\n")
