###############################################
# Mediation Analysis with Alcohol and Sleeplessness
# Author: Zhichao Liang
# Date: 2025-05-13
# Description:
# This script performs causal mediation analysis using
# the "mediation" package. It assesses whether alcohol
# consumption and sleeplessness mediate the relationship
# between intestinal diseases and Alzheimer’s disease (AD).
###############################################

# -------------------- Load Required Packages --------------------
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

if (!require("mediation")) install.packages("mediation")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")

library(dplyr)
library(tidyr)
library(mediation)

# -------------------- Load and Preprocess Data --------------------
data <- read.csv("follow_up_data.csv", stringsAsFactors = FALSE)

# Ensure necessary columns exist and drop missing values
data <- data %>%
  drop_na(Alcohol, Sleeplessness, Intestinal_disease, AD)

# -------------------- Mediation Analysis: Alcohol --------------------

# Step 1: Model for the mediator (Alcohol consumption)
medModel_alcohol <- glm(
  Alcohol ~ Intestinal_disease,
  family = gaussian(),
  data = data
)

# Step 2: Model for the outcome (AD) including the mediator
outModel_alcohol <- glm(
  AD ~ Intestinal_disease * Alcohol,
  family = binomial(link = "logit"),
  data = data
)

# Step 3: Mediation analysis
med_alcohol <- mediate(
  model.m = medModel_alcohol,
  model.y = outModel_alcohol,
  treat = "Intestinal_disease",
  mediator = "Alcohol",
  boot = TRUE,
  sims = 1000
)

# Extract results for alcohol mediation
alcohol_results <- data.frame(
  Estimate = c(med_alcohol$d.avg, med_alcohol$z.avg, med_alcohol$n.avg),
  `95% CI Lower` = c(med_alcohol$d.avg.ci[1], med_alcohol$z.avg.ci[1], med_alcohol$n.avg.ci[1]),
  `95% CI Upper` = c(med_alcohol$d.avg.ci[2], med_alcohol$z.avg.ci[2], med_alcohol$n.avg.ci[2]),
  p.value = c(med_alcohol$d.avg.p, med_alcohol$z.avg.p, med_alcohol$n.avg.p),
  Component = c("ACME (Average)", "ADE (Average)", "Proportion Mediated (Average)")
)

# -------------------- Mediation Analysis: Sleeplessness --------------------

# Step 1: Model for the mediator (Sleeplessness)
medModel_sleeplessness <- glm(
  Sleeplessness ~ Intestinal_disease,
  family = gaussian(),
  data = data
)

# Step 2: Model for the outcome (AD) including the mediator
outModel_sleeplessness <- glm(
  AD ~ Intestinal_disease * Sleeplessness,
  family = binomial(link = "logit"),
  data = data
)

# Step 3: Mediation analysis
med_sleeplessness <- mediate(
  model.m = medModel_sleeplessness,
  model.y = outModel_sleeplessness,
  treat = "Intestinal_disease",
  mediator = "Sleeplessness",
  boot = TRUE,
  sims = 1000
)

# Extract results for sleeplessness mediation
sleeplessness_results <- data.frame(
  Estimate = c(med_sleeplessness$d.avg, med_sleeplessness$z.avg, med_sleeplessness$n.avg),
  `95% CI Lower` = c(med_sleeplessness$d.avg.ci[1], med_sleeplessness$z.avg.ci[1], med_sleeplessness$n.avg.ci[1]),
  `95% CI Upper` = c(med_sleeplessness$d.avg.ci[2], med_sleeplessness$z.avg.ci[2], med_sleeplessness$n.avg.ci[2]),
  p.value = c(med_sleeplessness$d.avg.p, med_sleeplessness$z.avg.p, med_sleeplessness$n.avg.p),
  Component = c("ACME (Average)", "ADE (Average)", "Proportion Mediated (Average)")
)

# -------------------- Combine and Export Results --------------------
combined_results <- bind_rows(
  Alcohol = alcohol_results,
  Sleeplessness = sleeplessness_results,
  .id = "Mediator"
)

# Save results to file
write.table(combined_results, file = "mediation_analysis_results.txt", sep = "\t", row.names = FALSE)

# -------------------- Print Summary --------------------
cat("\n--- Mediation Model for Alcohol Frequency ---\n")
summary(medModel_alcohol)
cat("\n--- Outcome Model for Alcohol Frequency ---\n")
summary(outModel_alcohol)

cat("\n--- Mediation Model for Sleeplessness ---\n")
summary(medModel_sleeplessness)
cat("\n--- Outcome Model for Sleeplessness ---\n")
summary(outModel_sleeplessness)

cat("\n--- Mediation analysis results saved to 'mediation_analysis_results.txt' ---\n")
