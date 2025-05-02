#===============================
# Main Causal Mediation Analysis
#===============================

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Load required packages
if (!require("mediation")) install.packages("mediation")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
library(mediation)
library(dplyr)
library(tidyr)

# Load dataset
data <- read.csv("filtered_patient_info_long_balanced.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# Create a composite binary indicator for any intestinal disease
data <- data %>%
  mutate(Intestinal_disease = ifelse(
    Abdominal_Cancer == 1 | Intestinal_Inflammation == 1 |
    Irritable_Bowel_Syndrome == 1 | Constipation == 1 | Malabsorption == 1, 1, 0
  ))

# Drop rows with missing key variables
data <- data %>% drop_na(Alcohol_frequency, Sleeplessness, Intestinal_disease, AD)

#--------------------------
# Mediation via Alcohol
#--------------------------
medModel_alcohol <- glm(
  Alcohol_frequency ~ Intestinal_disease + Sex + Ethnic + Smoking + Cancer + Diabetes + Weight + Income_Category + Sleep_duration,
  family = gaussian(), data = data
)
outModel_alcohol <- glm(
  AD ~ Intestinal_disease * Alcohol_frequency + Sex + Ethnic + Smoking + Cancer + Diabetes + Weight + Income_Category + Sleep_duration,
  family = binomial(link = "logit"), data = data
)
med_alcohol <- mediate(medModel_alcohol, outModel_alcohol, treat = "Intestinal_disease", mediator = "Alcohol_frequency", boot = TRUE, sims = 1000)

alcohol_results <- data.frame(
  Estimate = c(med_alcohol$d.avg, med_alcohol$z.avg, med_alcohol$n.avg),
  `95% CI Lower` = c(med_alcohol$d.avg.ci[1], med_alcohol$z.avg.ci[1], med_alcohol$n.avg.ci[1]),
  `95% CI Upper` = c(med_alcohol$d.avg.ci[2], med_alcohol$z.avg.ci[2], med_alcohol$n.avg.ci[2]),
  p.value = c(med_alcohol$d.avg.p, med_alcohol$z.avg.p, med_alcohol$n.avg.p),
  Component = c("ACME", "ADE", "Prop Mediated")
)

#--------------------------
# Mediation via Sleeplessness
#--------------------------
medModel_sleep <- glm(
  Sleeplessness ~ Intestinal_disease + Sex + Ethnic + Smoking + Cancer + Diabetes + Weight + Income_Category + Sleep_duration,
  family = gaussian(), data = data
)
outModel_sleep <- glm(
  AD ~ Intestinal_disease * Sleeplessness + Sex + Ethnic + Smoking + Cancer + Diabetes + Weight + Income_Category + Sleep_duration,
  family = binomial(link = "logit"), data = data
)
med_sleep <- mediate(medModel_sleep, outModel_sleep, treat = "Intestinal_disease", mediator = "Sleeplessness", boot = TRUE, sims = 1000)

sleeplessness_results <- data.frame(
  Estimate = c(med_sleep$d.avg, med_sleep$z.avg, med_sleep$n.avg),
  `95% CI Lower` = c(med_sleep$d.avg.ci[1], med_sleep$z.avg.ci[1], med_sleep$n.avg.ci[1]),
  `95% CI Upper` = c(med_sleep$d.avg.ci[2], med_sleep$z.avg.ci[2], med_sleep$n.avg.ci[2]),
  p.value = c(med_sleep$d.avg.p, med_sleep$z.avg.p, med_sleep$n.avg.p),
  Component = c("ACME", "ADE", "Prop Mediated")
)

# Combine and export results
combined <- bind_rows(
  Alcohol = alcohol_results,
  Sleeplessness = sleeplessness_results,
  .id = "Mediator"
)

write.table(combined, file = "Mediation_Analysis_Results.txt", sep = "\t", row.names = FALSE)



#===========================
# Bootstrapping
#===========================

library(mediation)
library(dplyr)
library(tidyr)

data <- read.csv("filtered_patient_info_long_balanced.csv", fileEncoding = "UTF-8")

# Define composite disease status
data <- data %>%
  mutate(Intestinal_disease = ifelse(
    Abdominal_Cancer == 1 | Intestinal_Inflammation == 1 |
    Irritable_Bowel_Syndrome == 1 | Constipation == 1 | Malabsorption == 1, 1, 0
  )) %>%
  drop_na(Alcohol_frequency, Sleeplessness, Intestinal_disease, AD)

# Convert categorical variables to factors
data <- data %>%
  mutate(across(c(Education, Diabetes, Hypertension, Heart_Disease, Cancer, Cerebrovascular_Disease), as.factor))

# Function for bootstrap mediation
run_bootstrap <- function(mediator_var, output_file) {
  ACME <- ADE <- Prop <- data.frame(Estimate = rep(NA, 100), CI_lower = NA, CI_upper = NA)

  for (i in 1:100) {
    set <- data[sample(nrow(data), size = floor(0.5 * nrow(data))), ]
    medModel <- glm(reformulate(c("Intestinal_disease", "Sex", "Age", "Ethnic", "Education", "Smoking", "Cancer", "Diabetes", "Weight", "Income_Category"), response = mediator_var),
                    family = gaussian(), data = set)
    outModel <- glm(AD ~ Intestinal_disease + get(mediator_var) + Sex + Age + Ethnic + Education + Smoking + Cancer + Diabetes + Weight + Income_Category,
                    family = binomial(link = "probit"), data = set)
    med <- mediate(medModel, outModel, treat = "Intestinal_disease", mediator = mediator_var, boot = TRUE, sims = 100)
    s <- summary(med)
    ACME[i, ] <- c(s$d.avg, s$d.avg.ci)
    ADE[i, ] <- c(s$z.avg, s$z.avg.ci)
    Prop[i, ] <- c(s$n.avg, s$n.avg.ci)
  }

  write.csv(cbind.data.frame(ACME, ADE, Prop), file = output_file, row.names = FALSE)
}

# Run bootstrap analyses
run_bootstrap("Alcohol_frequency", "Intestinal_to_AD_via_Alcohol.csv")
run_bootstrap("Sleeplessness", "Intestinal_to_AD_via_Sleeplessness.csv")


#===========================
# Summary of Bootstrap Results
#===========================

summarize_mediation <- function(file_path, label) {
  data <- read.csv(file_path)
  colnames(data) <- c("ACME_Est", "ACME_Low", "ACME_High",
                      "ADE_Est", "ADE_Low", "ADE_High",
                      "Prop_Est", "Prop_Low", "Prop_High")

  median_prop <- median(data$Prop_Est)
  ci_prop <- quantile(data$Prop_Est, c(0.025, 0.975))
  p_val <- 2 * min(mean(data$Prop_Est > 0), mean(data$Prop_Est < 0))

  cat(paste0("\n--- ", label, " ---\n"))
  cat("Median Proportion Mediated: ", median_prop, "\n")
  cat("95% CI: ", ci_prop[1], " - ", ci_prop[2], "\n")
  cat("P-value: ", p_val, "\n")
}

summarize_mediation("Intestinal_to_AD_via_Alcohol.csv", "Alcohol")
summarize_mediation("Intestinal_to_AD_via_Sleeplessness.csv", "Sleeplessness")