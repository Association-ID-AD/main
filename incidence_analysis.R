#############################################
# Alzheimer's Disease (AD) Incidence Analysis
# Author: Zhichao Liang
# Description:
# This script analyzes AD onset proportion and cumulative
# incidence in relation to intestinal diseases using
# follow-up data. It includes both Kaplan-Meier survival
# analysis and proportion plots by age at AD onset.
#############################################

# ----------- Load Required Packages -----------
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("lubridate")) install.packages("lubridate")
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
if (!require("ggplot2")) install.packages("ggplot2")

library(dplyr)
library(readr)
library(lubridate)
library(survival)
library(survminer)
library(ggplot2)

# ----------- Load and Preprocess Data -----------
df <- read_csv("incidence_data.csv")

# Convert date columns to Date format
df <- df %>%
  mutate(
    start_time = as.Date(start_time, format = "%Y/%m/%d"),
    end_time = as.Date(end_time, format = "%Y/%m/%d")
  )

# Calculate survival time and event indicator
df <- df %>%
  mutate(
    start_year = year(start_time),
    end_year = year(end_time),
    Survival_Time = end_year - start_year,
    Event = ifelse(AD == 1, 1, 0)
  ) %>%
  filter(Survival_Time > 0)  # Remove invalid entries

# Create binary factor for intestinal disease
df <- df %>%
  mutate(intestinal_before = factor(Intestinal_disease, levels = c(0, 1), labels = c("0", "1")))

# ----------- Cox Proportional Hazards Model -----------
cox_model <- coxph(Surv(Survival_Time, Event) ~ intestinal_before + Sex + Age + BMI + Weight + Health + Education + 
                     Racitizen + Cancer + Diabetes + Hypertension + Heart_attack + Cerebrovascular, data = df)

# Print summary of model
summary(cox_model)

# ----------- Kaplan-Meier Cumulative Incidence Plot -----------
fit <- survfit(Surv(Survival_Time, Event) ~ intestinal_before, data = df, type = "kaplan-meier")

# Save KM plot
pdf("ad_cumulative_incidence_kmplot.pdf", width = 7, height = 5)
ggsurvplot(
  fit,
  data = df,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  censor = FALSE,
  risk.table = FALSE,
  fun = "event",  # Show cumulative incidence (1 - survival)
  pval = TRUE,
  xlab = "Follow-up Time (years)",
  ylab = "Cumulative Incidence of AD",
  legend.title = "Intestinal disease",
  legend.labs = c("0", "1"),
  surv.median.line = "hv",
  palette = c("#B7A7C5", "#E7B800"),
  xlim = c(0, 25),
  ylim = c(0, 0.15)
)
dev.off()

# ----------- Age Distribution of AD Onset by Group -----------
df_ad <- df %>%
  filter(AD == 1) %>%
  mutate(group = factor(Intestinal_disease, levels = c(0, 1), labels = c("0", "1")))

group_totals <- df %>%
  group_by(Intestinal_disease) %>%
  summarise(group_total = n()) %>%
  mutate(group = factor(Intestinal_disease, levels = c(0, 1), labels = c("0", "1")))

df_percent <- df_ad %>%
  group_by(Age, group) %>%
  summarise(AD_count = n(), .groups = "drop") %>%
  left_join(group_totals, by = "group") %>%
  mutate(Percentage = AD_count / group_total * 100)

# Save age-onset percentage plot
pdf("ad_onset_percentage.pdf", width = 7, height = 5)
ggplot(df_percent, aes(x = Age, y = Percentage, color = group, shape = group, fill = group)) +
  geom_point(size = 2.5, alpha = 0.8, stroke = 1.2) +
  geom_smooth(method = "lm", se = TRUE, size = 1.1, fullrange = TRUE, alpha = 0.3) +
  scale_shape_manual(values = c(1, 2)) +
  scale_color_manual(values = c("#B7A7C5", "#E7B800")) +
  scale_fill_manual(values = c("#B7A7C5", "#E7B800")) +
  labs(
    x = "Age at AD Onset (years)",
    y = "Proportion of AD Onset (%)",
    color = "Intestinal disease",
    shape = "Intestinal disease",
    fill = "Intestinal disease"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top",
    text = element_text(size = 13),
    axis.title = element_text(face = "plain"),
    legend.title = element_text(face = "plain"),
    legend.text = element_text(face = "plain")
  )
dev.off()

cat("\n--- AD incidence analysis completed. Output PDF files saved. ---\n")
