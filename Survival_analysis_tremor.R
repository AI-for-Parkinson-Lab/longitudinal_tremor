# Load libraries
library(survival)
library(survminer)
library(broom) # for tidy output 
library(ggplot2) # for graphing 
library(splines)
library(dplyr)
library(tidyr)

# Create dataframe
start_medication = read.csv("C:\\Users\\z835211\\OneDrive - Radboudumc\\Documents\\Tremor progression paper\\Matlab_results\\Cox_data_table.csv")
start_medication$ID <- as.factor(start_medication$ID)
start_medication$t_start <- start_medication$Time
start_medication$Time <- NULL
start_medication$t_stop <- start_medication$t_start +2
start_medication$Event <- 1 - start_medication$Event

# Fit cox model
start_medication.cox <- coxph(Surv(t_start,t_stop, Event) ~ X_MA_Lagged1, data = start_medication)
summary(start_medication.cox) 

# Check proportional hazards assumption
cox.zph(start_medication.cox) 
plot(cox.zph(start_medication.cox))

# Survival probability curve based on mean covariate values
surv.at.means <- survfit(start_medication.cox)
tidy(surv.at.means)
plot(surv.at.means, xlab="weeks", ylab="survival probability")

surv.probabilities <- survfit(start_medication.cox, newdata= start_medication, id = ID)
surv.probabilities$time <- surv.probabilities$time+4 # correct for lag
surv_probabilities <- data.frame(c(surv.probabilities[['time']]),c(surv.probabilities[['surv']]))
colnames(surv_probabilities)<- c("time","surv")
surv_probabilities <- surv_probabilities %>%
  mutate(ID = cumsum(time == 6))

# Convert to wide format
df_wide <- surv_probabilities %>%
  pivot_wider(
    id_cols = ID,
    names_from = time,
    values_from = surv,
    names_prefix = "time_"
  )

write.csv(df_wide, "C:\\Users\\z835211\\OneDrive - Radboudumc\\Documents\\Tremor progression paper\\Matlab_results\\survival_probabilities.csv")
