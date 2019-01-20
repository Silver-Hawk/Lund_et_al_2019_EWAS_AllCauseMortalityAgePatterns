#
# Clear working space and set working directory
#
rm(list=ls())
setwd("YOUR WORKING DIRECTORY")

#
# Load libraries (most libraries can be installed using install.packages)
#
library(survival)
library(foreach)
library(doMC)

#
# Read files
#
load(file = "Combat_adjusted_cohorts.Rda")
cpg_matrix <- combat_dataframe

clinical_information = read.csv("CLINICAL_TRAITS") #


#
# Register multicores for  
#
registerDoMC(6)


#
# Perform CpG-wise cox regressions
#
finalRes <- foreach (i = colnames(cpg_matrix), .export = c("i"), .combine = rbind) %dopar% 
{
  cox <- coxph(Surv(time = delta_time, event = status) ~ 
                      cpg_matrix[,i] + sex + age_in_yrs + 
                      CD8T + CD4T + NK + Bcell + Mono + Gran, 
                      data = clinical_information, ties = c('efron'))
  
  sum <- summary(cox)
  beta = coef(cox)
  se <- sqrt(diag(cox$var))
  p <- 1 - pchisq((beta/se)^2, 1)
  CI <- round(confint(cox), 3)
  
  new_row = data.frame(beta=beta[1], se=se[1], z=beta[1]/se[1], p=p[1], CI=t(CI[1,]), row.names = c(i))
  
  new_row 
}

#
# Save results from analysis to a .csv file
#
write.csv(x=finalRes, file="cox_results_after_combat.csv")
