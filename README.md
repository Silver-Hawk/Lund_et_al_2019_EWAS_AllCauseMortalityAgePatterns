# DNA methylome profiling of all-cause mortality in comparison with age-associated methylation patterns (Lund et al. 2019)
### Jesper Beltoft Lund, Shuxia Li, Jan Baumbach, Anne Marie Svane, Jacob Hjelmborg, Lene Christiansen, Kaare Christensen, Paul Redmond, Riccardo E. Marioni, Ian J. Deary, Qihua Tan.
###### (Github/Lund_et_al_2019_EWAS_AllCauseMortalityAgePatterns)
Set of scripts and instructions for the scientific article "DNA methylome profiling of all-cause mortality in comparison with age-associated methylation patterns" (Lund et al. 2019). 
*The paper is not yet published*.
The study uses the Lothian Birth Cohort (LBC) data which are accessible through the European Genome-phenome Archive (https://www.ebi.ac.uk/ega/home) with accession number EGAS00001000910, upon request.
### Requirements
In order to use these scripts, the following data, structure of data and software is assumed to be available:
+ Methylation data from the (450K) Illumina BeadChip array (CpGs as columns, participants as rows).
+ Annotation data with: age, sex + cell type compositions (assuming the data is based on whole blood samples), values as log2 based M-values. (Traits as columns, participants as rows)
+ The Illumina (450K) annotation file (this can be found at https://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit/downloads.html)
+ All scripts are written in the R statistical programming language (https://www.r-project.org/).
### Preprocessing (*Preprocessing.R*)
You should remove probes with European allele frequency above 1%, cross-reactive probes, and CpGs with more than 5% missing values and detection P values > 0.05 from your methylation dataset. See Chen et al. (2013) ”Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray” and Aryee et al. (2014) ”Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays.” for more info.
```R
#
# Clear working space and set working directory
#
rm(list=ls())
setwd("YOUR WORKING DIRECTORY")

#
# Load libraries (most libraries can be installed using install.packages)
#
library(sva)
library(zoo)
library(ggplot2)
library(gridExtra)

#
# Read files
#
pca <- read.csv("PCA_RESULTS.csv") # PCA was performed and analyized beforehand

# Convert to data.frame
pca <- data.frame(pc1 = pca$PC1, pc2 = pca$PC2, is_batched = pca$is_batched)

clinical_information = read.csv("CLINICAL_TRAITS") # 
load(file = "METHYLATION_MATRIX.RData") # Methylation matrix -> cpg_matrix
clinical_information$pca <- pca$is_batched

#
# Perform PCA
#

modcombat = model.matrix(~1, data=clinical_information)

combat_edata = ComBat(dat=t(cpg_matrix), 
                      batch=as.factor(clinical_information$pca), 
                      mod=modcombat,
                      ref.batch = 1)

#
# Save PCA results
#
write.csv(x=combat_edata, file="ComBat_Results_Matrix.csv")

combat_dataframe <- as.data.frame(t(combat_edata))

pr <- prcomp(combat_dataframe, rank. = 2)
pr_x <- pr$x

#
# Save new PCA results
#

write.csv(x=pr_x, file="PCA_after_combat.csv")

plot_data <- cbind(pr_x, clinical_information)

plot1 <- ggplot(aes(x=PC1, y=PC2, color=as.factor(status)), data=plot_data) + 
  geom_point(size=1) + labs(colour="Dead/Censor")  + theme_bw() +
  xlim(-400,400) + ylim(-400, 400) + coord_fixed(ratio=1)

plot2 <- ggplot(aes(x=PC1, y=PC2, color=as.factor(lbc)), data=plot_data) + 
  geom_point(size=1) + labs(colour="LBC Cohort") + theme_bw() + ggtitle("Status") +
  xlim(-400,400) + ylim(-400,400) + coord_fixed(ratio=1)

plot3 <- ggplot(aes(x=PC1, y=PC2, color=age_in_yrs), data=plot_data) + 
  geom_point(size=1)  + labs(colour="Age in years") + theme_bw() +
  xlim(-400, 400) + ylim(-400, 400) + coord_fixed(ratio=1)

plot4 <- ggplot(aes(x=PC1, y=PC2, color=as.factor(sex)), data=plot_data) + 
  geom_point(size=1)  + labs(colour="Gender") + theme_bw() +
  xlim(-400, 400) + ylim(-400,400) + coord_fixed(ratio=1)

plot5 <- ggplot(aes(x=PC1, y=PC2, color=as.factor(Sentrix.Barcode)), data=plot_data) + 
  geom_point(size=1)  + labs(colour="Barcode") + theme_bw() + theme(legend.position = "none") + ggtitle("Chip Barcodes") +
  xlim(-400,400) + ylim(-400,400) + coord_fixed(ratio=1)
  
plot6 <- ggplot(aes(x=PC1, y=PC2, colour=as.factor(pca)), data=plot_data) + 
  geom_point(size=1) + labs(colour="Old Batch effect")  + theme_bw() +
  xlim(-400,400) + ylim(-400,400) + coord_fixed(ratio=1)

#
# Draw plot
#

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 2, nrow = 3)


new_clinical_information <- cbind(clinical_information, pr_x)
merged <- cbind(new_clinical_information, combat_dataframe)

#
# Save adjusted data for further use
#

save(merged, file="Combat_adjusted_cohorts.Rda")
```
### Survival analysis (*Survival_analysis.R*)
Cox proportional hazard regressions for single CpG survival analysis.  
```R
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
```

### Figures
#### Figure 1 
Created using *CMplot* R-package (https://github.com/YinLiLin/R-CMplot) and *ggplot2* R-package (https://ggplot2.tidyverse.org/). 
See [Figure 1](Figure1.pdf).

#### Figure 2
Created using *fsmb* R-package (https://cran.r-project.org/web/packages/fmsb/index.html).
See [Figure 2](Figure2.pdf).

#### Figure 3
Created using *ggplot2* R-package (https://ggplot2.tidyverse.org/).
See [Figure 3](Figure3.pdf).




