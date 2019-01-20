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

#
# Save adjusted data for further use
#

save(combat_dataframe, file="Combat_adjusted_cohorts.Rda")