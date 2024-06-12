# Install if not installed
devtools::install_github("Mthrun/dbt.DataIO")
install.packages("AdaptGauss",dependencies = T)

# Load the libraries
library(dbt.DataIO)
library(AdaptGauss)
library(Metrics)
# Read the dataset
data_path <- "Flowcytometry.lrn"
data <- ReadLRN(data_path)

# Check the structure of the loaded data
print(head(data))
print(summary(data))

# Fitting a Finite Mixture Model (FMM) for each feature

# Check the data
print(head(data))

features <- colnames(data$Data)
###### Feature SS ###########
cat("Fitting model for feature:", features[1], "\n")
data_feature = data$Data[, "SS"]
hist(data_feature)
gmm_ss <-  AdaptGauss(data_feature)
# Otained (adapted in shiny) best Params with RMSD of 0.17:
#means_best = c(2.588828, 3.75499, 5.508684 )
#sd_best = c(0.336464, 0.393366, 0.314198)
#weight_best = c(0.2719, 0.3053, 0.3943)
gmm_ss <-  AdaptGauss(data_feature,
                      Means = c(2.588828, 3.75499, 5.508684 ),
                      SD=c(0.336464, 0.393366, 0.314198),
                      Weights = c(0.2719, 0.3053, 0.3943))
rms_ss <- gmm_ss$RMS
print(rms_ss)
qq <- QQplotGMM(data_feature, gmm_ss$Means, gmm_ss$SDs, gmm_ss$Weights)
# We can observe two deviations from the original distribution which
# result in higher RMSD...
Chi2testMixtures(data_feature, gmm_ss$Means, gmm_ss$SDs, gmm_ss$Weights)


###### Feature FS ###########
cat("Fitting model for feature:", features[2], "\n")
data_feature = data$Data[, "FS"]
hist(data_feature)
gmm_ss <-  AdaptGauss(data_feature)
# Obtained (adapted in shiny) best Params with RMSD of 0.16 (params already set):
#means_best = c(2.13697, 2.9648, 4.81905, 5.75355)
#sd_best = c(0.46847, 1.390585, 0.30243, 0.257955)
#weight_best = c(0.2026, 0.3441, 0.2215, 0.3649)
gmm_fs <-  AdaptGauss(data_feature,
                      Means = c(2.13697, 2.9648, 4.81905, 5.75355),
                      SD = c(0.46847, 1.390585, 0.30243, 0.257955),
                      Weights = c(0.2026, 0.3441, 0.2215, 0.3649))
rms_fs <- gmm_fs$RMS
print(rms_fc)
qq <- QQplotGMM(data_feature, gmm_fs$Means, gmm_fs$SDs, gmm_fs$Weights)
# We can observe slight but continuous deviations from the original distribution.
Chi2testMixtures(data_feature, gmm_fs$Means, gmm_fs$SDs, gmm_fs$Weights)

###### Feature HLA_DR ###########
cat("Fitting model for feature:", features[3], "\n")
data_feature = data$Data[, "HLA_DR"]
hist(data_feature)
gmm_hla_dr <-  AdaptGauss(data_feature)
# Obtained (adapted in shiny) best Params with RMSD of 0.18 (params already set):
#means_best = c(2.66528, 3.33126, 3.76365)
#sd_best = c(0.266656, 0.139738, 0.739907)
#weight_best = c(0.4771, 0.3169, 0.2148)
gmm_hla_dr <-  AdaptGauss(data_feature,
                          Means = c(2.66528, 3.33126, 3.76365),
                          SD = c(0.266656, 0.139738, 0.739907),
                          Weights = c(0.4771, 0.3169, 0.2148))
gmm_hla_dr_rms <- gmm_hla_dr$RMS
print(gmm_hla_dr_rms)
qq <- QQplotGMM(data_feature, gmm_hla_dr$Means, gmm_hla_dr$SDs, gmm_hla_dr$Weights)
# We can observe slight but continuous deviations from the original distribution.
Chi2testMixtures(data_feature, gmm_hla_dr$Means, gmm_hla_dr$SDs, gmm_hla_dr$Weights)

###### Feature CD7 ###########
cat("Fitting model for feature:", features[3], "\n")
data_feature = data$Data[, "CD7"]
hist(data_feature)
gmm_cd7 <-  AdaptGauss(data_feature)
# Obtained (adapted in shiny) best Params with RMSD of 0.18 (params already set):
#means_best = c(2.911157, 3.33126, 3.76365)
#sd_best = c(0.266656, 0.139738, 0.739907)
#weight_best = c(0.4771, 0.3169, 0.2148)
gmm_cd7 <-  AdaptGauss(data_feature,
                       Means = c(2.91962, 3.418674, 4.686419),
                       SD = c(0.482066, 0.124988, 0.312796),
                       Weights = c(0.7945, 0.0769, 0.1276))
gmm_cd7_rms <- gmm_cd7$RMS
print(gmm_cd7_rms)
qq <- QQplotGMM(data_feature, gmm_cd7$Means, gmm_cd7$SDs, gmm_cd7$Weights)
# We can observe slight but continuous deviations from the original distribution.
Chi2testMixtures(data_feature, gmm_cd7$Means, gmm_cd7$SDs, gmm_cd7$Weights)


# Output RMSD for each feature model
for (feature in features) {
  cat("RMSD for feature", feature, ":", models[[feature]]$RMSD, "\n")
}
