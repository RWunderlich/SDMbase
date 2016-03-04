# setwd("Your MaxEnt output directory")
# Following the Tutorial on MaxEnt by Steven Phillips

# Loading required packages
require(ROCR)
require(vcd)
require(boot)

# Reqading in and organizing the data
presence <- read.csv("YOUR_SPECIES_samplePredictions.csv")
background <- read.csv("YOUR_SPECIES_backgroundPredictions.csv")
pp <- presence$Logistic.prediction                   # get the column of predictions
testpp <- pp[presence$Test.or.train=="test"]         # select only test points
trainpp <- pp[presence$Test.or.train=="train"]       # select only test points
bb <- background$logistic

# Preparing input for ROCR
combined <- c(testpp, bb)                            # combine into a single vector
label <- c(rep(1,length(testpp)),rep(0,length(bb)))  # labels: 1=present, 0=random
pred <- prediction(combined, label)                  # labeled predictions, if there is an error make sure test % is not 0.
perf <- performance(pred, "tpr", "fpr")              # True / false positives, for ROC curve
plot(perf, colorize=FALSE)                            # Show the ROC curve
performance(pred, "auc")@y.values[[1]]               # Calculate the AUC
 
# Bootstrap to generate standard error and confidence interval for the AUC
# Function AUC
AUC <- function(p,ind) {
  pres <- p[ind]
  combined <- c(pres, bb)
  label <- c(rep(1,length(pres)),rep(0,length(bb)))
  predic <- prediction(combined, label)
  return(performance(predic, "auc")@y.values[[1]])
}

# Running the bootstrap with 100 AUC calculations
b1 <- boot(testpp, AUC, 100)
b1                                                   # gives estimates of standard error and bias of AUC
boot.ci(b1, conf = 0.95, type = "basic")             # confidence interval of AUC
