# Requiring external packages
require(raster)
require(sp)
require(localgauss)
# Calculate linear monotonic, linear non-monotonic and Gaussian non-linear correlations between predictors
# Reading in the continuous rasters in alphabetical order beware of bioclim variable naming omitting leading 0s
# Continuous layers's directory is located on the Desktop
fnames <- list.files(path = '~/Desktop/Layers/', pattern='asc', full.names=TRUE)
predictors <- stack(fnames[1:length(fnames)])
#### 1. pearson correlation ####
tempdat1<-cor(as.data.frame(predictors), use="pairwise", method = "pearson")
sqrt(mean(tempdat1^2))
write.table(as.data.frame(tempdat1),file="~/Desktop/pears.csv",row.names=TRUE,col.names=TRUE,sep=",")
#### 2. spearman correlation ####
tempdat2<-cor(as.data.frame(predictors), use="pairwise", method = "spearman")
sqrt(mean(tempdat2^2))
write.table(as.data.frame(tempdat2),file="~/Desktop/spear.csv",row.names=TRUE,col.names=TRUE,sep=",")
### localgauss ###
# identify and remove NA from all layers
predictorsnona <- complete.cases(predictors[1:length(predictors)])
predictorsnona <- predictors[predictorsnona]
# Computing Gauss correlation.. smart way to do it for all pairwise combinations?
x=predictorsnona[, 1]
y = predictorsnona[, 2]
lg.out <- localgauss(x = x, y = y, b1 = 1, b2 = 1)
# hthresh should be set to include most points but omit more remote points
# if mean(lg.out$eflag) =! 0 then NA
# else rhos<-(lg.out$par.est[,5])
