
#This script uses the np package to determine the conditional density bandiwdths by Max likelihood cross validation

rm(list = ls())
graphics.off()

library(MASS)
library(np)
library(mvtnorm)
library(R.matlab)

data = readMat("trainresults_CS1_PM_forR.mat")  

alldata_tr = data$alldata
nres = 1;
ncov = 1; #Max no. of covariates, for 2 to 8 vars.

#Define the index of covariates for each response variable:
covindex = matrix(data=NA, nrow=nres, ncol = ncov)
covindex[1,1] = 2
#for (i in 1:nres) {  
#  covindex[i,1] = i+8
#}
#covindex[,2] = c(3, 4, 1, 2, 7, 4, 5, 2)+8   #no longer having second x(t-1) covariate 
#covindex[2:8,2] = c(1, 1, 2, 3, 4,5,2)


#now cycle through each response variable and do the conditional bw estimation:
bw = matrix(data =NA, nrow = nres, ncol = ncov+1)
bwm = matrix(data = NA, nrow = nres, ncol = ncov)

for (i in 1:nres) {
  rm('x_train', 'y_train', 'bw_train')
  print(i)
  y_train = data.frame(alldata_tr[,i])  # residuals
  
#  if (i == 1){
#    ncovtemp = ncov - 1
#  } else {
    ncovtemp = ncov
#    }
  
  x_train = data.frame(alldata_tr[,covindex[i,1:ncovtemp]])  # covariates  
  bw_train = npcdensbw(xdat = x_train, ydat=y_train, bwtype= "fixed")
  bwm[i,1:ncovtemp] = bw_train$xbw
  bw[i,1] = bw_train$ybw
  bw[i,2:(ncovtemp+1)] = bw_train$xbw  
  
}

datatrain = data.matrix(alldata_tr,rownames.force = NA)
covindex = data.matrix(covindex,rownames.force = NA)

writeMat("bwresults_CS1_PM.mat", bw = bw, bwm = bwm, datatrain = datatrain, covindex = covindex)

