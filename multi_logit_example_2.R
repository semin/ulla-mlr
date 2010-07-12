
# New program using new data input structure

setwd("d:/blosum/logit"); # working directory
rm(list = ls()) # clear workingspace if necessary.
library(nnet) # load library

# Data read

data <- read.csv("ulla-env2-freq-toccata.csv")
szdata = dim(data);
attach(data);

b = data[, -(1:2)];
b = as.matrix(b)+1; # add 1 pseudo count
SA = as.factor(SA);
SSE = as.factor(SSE);

# data <- read.csv("ulla-env5-freq-toccata.csv")
# szdata = dim(data);
# attach(data);
# b = data[, -(1:5)];
# b = as.matrix(b)+1; # add 1 pseudo count
# SA = as.factor(SA);
# SSE = as.factor(SSE);
# HBOSH = as.factor(HBOSH);
# HBMCO = as.factor(HBMCO);
# HBMNH = as.factor(HBMNH);

options(contrasts=c("contr.treatment", "contr.poly"))
out <- multinom2(b ~ SA * SSE, Hess = FALSE)
# out <- multinom2(b ~ SA + SSE + HBOSH + HBMCO + HBMNH);

JP=list();
EP = list();
logo =  list();
logoE = list();

for(j in 1:szdata[1]){
	bb = b[j,];
	bbb = matrix(bb, 21,21, byrow=TRUE);
	ep =(bbb)/(sum(bbb));
        EP[[j]] = ep; 

	tmp = out$fitted.values[j,];
	JP[[j]] = matrix(tmp, 21,21, byrow=TRUE);

	logo[[j]] = BLOSUM(JP[[j]]);
	logoE[[j]] = BLOSUM(EP[[j]]);
}

