# New program using new data input structure
rm(list=ls()) # clear workingspace if necessary.
library(nnet) # load library

### Function that return PRB (probability) from frequency matrix
PRB <- function(FRQ){
  dd    <- dim(FRQ);
  prob  <- mat.or.vec(dd[1],dd[2]);
  tsum  <- sum(FRQ);
  for (i in 1:dd[1]) {
    for (j in 1:dd[2]) {
      prob[i,j] <- FRQ[i,j]/tsum;
    }
  }
  return(prob);
}

### Function that returns LOD (logg odds-ratio) from probability matrix
LOD <- function(JP){
  dd    <- dim(JP);
  logo  <- mat.or.vec(dd[1],dd[2]);
  csum  <- colSums(JP);
  rsum  <- rowSums(JP);
  fact  <- 3/log(2);
  for (i in 1:dd[1]) {
    for (j in 1:dd[2]) {
      logo[i,j] <- fact*log(JP[i,j]/(rsum[i]*csum[j]));
    }
  }
  # add 'U' ('C'+'J') residue type to the bottom of log odds ratio matrices
  aas <- strsplit("ACDEFGHIKLMNPQRSTVWYJ", "")[[1]];
  ci  <- which(aas == "C");
  ji  <- which(aas == "J");
  uJP <- JP[ci,]+JP[ji,];
  uLD <- c();
  for (j in 1:dd[2]) {
    uLD[j] <- fact*log(uJP[j]/(sum(uJP)*csum[j]));
  }
  ulogo <- rbind(logo, as.vector(uLD));
  return(ulogo);
}

### Data read
freq    <- read.csv("ulla-freq-toccata-maskA.csv",header=T) # usual case
szdata  <- dim(freq);
attach(freq);

b       <- freq[, -(1:5)]; # usual case
b       <- as.matrix(b) + 1; # pseudocount
SSE     <- as.factor(SSE);
SA      <- as.factor(SA);
HBOSH   <- as.factor(HBOSH);
HBMCO   <- as.factor(HBMCO);
HBMNH   <- as.factor(HBMNH);

options(contrasts=c("contr.treatment","contr.poly"))
PS <- list();
for(i in 1:21) {
  st      <- 21*i-20;
  ed      <- 21*i;
  ik      <- seq(st,ed);
  resy    <- b[,ik]; # add 1 for the case where no observations for entire row
  #out <- multinom(resy~SSE+SA);  
  #out <- multinom(resy~SSE*SA,Hess=F); # when considering reciprocal action
  out     <- multinom(resy~SSE+SA+HBOSH+HBMCO+HBMNH);
  out.sum <- summary(out);
  PS[[i]] <- out.sum$fitted.values;
}

### Joint probablity distn.
## This is only conditioned on X = A, hence to have the joint distribution of
## X & Y given the factors, we need to do the adjustments
## P(X, Y) = P(Y | X) P(X)
## to get the marginal P(X), just calculate empirical prob.

### First, convert PS into matrix with 21*21 matrix
PM  <- list();
tmp <- mat.or.vec(21,21);

for(j in 1:szdata[1]) {
  for(i in 1:length(PS)) {
    tmp[i,] <- PS[[i]][j,];
  }
  PM[[j]] <- tmp;
}

EP      <- list(); # empirical probability estimator
W       <- list();
JP      <- list(); # joint probability
logo    <- list(); # log odds ratio from JP
logoE   <- list(); # log odds ratio from EP

for(j in 1:szdata[1]) {
  bb          <- b[j,];
  bbb         <- matrix(bb,21,21,byrow=TRUE);
  wt          <- rowSums(bbb)/sum(bbb);
  W[[j]]      <- wt;
  JP[[j]]     <- PM[[j]]*wt;
  ep          <- (bbb)/(sum(bbb));
  EP[[j]]     <- ep;
  logo[[j]]   <- LOD(JP[[j]]);
  logoE[[j]]  <- LOD(EP[[j]]);
}

### calculate BLOSUM style log odds ratio matrices
#for(j in 1:length(logo)) {
  #logo[[j]]   <- LOD(JP[[j]]); # from multiple logistic regression
  #logoE[[j]]  <- LOD(EP[[j]]); # from emprical distribution
#}

### calculate a log odds ratio matrix for the total freq matrix
tot_freq <- colSums(freq[,6:dim(freq)[2]]);
tot_freq <- matrix(tot_freq,nrow=21);
tot_prob <- PRB(tot_freq);
tot_logo <- LOD(tot_prob);

### write results in files
## type 2
aas1      <- strsplit("ACDEFGHIKLMNPQRSTVWYJ", "")[[1]];
aas2      <- strsplit("ACDEFGHIKLMNPQRSTVWYJU", "")[[1]];
coln      <- c();
for (aa1 in aas1) {
  for (aa2 in aas2) {
    coln <- append(coln, paste(aa1,aa2,sep=""));
  }
}
new_nrow      <- dim(freq)[1]+1; # add 'total' table at the tail
new_ncol      <- dim(freq)[2]+21-4; # concatenate environmental class labels into one
                                    # and added 'U' increse the number of combination
csv           <- mat.or.vec(new_nrow,new_ncol);
colnames(csv) <- append("ENV",as.vector(coln));
for (k in 1:length(logo)) {
  csv[k, 1]           <- gsub("\\s", "", paste(as.matrix(freq[k, 1:5]), collapse="",sep=""));
  csv[k, 2:new_ncol]  <- as.vector(round(as.matrix(logo[[k]])));
}
csv[new_nrow,1]           <- "total";
csv[new_nrow,2:new_ncol]  <- as.vector(round(tot_logo));

file_name <- "ulla-mlr-logo-toccata-maskA.csv";
write.csv(csv, file_name);
