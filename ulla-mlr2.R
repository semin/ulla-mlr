# New program using new data input structure
rm(list=ls()) # clear workingspace if necessary.
library(nnet) # load library

# Slight modification of multinom function
# Add option for nnet.default() with max iterations
# and max weights. Max weights is related to the number of 
# parameter estimates. There are four places to change this.
# For more factors, you can increase it arbitraray.

multinom2 <- function (formula, data, weights, subset, na.action, contrasts = NULL, 
    Hess = FALSE, summ = 0, censored = FALSE, model = FALSE, 
    ...) 
{
    class.ind <- function(cl) {
        n <- length(cl)
        x <- matrix(0, n, length(levels(cl)))
        x[(1L:n) + n * (as.vector(unclass(cl)) - 1L)] <- 1
        dimnames(x) <- list(names(cl), levels(cl))
        x
    }
    summ2 <- function(X, Y) {
        X <- as.matrix(X)
        Y <- as.matrix(Y)
        n <- nrow(X)
        p <- ncol(X)
        q <- ncol(Y)
        Z <- t(cbind(X, Y))
        storage.mode(Z) <- "double"
        z <- .C(VR_summ2, as.integer(n), as.integer(p), as.integer(q), 
            Z = Z, na = integer(1L))
        Za <- t(z$Z[, 1L:z$na, drop = FALSE])
        list(X = Za[, 1L:p, drop = FALSE], Y = Za[, p + 1L:q])
    }
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$summ <- m$Hess <- m$contrasts <- m$censored <- m$model <- m$... <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    X <- model.matrix(Terms, m, contrasts)
    cons <- attr(X, "contrasts")
    Xr <- qr(X)$rank
    Y <- model.response(m)
    if (!is.matrix(Y)) 
        Y <- as.factor(Y)
    w <- model.weights(m)
    if (length(w) == 0L) 
        if (is.matrix(Y)) 
            w <- rep(1, dim(Y)[1L])
        else w <- rep(1, length(Y))
    lev <- levels(Y)
    if (is.factor(Y)) {
        counts <- table(Y)
        if (any(counts == 0L)) {
            empty <- lev[counts == 0L]
            warning(sprintf(ngettext(length(empty), "group %s is empty", 
                "groups %s are empty"), paste(sQuote(empty), 
                collapse = " ")), domain = NA)
            Y <- factor(Y, levels = lev[counts > 0L])
            lev <- lev[counts > 0L]
        }
        if (length(lev) < 2L) 
            stop("need two or more classes to fit a multinom model")
        if (length(lev) == 2L) 
            Y <- as.vector(unclass(Y)) - 1
        else Y <- class.ind(Y)
    }
    if (summ == 1) {
        Z <- cbind(X, Y)
        z1 <- cumprod(apply(Z, 2L, max) + 1)
        Z1 <- apply(Z, 1L, function(x) sum(z1 * x))
        oZ <- order(Z1)
        Z2 <- !duplicated(Z1[oZ])
        oX <- (seq_along(Z1)[oZ])[Z2]
        X <- X[oX, , drop = FALSE]
        Y <- if (is.matrix(Y)) 
            Y[oX, , drop = FALSE]
        else Y[oX]
        w <- diff(c(0, cumsum(w))[c(Z2, TRUE)])
        print(dim(X))
    }
    if (summ == 2) {
        Z <- summ2(cbind(X, Y), w)
        X <- Z$X[, 1L:ncol(X)]
        Y <- Z$X[, ncol(X) + 1L:ncol(Y), drop = FALSE]
        w <- Z$Y
        print(dim(X))
    }
    if (summ == 3) {
        Z <- summ2(X, Y * w)
        X <- Z$X
        Y <- Z$Y[, 1L:ncol(Y), drop = FALSE]
        w <- rep(1, nrow(X))
        print(dim(X))
    }
    offset <- model.offset(m)
    r <- ncol(X)
    if (is.matrix(Y)) {
        p <- ncol(Y)
        sY <- Y %*% rep(1, p)
        if (any(sY == 0)) 
            stop("some case has no observations")
        if (!censored) {
            Y <- Y/matrix(sY, nrow(Y), p)
            w <- w * sY
        }
        if (length(offset) > 1L) {
            if (ncol(offset) != p) 
                stop("ncol(offset) is wrong")
            mask <- c(rep(FALSE, r + 1L + p), rep(c(FALSE, rep(TRUE, 
                r), rep(FALSE, p)), p - 1L))
            X <- cbind(X, offset)
            Wts <- as.vector(rbind(matrix(0, r + 1L, p), diag(p)))
            fit <- nnet.default(X, Y, w, Wts = Wts, mask = mask, 
                size = 0, skip = TRUE, softmax = TRUE, censored = censored, 
                rang = 0, maxit = 3000, MaxNWts =15000, ...)
        }
        else {
            mask <- c(rep(FALSE, r + 1L), rep(c(FALSE, rep(TRUE, 
                r)), p - 1L))
            fit <- nnet.default(X, Y, w, mask = mask, size = 0, 
                skip = TRUE, softmax = TRUE, censored = censored, 
                rang = 0, maxit = 3000, MaxNWts = 15000, ...)
        }
    }
    else {
        if (length(offset) <= 1L) {
            mask <- c(FALSE, rep(TRUE, r))
            fit <- nnet.default(X, Y, w, mask = mask, size = 0, 
                skip = TRUE, entropy = TRUE, rang = 0, maxit = 3000, MaxNWts = 15000, ...)
        }
        else {
            mask <- c(FALSE, rep(TRUE, r), FALSE)
            Wts <- c(rep(0, r + 1L), 1)
            X <- cbind(X, offset)
            fit <- nnet.default(X, Y, w, Wts = Wts, mask = mask, 
                size = 0, skip = TRUE, entropy = TRUE, rang = 0, maxit = 3000, MaxNWts = 15000, 
                ...)
        }
    }
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$terms <- Terms
    fit$call <- call
    fit$weights <- w
    fit$lev <- lev
    fit$deviance <- 2 * fit$value
    fit$rank <- Xr
    edf <- ifelse(length(lev) == 2L, 1, length(lev) - 1) * Xr
    if (is.matrix(Y)) {
        edf <- (ncol(Y) - 1) * Xr
        if (length(dn <- colnames(Y)) > 0) 
            fit$lab <- dn
        else fit$lab <- 1L:ncol(Y)
    }
    fit$coefnames <- colnames(X)
    fit$vcoefnames <- fit$coefnames[1L:r]
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    fit$edf <- edf
    fit$AIC <- fit$deviance + 2 * edf
    if (model) 
        fit$model <- m
    class(fit) <- c("multinom", "nnet")
    if (Hess) 
        fit$Hessian <- multinomHess(fit, X)
    fit
}
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

## try 1
## estimation for comparted data
#PS <- list();
#for(i in 1:21) {
  #st      <- 21*i-20;
  #ed      <- 21*i;
  #ik      <- seq(st,ed);
  #resy    <- b[,ik]; # add 1 for the case where no observations for entire row
  ##out <- multinom(resy~SSE+SA);  
  ##out <- multinom(resy~SSE*SA,Hess=F); # when considering reciprocal action
  #out     <- multinom(resy~SSE+SA+HBOSH+HBMCO+HBMNH);
  #out.sum <- summary(out);
  #PS[[i]] <- out.sum$fitted.values;
#}

### Joint probablity distn.
## This is only conditioned on X = A, hence to have the joint distribution of
## X & Y given the factors, we need to do the adjustments
## P(X, Y) = P(Y | X) P(X)
## to get the marginal P(X), just calculate empirical prob.

### First, convert PS into matrix with 21*21 matrix
#PM  <- list();
#tmp <- mat.or.vec(21,21);

#for(j in 1:szdata[1]) {
  #for(i in 1:length(PS)) {
    #tmp[i,] <- PS[[i]][j,];
  #}
  #PM[[j]] <- tmp;
#}

## Try 2
## Estimation for whole data 
out     <- multinom2(b ~ SA + SSE + HBOSH + HBMCO + HBMNH);

#W       <- list();
EP      <- list(); # empirical probability estimator
JP      <- list(); # joint probability
logo    <- list(); # log odds ratio from JP
logoE   <- list(); # log odds ratio from EP

for(j in 1:szdata[1]) {
  bb          <- b[j,];
  bbb         <- matrix(bb,21,21,byrow=TRUE);
  #wt          <- rowSums(bbb)/sum(bbb);
  #W[[j]]      <- wt;
  #JP[[j]]     <- PM[[j]]*wt;
  ep          <- (bbb)/(sum(bbb));
  EP[[j]]     <- ep;
  tmp         <- out$fitted.values[j,];
  JP[[j]]     <- matrix(tmp,21,21,byrow=TRUE);
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
