##--------------------------------------------
## rfSubjectSpecific.R
## This script runs RF k-fold cross validation on the subject level
## Subject-specific implementation written by Pawel Gajer, updated 7/13/16, and adapted and updated by Steven Smith, 4/4/17
##--------------------------------------------

library(randomForest)
library(parallel)
library(rfPermute)

rfSubjectSpecific <- function(X, y, subjID, nfolds=10, verbose=FALSE,nrep=100, permute=TRUE,... )
  
  ## Arguments:
  ## X      - predictors; data frame or matrix
  ## y      - response; vector of the samle length as nrow(X)
  ## subjID - vector of length nrow(X) with subject assignment to each sample (used in )
  ## nfolds - number of folds
  ## nrep   - number of permutations (ignored if permute=FALSE)
  ## permute- whether rfPermute() should be used instead of randomForest()
  ## verbose -print progress/output
  ## ... - parameters to pass to randomForest() or rfPermute()

## Values:
## - error:  list of differences between prediction and true values for
##   regression and pred==true logical vectors for classification (one for each
##   fold)
## - rmse: root mean squared error
## - nmse: normalized MSE = MSE/MSE(mean(y) as predictor) for regression and
##   MSE/MSE(the highest frequency class as predictor ) for classification
## - mae: mean absolute error
## - nmae: normalized MAE
## - cl.err: classification error = accuracy = sum(pred!=y)/length(y)

{
  if(!is.numeric(y)){ ## Classification has 2 more importance metrics that regression does, plus one for each class
    nMetrics<-2+length(unique(y))
  }else{ ## regression only has 2 RF metrics: IncNodPurity & IncMSE
    nMetrics<-2
  }
  n <- length(y)
  
  if ( !is.numeric(y) && !is.factor(y) )
    stop("y should be either numeric or factor")
  
  if ( !is.data.frame(X) && !is.matrix(X) )
    stop("X is neither data frame nor matrix")
  
  if ( n!= nrow(X) )
    stop("n!= nrow(X)")
  
  if ( n!= length(subjID) )
    stop("n!= length(subjID)")
  
  if ( !is.numeric(nfolds) )
    stop("nfolds is not numeric")
  
  if ( nfolds < 1 )
    stop("nfolds < 1")
  
  if ( nfolds > n )
    stop("nfolds > n")
  
  if ( length(y[is.na(y)]) > 0 )
  {
    warning("y has some NAs; removing them and the corresponding rows of X and subjID")
    idx <- !is.na(y)
    y <- y[idx]
    X <- X[idx,]
    subjID <- subjID[idx]
  }
  
  subjID <- as.character(subjID)
  uqSubjIDs <- unique(subjID)
  nSubj <- length(uqSubjIDs)
  
  if ( nfolds > nSubj )
  {
    warning(paste("nfolds needs to be not greater than number of subjects. Changing it to number of subjects: ",nSubj, sep=""))
    nfolds <- nSubj;
  }
  
  sp <- split(1:n, list(factor(subjID))) # list with each entry being a
  # vector of indices of y that
  # correspond to the same subjID
  # (used as the label of the list
  # element)
  y.mean <- 0
  y.mostFr <- 0
  if ( is.numeric(y) )
  {
    y.mean <- mean(y)
  }
  
  ##     # The splitting of data is done on the subject level. In particular, if
  ##     # nfolds is equal to the number of subjects, we get a jack-knife
  ##     # leave-one-out CV.
  ##     # Note that y does not have to be constant over subjects.
  
  s0 <- split(sample(nSubj),rep(1:nfolds,length=nSubj)) # nfolds split of all subjects
  ## Turning each element of s0 from vector of subject indices to vector of
  ## sample indices corresponding to the given subjects
  s <- list()
  for ( i in seq(s0) )
  {
    v <- uqSubjIDs[s0[[i]]]
    s[[i]] <- as.vector(unlist(sp[v]))
  }
  
  x.null <- rep(y.mostFr, n)
  

  error <- list()      # list of prediction errors: prediction - y
  error.null <- list() # list of null model prediction errors: mean(y) - y;
  imp_permute<-list() ## container for importance metrics for rfPermute
  mdl<-list()
  # mean(y) is predicted value for each coordinate - null
  # model; I am returning error.null so we can test if the
  # current model is significantly better than the null
  # model - that is if the mean(abs(errors)) is
  # significantly different from the
  # mean(abs(null.errors))
  
  sampleIdx <- c()     # vector of sample indices from each run of cross validation, so we can match errors with samples
  r2.loc <- numeric(nfolds)
  r2.pearson.loc <- numeric(nfolds)
  r2.spearman.loc <- numeric(nfolds)
  rmse.loc <- numeric(nfolds)
  nmse.loc <- numeric(nfolds)
  mae.loc <- numeric(nfolds)
  nmae.loc <- numeric(nfolds)
  cl.err.loc <- numeric(nfolds)
  ncl.err.loc <- numeric(nfolds) # normalized classification error
  gError <- numeric(n) # "global" error array whose i-th entry is 1 is in the 10 fold CV the prediction of y[i] was correct
  imp <- matrix(0, nrow=ncol(X), ncol=nMetrics)
  importance_w_pval <- matrix(0, nrow=ncol(X), ncol=2*nMetrics) ##stores P vals and metrics, used in permutation
  for (i in seq(nfolds))
  {
    #i<-1
    if ( verbose )
      print(paste(" i=",i, sep=""))
    sampleIdx <- c(sampleIdx, s[[i]])
    trIdx <- setdiff(1:n, s[[i]])
            if(permute){
              print(paste0("Running CV fold ",i," out of ",nfolds ," using rfPermute"))
              m.rf <- rfPermute( X[trIdx,], y[trIdx], importance=TRUE,nrep = nrep, ... ) ## returns rfPErmute object, which is radnomForest object with additional results
              importance.i<-rp.importance(m.rf) ## includes p values calculated from rermutation model
            }else{
              print(paste0("Running CV fold ",i," out of ",nfolds ," using randomForest (no NULL distirbutions will be generated)"))
              m.rf <- randomForest( X[trIdx,], y[trIdx], importance=TRUE, ... )
              importance.i<-NULL
            }
            model.i<-m.rf ## saves the rf model object for iteration i
    
    ##m.rf <- randomForest( X, y, subset=setdiff(1:n, s[[i]]), importance=TRUE, ... )
    ##m.rf <- randomForest( X, y, subset=setdiff(1:n, s[[i]]), importance=TRUE)
    x <- predict(m.rf, newdata=X[s[[i]],], type="response")
    
    if ( is.numeric(y) ) ## If regression
    {
      error[[i]] <- x - y[s[[i]]]
      error.null[[i]] <- y.mean - y[s[[i]]]
      rmse.loc[i] <- sqrt(mean( error[[i]]^2 ))
      r2.loc[i] <- 100 * ( 1 - sum( error[[i]]^2 ) / sum( (x - y.mean)^2 ) ) # percentage of variance explained
      r2.pearson.loc[i] <- 100*cor(x, y[s[[i]]])^2
      r2.spearman.loc[i] <- 100*cor(x, y[s[[i]]], method="spearman")^2
      nmse.loc[i] <- mean( error[[i]]^2 ) / mean( (x - y.mean)^2 )
      mae.loc[i] <- mean( abs( error[[i]] ) )
      nmae.loc[i] <- mean( abs( error[[i]] ) ) / mean( abs(x - y.mean) )
      imp <- imp + randomForest::importance(m.rf)

    } else { ## If classification
      m <- length(s[[i]])
      error[[i]] <- as.character(x) != as.character(y[s[[i]]])
      error.null[[i]] <- as.character(x.null[1:m]) != as.character(y[s[[i]]])
      cl.err.loc[i] <- sum(error[[i]]) / m
      ncl.err.loc[i] <- cl.err.loc[i] / ( sum(error.null[[i]]) / m ) # NOTE that this will be NaN when the denominator is 0 (null model has no errors for the given y[s[[i]]]
      ##print(cbind(as.character(x),as.character(y[s[[i]]])))
      imp <- imp + randomForest::importance(m.rf)#[,3:4] ##why were only the last 2 being used?
      #head(importance(m.rf))
    }
    gError[s[[i]]] <- as.integer(error[[i]])
    imp_permute[[i]]<-importance.i
    mdl[[i]]<-model.i
    if(permute){
      importance_w_pval<-importance_w_pval+rp.importance(m.rf)
      head(importance_w_pval)
      }

}

    list(    #imp_permute=imp_permute,
             importance_w_pval=importance_w_pval/nfolds,
                  mdl=mdl,
             error=error,
       error.null=error.null,
       sampleIdx=sampleIdx,
       imp=imp/nfolds,
       rmse=mean(rmse.loc), nmse=mean(nmse.loc),
       mae=mean(mae.loc), nmae=mean(nmae.loc),
       r2=mean(r2.loc),
       r2.pearson=mean(r2.pearson.loc),
       r2.spearman=mean(r2.spearman.loc),
       gError=gError,
       cl.err=cl.err.loc,
       mean.cl.err=mean(cl.err.loc),
       ncl.err=ncl.err.loc,
       mean.ncl.err=mean(ncl.err.loc)
  )
}


## normalized accuracy and classification error for a two-class classification
## problem

## predVals <- pls.sPTB.v2.predict
## trueVals <- sptb2.char.f

normClErr <- function(predVals, trueVals)
  ## predVals - predicted values
  ## trueVals - true values
{
  if ( length(trueVals) != length(predVals) )
  {
    stop("ERROR: length(trueVals) != length(predVals)")
  }
  
  if ( !is.factor(trueVals) )
  {
    trueVals <- factor(trueVals)
  }
  
  predVals <- factor(predVals, levels=levels(trueVals))
  
  nElts <- length(trueVals)
  
  ## confusion matrix
  cm <- table(predVals, trueVals)
  
  ## accuracy
  acc <- (cm[1,1] + cm[2,2])/ nElts
  
  ## classification error
  clErr <- (cm[1,2] + cm[2,1])/ nElts
  
  
  ## accuracy of the naive classifier
  tt <- table(trueVals)
  i.mostFr <- which.max(tt)[[1]]
  true.mostFr <- names(tt)[i.mostFr]
  
  pred.null <- rep(true.mostFr, nElts)
  
  ## confusion matrix for the null model
  cm.null <- table(pred.null, trueVals)
  
  ## accuracy of the null model
  acc.null <- cm.null[1, i.mostFr] / nElts
  
  ## classification error of the null model
  if ( i.mostFr == 1 )
  {
    clErr.null <- cm.null[1,2]/ nElts
  } else {
    clErr.null <- cm.null[1,1]/ nElts
  }
  
  ## normalized accuracy
  norm.acc <- acc / acc.null
  
  ## normalized classification error
  norm.clErr <- clErr / clErr.null
  
  print(paste("Accuracy:", acc))
  print(paste("Accuracy of the naive classifier:", acc.null))
  print(paste("Relative Accuracy:", norm.acc, " ## Should be way greater than 1"))
  print(paste("Classification error:", clErr))
  print(paste("Classification error of the naive classifier:", clErr.null))
  print(paste("Relative Classification Error:", norm.clErr, " ## Should be very close to 0"))
  
  invisible(list(cm=cm, cm.null=cm.null, acc=acc, acc.null=acc.null, norm.acc=norm.acc, clErr=clErr, clErr.null=clErr.null, norm.clErr=norm.clErr))
}


