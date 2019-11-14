one_standard_error_rule <- function(complexityParam,cvResults){
  #
  # Input:
  #   complexityParam = the sampled complexity parameter (duplicate cross vaidation results expected)
  #   cvResults = matrix where
  #               each column is a different value of the complexity parameter (with complexity increasing from left to right) ...
  #               each row    is the square prediction error (SPE) of a different cross validation sample i.e. ( y_i - \hat{y}_i )^2 
  #
  # Output:
  #   optCP = the optimal complexity parameter
  #
  # Written by:
  # -- 
  # John L. Weatherwax                2009-04-21
  # 
  # email: wax@alum.mit.edu
  # 
  # Please send comments and especially bug reports to the
  # above email address.
  #
  #-----
  
  # this code is modeled after the code in "glmnet/cvelnet.R"
  #
  N = dim(cvResults)[1]
  
  # compute the complexity based mean and standard deviations:
  means = apply(cvResults,2,mean)
  stds  = sqrt( apply(cvResults,2,var)/N ) 
  
  # find the smallest espe:
  minIndex = which.min( means )
  
  # compute the confidence interval around this point:
  #cip = 1 - 2*pnorm(-1)
  #ciw   <- qt(cip/2, n) * stdev / sqrt(n)
  ciw = stds[minIndex] 
  
  # add a width of one std to the min point:
  maxUncertInMin = means[minIndex] + 0.5*ciw
  
  # find the mean that is nearest to this value and that is a SIMPLER model than the minimum mean model: 
  complexityIndex = which.min( abs( means[1:minIndex] - maxUncertInMin ) ) 
  complexityValue = complexityParam[complexityIndex]
  
  # package everything to send out:
  res = list(complexityValue,complexityIndex,maxUncertInMin, means,stds)
  
  return(res)
}




