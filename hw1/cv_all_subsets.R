cv_all_subsets <- function(k,D,numberOfCV=10){
  #
  # Sets up the predictor models to be used with the R command lm and all subsets cross validation.
  #
  # Input:
  #
  #   k = the number of predictors to include in the search 
  #   D = training data frames where the last column is the response variable :
  #   numberOfCV = number of cross validations to do
  #
  # Output:
  #   cvraw = ( y - y_hat )^2 i.e. the residual squared.
  #           Blocks of y_hat are predicted based on all the other data.
  #           Doing this repeatedly is the cross-validation part.
  #   cvm = mean of the residual squares vector.
  #   cvsd = stderr of the mean residual square error.
  
  p = dim( D )[2] - 1 # subtract the response which is the last column in this data frame
  stopifnot( (0 <= k) && (k <= p) ) # check the input argument k 
  responseName = names(D)[p+1] # the the name of the response 
  
  if( k==0 ){ 
    # do the k=0 (no features) subset as a special case: 
    #
    form                = paste( responseName, " ~ +1" ) # this is the only possible formula in this case 
    res                 = cvLMTrainNTest( form, D, numberOfCV )
    bestErrors          = res$cvraw
    bestFeaturesIndices = NULL 
    bestFormula         = form
    
  }else{
    # do all remaining k>0 subsets: 
    #
    allPosSubsetsSizeK = combn(p,k) # get all possible subsets of size k 
    numOfSubsets = dim(allPosSubsetsSizeK)[2]
    
    for( si in 1:numOfSubsets ){
      print(sprintf("si=%10d; subsetSize=%10d; numOfSubsets=%10d; percentDone=%10.6f",si,k,numOfSubsets,si/numOfSubsets))
      featIndices = allPosSubsetsSizeK[,si]
      featNames   = as.vector(names(D))[featIndices]
      
      # construct a formula needed for the linear regression:
      # 
      form = paste( responseName, " ~ " )
      for ( ki in 1:k ){
        if( ki==1 ){ 
          form = paste( form, featNames[ki], sep=" " )
        }else{
          form = paste( form, featNames[ki], sep="+" )
        }
      }
      
      res  = cvLMTrainNTest( form, D, numberOfCV )
      espe = res$cvraw 
      
      if(si==1){ # initial subset of size k becomes current best 
        bestErrors          = espe
        bestFeaturesIndices = featIndices
        bestFormula         = form
      }else{     # all subsequent subsets of size k must beat the current best 
        if( mean(espe) < mean(bestErrors) ){
          bestErrors          = espe
          bestFeaturesIndices = featIndices
          bestFormula         = form
        }
      }
    } # endfor subsets of size k loop we have an estimate of the best formula predicted on the training data
  } # endif else k==0 special case 
  
  # package a structure to send out:
  #
  res = list(bestErrors,bestFeaturesIndices,bestFormula) 
  
  return(res) 
}


cvLMTrainNTest <- function( form, D, numberOfCV ){
  #
  # Does Cross Validation of the Linear Model specified by the formula "form"
  #
  # Input:
  #   form = a string representation of the formula to use in the R call to "lm" 
  #   D = training data frame where the last column is the response variable
  #   numberOfCV = number of cross validations to do
  #
  # Output:
  #   MSPE = vector of length numberOfCV with components mean square prediction errors
  #          extracted from each of the numberOfCV test data sets
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
  # Checked: 2010-01-25 (no errors found)
  #
  #-----
  
  nSamples = dim( D )[1]
  p = dim( D )[2] - 1 # subtract the response which is the last column in this data frame 
  responseName = names(D)[p+1] # the the name of the response assumed to be the last column in the dataframe
  
  # needed for the cross vaildation (CV) loops:
  # 
  nCVTest    = round( nSamples*(1/numberOfCV) ) # each CV run will have this many test points
  nCVTrain   = nSamples-nCVTest # each CV run will have this many training points 
  
  for (cvi in 1:numberOfCV) {
    
    testInds = (cvi-1)*nCVTest + 1:nCVTest # this may attempt to access sample indexes > nSamples 
    testInds = intersect( testInds, 1:nSamples ) # so we restrict this list here
    DCVTest  = D[testInds,1:(p+1)] # get the predictor + response testing data
    
    # select this cross validation's section of data from all the training data:
    # 
    trainInds = setdiff( 1:nSamples, testInds ) 
    DCV       = D[trainInds,1:(p+1)] # get the predictor + response training data
    
    response            = DCV[,p+1]
    DCV[[responseName]] = NULL
    
    # Standardize the predictors and demean the response
    #
    # Note that one can get the centering value (the mean) with the command attr(DCV,'scaled:center')
    #                   and the scaling value (the sample standard deviation) with the command attr(DCV,'scaled:scale')
    # 
    DCV              = scale( DCV )                                                                  
    responseMean     = mean( response )                                                              
    response         = response - responseMean                                                       
    DCVb             = cbind( DCV, response ) # append back on the response                          
    DCVf             = data.frame( DCVb ) # a data frame containing all scaled variables of interest 
    names(DCVf)[p+1] = responseName # fix the name of the response
    
    # extract the centering and scaling information:
    #
    means = attr(DCV,"scaled:center")
    stds  = attr(DCV,"scaled:scale")
    
    # apply the computed scaling based on the training data to the testing data:
    #
    responseTest            = DCVTest[,p+1] # in physical units (not mean adjusted)
    DCVTest[[responseName]] = NULL 
    DCVTest                 = t( apply( DCVTest, 1, '-', means ) )                                             
    DCVTest                 = t( apply( DCVTest, 1, '/', stds ) )                                              
    DCVTestb                = cbind( DCVTest, responseTest - responseMean ) # append back on the response      
    DCVTestf                = data.frame( DCVTestb ) # a data frame containing all scaled variables of interest
    names(DCVTestf)[p+1]    = responseName # fix the name of the response
    
    # fit this linear model and compute the expected prediction error (EPE) using these features (this just estimates the mean in this case):
    # 
    mk = lm( formula = form, data=DCVf )
    
    pdt  = predict( mk, newdata=DCVTestf, interval="prediction" )[,1] # get predictions on the test set
    pdt  = pdt + responseMean # add back in the mean.  Now in physical units (not mean adjusted)
    
    if( cvi==1 ){
      predmat = pdt
      y       = responseTest
    }else{
      predmat = c( predmat, pdt )
      y       = c( y, responseTest )
    }
    
  } # endfor cvi loop
  
  # this code is modeled after the code in "glmnet/cvelnet.R"
  #
  N = length(y)
  cvraw = ( y - predmat )^2
  cvm   = mean(cvraw)
  cvsd  = sqrt( var(cvraw)/N )
  l = list( cvraw=cvraw, cvm=cvm, cvsd=cvsd, name="Mean Squared Error" )
  
  return(l) 
}

