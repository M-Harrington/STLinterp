require("stlplus")
require("tidyverse")

#######
# Subroutine for STLinterp to perform monte carlo cross validation
# mc_cross_val method takes time series (x), number of sets to average over (k),
# proportion of datapoints to include in each set (p), single set parameters(grid_row) for
# stlplus (as a vector). Returns the CV score from Monte Carlo method
# Note p should be high and k large for higher accuracy

mc_cross_val<- function(x, grid_row, k, p){
  
  #take grid_row, translate appropriate elemnts from NA to NULL
  grid_rowl <- as.list(grid_row)
  for (i in 1:length(grid_row)){ 
    if (is.na(grid_rowl[[i]])){grid_rowl[i]<- list(NULL)}
  }
  
  #init
  CVindicies <- vector(mode = "list",length=k)
  stlobjs <- vector(mode = "list",length=k)
  clone_holdout <- vector(mode = "list",length=k)
  score <- vector(mode = "numeric",length=k)
  
  #create indicies, prevent NA's from being selected
  a<- which(is.na(x)) %>% intersect(1:length(x))
  sampleindicies <- 1:length(x) %>% setdiff(a)
  
  for (i in 1:k){
    fold_indicies<- sample(sampleindicies,
                           size = round(length(sampleindicies)*p),
                           replace = F) %>% sort
    CVindicies[[i]]<- setdiff(sampleindicies, fold_indicies)}
  
  #Set indicies to missing, calculate stlplus with given params
  for (i in 1:k){
    clone_holdout[[i]]<- x
    clone_holdout[[i]][sort(CVindicies[[i]])] <- NA
    
    #create STL object
    stlobjs[[i]] <- stlplus(clone_holdout[[i]], s.window = grid_rowl[[1]], s.degree = grid_rowl[[2]],
                            t.window = grid_rowl[[3]], t.degree = grid_rowl[[4]],
                            fc.window = grid_rowl[[5]], fc.degree = grid_rowl[[6]])
    #reconstruct, compare to actual
    reconstructed <- seasonal(stlobjs[[i]])+trend(stlobjs[[i]])
    score[i] <-sum((x[sort(CVindicies[[i]])] - reconstructed[sort(CVindicies[[i]])] )^2)
  }
  return(mean(score))
}



#######
# Subroutine for STLinterp to perform k-fold cross validation
# kfold_cross_val method takes time series (x), number of folds(k), and a
# single set parameters(grid_row) for stlplus (as a vector). Returns the 
# CV score from k-fold method
kfold_cross_val<- function(x, grid_row, k){
  
  #take grid_row, translate appropriate elemnts from NA to NULL
  grid_rowl <- as.list(grid_row)
  for (i in 1:length(grid_row)){ 
    if (is.na(grid_rowl[[i]])){grid_rowl[i]<- list(NULL)}
  }
  
  #init
  CVindicies <- vector(mode = "list",length=k)    
  score <- vector(mode = "numeric",length=k)      
  
  #create indicies, prevent NA's from being selected
  a<- which(is.na(x)) %>% intersect(1:length(x))
  sampleindicies <- 1:length(x) %>% setdiff(a)
  random_ind <- sample(sampleindicies)
  
  
  if (k>length(random_ind)){stop("k is greater than number of non-missing values")}
  
  #break up random_indicies into semi-equal groups
  breaks <- sort(random_ind%%k)
  
  #gen held out indicies
  for (i in 1:k){
    CVindicies[[i]]<- random_ind[breaks==(i-1)]}
  
  
  #Set indicies to missing, calculate stlplus with given params
  for (i in 1:k){
    clone_holdout<- x
    clone_holdout[sort(CVindicies[[i]])] <- NA
    
    if (sum(is.na(clone_holdout))<1){stop("Time series too small")}
    
    #create STL object
    stlobj <- stlplus(clone_holdout, s.window = grid_rowl[[1]], s.degree = grid_rowl[[2]],
                      t.window = grid_rowl[[3]], t.degree = grid_rowl[[4]],
                      fc.window = grid_rowl[[5]], fc.degree = grid_rowl[[6]])
    #reconstruct, compare to actual
    reconstructed <- seasonal(stlobj)+trend(stlobj)
    score[i] <-sum((x[sort(CVindicies[[i]])] - reconstructed[sort(CVindicies[[i]])] )^2)
  }
  return(mean(score))
}




####STLinterp
# Takes a time series ts() obj (x), and vectors of hyper parameter 
# values for cross validation grid search.
# Type either takes on "kfold" or "mc" for kfold or monte carlo CV respectively 
# The function estimates the best parameters and return the estimated missing 
# values using the STL decompisition.
# Note that the parametrs can be single values or vector of values for grid
# Missing values must be given as NAs.
# Value can be either "best" for returning the best set of hyperparameter grid
# Or "grid" and returns the entire grid
#
# For information on the specific hyperparameters see the stlplus 
# documentation for more info 

#NOTE: k should be a high proportion of length(ts), otherwise might get errors

STLinterp <- function(x, s.window,  s.degree = 1, t.window = NA,
                      t.degree = 1, fc.window = NA, fc.degree = NA, n=NA, k=NA,
                      p=0.95, type ="kfold", value="best"){
  if (!is.na(n)) {if (n >length(x)/2){stop("n is too large for given time series")}}
  
  #init
  df <- list(s.window=s.window, s.degree=s.degree, t.window=t.window, 
             t.degree=t.degree, fc.window=fc.window, fc.degree=fc.degree)
  grid <- expand.grid(df)
  grid <- grid %>% mutate(CVscore = NA) 
  
  # apply correct cv function to parameter list, combine result with original grid
  #Try catch for errors in underlying STL algorithm
  if (type == "mc"){
    for (i in 1:nrow(grid)){
      res_err <- tryCatch(
        expr = {
          grid$CVscore[i]<-mc_cross_val(x, grid[i,1:(ncol(grid)-1)], k=k,p=p)        },
        error = function(e){
          message(paste("Problem with: ",paste(grid[i,1:(ncol(grid)-1)], collapse=" ")))
          grid$CVscore[i]<-NA
          
        }
      )
    }
    
  } else if (type=="kfold"){
    for (i in 1:nrow(grid)){
      res_err <- tryCatch(  
        expr = {
          grid$CVscore[i]<-kfold_cross_val(x=x, grid[i,1:(ncol(grid)-1)], k=k)        },
        error = function(e){ 
          message(paste("Problem with: ", paste(grid[i,1:(ncol(grid)-1)],collapse=" ")))
          grid$CVscore[i]<-NA
          
        }
      )
      
      
    }
  } else {stop("Type not supported")}
  
  
  # Return either best or grid
  if (value == "best"){
    best_params <- grid[grid$CVscore == min(grid$CVscore ,na.rm=TRUE)& !is.na(grid$CVscore),]
      return(best_params)
  } else if (value == "grid"){
      return(grid)
  } else {stop("Value not supported")}
  
  
}


### Example implementation
test.ts <- ts(c(3,4,3,4,5,NA,6,3,3,3,4,5,NA,5,15,3,4,5,NA,1,7,4,8), frequency=4)

# K-fold
best_param<-STLinterp(test.ts, s.window=c(5,7), t.window = c(NA,5),
                      t.degree = c(1,2), s.degree=c(1,2), fc.window = c(NA, 5), fc.degree = NA,
                      k =14, type="kfold", value ="best")
best_param

# Monte Carlo
grid<-STLinterp(test.ts, s.window=c(5,7), t.window = c(NA,5),
                t.degree = c(1,2), s.degree=c(1,2), fc.window = c(NA, 5), fc.degree = NA,
                k =10, type="mc", p=0.95, value ="grid")
grid