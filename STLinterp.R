require("stlplus")
require("tidyverse")

#######
# Subroutine for STLinterp to perform monte carlo cross validation
# mc_cross_val method takes time series (x), number of draws(n), single set parameters(grid_row) for
# stlplus (as a vector). Returns the CV score from Monte Carlo method
mc_cross_val<- function(x, grid_row, n){

  #take grid_row, translate appropriate elemnts from NA to NULL
  grid_rowl <- as.list(grid_row)
  for (i in 1:length(grid_row)){ 
    if (is.na(grid_rowl[[i]])){grid_rowl[i]<- list(NULL)}
    }

  #init
  CVindicies <- vector(mode = "list",length=n)
  stlobjs <- vector(mode = "list",length=n)
  clone_holdout <- vector(mode = "list",length=n)
  score <- vector(mode = "numeric",length=n)

  #create indicies, prevent NA's from being selected
  a<- which(is.na(x)) %>% intersect(1:length(x))
  sampleindicies <- 1:length(x) %>% setdiff(a)

  for (i in 1:n){
    CVindicies[[i]]<- sample(sampleindicies,
                             size = floor(length(sampleindicies)/17),
                             replace = T)}

  #Set indicies to missing, calculate stlplus with given params
  for (i in 1:n){
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
  CVindicies <- vector(mode = "list",length=n)
  stlobjs <- vector(mode = "list",length=n)
  clone_holdout <- vector(mode = "list",length=n)
  score <- vector(mode = "numeric",length=n)

  #create indicies, prevent NA's from being selected
  a<- which(is.na(x)) %>% intersect(1:length(x))
  sampleindicies <- 1:length(x) %>% setdiff(a)
  random_ind <- sample(sampleindicies)

  #break up random_indicies into semi-equal groups
  breaks <- sort(random_ind%%n)

  #gen held out indicies
  for (i in 1:n){
    CVindicies[[i]]<- random_ind[breaks==(i-1)]}


  #Set indicies to missing, calculate stlplus with given params
  for (i in 1:n){
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




####STLinterp
# Takes a time series ts() obj (x), and vectors of hyper parameter 
# values for cross validation grid search.
# Type either takes on "kfold" or "mc" for kfold or monte carlo CV respectively 
# The function estimates the best parameters and return the estimated missing 
# values using the STL decompisition.
# Note that the parametrs can be single values or vector of values for grid
# Missing values must be given as NAs.
#
# For information on the specific hyperparameters see the stlplus 
# documentation for more info 

STLinterp <- function(x, s.window,  s.degree = 1, t.window = NA,
                      t.degree = 1, fc.window = NA, fc.degree = NA, n=10, type ="kfold"){
  if (n >length(x)/2){stop("n is larger than given time series")}
  
  #init
  df <- list(s.window=s.window, s.degree=s.degree, t.window=t.window, 
             t.degree=t.degree, fc.window=fc.window, fc.degree=fc.degree)
  grid <- expand.grid(df)
  grid <- grid %>% mutate(CVscore = NA) 
  
  # apply correct cv function to parameter list, combine result with original grid
  if (type == "mc"){
    for (i in 1:nrow(grid)){
      grid$CVscore[i]<-mc_cross_val(x, grid[i,1:(ncol(grid)-1)], n=n) 
    }
  } else if (type=="kfold"){
    for (i in 1:nrow(grid)){
      grid$CVscore[i]<-kfold_cross_val(x=x, grid[i,1:(ncol(grid)-1)], k=n) 
    }
  } else {stop("Type not supported")}
  
  
  # Potential addition: return best param set
  best_params <- grid[grid$CVscore == min(grid$CVscore),]
  return(best_params)

  # return(grid)
}
  
  
  
# # Example implementation    
# test.ts <- ts(c(3,4,5,6,3,4,5,6,3,4,5,6,3,NA,5,6), frequency=4) 
# grid<-STLinterp(test.ts, s.window=c(2,4,5), t.window = c(NA,4),
#           t.degree = c(1,2), s.degree=c(1,2), fc.window = c(NA, 5), fc.degree = NA, n =3, type="kfold")
