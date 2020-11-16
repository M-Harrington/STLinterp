#' @title Interpolation and Cross Validation for Seasonal Trend Loess
#' @details The function performs cross validation using either k-fold CV or Monte Carlo CV to estimate the best hyperparameter set for missing value interpolation. This function relies on the "Seasonal Trend Loess" method suggested in (Cleveland et al, 1990).
#'
#' @param x a time series object of type \code{ts}
#' @param s.window a value or vector of either the character string "periodic" or the span (in lags) of the loess window for seasonal extraction, which should be odd. This has no default.
#' @param s.degree a value or vector of the degree of locally-fitted polynomial in seasonal extraction. Should be 0, 1, or 2.
#' @param t.window a value or vector of the span (in lags) of the loess window for trend extraction, which should be odd. If \code{NULL}, the default, \code{nextodd(ceiling((1.5*period) / (1-(1.5/s.window))))}, is taken.
#' @param t.degree a value or vector of the degree of locally-fitted polynomial in trend extraction. Should be 0, 1, or 2.
#' @param fc.window a value or vector of the vector lengths of windows for loess smoothings for other trend frequency components after the original STL decomposition has been obtained. See [stlplus::stlplus()] for more info.
#' @param fc.degree a value or vector of the degrees of locally-fitted polynomial in the loess smoothings for the frequency components specified in fc.window. Values of 0, 1 and 2 are allowed.
#' @param k is either the number of splits in k-fold, or number of runs in Monte Carlo.
#' @param p only used if \code{type = "mc"}, defines the proportion of data points withheld.
#' @param type either \code{"kfold"} or \code{"mc"} for k-fold and Monte Carlo cross validation respectively.
#' @param value return type of the hyperparameter sets, either \code{"best"} or \code{"grid"}.
#' @param ... any additional parameters passed to \code{stlplus}.
#'
#' @return Either a \code{data.frame} or a single row depending on \code{value}
#' @export
#'
#' @examples
#' test.ts <- ts(c(3,4,3,4,5,NA,6,3,3,3,4,5,NA,5,15,3,4,5,NA,1,7,4,8), frequency=4)
#'
#' # K-fold
#' best_param<-STLinterp(test.ts, s.window=c(5,7), t.window = c(NA,5),
#'                      t.degree = 2, s.degree=c(1,2), fc.window = c(NA, 5), fc.degree = NA,
#'                      k =14, type="kfold", value ="best")
#' best_param
#'
#' # Monte Carlo
#' grid<-STLinterp(test.ts, s.window=c(5,7), t.window = c(NA,5),
#'                t.degree = c(1,2), s.degree=c(1,2), fc.window = c(NA, 5), fc.degree = NA,
#'                k =10, type="mc", p=0.95, value ="grid")
#' grid
#'
#'
#' @author Matt R Harrington
#' @note Occasionally you might get errors from the underlying algorithm. Frequently this is an issue of \code{s.window} or \code{t.window} being missing for many values in a row. You can try increasing \code{k} or \code{p} (depending on the CV method), or else increasing \code{s.window} or \code{t.window}
#' @seealso \code{stlplus::stlplus()} for more parameter options and details
#' @importFrom stlplus tidyverse
STLinterp <- function(x, s.window,  s.degree = 1, t.window = NA,
                      t.degree = 1, fc.window = NA, fc.degree = NA, k=NA,
                      p=0.95, type ="kfold", value="best", ...){
  if (is.na(k)) stop("K must be user-defined, either k-folds or k runs of Monte Carlo")

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
          message(paste("Grid not evaluated: ",paste(grid[i,1:(ncol(grid)-1)],
                                                     " Try increasing k in kfold or p in montecarlo or increasing window parameters", collapse=" ")))
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
          message(paste("Grid not evaluated: ",paste(grid[i,1:(ncol(grid)-1)],
                                                     " Try increasing k in kfold or p in montecarlo or increasing window parameters", collapse=" ")))
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

#######
# MC and Kfold Subroutines for STLinterp
# These cross_val methods take time series (x), number of sets/folds to average over (k),
# proportion of datapoints to include in each set (p), single set of parameters(grid_row) for
# stlplus (as a vector). Returns the CV score for that hyperparameter's set sof runs.

#' @export
mc_cross_val<- function(x, grid_row, k, p){

  #take grid_row, translate appropriate elements from NA to NULL
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
                            fc.window = grid_rowl[[5]], fc.degree = grid_rowl[[6]], ...)
    #reconstruct, compare to actual
    reconstructed <- seasonal(stlobjs[[i]])+trend(stlobjs[[i]])
    score[i] <-mean((x[sort(CVindicies[[i]])] - reconstructed[sort(CVindicies[[i]])] )^2)
  }
  return(mean(score))
}



#' @export
kfold_cross_val<- function(x, grid_row, k){

  #take grid_row, translate appropriate elements from NA to NULL
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
                      fc.window = grid_rowl[[5]], fc.degree = grid_rowl[[6]], ...)
    #reconstruct, compare to actual
    reconstructed <- seasonal(stlobj)+trend(stlobj)
    score[i] <-mean((x[sort(CVindicies[[i]])] - reconstructed[sort(CVindicies[[i]])] )^2)
  }
  return(mean(score))
}
