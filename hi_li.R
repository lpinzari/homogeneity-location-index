
# Title       : hi_li.R
# Description : This file contains the implementation of the Homogeneity index (HI) and Location Index (LI)
# Author      : Ludovico Pinzari
# Usage       : Run all the functions included in this file and use the function uni.hom and uni.loc (see notes)
# R version   : 3.3 or later version
# Notes       : see below

################################################################################################
#                                                                                              #
#                               HOMOGENEITY INDEX  (HI)                                        #
#                                                                                              #
################################################################################################
# The HI function [uni.hom], consists of two main functions:                                   #
#                                                                                              #
#  - [uni.concI]    : computes the concentration index (CI) for a vector of observations       #
#                     (e.g. the number of people in each quantile of the socioeconomic index). #
#                                                                                              #
#  - [uni.div]      : computes the divergence index (DI) for a probability distribution vector #
#                  (e.g. the proportion of people in each quantile of the socioeconomic index),#
#                     by calling two support functions:                                        #
#                                                                                              #
#                     Divergence Index supporting functions                                    #
#                                                                                              #
#                       - [uni.conv] : computes the convolution of two vectors                 #
#                       - [uni.corr] : computes the autocorrelation of a vector                #
#                                                                                              #
#  - [uni.divConst] : computes the DI of a bimodal distribution with maximum variance          #
#                                                                                              #
################################################################################################

################################################################################################
#                                                                                              #
#                               LOCATION INDEX  (LI)                                           #
#                                                                                              #
################################################################################################
# The LI function [uni.loc],returns the position of the bins with the maximum concentration    #
# by calling a single support function:                                                        #
#                                                                                              #
#   - [uni.locVec] : returns the bin concentration vector for a probability distribution.      #
#                    The bin concentration vector provides the concentration value for each    #
#                    bin in the distribution.                                                  #                                                            #
#                                                                                              #
################################################################################################

#===============================================================================================#



#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                             CONCENTRATION INDEX (CI)                                        #
#---------------------------------------------------------------------------------------------#
#' @function uni.concI                                                                        #
#' @description function to compute the Concentration Index for a vector                      #
#' @param p numeric vector with non-negative integer values. (1 x n)                          #
#' @return the Concentration Index (CI) of \code{p}: real number  [0 1]                       #
#' @usage uni.concI(\code{p})                                                                 #
#' @examples                                                                                  #
#' uniform distribution                                                                       #
#' p <- c(10,10,10,10,10,10,10,10,10,10)                                                      #
#' r <- uni.concI(p)      r: 0                                                                #
#' p <- c(50,0,0,0,0,0,0,0,0,0,50)                                                            #   
#' r <- uni.concI(p)      r: 0.88                                                             #
#' p <- c(100,0,0,0,0,0,0,0,0,0)                                                              #
#' r <- uni.concI(p)      r: 1                                                                #
#' -------------------------------------------------------------------------------------------#

uni.concI <- function(p){
  
  d <- length(p)
  pop <- sum(p)
  i <- 1
  
  ## compute the pdf
  
  while(i <= d)
  {
    p[i] <- p[i]/pop
    i <- i+1
  }
  
  p<-sort(p) ## sort the pdf vector in ascending order
  
  ## Compute the cumulative frequencies Lorenz Curve
  
  lc <- rep(0,d+1)
  i <- 1
  
  while(i <= d)
  {
    if(i==1)
      lc[2]<-lc[1]+p[1]
    else
      lc[i+1]<-lc[i]+p[i]
    i<- i+1
  }
  
  ## Compute the Area vector Under the Lorenz Curve 
  
  b <- 1/d ## base of the trapezoid
  i <- 1
  imp_area <- 1-b ## to normalize the result
  la <- rep(0,d)
  
  while( i <= d)
  {
    la[i] <- b*(lc[i]+lc[i+1])/(2*imp_area)
    i <- i+1
  }
  
  ## Compute the Lorenz Curve of the unifotm distr.
  
  pu <- rep(1/d,d)
  lcu <- rep(0,d+1)
  i <- 1
  
  while(i <= d)
  {
    if(i==1)
      lcu[2] <- lcu[1]+pu[1]
    else
      lcu[i+1] <- lcu[i]+pu[i]
    i<- i+1
  }
  
  
  ## Compute the area under the Lorenz Curve of the unifoirm distribution
  
  lau <- rep(0,d)
  
  i<-1
  
  while( i <= d)
  {
    lau[i] <- b*(lcu[i]+lcu[i+1])/(2*imp_area)
    i <- i+1
  }
  
  ## Compute the Zonoid Area (Area between the lorenz Curve and Uniform distribution)
  
  Area <- 0
  i<- 1
  
  while( i <= d)
  {
    Area <- Area + lau[i]-la[i]
    i <- i+1
  }
  
  ## rounding the value
  Area <-round(Area,4) 
  
  return(2*Area)
}

#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                             DIVERGENCE INDEX SUPPORTING FUNCTIONS                           #
#---------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------#
#                                                                                            #
#                             UNIDIMENSIONAL CONVOLUTION                                     #
#--------------------------------------------------------------------------------------------#
#' @function uni.conv                                                                        #
#' @description function to compute the convolution of two vectors                           #
#' @param x numeric vector representing polynomial coefficient (1 x m)                       #
#' @param y numeric vector representing polynomial coefficient (1 x n)                       #
#' @return the coefficient vector resulting from multiplying the polynomial represented      #
#' by x by the polynomial represented by y (1 x m+n-1)                                       #
#' @usage uni_concv(\code{x},\code{y})                                                       #
#' @examples                                                                                 #
#' x <- c(0.25,0.25,0.25,0.25)                                                               #
#' y <- c(1,1,1,1)                                                                           #
#' r <- uni.conv(x,y)      r: num [1:7] 0.25 0.50 0.75 1.00 0.75 0.50 0.25                   #
#' x <- c(1,0,0,0)                                                                           #
#' y <- c(1,1,1,1)                                                                           #
#' r <- uni.conv(x,y)      r: num [1:7] 1 1 1 1 0 0 0                                        #
#' ------------------------------------------------------------------------------------------#

uni.conv <- function(x,y){
  
  m <- length(x)
  n <- length(y)
  z <- numeric(m+n-1)
  
  for(j in 1:m){
    for(k in 1:n){
      z[j+k-1] = z[j+k-1]+x[j]*y[k]
    }
  }
  return(z)
}

#--------------------------------------------------------------------------------------------#
#                                                                                            #
#                             UNIDIMENSIONAL AUTOCORRELATION                                 #
#--------------------------------------------------------------------------------------------#
#' @function uni.corr                                                                        #
#' @description \code{uni.corr} uses \code{uni.conv} to compute the autocorrelation of       #
#'               a vector                                                                    #
#' @param x numeric vector representing polynomial coefficient (1 x n)                       #
#' @return the coefficient vector resulting from the autocorrelation (1 x 2n-1)              #
#' @usage uni.corr(\code{x})                                                                 #
#' @seealso \code{uni.conv}                                                                  #
#' @examples                                                                                 #
#' x <- c(0.25,0.25,0.25,0.25)                                                               #
#' r <- uni.corr(x)      r: num [1:7] 0.0625 0.1250 0.1875 0.2500 0.1875 0.1250 0.0625       #
#' ------------------------------------------------------------------------------------------#

uni.corr <- function(x){
  
  R <- uni.conv(x,rev(x))
  return(R)
  
}

#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                               DIVERGENCE INDEX                                              #
#---------------------------------------------------------------------------------------------#
#' @function uni.div                                                                          #
#' @description \code{uni.div} uses \code{uni.conv} and \code{uni.corr} to compute the        #
#'              polarization divergence of a probability vector.                              #
#' @param x numeric vector representing polynomial coefficient distribution (pdf)             #
#' @return the Divergence index. (1 x 1) real number [0 1)                                    #
#' @usage uni.div(\code{x})                                                                   #
#' @seealso \code{uni.conv}, \code{uni.corr}                                                  #
#' @examples                                                                                  #
#' x <- c(0.5,0,0,0,0,0,0,0,0,0.5)                                                            #
#' r <- uni.div(x)      r: 0.2973122                                                          #
#' x <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1) Uniform distribution                       #
#' r <- uni.div(x)      r: 0.1229451                                                          #
#' x <- c(0,0,0,1,0,0,0,0,0,0)                                                                #
#' r <- uni.div(x)      r: 0                                                                  #
#' -------------------------------------------------------------------------------------------#

uni.div <- function(x){
  
  d <- length(x)
  comb <- rep(1,d)
  
  # compute the Bilateral Cumulative Distributive function and its autocorrelation (spectra)
  
  bcdf <- uni.conv(x,comb)
  sbcdf <- uni.corr(bcdf)
  
  # compute the Singleton autocorrelation function
  
  imp <- rep(0,d)
  imp[1] <- 1
  
  bcdfi <- uni.conv(imp,comb)
  sbcdfi <- uni.corr(bcdfi)
  
  # normalized energy (pdf of the signal)
  
  E <- sum(sbcdfi)
  
  sbcdf <- sbcdf/E
  sbcdfi <- sbcdfi/E
  
  # compute the binary logarithm of the signals
  
  l <- length(sbcdf)
  lgs1 <- rep(0,l)
  lgsi <- rep(0,l)
  
  i <- 1
  
  while (i <= l){
    
    if(sbcdf[i] > 0)
      lgs1[i] <- log(sbcdf[i],2)
    
    if(sbcdfi[i] > 0)
      lgsi[i] <- log(sbcdfi[i],2)
    
    i <- i+1
    
  }
  
  # compute the -log M(X)
  
  M <- rep(0,l)
  i <- 1
  
  while (i <= l){
    
    M[i] <- (sbcdf[i]+sbcdfi[i])/2
    if (M[i] > 0)
      M[i] <- -log(M[i],2)
    
    i <- i+1
    
  }
  
  # compute D(I,M) D(I,S)
  
  i <- 1
  div <- 0
  
  while (i <= l){
    
    im <- sbcdfi[i]*(lgsi[i]+M[i])
    is <- sbcdf[i]*(lgs1[i]+M[i])
    div <- div + im + is
    
    i <- i+1
  }
  
  return(div)
  
}

#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                            DIVERGENCE INDEX CONSTANT                                        #
#---------------------------------------------------------------------------------------------#
#' @function uni.divConst                                                                     #
#' @description \code{uni.divConst} uses \code{uni.div} to compute the distribution           #
#'               divergence constant for the maximum variance.                                #
#' @param n the number of bins in the pdf (pdf)                                               #
#' @return the Divergence index constant. (1 x 1) real number [0 1)                           #
#' @usage uni.divConst(\code{x})                                                              #
#' @seealso \code{uni.div}                                                                    #
#' @examples                                                                                  #
#' x <- c(0.5,0,0,0,0,0,0,0,0,0.5)                                                            #
#' r <- uni.div(x)       r: 0.2973122                                                         #
#' r <- uni.divConst(10) r: 0.2973122                                                         #
#' -------------------------------------------------------------------------------------------#

uni.divConst <- function(n){
  
  # create a bimodal distribution
  
  p <- rep(0,n) # pdf
  p[1] <- 0.5
  p[n] <- 0.5
  
  # compute the divergence index of the bimodal pdf
  
  const <- uni.div(p)
  
  return(const)
  
}

#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                               HOMOGENEITY INDEX                                             #
#---------------------------------------------------------------------------------------------#
#' @function uni.hom                                                                          #
#' @description \code{uni.hom} uses \code{uni.concI} and \code{uni.div} to compute            #
#'               the Homogeneity Index for a vector of observations.                          #
#' @param p numeric vector with non-negative integer values. (1 x n)                          #
#' @return the Homogeneity Index (HI) of \code{p}: real number  [0 1]                         #
#' @usage uni.hom(\code{p})                                                                   #
#' @seealso \code{uni.concI}, \code{uni.div}                                                  #
#' @examples                                                                                  #
#' p <- c(10,10,10,10,10)                                                                     #
#' r <- uni.hom(p)      r: 0                                                                  #
#' p <- c(0,0,50,0,0)                                                                         #
#' r <- uni.hom(p)      r: 1                                                                  #
#' -------------------------------------------------------------------------------------------#

uni.hom <- function(p){
  
  # compute the concentration index
  
  conc <- uni.concI(p)
  
  # compute the pdf of p
  
  d <- length(p)
  pop <- sum(p)
  i <- 1
  
  ## compute the pdf
  
  while(i <= d)
  {
    p[i] <- p[i]/pop
    i <- i+1
  }
  
  # compute the Divergence Index of the distribution
  
  div <- uni.div(p)
  
  # compute the Divergence Index for the Uniform distribution
  
  uniform <- rep(1/d,d)
  divUnif <- uni.div(uniform)
  
  # compute the Homogeneity Index
  
  H <- (conc + divUnif - div)/(1 + divUnif)
  
  return(H)
  
}

#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                          HOMOGENEITY INDEX - TRUE DIVERSITY                                 #
#---------------------------------------------------------------------------------------------#
#' @function uni.hom_mn                                                                       #
#' @description \code{uni.hom_mn} uses \code{uni.hom} to compute the Homogenity Index for a   #
#'               distribution of m bins and n equally abundant contiguous categories: pdf_mn  #
#' @param m integer indicating the number of bins in the distribution.                        #
#' @param n integer indicating the number of contiguous bins with equal number of observations#
#' @return the Homogeneity Index (HI) of pdf_mn: real number  [0 1]                           #
#' @usage uni.hom_mn(\code{m},\code{n})                                                       #
#' @seealso \code{uni.hom}                                                                    #
#' @examples                                                                                  #
#' r <- uni.hom_mn(5,5)  r: 0                                                                 #
#' r <- uni.hom_mn(5,1)  r: 1                                                                 #
#' p <- c(0.5,0.5,0,0,0)                                                                      #
#' r <- uni.hom(p)       r: 0.76                                                              #
#' r <- uni.hom(5,2)     r: 0.76                                                              #
#' -------------------------------------------------------------------------------------------#

uni.hom_mn <- function(m,n) {
  
  ## create a pdf_mn
  
  p <- rep(0,m)
  for(j in 1:n){
    p[j] <- 1/n
  }
  
  h <- uni.hom(p)
  
  return (h)
  
}

#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                          HOMOGENEITY INDEX - CLASSIFICATION                                 #
#---------------------------------------------------------------------------------------------#
#' @function uni.hom_class                                                                    #
#' @description \code{uni.hom_mn} uses \code{uni.hom} to compute the Homogenity class for a   #
#'              distribution                                                                  #
#' @param a real number indicating the lower bound of the first class - A                     #
#' @param b real number indicating the lower bound of the second class - B                    #
#' @param c real number indicating the lower bound of the third class - C                     #
#' @return String: Homogenity Class                                                           #
#' @usage uni.hom_class(\code{a},\code{b},\code{c},\{p})                                      #
#' @seealso \code{uni.hom}                                                                    #
#' -------------------------------------------------------------------------------------------#




uni.hom_class <- function(a,b,c,p) {
  
  h <- uni.hom(p)
  
  if (h >= a) {
    class <- 'A'
  } else if (h >= b) {
    class <- 'B'
  } else if (h >= c) {
    class <- 'C'
  } else {
    class <- 'D'
  }
  
  return (class)
}



#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                             LOCATION INDEX SUPPORTING FUNCTION                              #
#---------------------------------------------------------------------------------------------#
#' @function uni.locVec                                                                       #
#' @description function to compute the Concentration Location Index vector of a pdf          #
#' @param x probability density function vector. (1 x n)                                      #
#' @return the Location Index vector score (LIS) of \code{x}                                  #
#' @usage uni.locVec(\code{x})                                                                #
#' @examples                                                                                  #
#' x <- c(0.5,0,0,0,0.5)                                                                      #
#' v <- uni.locVec(x)   r: num[1:5] 0.6 0.6 0.6 0.6 0.6                                       #
#' x <- c(1,0,0,0,0)                                                                          #
#' x <- uni.locVec(x)   r: num[1:5] 1 0.8 0.6 0.4 0.2                                         #
#---------------------------------------------------------------------------------------------#

uni.locVec <- function(x){
  
  n <- length(x)
  i <- 1           ## iterator for the bins
  s <- 0           ## current interval score
  vs <- rep(0,n)   ## cumulative score
  k <- rep(0,n)    ## bins scores
  
  while(i <= n){
    
    j <- 0         ## iterator for the nested intervals (width)
    s <- 0
    vs <- rep(0,n)
    
    while(j < n){
      
      
      if(j == 0){
        
        s <- x[i]  ## interval width zero (initial point)
        
      }else{
        
        fw <- i+j  ## interval border right
        bk <- i-j  ## interval border left
        
        
        if (bk >= 1 && fw <= n){
          s <- s+x[bk]+x[fw]
        }else if (bk >= 1){
          s <- s+x[bk]
        }else if (fw <= n){
          s <- s+x[fw]
        }
        
      }
      vs[j+1] <- s
      j <- j+1
      
    }
    
    k[i] <- sum(vs)
    i<- i+1
    
  }
  
  ## Normalized score (singleton is n)
  
  k <- (1/n)*k
  
  
  return(k)
  
}

#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                                   LOCATION INDEX                                            #
#---------------------------------------------------------------------------------------------#
#' @function uni.loc                                                                          #
#' @description function to compute the Location Index (LI) and Compactness(C).               #
#' LI gives the minimum and maximum position of the bins with the maximum concentration       #
#' @param x probability density function vector. (1 x n)                                      #
#' @return the Location Index (LI) and Compactness of \code{x} (1x2)                          #
#' @usage uni.loc(\code{x})                                                                   #
#' @seealso \code{uni.locVec}                                                                 #
#' @examples                                                                                  #
#' x <- c(0.5,0,0,0,0.5)                                                                      #
#' v <- uni.loc(x)      r: 1 5                                                                #
#' x <- c(1,0,0,0,0)                                                                          #
#' x <- uni.loc(x)   r: 1 1                                                                   #
#' -------------------------------------------------------------------------------------------#

uni.loc <- function(x){
  
  n <- length(x)
  loc1 <- 0        ## minimum position of the bin with maximum score
  loc2 <- 0        ## maximum position of the bin with maximum score
  
  ## compute the Location Index vector and the maximum score
  
  v <- uni.locVec(x)
  m <- max(v)
  
  for (j in 1:n){
    if (v[j] == m){

      if (loc1 > 0)
          loc2 <- j
      else{
          loc1 <- j
          loc2 <- j
      }
    }
  }
  
  cl <- c(loc1,loc2)
  
  return (cl)
  
}


