# Title       : create_Sa3.R
# Description : This file contains the script to generate the Sa3 table of the publication
#               " A framework for the classifcation and identification of homogenous socioeconomic
#                 areas in the analysis of health care variation " IJHG
# Author      : Ludovico Pinzari
# Usage       : Run the uni.homTable function and the code below the main section
# R version   : 3.3 or later version
# Notes       : For the execution of this script you must run the hi_li.R first
#==================================================================================================#


#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                             GENERATE TABLE                                                  #
#---------------------------------------------------------------------------------------------#
#' @function uni.homTable                                                                     #
#' @description function to compute the HI, CI, DI and LI for a set of  decile distributions  #
#' @param x a data frame (table) n x 15. Columns:                                             #
#' [id, SA3_code, SA3_name, State_code, State_name, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10]  #
#' @return data frame df n x 20                                                               #
#' @usage uni.homTable(\code{x})                                                              #
#' @example                                                                                   #
#' t <- uni.homTable(x)                                                                       #
#' -------------------------------------------------------------------------------------------#

uni.homTable <- function(x){
  
  # data frame initialization
  
  df <- data.frame(id = integer(),
                   SA3_code = character(),
                   SA3_name = character(),
                   State_code = character(),
                   State_name = character(),
                   d1 = integer(),
                   d2 = integer(),
                   d3 = integer(),
                   d4 = integer(),
                   d5 = integer(),
                   d6 = integer(),
                   d7 = integer(),
                   d8 = integer(),
                   d9 = integer(),
                   d10 = integer(),
                   Hom = double(),
                   CI = double(),
                   DI = double(),
                   LI = integer(),
                   CL = character(),
                   stringsAsFactors = FALSE
  )
  
  # definition of Homogeneity class lower bound
  
  cl_a <- uni.hom_mn(10,4)
  cl_b <- uni.hom_mn(10,5)
  cl_c <- uni.hom_mn(10,6)
  
  # compute the number of rows
  
  rows <- nrow(x)
  
  
  for(i in 1:rows){
    pop <- as.numeric(x[i,6:15])
    
    ## compute the pdf
    
    tot <- sum(pop)
    pdf <- pop/tot
    
    ## copy the first fifteen columns
    
    df[i,1] <- as.numeric(x[i,1])   ## Sa3 id
    df[i,2] <- as.character(x[i,2]) ## Sa3 code
    df[i,3] <- as.character(x[i,3]) ## SA3 name
    df[i,4] <- as.character(x[i,4]) ## State code
    df[i,5] <- as.character(x[i,5]) ## State name
    df[i,6] <- as.numeric(x[i,6])   ## decile 1
    df[i,7] <- as.numeric(x[i,7])   ## decile 2
    df[i,8] <- as.numeric(x[i,8])   ## decile 3
    df[i,9] <- as.numeric(x[i,9])   ## decile 4
    df[i,10] <- as.numeric(x[i,10]) ## decile 5
    df[i,11] <- as.numeric(x[i,11]) ## decile 6
    df[i,12] <- as.numeric(x[i,12]) ## decile 7
    df[i,13] <- as.numeric(x[i,13]) ## decile 8
    df[i,14] <- as.numeric(x[i,14]) ## decile 9
    df[i,15] <- as.numeric(x[i,15]) ## decile 10
    
    ## compute statistical indices
    
    df[i,16] <- uni.hom(pop)        ## compute the Homogenetiy index
    df[i,17] <- uni.concI(pop)      ## compute the Concentration Index
    div <- uni.div(pdf)             ## compute the Divergence Index
    const <- uni.divConst(10)       ## ratio constant
    df[i,18] <- div/const
    
    loc <- uni.loc(pdf)             ## compute the location Index
    
    ## loc is a vector 1 x 2
    ## since we do not have bimodal pdf the two values are equal
    
    df[i,19] <- loc[1]              ## Location Index
    
    ## Homogeneity classifcation
    
    df[i,20] <- uni.hom_class(cl_a,cl_b,cl_c,pdf) ## Homogeneity class
    
    
  }
  
  return (df)
  
}

#---------------------------------------------------------------------------------------------#
#                                                                                             #
#                                       MAIN                                                  #
#---------------------------------------------------------------------------------------------#


## load data and libraries

library(data.table)
input <- read.csv("dataSa3.csv")
x <- data.frame(lapply(input,as.character),stringsAsFactors=FALSE)

## create the table with the HI, CI, DI,LI and HI class.

t <- uni.homTable(x)

## save the table as csv file

write.csv(t, file = "SA3db")






