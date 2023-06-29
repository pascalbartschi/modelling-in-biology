# Molecular Epidemiology
# prepared by Jessy Duran Ramirez
# 14.12.22


# R preparations ---------------------------------------------------------------------------------
# Remove current working environment and load packages

rm(list = ls())

if(!requireNamespace("pacman", quietly = T))   install.packages("pacman")
pacman::p_load("ape", "phangorn", "tidyverse")

require(ape)
# setwd("your path here")
primates <- read.dna(file = "sequences_primates.fasta",
                     format = "fasta", as.character = T,
                     as.matrix = F)


# Exercise 1: Measure the similarity between sequences  -------------------------------------------

## a) Write an R function to calculate the Hamming distance

hamming_2seq <- function(seq1, seq2){
  len <- length(seq1)
  dist <- rep(NA, len)
  
  for (i in 1:len){
    ifelse(seq1[i] == seq2[i], dist[i] <- 0, dist[i] <- 1)
  }
  
  return(list("vector" = dist, "dist" = sum(dist) / length(dist)))
}

# hamming_alignment <- function(alignment){
#   n <- length(names(alignment))
#   m <- matric(nrow = n, ncol = length(alignment[[1]]))
#   for (i in 1:n){
#     for (j in 1:m)
#   }
# }

## b) Implement the JC69 model, write an R function

hamming_matrix <- function(li){
  n <- length(names(li))
  m <- matrix(nrow = n, ncol = n)
  for (i in 1:n){
    for (j in i:n){
      m[j, i] <-  m[i, j] <- sum(hamming_2seq(li[[i]], li[[j]])$vector)
    }
  }
  colnames(m) <- rownames(m) <- names(li)
  
  return(m)
}

length_checker <- function(li){
  n <- length(li)
  len <- rep(NA, n)
  for (i in 1:n){
    len[i] <- length(li[[i]])
  }
  ifelse(mean(len) %in% len, NA, stop("Read in fasta files with sequences all of same length!"))
}


p_hat_matrix <- function(li){
  
  length_checker(li)
  
  n <- length(names(li))
  m <- matrix(nrow = n, ncol = n)
  for (i in 1:n){
    for (j in i:n){
      d <- hamming_2seq(li[[i]], li[[j]])
      p <- d$dist / length(d$vector)
      m[j, i] <-  m[i, j] <- (-3 / 4) * log(1 - (4 / 3) * p)
    }
  }
  colnames(m) <- rownames(m) <- names(li)
  
  return(round(m, 2))
}

## c) TASK: add an additional line to your code to check for the same length

primates_diff_len <- read.dna(file = "sequences_primates_unequal_length.fasta",
                     format = "fasta", as.character = T,
                     as.matrix = F)

hamming_matrix(primates)
p_hat_matrix(primates)
p_hat_matrix(primates_diff_len)

# stop is for stop and warning is for warning but continoues

## d)  calculating the distance 

ape::dist.dna(as.DNAbin(primates), model = "raw")
as.data.frame(as.matrix(ape::dist.dna(as.DNAbin(primates), model = "JC69")))
round(ape::dist.dna(as.DNAbin(primates), model = "K80"), 2)

distance_matrix <- hamming_matrix(primates)

## YOUR CODE HERE ##



## e)  Compare the three different distance matrices


# Exercise 2: A distance based method to build a phylogenetic tree (UPGMA)  ---------------------------

## a) TASK: Understand the Pseudocode
## DISCUSS ## 



## b) TASK: Manually apply the hamming distance

# check your result:

UPGMA_step <- function(d){
  n <- length(colnames(d))
  
  min_ <- max(d)
  min_ij <- c(NA, NA)
  
  for (i in 1:n){
    for (j in i:n){
      if ((d[i, j] < min_) & (d[i, j] != 0)){
        min_ <- d[i, j]
        min_ij <- c(i, j)
      }
      # print(d[i, j])
    }
  }
  i <- min_ij[1]
  j <- min_ij[2]
  new_ <-  (d[, i] + d[, j]) / 2
  new_names <- names(new_)[c(i, j)]
  new_name <- paste0(new_names[1], "/", new_names[2])
  new_ <- new_[-c(i, j)]
  new_[new_name] <- 0
  
  
  d <- d[-c(i, j), -c(i, j)]
  d <- rbind(d, new_[-c(length(new_))])
  d <- cbind(d, new_)
  
  colnames(d) <- rownames(d) <- names(new_)
  # d <- cbind()
  
  return(d)#list(new_, d))
}


### c) 1. BONUS TASK, R skeleton provided
# please only do this after having completed all other exercises
## YOUR CODE HERE ##


# Exercise 3: Application of phylogenetics: the origin of HIV ----------------------------------------
## TASK: Have a look at the HIV_SIV_alignment.fasta

## YOUR CODE HERE ##
HIV_SIV<- read.dna()


# Exercise 4 - Application of phylogenetics: Phylogenetics in court ----------------------------------
## TASK: READ & DISCUSS ## 

## TASK: Have a look at the criminal_case_aligned.fasta

## YOUR CODE HERE ##
CRIME <- read.dna()


# 2. Bonus ----------------------------------
## if you want, can make prettier plots using ggtree::ggtree(), 
## already used it yesterday

## YOUR CODE HERE ##
