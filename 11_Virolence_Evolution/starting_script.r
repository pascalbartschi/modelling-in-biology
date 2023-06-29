# Setup ----


setwd("/Users/bp/Uni/Computational/HS22/BIO445/BIO445_day11")
rm(list = ls())
if(!requireNamespace("pacman", quietly = T))   install.packages("pacman")
pacman::p_load("gganimate", "ggforce", "tidyverse", "av", "ape", "formatR")

# R toolbox ----
example <- list(1:5, rep(5,10), 10:5, c(3,6))
help(lapply)
help(sapply)


sapply(example,length)

sapply(example,mean)

lapply(example, sum)

lapply(example,function(x){
  return (x + 1)
})


help(match)
letters # the whole alphabet as a character vector
letters[1:3] # the letters a, b, and c

match("s",letters)

match("s",sample(letters))


mean_pos <- function(n,letters){
  all_letters <- matrix(rep(letters, n),ncol=n)
  all_letters <- apply(all_letters, MARGIN=2,function(x){
    return(sample(x))
  })
  all_position <- apply(all_letters, MARGIN=2,function(x){
    return(match("s",x))
  })
  mean_pos <- mean(all_position)
  return (letters)
}







pos = c()
for (i in 1:1000){
  p = match("s",sample(letters))
  pos = c(pos,p)
}

mean(pos)

# is package ape installed?
require(ape)



# Problem 1 ----
## 1.1
library(readr)
vl_tr <- read_csv("vl_tr_data.csv", 
                  col_names = FALSE)

colnames(vl_tr) <- c("vl", "tr")

vl_tr$logvl <- log10(vl_tr$vl)


vl_tr_fun <- approxfun(vl_tr$logvl, vl_tr$tr, method = "linear", rule = 2)

plot(seq(1, 7, by = 0.1), sapply(seq(1, 7, by = 0.1), vl_tr_fun), 
     type = "l", lty = 2, xlim = c(3, 6),
     xlab = "log10 viral load", ylab = "transmission rate")


## 1.2
vl_aids_data <- read_csv("vl_aids_data.csv", 
                         col_names = FALSE)

colnames(vl_aids_data) <- c("vl", "time_to_aids")

vl_aids_data$logvl <- log10(vl_aids_data$vl)

vl_aids_fun <- approxfun(vl_aids_data$logvl, vl_aids_data$time_to_aids, method = "linear", rule = 2)


plot(seq(1, 7, by = 0.1), sapply(seq(1, 7, by = 0.1), vl_aids_fun), 
     type = "l", lty = 2, xlim = c(3, 6),
     xlab = "log10 viral load", ylab = "time to aids")


## 1.3

fitness <- sapply(seq(1, 7, by = 0.1), vl_aids_fun) * sapply(seq(1, 7, by = 0.1), vl_tr_fun)
plot(seq(1, 7, by = 0.1), fitness, 
     type = "l", lty = 2, xlim = c(3, 6),
     xlab = "log10 viral load", ylab = "fitness")



## 1.4
population_size <- 1000
number_generations <- 500

sigma_environmental <- 0.8
sigma_mutation <- 0.1

population_genetic <- matrix(data = rep(0,population_size * number_generations), ncol = population_size, nrow = number_generations)
population_realized <- matrix(data = rep(0,population_size * number_generations), ncol = population_size, nrow = number_generations)
  
initial_vl <- 3
population_genetic[, 1] <- log10(3)
  
#for loop over the number of generations
for(k in 1:(number_generations - 1)) { 
  #extract the genetic component of virus load in the kth generation of the matrix
  genetic_vl <- population_genetic[,k]
    
  #add normal distributed environmental term to get the realized virus load
  population_realized[, k] <- 
        
  #get the fitness of the kth population
  fitness <- sapply(??, ??)
        
  #get the number of offspring by using a multinomial distribution
  #the population size remains constant but fitter individuals are
  #more likely to have offspring
  nr_offspring <- rmultinom(1, population_size, fitness)
        
  #genetic component is passed on to the offspring
  genetic_vl_unmut <- unlist(sapply(??, FUN = function(x) { rep(genetic_vl[x], nr_offspring[x]) } ))
        
  #add mutations to the viral component, approximated by normal distribution,
  #save it as genetic components of next generation
  population_genetic[, k + 1] <- ??
}

## 1.5
#plot from 1.3 (viral load vs. viral fitness)
plot(??)

#nice colors: start with blue colors and go more and more towards coral colored ones
palette <- colorRampPalette(c("darkblue", "coral"))

#which generation do we want to show?
seq <- seq(10, 80, by = 10)
palette <- palette(length(seq))

for(k in seq){ 
  lines(density(??), col = palette[k/10])
}

legend("topright", col = palette, legend = paste("gen", seq), lwd = 1)



# Problem 2 ----
## 2.1


## 2.2
tree_spvl <- read.tree(file = "tree_spvl")

## 2.3
#This function gives the label of a node pointing to the tip
find_node <- function(??, ??) {
  tip <- match(id, tree$tip.label)
  node <- tree$edge[match(as.numeric(??), tree$edge[, 2]), ][1]
  
  if(is.na(node) == "FALSE") {
    return(as.numeric(node))
  }
}

#This function returns c(NA,NA) if this patient (tiplabel) is not in a transmission
#pair and the name of the tiplabels if he/she is in a transmission pair
#(so the tiplabel of the patient and the other member of the pair)
find_pair <- function(??, ??) { 
  clade <- extract.clade(??)
  
  if (length(clade$tip.label) != ??) {
    return (c(NA, NA)) 
  } else { 
    return (??) 
  }
}

#This function returns all transmission pairs
transmission_pairs <- function(??) { 
  pair <- matrix(NA, length(tree$tip.label), 2)
  pair <- t(sapply(tree$tip.label, FUN = function(a) { find_pair(??) } ))
  pair <- pair[which(!is.na(??)), ]
  order <- function(row) { paste0(min(row), "_", max(row)) } 
  pair <- unique(sapply(pair, order))
  pair <- lapply(pair, FUN = function(a) { unlist(strsplit(a, "_")) } )
  return(??)
}

pairs <- ??


## 2.4
parent <- sapply(pairs, FUN = function(a) { unlist(a)[??] } )
offspring <- sapply(pairs, FUN = function(a) { unlist(a)[??] } )
??



# Functions ----
find_node <- function(id, tree) {
  tip <- match(id, tree$tip.label)
  first_node <- tree$edge[match(as.numeric(tip), tree$edge[, 2]), ][1]
  
  if(is.na(first_node) == "FALSE") {
    return(as.numeric(first_node))
  }
}

find_pair <- function(tree, id){
  clade <- extract.clade(tree, find_node(id, tree))
  
  if (length(clade$tip.label) != 2) { 
    return (c(NA, NA)) 
  } else { 
    return (clade$tip.label) 
  }
}

transmission_pair<-function(tree){
  pair <- lapply(tree$tip.label, FUN = function(a) {
    find_pair(tree, a) } )
  
  pair <- pair[sapply(pair, function(a) { !is.na(a[1]) } )]
  pair <- unique(lapply(pair, function(row) {paste0(min(as.numeric(row)), "_", max(as.numeric(row))) } ))
  pair <- lapply(pair, FUN = function(a) {unlist(strsplit(a, "_")) } )
  return(pair)
}
