##### STARTING SCRIPT - ECOLOGICAL NETWORKS

# SETUP #
rm(list=ls())

#libraries
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rjson, data.table, bipartite, igraph, Matrix)



# look at one specific community

##### 1. Network Science -----------------------------------------------

# A) ----------
download_and_write_csv_single_network <- function(network_name){
  networkName <- network_name
  speciesName <- "yes"
  url <- paste("http://www.web-of-life.es/download/", networkName, "_", speciesName, ".csv", sep = "")
  data <- fread(url)
  write.csv(data, file = paste(trimws(networkName), ".csv", sep = ""))
  return(data)
}

data <- download_and_write_csv_single_network( "M_PL_052")


# B) ----------
n_row <- 15 # plants
n_col <- 39 # pollinators
matrix <- as.matrix(data, rownames = "V1")
heatmap(matrix, scale = "none")

poll_edges <- rowSums(matrix) # sums of all pollinators
plant_edges <- colSums(matrix) # sums of all plants

max_name <- function(edges){
  name <- edges[edges == max(edges)]
  
  return(names(name))
}
max_poll <- max_name(poll_edges)
max_plant <- max_name(plant_edges)

# C) ----------
full_adjm <- function(data, n_row,n_col) {
  full_adjm <- matrix(0, n_row + n_col, n_row + n_col) # square zero matrix
  full_adjm[(n_row + 1):(n_row + n_col), 1:n_row] = t(as.matrix(data > 0))
  full_adjm[1:n_row, (n_row + 1):(n_row + n_col)] = as.matrix(data > 0)
  return(full_adjm)
}

rownames(data) <- data[["V1"]]
adjm <- full_adjm(data, n_row, n_col)
net <- graph_from_adjacency_matrix(adjm, mode = "undirected")
colrs <- c("green", "gold")
V(net)$color <- colrs[as.numeric(V(net) > n_row) + 1]
V(net)$label <- append(row.names(data), names(data))
plot(net)
plot(net, layout = layout_in_circle)



# D) ----------
# histograms of poll_edges and row_egdes


poll_edges <- rowSums(matrix) # sums of all pollinators
plant_edges <- colSums(matrix) # sums of all plants
hist(as.vector(plant_edges), main = "Hist Plant Edges", breaks = c(0:10))
hist(as.vector(poll_edges), main = "Hist poll Edges", breaks = c(0:20))



# E) ----------
# Hint: Density is the ratio of observed interactions to all possible interactions in a network

density1 <- function(matrix){
  result <- matrix / sum(matrix)
  return(result)
}

density_matrix <- density1(matrix)
poll_edges <- rowSums(density_matrix) # sums of all pollinators
plant_edges <- colSums(density_matrix) # sums of all plants
hist(as.vector(plant_edges), main = "Hist Plant Edges", breaks = seq(0,0.2, 0.001))
hist(as.vector(poll_edges), main = "Hist poll Edges", breaks = seq(0,0.3, 0.001))

# F) ----------
# Hint: write a function that returns the nestedness

nestedness_binmatnest2 <- function(matrix) {
  # The following function calculates the nestedness temperature, which ranges between 0 (fully nested) and 100 (not nested)
  value <- as.numeric(nested(matrix, "binmatnest"))
  value <- (100 - value) / 100
  # modify value here as indicated in the exercise sheet
  return(value)
}

nestedness <- nestedness_binmatnest2(matrix)

# G) ----------
# Hint: write a function that takes a matrix as input. 
#       First generate an empty matrix of the same size.
#       Next read in Prof. Bascompte´s paper, how the empty matrix will be filled in, depending on the input matrix

nullmodel1 <- function(matrix){
  connections_row <- 0
  connections_col <- 0
  while (connections_row == 0 | connections_col == 0){
    prob <- sum(matrix) / (n_row * n_col)
    samples_ <- rbinom((n_row * n_col), 1, prob = prob)
    zero_matrix <- matrix(samples_, nrow = n_row, ncol = n_col)
    connections_row <- min(rowSums(zero_matrix))
    connections_col <- min(colSums(zero_matrix))
  }
  
  return(zero_matrix)
  
}

nullmodel2 <- function(matrix){
  connections_row <- 0
  connections_col <- 0
  prob_rows_matrix <- matrix(rep(as.vector(rowSums(matrix) / n_col), n_col), nrow = n_row, ncol = n_col)
  prob_cols_matrix <- matrix(rep(as.vector(colSums(matrix) / n_row), n_row), nrow = n_col, ncol = n_row)
  prob_matrix <- (prob_rows_matrix + prob_rows_matrix) / 2
  while (connections_col == 0 | connections_row == 0){
    zero_matrix <- matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix))
    for (i in 1:nrow(matrix)){
      for (j in 1:ncol(matrix)){
        proba <- prob_matrix[i, j]
        # print(proba)
        zero_matrix[i, j] <- sample(c(0, 1), size = 1, prob = c(1 - proba, proba))
        
      }
    }
    connections_row <- min(rowSums(zero_matrix))
    connections_col <- min(colSums(zero_matrix))
  }
  return(zero_matrix)
}




# H) ----------
# Hint: you have a function to calculate nestedness, and you have a function that generates you a random nullmodel. 

null_simulation1 <- function(N_Steps, matrix){
  nestedness0 <- rep(NA, N_Steps) 
  for (i in 1:N_Steps){
    zero_m <- nullmodel1(matrix)
    nestedness0[i] <- nestedness_binmatnest2(zero_m)
  }
  return(nestedness0)
}

null_simulation2 <- function(N_Steps, matrix){
  nestedness0 <- rep(NA, N_Steps) 
  for (i in 1:N_Steps){
    zero_m <- nullmodel2(matrix)
    nestedness0[i] <- nestedness_binmatnest2(zero_m)
  }
  return(nestedness0)
}

nestedness01 <- null_simulation1(1000, matrix)
nestedness02 <- null_simulation2(1000, matrix)
our_nestedness <- nestedness_binmatnest2(matrix)

# p_value <- length(nestedness0[nestedness0 > our_nestedness]) / length(nestedness0)


# plotting 1

pnorm(our_nestedness, mean(nestedness01), sd(nestedness01), lower.tail = F) 

normal_distr <- rnorm(1000, mean = mean(nestedness01), sd=sd(nestedness01))

library(ggplot2)
nested_dd_null <- as.data.frame(nestedness01)
ggplot(nested_dd_null, aes(nestedness01))+
  #geom_histogram(bins=20, linetype=8)
  geom_density(kernel = "gaussian" , aes(normal_distr, color = "gaussian"), col="blue", linewidth=1, ,linetype="dashed")+
  geom_density(kernel = "gaussian", outline.type= "upper", col="green", linewidth=2)+
  xlim(0,1)+
  geom_vline(aes(xintercept=mean(our_nestedness)), col="red", linewidth=2)+
  #ggtitle("Simulated distribution of nestedness") + 
  labs(x="Temperature of nestedness", y="Frequency")+
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        title = element_text(size = 28))

ggsave("plot_model1.png", height = 15, width = 20, units = "cm")

# plotting 2

pnorm(our_nestedness, mean(nestedness02), sd(nestedness02), lower.tail = F) 

normal_distr <- rnorm(1000, mean = mean(nestedness02), sd=sd(nestedness02))

library(ggplot2)
nested_dd_null <- as.data.frame(nestedness02)
ggplot(nested_dd_null, aes(nestedness02))+
  #geom_histogram(bins=20, linetype=8)
  geom_density(kernel = "gaussian" , aes(normal_distr), col="blue", size=1, ,linetype="dashed" )+
  geom_density(kernel = "gaussian", outline.type= "upper", col="green", size=2)+
  xlim(0,1)+
  geom_vline(aes(xintercept=mean(our_nestedness)), col="red", size=2)+
  # ggtitle("Simulated distribution of nestedness") + 
  labs(x="Temperature of nestedness", y="Frequency")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        title = element_text(size = 28))

ggsave("plot_model2.png", height = 15, width = 20, units = "cm")



##### 2. Scale-free networks -----------------------------------------------

# B) ---------- 
# create a network with preferential attachement
edgelist <- c(1, 2) # list of the edges
degreelist <- c(1, 1) # list of the degrees of each node

# Hint: use a loop. You already have two connected nodes, you will start with the third
#       sample.int(number_of_existing_nodes, number_of_edges, replace = F, prob = degreelist) should help you. start off with allowing only 1 edge per new node

A <- make_graph(edgelist, directed = FALSE)
V(A)$label <- ""
V(A)$size <- degreelist + 2
plot(A)
