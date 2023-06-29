pacman::p_load("ape", "phangorn", "tidyverse")
library(ggtree)


distdna <- function(list, modelused){   #give model in " "
  dist <- dist.dna(as.DNAbin(list), model = modelused)
  return(dist)
}

env_plasma <- read.dna(file = "env_plasma.fasta" , format = "fasta",
                     as.character = TRUE, as.matrix = FALSE)
env_pbmc <- read.dna(file = "env_pbmc.fasta" , format = "fasta",
                       as.character = TRUE, as.matrix = FALSE)

pol_plasma <- read.dna(file = "pol_genomic.fasta" , format = "fasta",
                as.character = TRUE, as.matrix = FALSE)
pol_pbmc <- read.dna(file = "pol_pbmc.fasta" , format = "fasta",
                       as.character = TRUE, as.matrix = FALSE)


env_dist_plasma <- upgma(distdna(env_plasma, "JC69"))
env_dist_pbmc <- upgma(distdna(env_pbmc, "JC69"))

pol_dist_plasma <- upgma(distdna(pol_plasma, "JC69"))
pol_dist_pbmc <- upgma(distdna(pol_pbmc, "JC69"))




ournames <- c("RS0544", "RS0547", "RS0548")

pol_dist_pbmc$tip.label <- pol_dist_plasma$tip.label <- env_dist_pbmc$tip.label <- env_dist_plasma$tip.label <- ournames 

# ggtree(pol_dist_plasma, branch.length = 0.01) + 
#   geom_tiplab() + 
#   labs(title = "HIV-polymerase")
# 
# ggsave("polymerase_plot.png", height = 5, width = 10, units = "cm")

# plot(pol_dist_pbmc)
plot(pol_dist_plasma, main = "HIV polymerase")
# plot(env_dist_pbmc)
plot(env_dist_plasma, main = "HIV envelope")
