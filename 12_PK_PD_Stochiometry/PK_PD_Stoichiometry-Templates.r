######################################################################
### Problem 2
######################################################################

### a) paramters for KD adb receotir occupancy RO
setwd("~/BSc_UZH/UZH_22HS/BIO445_BK/12_PK_PD_Stochiometry")

E <- function(conc, E0, Emax, EC50, gamma) {
  
  tot_eff <- E0 + ((Emax * conc ** gamma) / (EC50 ** gamma + conc ** gamma))
    
  return(tot_eff)
}

### b) minimum efficacy and potency

Emax <- 90 # maximum treatment effect
E0 <- 10 # baseline effect
EC50 <- 0.3 # 50 % of drug effect
# Hill coefficient
gammaA <- 1
gammaB <- 4
gammaC <- 0.25
conc <- seq(0, 1.5, by = 0.015)

# drug_A <- drug_B <- drug_C <- rep(NA, 10)

# for (i in 1:length(conc)){
#   drug_A[i] <- E(conc[i], E0, Emax, EC50, gammaA)
#   drug_B[i] <- E(conc[i], E0, Emax, EC50, gammaB)
#   drug_C[i] <- E(conc[i], E0, Emax, EC50, gammaC)
# }
drug_A <- E(conc, E0, Emax, EC50, gammaA) # possible to directly plug in without for loop
drug_B <- E(conc, E0, Emax, EC50, gammaB)
drug_C <- E(conc, E0, Emax, EC50, gammaC)

drug_frame <- data.frame(values = c(drug_A, drug_B, drug_C),
                         conc = rep(conc, 3),
                         drug = c(rep("drug_A", length(drug_A)), 
                                  rep("drug_B", length(drug_B)),
                                  rep("drug_C", length(drug_C))))

library(ggplot2)
library(ggforce)
library(tidyverse)
           
p <- ggplot(data = drug_frame, mapping = aes(x = conc, y = values, color = drug)) + 
  geom_line() + 
  geom_vline(xintercept = EC50, linetype = "dashed") + 
  labs(x = "concentration", y = "effect size", title = "Sigmoid curves of different drugs") + 
  theme_bw()

# required: minimum therapeutic effect, c > too toxic

min_ther <- 20 # min therapeutic effect
max_c <- 0.25 # max concentration that can be taken

p + ggforce::facet_zoom(xlim = c(0, max_c), ylim = c(min_ther, 100)) +
  labs(title = "Drug C is best for (i)")

ggsave("2_b_i_drug_C.png", height = 10, width = 20, units = "cm")

ther_window <- c(0.45, 1.2) # wished therapeutic conc

p + ggforce::facet_zoom(xlim = ther_window, ylim = c(min_ther, 100)) + 
  labs(title = "Drug B is best for (ii)")

ggsave("2_b_i_drug_B.png", height = 10, width = 20, units = "cm")

### Problem 2: direct response and biophase distribution models

source("PD_scripts.r")

some_frame <- as.data.frame(conc_Blood_Brain())
exploratory_frame <- pivot_longer(data = some_frame,
                                  cols = c(Blood, Brain),
                                  names_to = "Type",
                                  values_to = "C")

ggplot(data = exploratory_frame, 
       mapping = aes(x = Time, y = C, color = Type)) + 
  geom_line() + 
  labs(x = "Time [min]", y = "Concentration", 
       title = "Drug concetration and blood brain barrier") + 
  theme_bw()
  
ggsave("3_Exploratory_concBloodBrain_standard_params.png", 
       height = 10, width = 20, units = "cm")

# for quick exploration
plot(conc_Blood_Brain())

#** Problem b)
E0 <- 0
Emax <- 100
EC50 <- 1
gamma <- 2


out_frame <- conc_Blood_Brain()

effect_brain <- E(out_frame[,"Brain"],E0, Emax, EC50, gamma)

# determines highest effect
max_EGG <- which.max(effect_brain)

ggplot(mapping = aes(x = time_line, y = effect_brain)) + 
  geom_vline(xintercept = conc_brain[max_EGG], color = "red") +
  geom_line() + 
  theme_bw()


exploratory_frame <- pivot_longer(data = out_frame,
                                  cols = c(Blood, Brain),
                                  names_to = "Type",
                                  values_to = "C")

ggplot(data = exploratory_frame, 
       mapping = aes(x = Time, y = C, color = Type)) + 
  geom_line() + 
  labs(x = "Time [min]", y = "Concentration", 
       title = "Drug concetration and blood brain barrier")
theme_bw()