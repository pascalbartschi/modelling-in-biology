# Modeling Cancer
# 9.12.22

# Remove current working environment
# rm(list = ls())

if(!requireNamespace("pacman", quietly = T))   install.packages("pacman")
pacman::p_load("deSolve", "ggforce", "tidyverse")



# Exercise 1 - Estimate tumor growth rate ------------------------------------------------------



## a) Load the database into R and explore it
db <- read.csv2("Nakamura_db_BT_start.csv")

# TASK: Explore and get familiarized with the dataset

db <- db %>%
  mutate(abs_growth = Latest_tumor_vol - Initial_tumor_vol) %>%
  mutate(abs_growth_rate = abs_growth / (Follow.up.time / 12)) %>%
  mutate(rel_growth_rate = ((Latest_tumor_vol / Initial_tumor_vol) ** (1 / (Follow.up.time / 12)) - 1) * 100)


quantile(db$rel_growth_rate, probs = c(0.2, 0.8))
# FOR YOU TO COMPARE: if your relative growth-rate is correctly calculated, 
# quantile(db$relative.growth.rate, probs = c(0.2, 0.8)) should return 3.07455 25.21822 



## c) Visualization
# Make a plot of the age of the patient versus absolute growth rate.
# Is there an apparent correlation?

ggplot(data = db, mapping = aes(x = abs_growth_rate, y = Age)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw()



## d) Linear regression
# Construct a linear regression using lm(), with age and absolute growth rate and plot the regression fit.

linmod <- lm(Age ~ abs_growth_rate, data = db)
db$fit.age_abs_grate <- db$abs_growth_rate * linmod$coefficients["abs_growth_rate"] + linmod$coefficients["(Intercept)"]

ggplot(data = db) + 
  geom_point(aes(x = abs_growth_rate, y = Age)) + 
  geom_line(aes(x = abs_growth_rate, y = fit.age_abs_grate), color = "blue") + 
  geom_abline(intercept = coef(linmod)[1], slope = coef(linmod)[2], color = "green") + 
  theme_bw()



## e) Calculate the tumor doubling
db$doubling.time =  (db$Follow.up.time / 12)* (log10(2) / (log10(db$Latest_tumor_vol / db$Initial_tumor_vol)))
quantile(db$doubling.time, probs = c(0.2, 0.8))
# FOR YOU TO COMPARE: if your relative growth-rate is correctly calculated, 
# quantile(db$relative.growth.rate, probs = c(0.2, 0.8)) should return 3.082191 22.889495 



### e.1) 
# Do tumors that were larger at the first screen also tend to have higher doubling time?
# How can you test this hypothesis?

linmod2 <- lm(doubling.time ~ Initial_tumor_vol,  data = db)
p_linmod2 <- summary(linmod2)$coefficients[, 4][[2]]


### e.2)
# Do calcified tumors have a higher doubling time?
# How can you test this hypothesis?
# Think of a way to visualize your conclusion.

# explore

ggplot(data = db) + 
  geom_boxplot(aes(x = Calcification, y = doubling.time)) + 
  geom_jitter(aes(x = Calcification, y= doubling.time), width = .2) # jitter is point

ggplot(data = db) + 
  geom_violin(aes(x = Calcification, y= doubling.time), scale = "width", alpha = .3) + # adjusts scaling
  geom_boxplot(aes(x = Calcification, y = doubling.time), fill=NA, outlier.shape = NA,
               width = .2)
   # jitter is point

# approximately normal distributed for t.test

ggplot(data = db) + 
  geom_histogram(aes(x = log(doubling.time))) + 
  facet_wrap(~Calcification)

# ttest

ttest1 <- t.test(log(db$doubling.time) ~ db$Calcification)
kalktest <- kruskal.test(db$doubling.time ~ db$Calcification) # rank based test

db$Calc_binary <- ifelse(db$Calcification == "yes", 1, 0) 
              

ggplot(data = db, aes(x=doubling.time, y=Calc_binary)) +
  geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = T) 

glm(Calc_binary ~ doubling.time, data = db, family = "binomial")


# Exercise 3 - Modeling with ODEs ----------------------------------------



## a)
# define the ODE 
Gompertz <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dN  <- r * N * log(K/N)
    return(list(c(dN)))
  })
}

pars  <- c(r = 0.015, K = 40)  # parameters   
yini  <- c(N = 0.001)          # initial values
times <- seq(0, 2000, by = 1)  # time in days



## b)
# calling the ODE solver
out   <- as.data.frame(ode(yini, times, Gompertz, pars))

# plot the tumor volume over time

ggplot(data = out, mapping = aes(x = time, y = N)) + 
  # geom_point() + 
  geom_line(linewidth = .5, color = "darkgreen") +
  theme_bw()

## c) Effect of therapy
# Now let's examine the effect of therapy (chemotherapy or immunotherapy) by adding a right hand side negative term



### c.1) with treatment
Gompertz_Treat1 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity), then introduce treatment at day 1095 (3 years)
    # HINT: overwrite alpha before time is 1095
    # alpha <- ifelse(Time < 1095, 0, alpha)
    dN  <- (r * N * log(K/N)) - (((Time > 1095) * alpha) * c * N)
    
    return(list(c(dN)))
  })
}
pars <- c(r = 0.015, K = 40, alpha = 0.7, c = 0.3)
yini <- c(N = 0.001)
times <- seq(0, 2000, by = 1)

out <- as.data.frame(ode(yini, times, Gompertz_Treat1, pars))

# plot the tumor volume considering treatment:
# HINT:
# if you're using ggplot, have a look at facet_zoom from the ggforce package

ggplot(data = out, mapping = aes(x = time, y = N)) + 
  geom_point() + 
  geom_line(linewidth = .5, color = "darkgreen") +
  ggforce::facet_zoom(xlim = c(0, 300)) + 
  geom_abline(intercept = max(out$N) / 2, slope = 0) +
  theme_bw()


# 50% of the tumor volume after treatment initiation
  
days2half <- out %>% 
  filter(time > 1095 & N > (min(out$N) + max(out$N)) / 2) %>%
  nrow() + 1



### c.2)
# Elimination: how many days does it take for tumor elimination?
# calculate and visualize
# the model is updated for c decline

treatment_start <- 1095

days2tumor_elim <- out %>%
  filter(time > treatment_start & N < 0.001) %>%
  nrow()

ggplot(data = out, mapping = aes(x = time, y = N)) + 
  geom_point() + 
  geom_line(linewidth = .5, color = "darkgreen") +
  ggforce::facet_zoom(xlim = c(treatment_start, treatment_start + days2tumor_elim),
                      ylim = c(0, 0.001)) + 
  geom_abline(intercept = max(out$N) / 2, slope = 0) +
  labs(title = "days to tumour elimination") + 
  theme_bw()


# trials ==============================================================================


### c.3) 
c_transform <- rep(1, 2000)
for (i in 2:2000){
  c_transform[i] <- c_transform[i - 1]  * 0.75
}

index <- 1
t <- 0
i <- 1

# drug concentration decreases with time
Gompertz_Treat2 <- function(time, State, Pars, list) {
  N = State[2]
  C = State[1]
  
  
  # list(r = 0.015, K = 40, alpha = 0.7, treatments = c(2, 1095, 1105, 1125))
  r = Pars[1]
  K = Pars[2]
  alpha = Pars[3]
  treatments = list
  
  if(time %in% treatments){dC <- C + 0.3} else {dC <- -0.25 * C} 
  dN  <- (r * N * log(K/N)) - (alpha * C * N)
  list(c(dC, dN))
}

# solve_function = function(times, y, params){
#   return(Gompertz_Treat2(times = times, State = y, Pars = params))
# }
Gompertz_Treat2
out = ode(method = "lsoda",
    y = c(C = 0, N = 0.001), 
    time = seq(0, 2000, by = 1), 
    func = Gompertz_Treat2, 
    parms = c(r = 0.015, K = 40, alpha = 0.7),
    list = c(2, 1095, 1105, 1125))

pars <- list(r = 0.015, K = 40, alpha = 0.7, treatments = c(2, 1095, 1105, 1125))
yini <- c(N = 0.001, C = 0)
times <- seq(0, 2000, by = 1)


out <- as.data.frame(ode(yini, times, Gompertz_Treat2, pars))

# Visualize Tumor dynamics considering the decreasing drug concentration

ggplot(data = out, mapping = aes(x = time, y = C)) + 
  geom_point() + 
  geom_line(linewidth = .5, color = "darkgreen") +
  ggforce::facet_zoom(xlim = c(1090, 1125)) +
  geom_abline(intercept = max(out$N) / 2, slope = 0) +
  theme_bw()


# 50% of the tumor volume after treatment initiation
## YOUR CODE HERE ##


# Elimination:
## YOUR CODE HERE ##



### c.4)
# increasing the drug effectivness
# increase the effectivness (alpha)
pars <- c(r = 0.015, K = 40, alpha = ## YOUR CODE HERE ## ) 
yini <- c(N = 0.001, C = 0.3)
times <- seq(0, 2000, by = 1)

# Plot the tumor rebound and calculate time to 50% of the tumor volume after treatment initiation

## YOUR CODE HERE ##


# Elimination:
## YOUR CODE HERE ##






# Exercise 4 - Resistance -------------------------------------------------------------

# Now let’s assume that 1% of the treated cells develop resistance to treatment.
# Modify the system accord- ingly (assume that c(t) is constant again).
Gompertz_Treat_Resistance1 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity),
    # then introduce treatment at day 1095 (3 years)
    # in this the drug concentration also has to descrease only
    # since 1095 (before it was not in the body..)
    # let cells leaving the N compartment become resistant with rate m (CAVE: gompertz, not exponential)
    
    ## YOUR CODE HERE ##
    
    return(list(c(dN, dC, dR)))
  })
}



## a.) 
# One month after the first treatment, what would be the volume of the resistant cells?

## YOUR CODE HERE ##



## b.)
# Plot the dynamics of both cells types. When will the resistant cells predominate?
# Plot both strains

## YOUR CODE HERE ##
  


## c.*)
# Now assume that resistant cells can revert back and become susceptible again,
# with a certain rate ”q”, modify the system accordingly and explore.
Gompertz_Treat_Resistance2 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity),
    # then introduce treatment at day 1095 (3 years)
    # in this the drug concentration also has to descrease only
    # since 1095 (before it was not in the body..)
    # let cells leaving the N compartment become resistant with rate m (CAVE: gompertz, not exponential)
    # let resistant cells revert with rate q 
    
    ## YOUR CODE HERE ##
    
    return(list(c(dN, dC, dR)))
  })
}

# Plot both strains

## YOUR CODE HERE ##
