# Model by Chrisitan Althaus
# Necessary libraries
# Necessary libraries
if(!requireNamespace("pacman", quietly = T))   install.packages("pacman")
pacman::p_load("deSolve", "chron", "tidyverse")


# Read the data)
ebola <- read.csv("Ebola_outbreak_West_Africa_data.csv")
# write.csv(ebola, file = "Ebola_outbreak_West_Africa_data_reformat.csv", quote = F, row.names = F)
ebola$Date <- as.Date(ebola$Date, format = ("%d %b %Y"))

# exploratory plots

ggplot(data = ebola, mapping = aes(x = Guinea_Cases, y = Date)) +
  geom_line() + 
  geom_point() + 
  labs(title = "Guinea", x = "Cases", y = "Timeline")

ggplot(data = ebola, mapping = aes(x = SierraLeone_Cases, y = Date)) +
  geom_line() + 
  geom_point() + 
  labs(title = "Sierra Leone", x = "Cases", y = "Timeline") + 
  theme_bw()

ggplot(data = ebola, mapping = aes(x = Liberia_Cases, y = Date)) +
  geom_line() + 
  geom_point() + 
  labs(title = "Liberia", x = "Cases", y = "Timeline")

# fit a linear regression model to the data
ebola1 <- ebola[1:16,]
ebola1 <- ebola1 %>% mutate(Timeline = Date - min(Date))
ebola1 <- ebola1 %>% mutate(logGuinea_Cases = log(Guinea_Cases))

model1 <- lm(logGuinea_Cases ~ Timeline, data = ebola1)

ggplot(data = ebola1, mapping = aes(x = Timeline, y = logGuinea_Cases)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Timeline", y = "ln(Cumulative Guinea Cases)") + 
  theme_bw()

delta <- model1$coefficients["Timeline"]
incubation = 5.3 # incubation 
infection = 5.61 # infectious
R0 <- (1 + delta * incubation) * (1 + delta * infection)
  
ggsave("Plot_Day2.png", height = 10, width = 20, units ="cm")

model# Definition of the SEIR model
SEIR <- function(t, x, parms) {
  with(as.list(c(parms,x)),{
    if(t < tau1) beta <- beta0
    else beta <- beta0*exp(-k*(t-tau1)) 
    # beta0 = per contact infection rate    
    #     k = decay induced by the control measures
    #  tau1 = time of introduction of control measures
    N <- S + E + I + R
    # N = Total, S = Susceptible, E = Exposed, I = infected, R = Recovered 
    dS <- - beta*S*I/N
    dE <- beta*S*I/N - sigma*E
    dI <- sigma*E - gamma*I
    # 1/sigma = average duration of incubation
    # 1/gamma = average duration of infectiousness 
    dR <- (1-f)*gamma*I
    #f = fatality (death) rate
    dD <- f*gamma*I
    #D = deaths
    dC <- sigma*E
    #C = Cumulative cases of infected 
    der <- c(dS,dE,dI,dR,dD,dC)
    list(der)
  })
}

# Cost function to calculate the sum of squared residuals (SSR) of predicted vs observed data
cost <- function(free,fixed,init,data) {
  pars <- c(free,fixed)
  pars <- trans(pars)
  times <- c(0,data$times+pars["tau0"])
  simulation <- as.data.frame(ode(init,times,SEIR,parms=pars))
  simulation <- simulation[-1,]
  cost <- sum((simulation$C-data$cases)^2+(simulation$D-data$deaths)^2)
  return(cost)
}

# Parameter transformation
trans <- function(pars) {
  pars["beta0"] <- exp(pars["beta0"])
  pars["k"] <- exp(pars["k"])
  pars["f"] <- plogis(pars["f"])
  pars["tau0"] <- exp(pars["tau0"])
  pars["tau1"] <- exp(pars["tau1"])
  return(pars)
}

# Fit the model to the data

#####################################################
# GUINEA: Prepare the data and set the initial values
#####################################################
data <- na.omit(ebola[c("Date","Guinea_Cases","Guinea_Death")])
names(data) <- c("times","cases","deaths")
data$times <- as.Date(data$times,format = "%d %b %Y")

begin <- as.Date("2 Dec 2013", format="%d %b %Y") 
delay <- as.numeric(data$times[1] - begin)

data$times <- as.numeric(data$times - data$times[1])
N <- 1e6		
init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 1)
fixed <- c(tau0 = log(delay), tau1 = -Inf, sigma = 1/5.3, gamma = 1/5.61)
free <- c(beta0 = log(0.2), k = log(0.001), f = 0)

# Fit the model to the data
fit <- optim(free, cost, gr=NULL, fixed, init, data, method="Nelder-Mead", hessian=TRUE)

# Plot the best-fit model to the data
# back transform the parameteres:
pars <- trans(c(fit$par,fixed))
end <- as.Date("1 Sep 2014", format="%d %b %Y") 

months <- as.Date(c("1 Dec 2013","1 Mar 2014","1 Jun 2014","1 Sep 2014"), format="%d %b %Y") 
timepoints <- seq(0,as.numeric(end-begin),1)

# here we run the model with the estmated parameters (in the "pars" vector)
simulation <- as.data.frame(ode(init,timepoints,SEIR,parms=pars))

# Hint: Ploting the predicted cases
plot(x = timepoints, y = simulation$C, type="l", col="red", ylim=c(0,700), xlab='Time', 
     ylab="Cumulative number of cases",frame=FALSE,axes=FALSE,main="Guinea")
axis(1,months-begin,months)
axis(2)

#Hint: While plotting the observed cases, do not forget to account for the delay in reporting.


simulation %>% 
  ggplot(aes(x = time, y = C, color = "Cases")) + 
  geom_line() + 
  geom_line(aes(y = D, color = "Deaths")) +
  geom_point(data = data, mapping = aes(x = times + pars["tau0"], y = cases, color = "Cases")) +
  geom_point(data = data, mapping = aes(x = times + pars["tau0"], y = deaths, color = "Deaths")) +
  labs(x="Time", y="Cumulative number of cases", color = "Legend") +
  scale_color_manual(values = colors <- c("Cases" = "blue", "Deaths" = "red")) +
  theme_classic() +
  scale_x_continuous(labels = months, breaks = as.numeric(months-begin))
