rm(list=ls())
library(deSolve)

TI.model <- function(t, x, parms) {
  Ta <- x[1] # pop 1
  I  <- x[2] # pop 2
  with(as.list(parms),{
    dT <- lambda- deltaT * Ta - beta * Ta * I
    dI <- beta * Ta * I - deltaI * I 
    res <- c(dT, dI)
    list(res)
  })
}

### by varying beta 
parmsTI  <- c(lambda = 1000, 
               deltaT = 0.1,
               deltaI = 0.1,
               beta = NA)
parmsTI["beta"] <- 5/(parmsTI["lambda"]/parmsTI["deltaT"]/parmsTI["deltaI"])
xstartTI <- c(Ta = as.numeric(parmsTI["lambda"]/parmsTI["deltaT"]),
               I = 1)
times  <- seq(0, 100, length=500)

#here rk4 is a function provided by R to solve differential equations which you do not have to understand in further detail
out.TI <- as.data.frame(rk4(xstartTI, times, TI.model, parmsTI))


library(ggplot2)

ggplot(data=out.TI) +
  geom_line(aes(x=time, y=I, color="I")) +
  geom_line(aes(x=time,y=Ta, color="Ta"))
