library(deSolve)
library(ggplot2)


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

TIV.model <- function(t, x, parms) {
  Ta <- x[1] # pop 1
  I  <- x[2] # pop 2
  V  <- x[3]
  with(as.list(parms),{
    dT <- lambda - deltaT * Ta - beta * Ta * V
    dI <- beta * Ta * V - deltaI * I
    dV <- p * I - deltaV * V
    res <- c(dT, dI, dV)
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

parmsTIV  <- c(lambda = 1000, 
              deltaT = 0.1,
              deltaI = 0.5,
              deltaV = 0.1,
              beta = NA, 
              p = 1000)

parmsTIV["beta"] <- 5/(parmsTIV["lambda"]/parmsTIV["deltaT"]/parmsTI["deltaI"]/parmsTI["p"]/parmsTI["deltaV"])
                                         
xstartTIV <- c(Ta = as.numeric(parmsTI["lambda"]/parmsTI["deltaT"]),
              I = 1, V = 1)

times  <- seq(0, 100, length=500)

out.TI <- as.data.frame(rk4(xstartTI, times, TI.model, parmsTI))
out.TIV <- as.data.frame(rk4(xstartTIV, times, TIV.model, parmsTIV))

ggplot(data=out.TI) +
  geom_line(aes(x=time, y=I, color="I")) +
  geom_line(aes(x=time,y=Ta, color="Ta"))

ggplot(data=out.TIV) +
  geom_line(aes(x=time, y=I, color="I")) +
  geom_line(aes(x=time,y=Ta, color="Ta")) + 
  geom_line(aes(x=time,y=V, color="V"))





