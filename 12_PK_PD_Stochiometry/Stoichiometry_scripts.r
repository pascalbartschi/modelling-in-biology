############################################################
#############                                  #############
#############       PK/PD & Stoichiometry      #############
#############                                  #############
############################################################

############# Stoichiometry scripts


# Sample a trimer ---------------------------------------------------------

sample_trimer <- function(f_M) {

  # "M": mutant Env protein
  # "W": wild-type Env protein
  proteins <- c("M", "W") # possible proteins to sample from

  # randomly sample 3 proteins with replacement and prob=f_M
  trimer <- sample(proteins,
                   size = 3,
                   replace = TRUE,
                   prob = c(f_M, 1-f_M))

    return(trimer)

}

# Sample a virion ---------------------------------------------------------

sample_virion <- function(s, f_M) {

  virion_trimers <- replicate(n = s,
                              expr = sample_trimer(f_M),
                              simplify = FALSE)
  return(virion_trimers)

}




# Stoichiometry estimator -------------------------------------------------

estimate_T <- function(data, trimer_number_sample) {

  ## Make sure that the data format is data.frame
  data <- as.data.frame(data)

  ## Fit trimer number distribution
  pdf_eta <- function(trimer_number_sample) {

    # trimer_number_sample = vector containing trimer numbers

    mu_s <- mean(trimer_number_sample)
    var_s <- 49/14*mu_s
    max_s <- 100

    #### definition of the trimer number distribution
    # correction for infinity
    eta_muv.b<-function(mu,v,s.max){
      # mu:		mean of trimer number distribution
      # v:		variance of trimer number distribution
      # s.max:	maximal number of trimers on virion surface (minimal number is 0)

      # calculation of parameters of B-distribtion

      mu<-mu/s.max
      v<-v/s.max^2

      p<-(mu^2-mu^3-mu*v)/v
      q<-(mu-2*mu^2+mu^3-v+v*mu)/v
      dos<-dbeta((0:s.max)/s.max,p,q)

      if(dos[1]==Inf){dos[1]<-ceiling(dos[2])
      }
      if(dos[s.max+1]==Inf){dos[s.max+1]<-ceiling(dos[s.max])
      }

      out<-dos/sum(dos)

    }

    distr <- eta_muv.b(mu_s, var_s, max_s)
    names(distr) <- seq(0, max_s, by=1)

    return(distr)

  }

  eta <- pdf_eta(trimer_number_sample)

  # RI function from the basic model
  entryRI.basic<-function(eta,f.M,TT){

    s.max<-length(eta)-1
    p3<-(1-f.M)^3

    if(f.M==0 | f.M==1){out<-1-f.M}

    else{
      out<-sum(eta[(TT:s.max)+1]*sapply(TT:s.max,function(s) sum(dbinom(TT:s,s,p3))))/sum(eta[(TT:s.max)+1])
    }

    out

  }

  # vectorization
  entryRI.basic.v<-function(eta,F.M,TT)sapply(F.M, entryRI.basic,eta=eta,TT=TT)

  # arcsin sqrt - transformation
  trafas<-function(x){
    if(!is.na(x)){
      if(x<0){x<-x}
      else{x<-asin(sqrt(x))}
    }
    x
  }
  arcsinsqrt<-function(x) sapply(x,trafas)

  # rss. matrix
  rss.matrix.entry.basic<-function(eta,data){
    T.max<-length(eta)-sum(eta==0)
    #	print(paste("T.max",T.max))
    out<-matrix(NA,nrow=T.max,ncol=2,dimnames=list(rep("",T.max),c("TT","rss")))
    for(TT in 1:T.max){

      #print(paste("T=",TT,sep=""))

      # no transformation
      # out[TT,]<-c(TT,sum(((data[,"RI"])-(entryRI.basic.v(eta,data[,"f.M"],TT)))^2))
      # arcsin sqrt- transformation
      out[TT,]<-c(TT,sum((arcsinsqrt(data[,"RI"])-arcsinsqrt(entryRI.basic.v(eta,data[,"f.M"],TT)))^2))
      # logit transformation
      # out[TT,]<-c(TT,sum((logit(data[,"RI"])-logit(entryRI.basic.v(eta,data[,"f.M"],TT)))^2))
    }
    out
  }

  # estimator
  T.estimator.entry.basic<-function(eta,data){
    W<-rss.matrix.entry.basic(eta,data)
    array(c(W[which(W[,"rss"]==min(W[,"rss"],na.rm=TRUE)),c("TT","rss")]),dim=c(1,2),dimnames=list("",c("TT","rss")))
  }

  T_est <- T.estimator.entry.basic(eta, data)
  return(T_est)
}
