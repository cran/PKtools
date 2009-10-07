%copyright (c) 2009 MS Blanchard
\name{RunWB}
\alias{RunWB}
\title{RunWB}
\description{RunWB uses WinBUGS to estimate parameters for a
single dose population PK model with hierarchical data. }
\usage{
RunWB(inputStructure, data, nameData, WBargs) 
}

\arguments{
  \item{inputStructure}{the inputStructure for WBoutput is .txt file that
  defines the model using standard WinBUGS code}
  \item{data}{list of data files following Gelman (2003)}
  \item{nameData}{list of names, including, covnames, yvarlab, xvarlab, coef, 
     params, tparams,reparams, varparams}
  \item{WBargs}{Additional arguements required to run WinBUGS, including
  parameters- the list of parameters to be sampled, initial values (inits),
  n.chain- number of chains, n.iter-number of iterations, n.burnin- length of the
  burnin, n.thin, and debug-T/F if T the run will stop at the end of the WinBUGS
  run to allow use of WinBUGS to study mixing and convergence.}
}

\details{ 
nameData is a list of the labels including the names of the covariates in the order
they are given in the covariate dataset, y and x variable, coef are the model coefficient 
names, params are PK parameter names including fixed PK parameters, reparams 
are the parameters in the params list that are "not" fixed. tparams are the labels for 
transformed parameters (in the Theo example the model parameters are on a 
log scale tparam=c("log(Ka)","log(V)","log(Cl)") and finally varparams are 
the names of the parameters in the full covariance matrix.
}

\value{
The output from this function is an WinBUGS object that includes the mean and sims.list values
as described by Gelman, and the input data set, nameData, model predictions, and a covariate
data set by id.
}

\references{ Lunn, D.J. and Best, N. and Thomas, A. and
Wakefield, J. and Spiegelhalter, D. (2002). "Bayesian analysis of
population PK/PD Models: General concepts and software. Journal of
Pharmacokinetics and Pharmacodynamics", 29 (3), 271-307.

Lunn, D.J. and Wakefield, J. and Thomas, A. and Best, N. and
Spiegelhalter,D. (1999). PKBugs User Guide (version 1.1). Imperial
College: London.

Spiegelhalter, D. and Thomas, A. and Best, N. (2001). "Winbugs
Version 1.4 User Manual.", Imperial College School of
Medicine:London. 
}


\author{M.S. Blanchard <sblanchard@coh.org>}

\seealso{\code{\link{bugs}}, \code{\link{RunNLME}}, \code{\link{RunNM}} }

\examples{ 
if (.Platform$OS.type =="windows"){

curwd=getwd()
if (file.exists("C:/bugsR")) {
setwd("C:/bugsR")
library(PKtools)
library(nlme)
data(Theoph)
Theoph<-Theoph[Theoph$Time!=0,]
id<-as.numeric(as.character(Theoph$Subject))
dose<-Theoph$Dose
time<-Theoph$Time
conc<-round(sqrt(Theoph$conc),4)
sid<-split(id,id)
hist<-sapply(sid,length)
n.ind<-12
off.data<-matrix(NA,n.ind+1,1)
off.data[1,1]<-1
for (i in 2:(n.ind+1)) off.data[i,1]<-off.data[i-1,1]+ hist[i-1]
off.data<-c(off.data)
mean <- c(.5, -.6, -3)
R<-structure(.Data=diag(rep(.1,3)))
prec<-structure(.Data=diag(rep(.000001,3)))
data<-list(n.ind=n.ind,off.data=off.data,dose=dose,conc=conc,
time=time,mean=mean,R=R,prec=prec)

inits<- function(){
    list(beta = structure(
        .Data = c(rep(.5,12),rep(-.6,12),rep(-3,12)),
        .Dim = c(12, 3)), 
        mu = c(.5, -.6,  -3), 
        tau = structure(.Data = c(0.1, 0, 0,
                                  0, 0.1, 0,
                                  0, 0, 0.1), .Dim = c(3, 3)),
        tauC = 20)

list(beta = structure(
        .Data = c(rep(-.5,12),rep(-.8,12),rep(-3.5,12)),
        .Dim = c(12, 3)), 
        mu = c(-.5, -.8,  -3.5), 
        tau = structure(.Data = c(0.1, 0, 0,
                                  0, 0.1, 0,
                                  0, 0, 0.1), .Dim = c(3, 3)),
        tauC = 20)

list(beta = structure(
        .Data = c(rep(1.5,12),rep(-.4,12),rep(-2.8,12)),
        .Dim = c(12, 3)), 
        mu = c(1.5, -.4,  -2.8), 
        tau = structure(.Data = c(0.1, 0, 0,
                                  0, 0.1, 0,
                                  0, 0, 0.1), .Dim = c(3, 3)),
        tauC = 20)

}

#covariates
wt.v<-Theoph$Wt

parameters <- c("sigma2","ka","cl","v","beta","mu","re","itau","ipredwb","ppredwb")

nameData<-list(covnames=c("wt"), 
               yvarlab="Sqrt(Theop. Conc.) Sqrt(mg/L)",
               xvarlab="Time since dose (hrs)",
               coef=c("Ka", "V", "Cl"),
               params=c("Ka", "V", "Cl"),
               reparams=c("Ka", "V", "Cl"),
               tparams=c("log(Ka)","log(V)","log(CL)"),
               varnames=c("D[1,1]","D[1,2]","D[1,3]",
               "D[2,1]","D[2,2]","D[2,3]",
               "D[3,1]","D[3,2]","D[3,3]")
               )
                                             
data<-list(data=data, cov=wt.v, id=id)

WBargs<-list(parameters=parameters, inits=inits, n.chains=3, 
             n.iter=12000, n.burnin=4000, n.thin=3, debug=TRUE)
      
WB<-RunWB(inputStructure="theosw.txt", data=data, nameData=nameData, WBargs=WBargs) 
print(WB)
setwd(curwd)
}
else{
    "You do not have the C:/BugsR directory."
    }
}
}

\keyword{ models }%
