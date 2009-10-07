\name{bugs}
\alias{bugs}
\title{bugs}
\description{
Interface from R to WinBUGS by A. Gelman
}
\usage{
bugs(data, inits, parameters.to.save, model.file="model.bug", 
n.chains=3, n.iter=2000, n.burnin=floor(n.iter/2), 
n.thin=max(1, floor(n.chains * (n.iter - n.burnin)/1000)), 
debug=FALSE, attach.sims=TRUE, print.summary=TRUE, plot.summary=TRUE, 
digits.summary=1, display.parallel=FALSE, DIC=TRUE, 
bugs.directory="c:/Program Files/WinBUGS14/", 
dos.location="c:/progra~1/winbug~1/winbug~1")
}

\arguments{
  \item{data}{ a list of the data for the Winbugs model, or a vector of the 
               names of the data objects used by the model}
  \item{inits}{  a list with n.chains elements; each element of the list is
             itself a list of starting values for the Winbugs model,
             or a function creating (possibly random) initial values}
  \item{parameters.to.save}{  vector of the names of the parameters to save}
  \item{model.file}{  location of the model.  (Default is "model.txt".)}
  \item{n.chains}{  number of chains.  Must be at least 2.  (Default is 3.)}
  \item{n.iter}{ number of iterations per chain.  (Default is 2000.)}
  \item{n.burnin}{  number of iterations to discard at the beginning.
                (Default is n.burnin=n.iter/2, that is, discarding the first
                half of the simulations.)}
  \item{n.thin}{ thinning rate.  Must be a positive integer.  Set n.thin>1 to
              save memory and computation time if n.iter is large.
              (Default is n.thin=max(1,floor(n.chains*(n.iter-n.burnin)/1000))
              which will only thin if there are at least 2000 simulations.)}
  \item{debug}{   option to not automatically quit out of WinBugs when the
             script has finished running, so that you can look at what's
             going on within WinBugs.  (Default is debug=F.)}
  \item{attach.sims}{  option to save all the parameters in parameters.to.save as
                   R objects, overwriting any existing variables with these
                   names.  (Default is attach.sims=T.)}
  \item{print.summary}{ option to print summary statistics and convergence
                     information.  (Default is print.summary=T.)}
  \item{plot.summary}{option to display summary statistics and convergence
                    information as a graph.  (Default is plot.summary=T.)}
  \item{digits.summary}{  rounding for tabular output on the console.
  (Default is to round to 1 decimal place.)}
  \item{display.parallel}{option to display parallel intervals in both halves of
                        the summary plots.  This is a convergence-monitoring
                        tool and is not necessary once you have approximate
                        convergence.  (Default is display.parallel=F.)}
  \item{DIC}{  option to compute deviance, pD, and DIC.  (Default is DIC=T.)}
  \item{bugs.directory}{ bugs.directory}
  \item{dos.location}{ dos.location}
}


\value{
  Output is a list (sims.array, sims.list, sims.matrix, summary):

     n.keep:  number of iterations kept per chain
              (equal to (n.iter-n.burnin)/n.thin)

     n.sims:  number of posterior simulations (equal to n.chains*n.keep)

     sims.array:  3-way array of simulation output, with dimensions
                  n.keep, n.chains, and length of combined parameter vector

     sims.list:  list of simulated parameters:
                 for each scalar parameter, a vector of length n.sims
                 for each vector parameter, a 2-way array of simulations,
                 for each matrix parameter, a 3-way array of simulations, etc.

     sims.matrix:  matrix of simulation output, with n.chains*n.keep rows and
                   one column for each element of each saved parameter.
             (For convenience, the n.keep*n.chains simulations in
             sims.array and sims.list have been randomly permuted.)

     summary:  summary statistics and convergence information for each
               element of each saved parameter.  Approximate convergence
               has been reached when R-hat < 1.2 for all parameters.

     mean:  a list of the estimated parameter means

     sd:  a list of the estimated parameter sd's

     median:  a list of the estimated parameter medians
       (The information in "mean", "sd", and "median" is already included in
       "summary"; it is included in list form for convenience in later
       analyses.)

     pD:  var(deviance)/2, an estimate of the effective number of parameters
       (The variance is computed as the average of the within-chain variances,
       which gives a more reasonable estimate when convergence has not been
       reached.)

     DIC:  mean(deviance) + pD

     last.values:  list of simulations from the most recent iteration.  They
       can be used as starting points if you wish to run Bugs for further
       iterations

     In addition, the simulated parameter values are automatically saved as R
     objects (in the same form as the elements of sims.list).
     And the summary elements are also saved as R objects.
     (For example, if "beta" is a 10 x 3 array in the model, then it will
     be saved as an array of dimensions n.sims x 10 x 3.)

}
\references{
Gelman, A. and Carlin, J.B. and Stern, H.S. and Rubin, D.B.  (2003).
"Bayesian Data Analysis (2nd edition)". Chapman \& Hall/CRC:New York.
  }
\author{A. Gelman}

\seealso{\code{\link{RunWB}}}

\examples{
if (.Platform$OS.type =="windows"){
library(PKtools)
library(nlme)
curwd=getwd()
if (file.exists("C:/bugsR")) {
  setwd("C:/bugsR")
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

  parameters <- c("sigma2","ka","cl","v","beta","mu","re","itau","ipredwb","ppredwb")

  theo.sim <- bugs(data = data, inits = inits,
        parameters.to.save = parameters, model.file = "theosw.txt",
        n.chains = 3, n.iter = 12000, debug = TRUE,
        n.burnin = 4000 , n.thin = 8 , print.summary = FALSE,
        plot.summary = TRUE)
  print(names(theo.sim)) 
  setwd(curwd)
  }
else{
    "You do not have a C:/BugsR directory."
    }
}
}

\keyword{ models }
