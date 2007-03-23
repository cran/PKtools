%copyright (c) 2005 MS Blanchard
\name{diagplot}
\alias{diagplot}
\alias{diagplot.NONMEM}
\alias{diagplot.PKNLME}
\alias{diagplot.WinBUGS}
\title{diagplot}
\description{
diagplot creates plots of observed versus predicted values and residuals (ordinary and standardized) 
  versus predicted values for both the population (marginal) and individual (conditional) 
  levels.
}
\usage{
diagplot(x,...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{object of class, NONMEM, PKNLME, WinBUGS}
  \item{...}{additional arguments to be passed to lower level functions}
}

\value{
  Plots of observed versus predicted values and residuals (ordinary and standardized) 
  versus predicted values for both the population (marginal) and individual (conditional) 
  levels.
}


\author{M.S. Blanchard <sblanchard@coh.org>}

\seealso{\code{\link{residplot}}, \code{\link{obvsprplot}}, \code{\link{identify}}}

\examples{
if (.Platform$OS.type =="windows"){
curwd=getwd()
if (file.exists("C:/nmv/run")) {
  setwd("C:/nmv/run")
  library(nlme) 
  library(PKtools)
  data(Theoph) 
  Theoph<-Theoph[Theoph$Time!=0,]
  id<-as.numeric(as.character(Theoph$Subject)) 
  dose<-Theoph$Dose
  time<-Theoph$Time 
  conc<-round(sqrt(Theoph$conc),4)
  Theo<-data.frame(cbind(id,dose,time,conc))
  names(Theo)<-c("id","dose","time","conc") 
  wt.v<-Theoph$Wt

  nameData<-list(covnames=c("wt"),
               yvarlab="Sqrt(Theop. Conc.) Sqrt(mg/L)",
               xvarlab="Time since dose (hrs)",
               reparams=c("Ka", "V", "Cl"),
               params=c("Ka", "V", "Cl"),
               tparams=c("log(Ka)", "log(V)", "log(Cl)"),
               varnames=c("D[1,1]","D[1,2]","D[2,2]","D[1,3]","D[2,3]","D[3,3]")
               )

  data<-list(pkvar=Theo, cov=wt.v)
  NM<-RunNM(inputStructure="control.model5", data=data, nameData=nameData)
  diagplot(NM)
  setwd(curwd)
  }
  else{
    "You do not have NONMEM."
    } 
}
}
\keyword{ hplot }% 
