%copyright (c) 2005 MS Blanchard
\name{obvsprplot}
\alias{obvsprplot}
\alias{obvsprplot.NONMEM}
\alias{obvsprplot.PKNLME}
\alias{obvsprplot.WinBUGS}
\title{obvsprplot}
\description{
obvsprplot creates individual observed vs predicted plots at the population (marginal) and 
individual (conditional) levels of the mixed model the can be used with the 
method identify to identify outliers.
}
\usage{
obvsprplot(x,...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{object of class, NONMEM, PKNLME, WinBUGS}
  \item{...}{additional arguments to be passed to lower level functions}
}
\details{
The method identify can be used with objects of class NONMEM, PKNLME, and WinBUGS 
by including the following code.
\itemize{
\item NONMEM:
\itemize{
\item population level: identify(NM\$pred\$PRED, NM\$pred\$CONC)
\item individual level: identify(NM\$pred\$IPRE, NM\$pred\$CONC)
}
\item PKNLME:
\itemize{
\item population level: identify(MM\$mm\$fitted[,1], MM\$pkdata\$conc)
\item individual level: identify(MM\$mm\$fitted[,2], MM\$pkdata\$conc)
}
\item WinBUGS:
\itemize{
\item population level: identify(WB\$pred\$ppred, WB\$pred\$conc)
\item individual level: identify(WB\$pred\$ipred, WB\$pred\$conc)
}
}
}
\value{
  plots of observed versus predicted values for both the 
  population (marginal) and individual (conditional) 
  levels.
}


\author{M.S. Blanchard <sblanchard@coh.org>}

\seealso{ \code{\link{identify}}, \code{\link{obvsprplot}}, \code{\link{diagplot}}}

\examples{
if (.Platform$OS.type =="windows"){
library(PKtools)
library(nlme) 

curwd=getwd()
if (file.exists("C:/nmv/run")) {
  setwd("C:/nmv/run")
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
  obvsprplot(NM, "p")
  setwd(curwd)
  }
  else{
      "You do not have NONMEM."
      }
}
}
\keyword{ hplot }% 
