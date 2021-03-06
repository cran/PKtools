%copyright (c) 2009 MS Blanchard
\name{residplot}
\alias{residplot}
\alias{residplot.PKNLME}
\alias{residplot.NONMEM} 
\alias{residplot.WinBUGS}
\title{residplot}
\description{
resid creates individual residual vs predicted plots at the population (marginal) and 
individual (conditional) levels of the mixed model the can be used with the 
method identify to identify outliers.
}
\usage{
residplot(x,level, ...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{object of class, NONMEM, PKNLME, or WinBUGS}
  \item{level}{level of mixed model ("p"-population, "i"-individual)}.
  \item{...}{additional arguments to be passed to lower level functions}
}

\details{
The method identify can be used with objects of class NONMEM, PKNLME, and WinBUGS 
by including the following code.
\itemize{
\item NONMEM:
\itemize{
\item population level: identify(NM\$pred\$PRED, NM\$pred\$WRES)
\item individual level: identify(NM\$pred\$IPRE, NM\$pred\$IWRE)
}
\item PKNLME:
\itemize{
\item population level: identify(MM\$mm\$fitted[,1], MM\$mm\$RES)
\item individual level: identify(MM\$mm\$fitted[,2], MM\$mm\$IRES)
}
\item WinBUGS:
\itemize{
\item population level: identify(WB\$pred\$ppred, WB\$pred\$presid)
\item individual level: identify(WB\$pred\$ipred, WB\$pred\$iresid)
}
}
}


\value{
  plots of residual versus predicted values for both the 
  population (marginal) and individual (conditional) 
  levels.
}


\author{M.S. Blanchard <sblanchard.coh.org>}


\seealso{\code{\link{identify}}}

\examples{

if (.Platform$OS.type =="windows"){
library(PKtools)
library(nlme) 
curwd=getwd()
if (file.exists("C:/nmv/run")){
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
               varnames=c("D[1,1]","D[1,2]","D[2,2]","D[1,3]","D[2,3]","D[3,3]"))

  data<-list(pkvar=Theo, cov=wt.v)
  NM<-RunNM(inputStructure="control.model5", data=data, nameData=nameData)
  residplot(NM, level="p")
  setwd(curwd)
  }
  else{
      "You do not have NONMEM."
      }
}
}
\keyword{ hplot }% 
