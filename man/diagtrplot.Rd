%copyright (c) 2005 MS Blanchard
\name{diagtrplot}
\alias{diagtrplot}
\alias{diagtrplot.NONMEM}
\alias{diagtrplot.PKNLME}
\alias{diagtrplot.WinBUGS}
\title{diagtrplot}
\description{
diagtrplot creates a trellis plot of the observed concentrations and predicted values vs time by subject.
}
\usage{
diagtrplot(x, level, xvarlab, yvarlab, pages,...)
}
\arguments{
  \item{x}{variable identifying the clustering variable}
  \item{level}{level of mixed model ("p"-population, "i"-individual)}
  \item{xvarlab}{label for x variable}
  \item{yvarlab}{label for y variable}
  \item{pages}{number of pages to print, 1 prints first page}
  \item{...}{additional arguments to be passed to lower level functions}
}

\value{
diagtrplot produces a trellis plot of observed concentrations and predicted values vs time by subject.
}


\author{M.S. Blanchard<sblanchard@coh.org>}


\seealso{\code{\link{trplot}}, \code{\link{diagplot}}, \code{\link{residplot}}, \code{\link{obvsprplot}}, \code{\link{tex}}, \code{\link{HTMLtools}}}

\examples{
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
data<-list(pkvar=Theo, cov=wt.v) 

nameData<-list(covnames=c("wt"),
               yvarlab="Sqrt(Theop. Conc.) (mg/L)",
               xvarlab="Time since dose (hrs)",
               reparams=c("Cl"),
               params=c("Ka","V", "Cl"),
               tparams=c("log(Ka)","log(V)","log(CL)"))

model.def<-list(fixed.model=lKa+lV+lCl~1,random.model=lCl~1,
start.lst=c(lKa=.3,lV=-.6,lCl=-3), form=conc~sonecpmt(dose, time,
lV, lKa, lCl), control=nlmeControl(returnObject=FALSE))

MM<-RunNLME(inputStructure=model.def,data=data, nameData=nameData)

diagtrplot(x=MM,level="p", xvarlab=nameData$xvarlab, 
yvarlab=nameData$xvarlab, pages=1)
}

\keyword{ hplot }% 
