%copyright (c) 2009 MS Blanchard
\name{trplot}
\alias{trplot} 
\alias{trplot.PKNLME}
\alias{trplot.NONMEM} 
\alias{trplot.WinBUGS}
\title{trplot}
\description{
trplot creates a trellis plot of concentration vs time by subject.
}
\usage{
trplot(x,xvarlab,yvarlab,pages,...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{object of class, NONMEM, PKNLME, WinBUGS}
  \item{xvarlab}{x variable label}
  \item{yvarlab}{y variable label}  
  \item{pages}{number of pages you want to print, pages=1 prints the first page}
  \item{...}{additional arguments to be passed to lower level functions}
}

\value{
  A trellis plot of concentration vs time by subject.
}


\author{M.S. Blanchard <sblanchard@coh.org>}
\examples{
#NLME example
library(nlme) 
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

nameData<-list(covnames=c("wt"),               xvarlab="Time since dose (hrs)",
               yvarlab="Sqrt(Theop. Conc.) (mg/L)",
               reparams=c("Cl"),
               params=c("Ka","V", "Cl"),
               tparams=c("log(Ka)","log(V)","log(CL)"))

model.def<-list(fixed.model=lKa+lV+lCl~1,random.model=lCl~1,
start.lst=c(lKa=.3,lV=-.6,lCl=-3), form=conc~sonecpmt(dose, time,
lV, lKa, lCl), control=nlmeControl(returnObject=FALSE))

MM<-RunNLME(inputStructure=model.def, data=data, nameData=nameData)
trplot(x=MM,xvarlab=nameData$xvarlab,yvarlab=nameData$yvarlab,pages=1)
}

\keyword{ models }% 
