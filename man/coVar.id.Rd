%copyright (c) 2009 MS Blanchard
\name{coVar.id}
\alias{coVar.id}
\title{coVar.id}
\description{
coVar.id creates a data set of the covariates one line per id with id as the first column.
}
\usage{
coVar.id(id, coVar, nameData) 
}

\arguments{
  \item{id}{cluster id}
  \item{coVar}{data set of the covariates with length equal to the full data set}
  \item{nameData}{list of names, including, covnames}
}

\value{coVar.id outputs a data set of the covariates one line per id with id as the first column.}

\author{M.S. Blanchard <sblanchard.coh.org>}
\seealso{RunNM, RunNLME, RunWB}
\examples{
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

nameData<-list(covnames=c("wt"))

descStructure<-list(pcts=c(.025,.05,.95,.975),nsig=4)


cov.id <- coVar.id(data$pkvar$id, data$cov, nameData)
cov.id
}
\keyword{ models }% 
