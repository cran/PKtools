%copyright (c) 2005 MS Blanchard
\name{pk}
\alias{pk}
\title{pk}
\description{
pk creates two data sets, a rectangular data set for R, and a NONMEM
ready data set. 
}
\usage{ pk(pkvar, covdata, covnames)
}

\arguments{
  \item{pkvar}{PK data set including; id, dose, conc, and time}
  \item{covdata}{matrix/vector of covariate data}
  \item{covnames}{vector of names of covarites in the cov matrix/vector}
}

\value{
pk creates a pk data file pkdat including: id, dose, time, conc, plus the covaraites, 
and also creates NMdata, a NONMEM ready data file.
}


\author{M.S. Blanchard <sblanchard@coh.org>}


\seealso{\code{\link{RunNM}}}

\examples{
library(PKtools)
library(nlme)
data(Theoph)
Theoph<-Theoph[Theoph$Time!=0,]
id<-as.numeric(as.character(Theoph$Subject))
dose<-Theoph$Dose
time<-Theoph$Time
concblq<-round(sqrt(Theoph$conc),4)
conc<-concblq
Theo<-data.frame(cbind(id,dose,time,conc))
names(Theo)<-c("id","dose","time","conc")
wt.v<-Theoph$Wt

data<-list(pkvar=Theo, cov=wt.v)

nameData<-list(covnames=c("wt"))

pk(pkvar=data$pkvar, covdata=data$cov, covnames=nameData$covnames)
}

\keyword{ models }% 
