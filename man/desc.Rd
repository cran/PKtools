%copyright (c) 2005 MS Blanchard
\name{desc}
\alias{desc}
\title{desc}
\description{
desc calculates select descriptive statistics for the variable X.
}
\usage{
desc(y, pcts, nsig)
}

\arguments{
  \item{y}{variable of interest}
  \item{pcts}{percentiles of interest, the default is c(0.025, 0.05, 0.95, 0.975)}
  \item{nsig}{number of significant figures, the default is 4}
}

\value{desc prints descriptive statistics including mean, median, standard deviation,
range, min, max, and select percentiles.
}


\author{M.S. Blanchard <sblanchard@coh.org>}

\seealso{\code{\link{tex}}}

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
desc(Theo$conc)
}

\keyword{ univar }% 
