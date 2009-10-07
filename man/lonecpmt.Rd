%copyright (c) 2009 MS Blanchard
\name{lonecpmt}
\alias{lonecpmt}
\title{lonecpmt}
\description{
A function for a one compartment model with first order absorption with a natural log transformation.
}
\usage{lonecpmt(dose,time,lVol,lKa, lCl)
}

\arguments{
  \item{dose}{numeric representation of the dose}
  \item{time}{a vector of the time measurements}
  \item{lVol}{the numeric parameter log of volume of distribution}
  \item{lKa}{the numeric parameter log of absorption}
  \item{lCl}{the numeric parameter log of clearance}
}

\value{
The value of the expression
  log(dose * exp(lKa) * (exp(-(exp(lCl)/exp(lVol)) * time) - 
        exp(-exp(lKa) * time))/(exp(lVol) * (exp(lKa) - (exp(lCl)/exp(lVol)))))

}


\author{M.S. Blanchard<sblanchard.coh.org}


\seealso{\code{\link{RunNM}}}

\examples{
library(nlme) 
data(Theoph) 
Theoph=Theoph[Theoph$Subject==2 & Theoph$Time>0,]
dose=Theoph$Dose
time=Theoph$Time
exp.lconc=lonecpmt(dose,time,lVol=-.6, lKa=.3, lCl=-3)
exp.lconc
}

\keyword{ models }% 
