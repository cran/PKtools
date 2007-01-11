%copyright (c) 2005 MS Blanchard
\name{RunNLME}
\alias{RunNLME}
\title{RunNLME}
\description{
RunNLME uses the NLME software to estimate parameters for
a single dose population PK model with hierarchical data.
}
\usage{
RunNLME(inputStructure,data, nameData)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{inputStructure}{NLME-model.def}
  \item{data}{list of data files including pk data and covariate data the length of the full dataset}
  \item{nameData}{list of names, including, covnames, yvarlab, xvarlab, paramester names}
}
\details{
model.def is a list of the definitions of the model form, fixed and random effects, 
the starting values and control argument from the nlme function. The following is 
an example.
\itemize{
\item model.def<-list(
fixed.model=lKa+lVol+lCl \eqn{\mbox{\textasciitilde}}{~}1,
random.model=lVol+lCl \eqn{\mbox{\textasciitilde}}{~}1, 
start.lst=c(lKa=.3,lVol=-.6,lCl=-3), 
form=conc \eqn{\mbox{\textasciitilde}}{~} sonecpmt(dose, time, lV, lKa, lCl),
control=nlmeControl(returnObject=FALSE)).
}

nameData is a list of the labels including the names of the covariates in the order
they are given in the covariate dataset, y and x variable, the random parameters
(reparams -should match the list for random.model in the model.def), fixed parameters (params -should
match the list for fixed.model in the model.def), label for transformed parameters (
in the Theo example the model parameters are on a log scale tparam=c("log(Ka)","log(V)","log(Cl)") 
and the names of the variance parameters are not required for NLMEoutput.
}

\value{
Output datasets include the input data, the parameter estimates, covariates, model
residuals at the population and individual levels,and model predicted values for the population
and individual levels.
}

\references{
Pinheiro, J.C. and Bates, D.M. (2000). "Mixed-Effects Models in S and SPLUS."
Springer: New York.
}

\author{M.S. Blanchard <sblanchard@coh.org>}


\seealso{ \code{\link{pk}}, \code{\link{coVar.id}}, \code{\link{RunNM}}, \code{\link{RunWB}}}

\examples{
#NLME example
library(PKtools)
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

nameData<-list(covnames=c("wt"),
               yvarlab="Sqrt(Theop. Conc.) (mg/L)",
               xvarlab="Time since dose (hrs)",
               reparams=c("Cl"),
               params=c("Ka","V", "Cl"),
               tparams=c("log(Ka)","log(V)","log(CL)"))

model.def<-list(fixed.model=lKa+lV+lCl~1,random.model=lCl~1,
start.lst=c(lKa=.3,lV=-.6,lCl=-3), form=conc~sonecpmt(dose, time,
lV, lKa, lCl), control=nlmeControl(returnObject=FALSE))

MM<-RunNLME(inputStructure=model.def, data=data, nameData=nameData)
MM

}

\keyword{ models }% 
