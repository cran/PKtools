%copyright (c) 2009 MS Blanchard
\name{PKtools.AIC}
\alias{PKtools.AIC}
\title{PKtools.AIC}
\description{ PKtools.AIC calculates the AIC and AICc. }
\usage{PKtools.AIC(loglike,n,K,...)}
\arguments{
  \item{loglike}{loglikelihood}
  \item{n}{total number of samples}
  \item{K}{number of fixed parameters including both mean and variance parameters}
  \item{...}{additional arguments to be passed to lower level functions}
}

\value{
This function outputs the AIC and and the small sample AIC, AICc, 
as well as the objective function (-2 \eqn{\times}{x} loglikelihood) and K.
}

\references{ Burnham, K.P. and Anderson,D.R., (2002). Model
Selection and Multimodel Inference: A Practical Information -
Theoretic Approach (2nd edition). Springer: New York. }

\author{M.S. Blanchard <sblanchard.coh.org>}

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

nameData<-list(covnames=c("wt"),
               yvarlab="Sqrt(Theop. Conc.) (mg/L)",
               xvarlab="Time since dose (hrs)",
               reparams=c("V","Cl"),
               params=c("Ka","V", "Cl"),
               tparams=c("log(Ka)","log(V)","log(CL)"))


model.def<-list(fixed.model=lKa+lV+lCl~1,random.model=lV+lCl~1,
start.lst=c(lKa=.3,lV=-.6,lCl=-3), form=conc~sonecpmt(dose, time,
lV, lKa, lCl), control=nlmeControl(returnObject=FALSE))

MM<-RunNLME(inputStructure=model.def, data=data, nameData=nameData)

K = attr(logLik(MM$mm), "df")
n<-nrow(MM$pkdata) 
AIC.table<-data.frame(PKtools.AIC(loglike=logLik(MM$mm),n=n,K=K), row.names="")
AIC.table 

}

\keyword{ models }%
