%copyright (c) 2004 MS Blanchard
\name{RunNM}
\alias{RunNM}
\title{RunNM}
\description{ RunNM runs the function pk to create
pharmacokinetics datasets for R and NONMEM, runs the system
commend to run NONMEM, and reads the NONMEM datasets. }
\usage{ RunNM(inputStructure, data, nameData) }
\arguments{
  \item{inputStructure}{the standard NONMEM control file}
  \item{data}{list of data files including pk data and covariate data with length of the full dataset}
  \item{nameData}{list of names, including, covnames, yvarlab, xvarlab, params}
}
\details{ 
nameData is a list of the labels including the names of the covariates in the order
there are given in the covariate dataset, y and x variable, the random parameters
(reparams -should match the list random effects defined in the control file), 
fixed parameters (params -should match the list for fixed effects in the control file), 
label for transformed parameters ( in the Theo example the model parameters are on a log scale 
tparam=c("log(Ka)","log(V)","log(Cl)") and the names of the variance parameters should
list the parameters for the upper triangle of variance covariance table.
}

\value{
The output from NMoutput are data tables of the results, including the  
objective function (ob), population parameters
(params), random effects (re), individual parameters (ip),
covariates (cov), predicted values (pred). If the objects of class NONMEM is called
NM, then the objective function can be accessed by typing NM\$ob, similarly the population parameters
can be accessed by  typing NM\$param.} 

\references{
Boeckmann, A.J.  and Sheiner, L.B. and Beal, S.L. (1994). "NONMEM
Users Guide- Part V, Introductory Guide". NONMEM Project
Group:UCSF. }

\author{M.S. Blanchard <sblanchard@coh.org>}


\seealso{\code{\link{pk}}, \code{\link{coVar.id}}, \code{\link{RunNLME}},\code{\link{RunWB}}}

\examples{ 
#NONMEM example
if (.Platform$OS.type =="windows"){
setwd("C:/nmv/run")
library(PKtools)
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

nameData<-list(covnames=c("wt"), 
               yvarlab="Sqrt(Theop. Conc.) Sqrt(mg/L)",
               xvarlab="Time since dose (hrs)",
               reparams=c("Ka", "V", "Cl"),
               params=c("Ka", "V", "Cl"),
               tparams=c("log(Ka)", "log(V)", "log(Cl)"),
               varnames=c("D[1,1]","D[1,2]","D[2,2]","D[1,3]","D[2,3]","D[3,3]")
               )

NM<-RunNM(inputStructure="control.model5", data=data, nameData=nameData) 
NM
}
}
\keyword{ models }%
