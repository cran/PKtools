%copyright (c) 2005 MS Blanchard
\name{AICcomp}
\alias{AICcomp}
\title{AICcomp}
\description{
AICcomp  calculates and or prints the AIC, AICc (small
sample AIC) and the loglikelihood from NONMEM and NLME for each
of any number of models
}
\usage{
AICcomp(PKNLMEobjects, NONMEMobjects)
}

\arguments{
  \item{PKNLMEobjects}{a list of PKNLME objects}
  \item{NONMEMobjects}{a list of NONMEM objects}
}
\details{
The lists of PKNLME objects and NONMEM objects must be in the same order and
must be of the same length.
}

\value{data frame of the AIC, AICc (small sample AIC), the loglikelihood and the 
K, number of population parameters including both means and variance parameters.
}

\references{ Burnham, K.P. and Anderson,D.R., (2002). Model
Selection and Multimodel Inference: A Practical Information -
Theoretic Approach (2nd edition). Springer: New York.}

\author{M.S. Blanchard <sblanchard@coh.org>}
\note{ }
\seealso{AIC}
\examples{
if (.Platform$OS.type =="windows"){
  library(PKtools)
  library(nlme)  
curwd=getwd()
if (file.exists("C:/nmv/run"))  {
  setwd("C:/nmv/run")
  #NLME code models 3 and 6
  #data definition for NLME and NONMEM
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

  #model 3
  nameData<-list(covnames=c("wt"), 
               yvarlab="Sqrt(Theop. Conc.) (mg/L)", 
               xvarlab="Time since dose (hrs)",
               reparams=c("Ka","Cl"),
               params=c("Ka","V int", "V slope", "Cl"),
               tparams=c("log(Ka)","log(V) int"," log(V) slope", "log(CL)"))
 
  model.def<-list(fixed.model=list(lKa~1,lV~wt,lCl~1), random.model=lKa+lCl~1,
  start.lst=c(.3,-.6,0,-3),form=conc~sonecpmt(dose, time, lV, lKa, lCl),
  control=nlmeControl(returnObject=FALSE))
  results.nlme3<-RunNLME(inputStructure=model.def, data=data, nameData=nameData)


  #model 6
  nameData<-list(covnames=c("wt"), 
               yvarlab="Sqrt(Theop. Conc.) (mg/L)", 
               xvarlab="Time since dose (hrs)",
               reparams=c("Ka","V","Cl"),
               params=c("Ka","V", "Cl"),
               tparams=c("log(Ka)","log(V)","log(CL)"))

  model.def<-list(fixed.model=c(lKa+lV+lCl~1),random.model=pdDiag(form=lKa+lV+lCl~1),
  start.lst=c(.5,-.6,-3), form=conc~sonecpmt(dose, time, lV, lKa, lCl),
  control=nlmeControl(returnObject=FALSE))
  results.nlme6<-RunNLME(inputStructure=model.def, data=data, nameData=nameData)


 
  #NONMEM code models 3 and 6 
  #note control files must be placed in the C:/nmv/run directory


  #model 3
  nameData<-list(covnames=c("wt"), 
               yvarlab="Sqrt(Theop. Conc.) Sqrt(mg/L)",
               xvarlab="Time since dose (hrs)",
               reparams=c("Ka", "Cl"),
               params=c("Ka", "V", "Cl", "V slope"),
               tparams=c("log(Ka)", "log(V)", "log(Cl)","log(V slope)"),
               varnames=c("D[1,1]","D[1,2]","D[2,2]")   
               )

  results3<-RunNM(inputStructure="control.model3", data=data, nameData=nameData)


  #model 6
  nameData<-list(covnames=c("wt"), 
               yvarlab="Sqrt(Theop. Conc.) Sqrt(mg/L)",
               xvarlab="Time since dose (hrs)",
               reparams=c("Ka", "V", "Cl"),
               params=c("Ka", "V", "Cl"),
               tparams=c("log(Ka)", "log(V)", "log(CL)"),
               varnames=c("D[1,1]","D[1,2]","D[2,2]","D[1,3]","D[2,3]","D[3,3]")
               )

  results6<-RunNM(inputStructure="control.model6", data=data, nameData=nameData)

  #Multimodel Code
  PKNLMEobjects=list(results.nlme3,results.nlme6)
  NONMEMobjects=list(results3,results6)
  print(AICcomp(PKNLMEobjects=PKNLMEobjects, NONMEMobjects=NONMEMobjects))
  setwd(curwd)
  }
  else{
    "You do not have NONMEM."
  }
}
}

\keyword{ models }% 
