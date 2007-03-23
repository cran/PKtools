%copyright (c) 2005 MS Blanchard
\name{indEst}
\alias{indEst}
\title{indEst}
\description{
indEst outputs the individual level parameter estimates from NONMEM, PKNLME and WinBUGS.
}
\usage{
indEst(PKNLMEobject, NMobject, WBobject, outputType) 
}

\arguments{
  \item{PKNLMEobject}{PKNLME object from RunNLME}
  \item{NMobject}{NONMEM object from RunNM}
  \item{WBobject}{WinBUGS object from RunWB}
  \item{outputType}{"tex" or "R" outputs are available}
}
\details{The PKNLME, NM and WB objects should all be from the same model}
}

\value{The output is a dataframe of the indiviudal parameter estimates.
}

\references{
Boeckmann, A.J.  and Sheiner, L.B. and Beal, S.L. (1994). "NONMEM
Users Guide- Part V, Introductory Guide". NONMEM Project
Group:UCSF.

Pinheiro, J.C. and Bates, D.M. (2000). "Mixed-Effects Models in S and SPLUS."
Springer: New York.

Spiegelhalter, D. and Thomas, A. and Best, N. and Lunn, D. (2001). "Winbugs
Version 1.4 User Manual.", Imperial College School of
Medicine:London. 
}

\author{M.S. Blanchard <sblanchard@coh.org>}
\seealso{RunNM, RunNLME, RunWB}
\examples{

if (.Platform$OS.type =="windows"){
library(PKtools)
library(nlme)
out<-0
curwd=getwd()
if (file.exists("C:/bugsR")) {
  setwd("C:/bugsR")
  data(Theoph) 
  Theoph<-Theoph[Theoph$Time!=0,]
  id<-as.numeric(as.character(Theoph$Subject)) 
  dose<-Theoph$Dose
  time<-Theoph$Time 
  conc<-round(sqrt(Theoph$conc),4)
  sid<-split(id,id) 
  hist<-sapply(sid,length) 
  n.ind<-12
  off.data<-matrix(NA,n.ind+1,1)
  off.data[1,1]<-1
  for (i in 2:(n.ind+1)) off.data[i,1]<-off.data[i-1,1]+ hist[i-1]
  off.data<-c(off.data)
  mean <- c(.5, -.6, -3)
  R<-structure(.Data=diag(rep(.01,3)))
  prec<-structure(.Data=diag(rep(.000001,3)))
  data<-list(n.ind=n.ind,off.data=off.data,dose=dose,conc=conc,
  time=time,mean=mean,R=R,prec=prec)

  inits<- function(){
    list(beta = structure(
        .Data = c(rep(.5,12),rep(-.6,12),rep(-3,12)),
        .Dim = c(12, 3)),
        mu = c(.5, -.6,  -3),
        tau = structure(.Data = c(0.1, 0, 0,
                                  0, 0.1, 0,
                                  0, 0, 0.1), .Dim = c(3, 3)),
        tauC = 20)

    list(beta = structure(
        .Data = c(rep(-.5,12),rep(-.8,12),rep(-3.5,12)),
        .Dim = c(12, 3)),
        mu = c(-.5, -.8,  -3.5),
        tau = structure(.Data = c(0.1, 0, 0,
                                  0, 0.1, 0,
                                  0, 0, 0.1), .Dim = c(3, 3)),
        tauC = 20)

    list(beta = structure(
        .Data = c(rep(1.5,12),rep(-.4,12),rep(-2.8,12)),
        .Dim = c(12, 3)),
        mu = c(1.5, -.4,  -2.8),
        tau = structure(.Data = c(0.1, 0, 0,
                                  0, 0.1, 0,
                                  0, 0, 0.1), .Dim = c(3, 3)),
        tauC = 20)

  }

  #covariates 
  wt.v<-Theoph$Wt

  parameters <-
  c("sigma2","ka","cl","v","beta","mu","itau","ipredwb","ppredwb")

  nameData<-list(covnames=c("wt"),
               yvarlab="Sqrt(Theop. Conc.) Sqrt(mg/L)",
               xvarlab="Time since dose (hrs)",
               params=c("Ka", "V", "Cl"),
               tparams=c("log(Ka)","log(V)","log(CL)"),
               varnames=c("D[1,1]","D[1,2]","D[1,3]",
               "D[2,1]","D[2,2]","D[2,3]",
               "D[3,1]","D[3,2]","D[3,3]")
               )


  data<-list(data=data, cov=wt.v, id=id)

  WBargs<-list(parameters=parameters, inits=inits, n.chains=3,
             n.iter=12000, n.burnin=4000, n.thin=3, debug=T)

  WB2<-RunWB(inputStructure="theosw.txt", data=data, nameData=nameData, WBargs=WBargs)
  setwd(curwd)
  }
else {
   print("You do not have C:/BugsR directory.")
   out<-1
   }
#NLME code model 5
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
               reparams=c("Ka","V","Cl"),
               params=c("Ka","V", "Cl"),
               tparams=c("log(Ka)","log(V)","log(CL)"))

#mat<-matrix(c(.5, 0, 0, 0,.03, 0, 0,0,.08),nrow=3)
model.def<-list(fixed.model=c(lKa+lV+lCl~1),random.model=lKa+lV+lCl~1,
start.lst=c(.5,-.6,-3), form=conc~sonecpmt(dose, time, lV, lKa, lCl),
control=nlmeControl(returnObject=TRUE, opt=c("nlm")))
results.nlme5<-RunNLME(inputStructure=model.def, data=data, nameData=nameData)

 
#NONMEM code model 5
curwd=getwd()
if (file.exists("C:/nmv/run")) {
  setwd("C:/nmv/run")
  nameData<-list(covnames=c("wt"), 
               yvarlab="Sqrt(Theop. Conc.) Sqrt(mg/L)",
               xvarlab="Time since dose (hrs)",
               reparams=c("Ka", "V", "Cl"),
               params=c("Ka", "V", "Cl"),
               tparams=c("log(Ka)", "log(V)", "log(Cl)"),
               varnames=c("D[1,1]","D[1,2]","D[1,3]","D[2,2]","D[2,3]","D[3,3]")
               )

  results5<-RunNM(inputStructure="control.model5", data=data, nameData=nameData)
  setwd(curwd)
  }
  else {
  print("You do not have NONMEM.")
  out<-1
    } 
if (out==0) print(try(indEst(PKNLMEobject=results.nlme5, NMobject=results5, WBobject=WB2, outputType="R"))) 
}
}


\keyword{ models }% 
