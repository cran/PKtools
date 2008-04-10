%copyright (c) 2005 MS Blanchard
\name{tex}
\alias{tex}
\alias{tex.NONMEM}
\alias{tex.PKNLME}
\alias{tex.WinBUGS}
\title{tex}
\description{
tex is a method that outputs a tex file of the parameter estimates
and diagnostic plots from an object of class NONMEM, PKNLME, or WinBUGS 
for a single dose population PK model with hierarchical data.
}
\usage{
tex(x, nameDir, nameFile, descStructure,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{an object from one of the following classes NONMEM, PKNLME, or WinBUGS}
  \item{nameDir}{the path and name of the directory where the HTML file will reside}
  \item{nameFile}{lists the name of the tex file and of the plots being output to the .tex file}
  \item{descStructure}{list of variables (pcts,nsig) for the function desc}
  \item{...}{additional arguments to be passed to lower level functions}
}
\details{
RunNM, RunNLME, and RunWB create the NONMEM, PKNLME and WinBUGS 
objects NM, MM, and WB,respectively, that can be read by this method.

nameFile lists the name of the tex file and the names of the plots being output to .tex file.
note the tex file should have a tex extension and the plots should have a ps extension.  Finally,
note that there is a file0 in tex.WinBUGS which includes the density plots of the model coefficients.

nameFile<-list(file="wb.tex", file0=hist.ps, file1="trplt.wb.ps",
               file2="diagplt.wb.ps", file3="qqploti.wb.ps",
               file4="qqnormre.wb.ps", file5="covre.wb.ps",
               file6="diagtrplti.wb.ps", file7="diagtrpltp.wb.ps")

}

\value{
A tex file of the results from the selected object.

The trellis plots including those from trplt, diagtrplt output the first page of plots to the
tex file and all pages to an accompanying postscript file. The covariate plot allows for up
to 16 covariates. printing the first page in the tex file and any additional plots to an  
accompanying postscript file.
}

\author{M.S. Blanchard <sblanchard@coh.org>}

\seealso{\code{\link{HTMLtools}}}

\examples{
#NLME example
setwd(tempdir())
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

nameFile<-list(file="nlme.tex", file1="trplt.nl.ps",
               file2="diagplt.nl.ps", file3="qqploti.nl.ps",
               file4="qqnormre.nl.ps", file5="covre.nl.ps",
               file6="diagtrplti.nl.ps", file7="diagtrpltp.nl.ps")

descStructure<-list(pcts=c(.025,.05,.95,.975),nsig=4)

model.def<-list(fixed.model=lKa+lV+lCl~1,random.model=lCl~1,
start.lst=c(lKa=.3,lV=-.6,lCl=-3), form=conc~sonecpmt(dose, time,
lV, lKa, lCl), control=nlmeControl(returnObject=FALSE))

MM<-RunNLME(inputStructure=model.def, data=data, 
nameData=nameData)

tex(MM, nameDir=tempdir(), nameFile = nameFile, descStructure = descStructure) 
}


\keyword{ models }% 
