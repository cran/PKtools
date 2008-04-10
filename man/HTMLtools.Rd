%copyright (c) 2005 MS Blanchard
\name{HTMLtools}
\alias{HTMLtools} 
\alias{HTMLtools.NONMEM} 
\alias{HTMLtools.PKNLME} 
\alias{HTMLtools.WinBUGS}
\title{HTMLtools}
\description{
HTMLtools is a method that outputs a HTML file of the parameter estimates
and diagnostic plots from an object of class NONMEM, PKNLME, or WinBUGS 
for a single dose population PK model with hierarchical data.
}
\usage{
HTMLtools(x, nameDir, nameFile, descStructure, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{an object from one of the following classes NONMEM, PKNLME, or WinBUGS}
  \item{nameDir}{the path and name of the directory where the HTML file will reside}
  \item{nameFile}{list of the names of the plots being output to the .html file}
  \item{descStructure}{list of arguments (pcts,nsig) for the function desc}
  \item{...}{additional arguments to be passed to lower level functions}
}
\details{
RunNM, RunNLME, and RunWB create the objects of the respective classes NONMEM, PKNLME and WinBUGS 
that can be read by this method.

nameFile lists the name of the output html file and the names of the plots being output to the html file.
Note the html file name should not have an html extension that will be added by the program and 
the plots should not have a png extension again that will be added by the program.  Finally,
note that there is a file0 in HTMLtools for the WinBUGS class to allow for inclusion of the density plots 
of the model coefficients.

\itemize{
\item nameFile<-list(file="wb", file0=hist, file1="trplt.wb",file2="diagplt.wb", 
file3="qqploti.wb", file4="qqnormre.wb", file5="covre.wb", file6="diagtrplti.wb", 
file7="diagtrpltp.wb")
}

Finally for the HTML file to be in color the correct path must be given in nameDir.
}


\value{
An HTML file of the results from the selected object.

The trellis plots including those from trplot, diagtrplot output only the first page of plots to the
HTML file and a png file of all pages is also created. The covariate plot allows for up
to 16 covariates again printing the first page in the HTML file. The covariates will be printed in order for
each random effect.
}


\author{M.S. Blanchard <sblanchard@coh.org>}

\seealso{RunNM, RunNLME, RunWB}
\examples{
   #NLME example
     setwd(tempdir())
     library(PKtools) 
     library(R2HTML)
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

     nameFile<-list(file="nlme.output", file1="trplt.nl",
                    file2="diagplt.nl", file3="qqploti.nl",
                    file4="qqnormre.nl", file5="covre.nl",
                    file6="diagtrplti.nl", file7="diagtrpltp.nl")

     descStructure<-list(pcts=c(.025,.05,.95,.975),nsig=4)

     model.def<-list(fixed.model=lKa+lV+lCl~1,random.model=lCl~1,
     start.lst=c(lKa=.3,lV=-.6,lCl=-3), form=conc~sonecpmt(dose, time,
     lV, lKa, lCl), control=nlmeControl(returnObject=FALSE))


MM<-RunNLME(inputStructure=model.def, data=data, nameData=nameData)


HTMLtools(x=MM, nameDir=tempdir(), nameFile = nameFile, 
       descStructure = descStructure) 
}

\keyword{ models }% 
