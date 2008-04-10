# copyright (c) 2005 M S Blanchard 
# with exception of excerpts noted below

.PKtoolsOptions <- list(NONMEMloc = "C:/nmv/run", bugsRloc = "C:/bugsR")

.packageName <- "PKtools"

require(xtable)
require(R2HTML)
require(lattice)
require(nlme)

sonecpmt <- function (dose, time, lVol, lKa, lCl) {
    sqrt(dose * exp(lKa) * (exp(-(exp(lCl)/exp(lVol)) * time) -
        exp(-exp(lKa) * time))/(exp(lVol) * (exp(lKa) - (exp(lCl)/exp(lVol)))))
}

lonecpmt <- function (dose, time, lVol, lKa, lCl) {
    log(dose * exp(lKa) * (exp(-(exp(lCl)/exp(lVol)) * time) -
        exp(-exp(lKa) * time))/(exp(lVol) * (exp(lKa) - (exp(lCl)/exp(lVol)))))
}


PKtools.AIC <- function (loglike, n, K,...)
{
    AIC <- round(-2 * loglike + 2 * K, 3)
    AICc <- round(AIC + 2 * K * (K + 1)/(n - K - 1),3)
    of <- round(-2 * loglike,3)
    tmp<- cbind(AIC, AICc, of, K)
    tmp   
}


AICcomp<-function(PKNLMEobjects, NONMEMobjects){
 n1 <- length(PKNLMEobjects)
 n2 <- length(NONMEMobjects)
 if (n1 != n2) stop("model lists of unequal length")
AIC.table <- matrix(NA,nr=n1,nc=7)
for (i in 1:n1){
n<-nrow(PKNLMEobjects[[i]]$mm$fitted)
loglike<-logLik(PKNLMEobjects[[i]]$mm)[1]
K<-attr(logLik(PKNLMEobjects[[i]]$mm),"df")

#PKNLME
AIC.NLME<-data.frame(PKtools.AIC(loglike=loglike,n=n,K=K), row.names="")
#print(PKtools.AIC)s

#NONMEM
loglike.nm<- -.5*NONMEMobjects[[i]]$of
K.nm<-sum(NONMEMobjects[[i]]$param[,1]!=0)
AIC.NM<-data.frame(PKtools.AIC(loglike=loglike.nm,n=n,K=K.nm))
#print(AIC.NM)

#AIC Table
AIC.table[i,]<-c(AIC.NM[,1],AIC.NLME[,1],AIC.NM[,2],AIC.NLME[,2],AIC.NM[,3],AIC.NLME[,3],AIC.NM[,4])
}
AIC.table<-data.frame(AIC.table)
names(AIC.table)<-c("NM AIC", "NLME AIC","NM AICc", "NLME AICc", "NM of","NLME of", "K")
return(AIC.table)
}

desc <-
function (y, pcts = c(0.025, 0.05, 0.95, 0.975), nsig = 4) {
    x <- y[!is.na(y)]
    rn <- range(x)
    bin <- ifelse(is.na(y), 1, 0)
    NAs <- sum(bin)
if (NAs == 0){
    tmp <- as.numeric(signif(cbind(N = length(x),
        mean = mean(x), med = median(x), s = sqrt(var(x)), t(quantile(x,
            pcts)), min = rn[1], max = rn[2]), nsig))
    names(tmp) <- c("N", "Mean", "Med", "S", paste(pcts),
        "Min", "Max")}
else{
   tmp <- as.numeric(signif(cbind(N = length(x), NAs = NAs,
        mean = mean(x), med = median(x), s = sqrt(var(x)), t(quantile(x,
            pcts)), min = rn[1], max = rn[2]), nsig))
    names(tmp) <- c("N", "NA", "Mean", "Med", "S", paste(pcts),
        "Min", "Max")
    }
    tmp
}


coVar.id<-function (id, coVar, nameData) 
{
    splitu <- function(x, y) {
        sx <- split(x, y)
        usx <- lapply(sx, unique)
        unusx <- unlist(usx)
    }
    id.id <- cbind(splitu(id, id))
    if (is.data.frame(coVar) == "FALSE") 
        cov.id <- cbind(splitu(coVar, id))
    if (is.data.frame(coVar) == "TRUE") {
        ncov <- ncol(coVar)
        nid <- length(unique(id))
        cov.id <- matrix(NA, nr = nid, nc = ncov)
        for (j in 1:ncov) {
            cov.id[, j] <- splitu(coVar[, j], id)
        }
        cov.id <- data.frame(cov.id)
    }
    cov.id <- data.frame(cbind(id.id, cov.id))
    names(cov.id) <- c("id", nameData$covnames)
    return(cov.id=cov.id)
}

indEst<-function (PKNLMEobject, NMobject, WBobject, outputType="R") {
library(xtable)
    num <- ncol(WBobject$mean$beta)
    for (i in 1:num) {
        param.table <- data.frame(cbind(log(NMobject$ip[, i +
            1]), coef(PKNLMEobject$mm)[, i], WBobject$mean$beta[,
            i]))
        names(param.table) <- c("NONMEM", "NLME", "WinBUGS")
if (outputType=="tex"){
        print(paste("Individual Estimates -", WBobject$nameData$reparam[i], sep = " "))
        print(xtable(param.table))
        }
   else { 
        print(paste("Individual Estimates -", WBobject$nameData$reparam[i], sep = " "))
        print(param.table)
        }
}
}

pk <- function (pkvar = pkvar, covdata = cov, covnames = covnames)
{
    splitu<-function(x,y){
       sx<-split(x,y)
       usx<-lapply(sx,unique)
       unusx<-unlist(usx)
       }

    if (is.null(covdata)){ 
    dat<-pkvar
    pkdata<-pkvar
    }
    else{
    newcov <- data.frame(covdata)
    names(newcov) <- covnames
    pkdata <- cbind(pkvar, newcov)
    dat <- pkdata
    }

    if (is.null(covdata)){
    unusid <- splitu(pkvar$id, pkvar$id)
    unusdose <- splitu(pkvar$dose, pkvar$id)
    luid <- length(unique(pkvar$id))
    concd <- rep(NA, luid)
    timed <- rep(0, luid)
    names <- cbind(t(names(pkvar)))
    dosedt <- data.frame(unusid, unusdose, timed, concd)
    names(dosedt) <- names
    } 
    else{
    unusid <- splitu(pkvar$id, pkvar$id)
    unusdose <- splitu(pkvar$dose, pkvar$id)
    luid <- length(unique(pkvar$id))
    concd <- rep(NA, luid)
    timed <- rep(0, luid)
    if (is.data.frame(covdata)=="FALSE") cov.id<-splitu(covdata,pkvar$id)
    if (is.data.frame(covdata)=="TRUE") {
    ncov<-ncol(covdata)
    nid<-length(unique(pkvar$id))
    cov.id <- matrix(NA, nr =nid , nc = ncov)
    for (j in 1:ncov){
       cov.id[,j]<-splitu(covdata[,j],pkvar$id)
       }
    cov.id<-data.frame(cov.id)
    }
    names <- cbind(t(names(pkvar)), t(covnames))
    dosedt <- data.frame(unusid, unusdose, timed, concd, cov.id)
    names(dosedt) <- names
    }
    pkdata$dose <- NA
    newpkdata <- rbind(dosedt, pkdata)
    names(newpkdata) <- names
    oid <- order(newpkdata$id)
    NMdata <- newpkdata[oid, ]
    NMdata$evid <- ifelse(is.na(NMdata$dose), 0, 1)
    NMdata$dose <- ifelse(is.na(NMdata$dose), ".", NMdata$dose)
    NMdata$conc <- ifelse(is.na(NMdata$conc), ".", NMdata$conc)
    NMdata[12:22, ]

    write(t(NMdata), file = "NMdata", ncolumns = ncol(NMdata))
    return(list(dat = dat))
}



pk.nlme <- function (pkvar = pkvar, covdata = cov, covnames = covnames)
{
    if (is.null(covdata)) pkdata <-pkvar
    else{
    newcov <- data.frame(covdata)
    names(newcov) <- covnames
    pkdata <- cbind(pkvar, newcov)
    }	
}


paramEst <- function(PKNLMEobject, NMobject, WBobject){

#NONMEM
NMparam<-data.frame(NMobject$param)  
names(NMparam)<-c("NMparams","NMse")

#NLME
r<-summary(PKNLMEobject$mm)
fp.nlme<-r$tTable[,1]
var.nlme<-getPsi(r)
n<-ncol(var.nlme)
Gv<-NA
 for (j in n:1){
 for (i in j:1){
   G<-var.nlme[i, j]
   Gv<-cbind(Gv,G)  
}
}     
Gv<-Gv[-c(1)]
m<-n*(n+1)/2
rseq<-c(seq(1:m))
o<-order(rseq,decreasing=TRUE)
sigma2<-r$sigma^2
NLmean<-as.vector(rbind(cbind(fp.nlme),cbind(Gv[o]),sigma2))
fpst.nlme<-cbind(r$tTable[,2])
col2<-cbind(rep(NA,m+1))
NLse<-as.vector(rbind(fpst.nlme,col2))
NLparam<-data.frame(cbind(NLmean,NLse))
names(NLparam)<-c("NLparams", "NLse")

#WinBUGS
    mu.table <- t(apply(WBobject$sims.list$mu[,], 2, desc))
    mu.table <- data.frame(mu.table[,,drop=FALSE])
    names(mu.table) <- gsub("X", " ", names(mu.table))

n<-ncol(WBobject$mean$itau)
itau.meanv<-NA
itau.sev<-NA
 for (j in n:1){
 for (i in j:1){
   itau.mean<-desc(WBobject$sims.list$itau[,i, j])[2]
   itau.meanv<-rbind(itau.meanv,itau.mean)
   itau.se<-desc(WBobject$sims.list$itau[,i, j])[4]
   itau.sev<-rbind(itau.sev,itau.se)
}
}
itau.meanv<-cbind(itau.meanv[-c(1)])
itau.sev<-cbind(itau.sev[-c(1)])

m<-n*(n+1)/2
rseq<-c(seq(1:m))
o<-order(rseq,decreasing=TRUE)
itau.df<-data.frame(cbind(rseq,itau.meanv[o], itau.sev[o]),row.names=NULL)

sigma.table<-as.matrix(t(desc(WBobject$sims.list$sigma)))
sigma.table<-data.frame(sigma.table[,,drop=FALSE])

WBmean<-rbind(cbind(mu.table[,2]),cbind(itau.df[,2]),sigma.table[,2])
WBse<-rbind(cbind(mu.table[,4]),cbind(itau.df[,3]),sigma.table[,4])
WBparam<-data.frame(cbind(WBmean,WBse))
names(WBparam)<-c("WBparams","WBse")

params<-data.frame(cbind(NMparam,NLparam,WBparam))
row.names(params)<-rbind(cbind(NMobject$nameData$tparams),cbind(NMobject$nameData$varnames),"sigma^2")
return(params)
}



RunNLME<-function(inputStructure=inputStructure, data=data, 
nameData=nameData){
    pkdata <- pk.nlme(pkvar = data$pkvar, covdata = data$cov, 
        covnames = nameData$covnames)
    form <- inputStructure$form
    mm<-nlme(form, data = pkdata, fixed = inputStructure$fixed.model, 
    control = inputStructure$control, random = inputStructure$random.model, 
    groups = ~id, start = inputStructure$start.lst)
    if (is.null(data$cov))    MM <- list(mm=mm, pkdata=pkdata, nameData = nameData)
    else{
    cov.id <- coVar.id(data$pkvar$id, data$cov, nameData)
    MM <- list(mm=mm, cov.id = cov.id, pkdata=pkdata, nameData = nameData)
    }
    class(MM) <- "PKNLME"
    MM
}


RunNM<-function (inputStructure = inputStructure, data=data, nameData = nameData) {
    options(width = 100)
    pkdat <- pk(pkvar = data$pkvar, cov = data$cov, covnames = nameData$covnames)
    unlink("MVOF.TBL")
    system(paste("nmfe5", inputStructure, "report.nonmem"),show.output.on.console=FALSE)
    of <- as.data.frame(read.table("MVOF.TBL"))[1, 1]
    param <- as.data.frame(read.table("PAR.TAB", header = TRUE, 
        row.names = "Parameter"))
    re <- as.data.frame(read.table("NMRE.TBL", header = TRUE, 
        skip = 1))
    ip <- as.data.frame(read.table("NMIP.TBL", header = TRUE, 
        skip = 1))
    if (!is.null(data$cov)) cov <- as.data.frame(read.table("NMCOV.TBL", header = TRUE, 
        skip = 1))
    NMpred <- as.data.frame(read.table("NMPR.TBL", header = TRUE, 
        skip = 1))
    pred <- NMpred[NMpred$EVID == 0, ]
    if (is.null(data$cov)) NM <- list(of = of, param = param, ip = ip, re = re, pred = pred, 
        pkdat = pkdat, nameData = nameData)
    else NM <- list(of = of, param = param, ip = ip, re = re, pred = pred, 
        cov = cov, pkdat = pkdat, nameData = nameData)
    class(NM) <- "NONMEM"
    NM
}

RunWB<-function (inputStructure = inputStructure, data = data, nameData = nameData, WBargs = WBargs) 
{
    out<-bugs(data = data$data, inits = WBargs$inits, parameters.to.save = WBargs$parameters, 
        model.file = inputStructure, n.chains = WBargs$n.chains, 
        n.iter = WBargs$n.iter, debug = WBargs$debug, n.burnin = WBargs$n.burnin, 
        n.thin = WBargs$n.thin, print.summary = FALSE, plot.summary = TRUE)
    if ("conc" %in% WBargs$parameters) data$data$conc[is.na(data$data$conc)]<-mean$conc
    id.v <- data$id
    mean$iresid <- data$data$conc - mean$ipred
    mean$presid <- data$data$conc - mean$ppred
    n<-length(data$data$conc)
    ipred.v <- t(structure(.Data = c(t(mean$ipred)), .Dim = c(1, n)))
    ppred.v <- t(structure(.Data = c(t(mean$ppred)), .Dim = c(1, n)))
    conc.v <- t(structure(.Data = c(t(data$data$conc)), .Dim = c(1, n)))
    time.v <- t(structure(.Data = c(t(data$data$time)), .Dim = c(1, n)))
    pred <- data.frame(cbind(id.v, time.v, conc.v, ppred.v, ipred.v,mean$presid, mean$iresid))
    names(pred) <- c("id", "time", "conc", "ppred", "ipred","presid","iresid")
    if (is.null(data$cov))  WB <- list(mean = mean, sims.list = out$sims.list, data=data, nameData = nameData, pred=pred)
    else{ 
    cov.id<-coVar.id(data$id, data$cov, nameData)
    WB <- list(mean = mean, sims.list = out$sims.list, data=data, nameData = nameData, pred=pred, cov.id=cov.id)
    } 
    class(WB) <- "WinBUGS"
    WB
}



print.PKNLME<-function (x,...){
    param<-x$mm
    print(param)
    invisible(NULL)
}

print.NONMEM <- function (x,...) {
    NM.param<-data.frame(x$param)
    names(NM.param)<-c("Estimate","Standard Error")
#   everything must come from x
    rownames(NM.param)<-rbind(cbind(x$nameData$tparams),cbind(x$nameData$varnames),"sigma^2")
    cat("object of class NONMEM\n")
    cat(paste("the objective function is:\n"))
    print(x$of)
    cat(paste("the population parameters are:\n"))
    print(NM.param)
    invisible(NULL)
}


print.WinBUGS <-function (x,...) 
{
    cat(paste("the population parameters are:\n"))
    cat(paste("\n"))
    mu <- x$mean$mu
    cat(paste("mu\n"))
    print(mu)
    cat(paste("\n"))
    D <- x$mean$itau
    cat(paste("D\n"))
    print(D)
    cat(paste("\n"))
    sigma2 <- x$mean$sigma2
    cat(paste("sigma2\n"))
    print(sigma2)
    invisible(NULL)
}





diagplot.PKNLME<- function (x,...) {
    CONC<-x$pkdata$conc
    PRED<-x$mm$fitted[,1]
    IPRE<-x$mm$fitted[,2]
    RES<-x$mm$resid[,1]
    IRES<-x$mm$resid[,2] 
    dat<-data.frame(cbind(CONC,PRED,IPRE,RES,IRES))
    par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))
    plot(dat$PRED, dat$CONC, xlab = "Predicted values", ylab = "Observed values",
        main = "Population Model")
    abline(c(0, 1))
    plot(dat$IPRE, dat$CONC, xlab = "Predicted values", ylab = "Observed values",
        main = "Individual Model")
    abline(c(0, 1))
    plot(dat$PRED, dat$RES, xlab = "Predicted values", ylab = "Residuals")
    abline(h = 0)
    plot(dat$IPRE, dat$IRES, xlab = "Predicted values", ylab = "Residuals")
    abline(h = 0)
}
diagplot.NONMEM <- function (x,...) {
    dat<-x$pred
    par(mfrow = c(3, 2), mar = c(5, 4, 4, 2))
    plot(dat$PRED, dat$CONC, xlab = "Predicted values", ylab = "Observed values",
        main = "Population Model")
    abline(c(0, 1))
    plot(dat$IPRE, dat$CONC, xlab = "Predicted values", ylab = "Observed values",
        main = "Individual Model")
    abline(c(0, 1))
    plot(dat$PRED, dat$RES, xlab = "Predicted values", ylab = "Residuals")
    abline(h = 0)
    plot(dat$IPRE, dat$IRES, xlab = "Predicted values", ylab = "Residuals")
    abline(h = 0)
    plot(dat$PRED, dat$WRES, xlab = "Predicted values", ylab = " Weighted Residuals")
    abline(h = 0)
    plot(dat$IPRE, dat$IWRE, xlab = "Predicted values", ylab = "Weighted Residuals")
    abline(h = 0)
}


diagplot.WinBUGS<-function (x,...) 
{
dat<-x$pred
    par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))
    plot(dat$ppred, dat$conc, xlab = "Predicted values", ylab = "Observed values", 
        main = "Population Model")
    abline(c(0, 1))
    plot(dat$ipred, dat$conc, xlab = "Predicted values", ylab = "Observed values", 
        main = "Individual Model")
    abline(c(0, 1))
    plot(dat$ppred, dat$presid, xlab = "Predicted values", ylab = "Residuals")
    abline(h = 0)
    plot(dat$ipred, dat$iresid, xlab = "Predicted values", ylab = "Residuals")
    abline(h = 0)
}



diagplot<-function(x, ...)
UseMethod("diagplot")

diagplot.default <- function(x, ...)
stop(paste("diagplot method does not yet exist for instances of class",
class(x)))



residplot.PKNLME <- function (x,level,...) {
    CONC<-x$pkdata$conc
    PRED<-x$mm$fitted[,1]
    IPRE<-x$mm$fitted[,2]
    RES<-x$mm$resid[,1]
    IRES<-x$mm$resid[,2] 
    dat<-data.frame(cbind(CONC,PRED,IPRE,RES,IRES))
    par(mfrow = c(1, 1))
    if (level == "p"){
    plot(dat$PRED, dat$RES, xlab = "Population Predicted values", 
    ylab = "Population Residual values")
    abline(h=0)
    }
    else
    plot(dat$IPRE, dat$IRES, xlab = "Individual Predicted values", 
    ylab = "Individual Residual values")
    abline(h=0)
}


residplot.NONMEM <- function (x,level,...) {
    par(mfrow = c(1, 1))
    dat<-x$pred
    if (level == "p"){
    plot(dat$PRED, dat$WRES, xlab = "Population Predicted values", 
    ylab = "Weighted Population Residual values")
    abline(h=0)
    }
    else
    plot(dat$IPRE, dat$IWRE, xlab = "Individual Predicted values", 
    ylab = "Weighted Individual Residual values")
    abline(h=0)
}



residplot.WinBUGS <- function (x,level,...) {
diagdata<-x$pred
    par(mfrow = c(1, 1))
    dat<-x$pred
    if (level == "p"){
    plot(diagdata$ppred, diagdata$presid, xlab = "Population Predicted values", 
    ylab = "Weighted Population Residual values")
    abline(h=0)
    }
    else
    plot(diagdata$ipred, diagdata$iresid, xlab = "Individual Predicted values", 
    ylab = "Weighted Individual Residual values")
    abline(h=0)
}



residplot<-function(x, ...)
UseMethod("residplot")

residplot.default <- function(x, ...)
stop(paste("residplot method does not yet exist for instances of class",
class(x)))



obvsprplot.PKNLME <- function (x,level,...) {
    CONC<-x$pkdata$conc
    PRED<-x$mm$fitted[,1]
    IPRE<-x$mm$fitted[,2]
    RES<-x$mm$resid[,1]
    IRES<-x$mm$resid[,2] 
    dat<-data.frame(cbind(CONC,PRED,IPRE,RES,IRES))
    par(mfrow = c(1, 1))
    if (level == "p"){
    plot(dat$PRED, dat$CONC, xlab = " Population Predicted values", ylab = "Observed values")
    abline(c(0,1))
    }
    else
    plot(dat$IPRE, dat$CONC, xlab = "Individual Predicted values", ylab = "Observed values")
    abline(c(0,1))
}



obvsprplot.NONMEM <- function (x,level,...) {
    par(mfrow = c(1, 1))
    dat<-x$pred
    if (level == "p"){
    plot(dat$PRED, dat$CONC, xlab = " Population Predicted values", ylab = "Observed values")
    abline(c(0,1))
    }
    else
    plot(dat$IPRE, dat$CONC, xlab = "Individual Predicted values", ylab = "Observed values")
    abline(c(0,1))
}

obvsprplot.WinBUGS <- function (x, level,...) {
diagdata<-x$pred
    par(mfrow = c(1, 1))
    dat<-x$pred
    if (level == "p"){
    plot(diagdata$ppred, diagdata$conc, xlab = " Population Predicted values", ylab = "Observed values")
    abline(c(0,1))
    }
    else
    plot(diagdata$ipred, diagdata$conc, xlab = "Individual Predicted values", ylab = "Observed values")
    abline(c(0,1))
}




obvsprplot<-function(x, ...)
UseMethod("obvsprplot")

obvsprplot.default <- function(x, ...)
stop(paste("obsvpplot method does not yet exist for instances of class",
class(x)))


trplot.PKNLME<- function (x,xvarlab,yvarlab,pages,...){
    id  <-x$pkdata$id
    CONC<-x$pkdata$conc
    TIME<-x$pkdata$time
    dat<-data.frame(cbind(id,CONC,TIME))
    library(lattice)
    ID <- as.factor(dat$id)
    xyplot(CONC ~ TIME | ID, data = dat, xlab = xvarlab, xlim = c(-1,
        max(TIME)), ylab = yvarlab, ylim = c(min(CONC), max(CONC)), 
        span = 2, layout = c(4, 4, pages), aspect = 1,
        panel = function(x,y, span) {
            panel.xyplot(x, y, col = "black", cex = 1.2, type="b")})
}



trplot.NONMEM<- function (x,xvarlab,yvarlab,pages,...){
dat<-x$pred
    library(lattice)
    Id <- as.factor(dat$ID)
    xyplot(CONC ~ TIME | Id, data = dat, xlab = xvarlab, xlim = c(-1,
        max(time)), ylab = yvarlab, ylim = c(min(dat$CONC), max(dat$CONC)),
        span = 2, layout = c(4, 4, pages), aspect = 1,panel = function(x,
            y, span) {
            panel.xyplot(x, y, col = "black", cex = 1.2, type="b")
        })
}


trplot.WinBUGS<- function (x, xvarlab,yvarlab,pages,...){
library(lattice)
dat<-x$pred
    ID <- as.factor(dat$id)
    xyplot(conc ~ time | ID, data = dat, xlab = xvarlab, xlim = c(-1,
        max(time)), ylab = yvarlab, ylim = c(min(dat$conc), max(dat$conc)),
        span = 2, layout = c(4, 4, pages), aspect = 1,panel = function(x,
            y, span) {
            panel.xyplot(x, y, col = "black", cex = 1.2, type="b")
        })
}



trplot<-function (x,xvarlab,yvarlab,pages,...)
UseMethod("trplot")

trplot.default <- function (x,xvarlab,yvarlab,pages,...)
stop(paste("trplot method does not yet exist for instances of class",
class(x)))

diagtrplot.PKNLME<-
function (x, level , xvarlab = xvarlab, yvarlab = yvarlab, pages = pages,...) 
{
    ID<-x$pkdata$id
    CONC<-x$pkdata$conc
    TIME<-x$pkdata$time
    PRED<-x$mm$fitted[,1]
    IPRE<-x$mm$fitted[,2]
    RES <-x$mm$resid[,1]
    IRES<-x$mm$resid[,2] 
    dat<-data.frame(cbind(ID,CONC,TIME,PRED,IPRE,RES,IRES))
    library(lattice)
    typedv <- rep(1, length(dat$ID))
    typepred <- rep(2, length(dat$ID))
    dvdata <- cbind(dat$ID, dat$TIME, dat$CONC, typedv)
    if (level =="p") preddata <- cbind(dat$ID, dat$TIME, dat$PRED, typepred)
    else preddata <- cbind(dat$ID, dat$TIME, dat$IPRE, typepred)
    plotdata <- as.data.frame(rbind(dvdata, preddata))
    names(plotdata) <- c("id", "time", "dv.pred", "typepred")
    newid <- factor(plotdata$id)
    sps <- trellis.par.get("superpose.symbol")
    sps$pch <- 1:7
    trellis.par.set("superpose.symbol", sps)
    print(xyplot(dv.pred ~ time | newid, data = plotdata, groups = typepred, 
        xlab = xvarlab, ylab = yvarlab, span = 2, layout = c(4, 
            4, pages), aspect = 1, type = "b", panel = "panel.superpose", 
        panel.group = "panel.xyplot", key = list(columns = 2, 
            text = list(paste(c("Values:  Observed", " Predicted"))), 
            points = Rows(sps, 1:2))))
}


diagtrplot.NONMEM<-
function (x, level , xvarlab = xvarlab, yvarlab = yvarlab, pages = pages,...) 
{
    dat<-x$pred
    library(lattice)
    typedv <- rep(1, length(dat$ID))
    typepred <- rep(2, length(dat$ID))
    dvdata <- cbind(dat$ID, dat$TIME, dat$CONC, typedv)
    if (level =="p") preddata <- cbind(dat$ID, dat$TIME, dat$PRED, typepred)
    else preddata <- cbind(dat$ID, dat$TIME, dat$IPRE, typepred)
    plotdata <- as.data.frame(rbind(dvdata, preddata))
    names(plotdata) <- c("id", "time", "dv.pred", "typepred")
    newid <- factor(plotdata$id)
    sps <- trellis.par.get("superpose.symbol")
    sps$pch <- 1:7
    trellis.par.set("superpose.symbol", sps)
    print(xyplot(dv.pred ~ time | newid, data = plotdata, groups = typepred, 
        xlab = xvarlab, ylab = yvarlab, span = 2, layout = c(4, 
            4, pages), aspect = 1, type = "b", panel = "panel.superpose", 
        panel.group = "panel.xyplot", key = list(columns = 2, 
            text = list(paste(c("Values:  Observed", " Predicted"))), 
            points = Rows(sps, 1:2))))
}


diagtrplot.WinBUGS<-function (x, level , xvarlab = xvarlab, yvarlab = yvarlab, pages = pages,...) 
{
trdata<-x$pred
    library(lattice)
    typedv <- rep(1, length(trdata$id))
    typepred <- rep(2, length(trdata$id))
    dvdata <- cbind(trdata$id, trdata$time, trdata$conc, typedv)
    if (level =="p") preddata <- cbind(trdata$id, trdata$time, trdata$ppred, typepred)
    else preddata <- cbind(trdata$id, trdata$time, trdata$ipred, typepred)
    plotdata <- as.data.frame(rbind(dvdata, preddata))
    names(plotdata) <- c("id", "time", "dv.pred", "typepred")
    newid <- factor(plotdata$id)
    sps <- trellis.par.get("superpose.symbol")
    sps$pch <- 1:7
    trellis.par.set("superpose.symbol", sps)
    print(xyplot(dv.pred ~ time | newid, data = plotdata, groups = typepred, 
        xlab = xvarlab, ylab = yvarlab, span = 2, layout = c(4, 
            4, pages), aspect = 1, type = "b", panel = "panel.superpose", 
        panel.group = "panel.xyplot", key = list(columns = 2, 
            text = list(paste(c("Values:  Observed", " Predicted"))), 
            points = Rows(sps, 1:2))))
}
      
            
            
diagtrplot <- function (x, level , xvarlab = xvarlab, yvarlab = yvarlab, pages = pages,...) 
UseMethod("diagtrplot")


tex.PKNLME <-function (x=x, nameDir=nameDir, 
    nameFile = nameFile, descStructure = descStructure,...) 
{ 


cat ("Results written to", nameDir, "\n")

if (!is.null(x$cov.id)){
#Covariates
    cov.desc <- t(apply(x$cov.id[, ], 2, desc, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    cov.desc <- data.frame(cov.desc[-c(1),, drop=FALSE])
    names(cov.desc) <- gsub("X", " ", names(cov.desc))
#    row.names(cov.desc) <- x$nameData$covnames
#    cov.desc <- as.matrix(cov.desc)
#    cov.desc <- formatC(cov.desc, format = "fg", digits = 3)
}

#NLME parameters
    r <- summary(x$mm)
    r.tTable <- as.matrix(r$tTable)
    r.tTable <- formatC(r.tTable, format = "fg", digits = 3)
    sigma.nlme <- data.frame(r$sigma^2)
    names(sigma.nlme) <- c("sigma2")
    sigma.nlme <- as.matrix(sigma.nlme)
    sigma.nlme <- data.frame(formatC(sigma.nlme, format = "fg", digits = 3), row.names="")
    var.nlme <- getPsi(r)
    var.nlme <- as.matrix(var.nlme)
    var.nlme <- formatC(var.nlme, format = "fg", digits = 3)
    n <- nrow(data$pkvar)

#AIC table
    AIC.table <- data.frame(PKtools.AIC(loglike = logLik(x$mm), n = n, 
        K = attr(logLik(x$mm), "df")),row.names="")
    colnames(AIC.table) <- c("AIC", "AICc", "OF", "K")

#individual PK parameters
    coef.table <- t(apply(coef(x$mm)[, ], 2, desc, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    coef.table <- data.frame(coef.table[,,drop=FALSE])
    names(coef.table) <- gsub("X", " ", names(coef.table))
    coef.table <- as.matrix(coef.table)
    coef.table <- formatC(coef.table, format = "fg", digits = 3)

    pages = round(length(unique(data$pkvar$id))/16)

    postscript(file = nameFile$file1, paper = "letter")
    library(lattice)
    print(trplot(x, xvarlab=x$nameData$xvarlab, yvarlab=x$nameData$yvarlab, pages=pages))
    dev.off()
    postscript(file = "trplt.ps", paper = "letter")

    print(trplot(x, xvarlab=x$nameData$xvarlab, yvarlab=x$nameData$yvarlab, pages=1))
    dev.off()
    postscript(file = nameFile$file2, paper = "letter")
    diagplot(x)
    dev.off()
    
    postscript(file = nameFile$file3, paper = "letter")
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 5))
    qqnorm(x$mm$residuals[, 1], main = "Residuals from population model")
    qqline(x$mm$residuals[, 1])
    qqnorm(x$mm$residuals[, 2], main = "Residuals from individual model")
    qqline(x$mm$residuals[, 2])
    dev.off()


    postscript(file = nameFile$file4, paper = "letter")
    re <- x$mm$coef$random$id
    nre <- ncol(re)
    par(mfrow = c(nre, 1), mar = c(4, 5, 3, 5))
    for (i in 1:nre) {
        qqnorm(re[, i], main = paste("Random effects for log(", 
            x$nameData$reparams[i], ")", sep = ""))
        qqline(re[, i])
    }
    dev.off()

if (!is.null(x$cov.id)){
    postscript(file = nameFile$file5, paper = "letter")
    par(mfrow = c(1, 1), mar = c(4, 5, 3, 5))
    ncov <- ncol(x$cov.id)
    for (j in 2:ncov) {
        for (i in 1:nre) {
            plot(x$cov.id[, j], re[, i], xlab = x$nameData$covnames[j - 
                1], ylab = paste("Random effects for log(", x$nameData$reparams[i], 
                ")", sep = ""))
            abline(h = 0)
        }
    }
    dev.off()
    postscript(file = "re1.ps", paper = "letter")
    par(mar=c(4, 5, 3, 5), bg="grey")
    if (ncov - 1 <= 16) 
        par(mfrow = c(ceiling((ncov - 1)/4), 4))
    if (ncov - 1 <= 12) 
        par(mfrow = c(ceiling((ncov - 1)/3), 3))
    if (ncov - 1 <= 8) 
        par(mfrow = c(ceiling((ncov - 1)/2), 2))
    if (ncov - 1 <= 4) 
        par(mfrow = c(ncov - 1, 1))
    for (j in 2:ncov) {
        plot(x$cov.id[, j], re[, 1], xlab = x$nameData$covnames[j - 
            1], ylab = paste("Random effects for log(", x$nameData$reparams[1], 
            ")", sep = ""))
        abline(h = 0)
    }
    dev.off()
}
    postscript(file = nameFile$file6, paper = "letter")
    diagtrplot(x, "i", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = pages) 
    dev.off()
    postscript(file = "diagtrplti.ps", paper = "letter")
    diagtrplot(x, "i", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = 1) 
    dev.off()
    postscript(file = nameFile$file7, paper = "letter")
    diagtrplot(x, "p", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = pages) 
    dev.off()
    postscript(file = "diagtrpltp.ps", paper = "letter")
    diagtrplot(x, "p", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = 1) 
    dev.off()
    library(xtable)
    format.NLME <- function(file = "mynlme.tex", file1 = "trplt.nl.ps", 
        file2 = "diagplt.nl.ps", file3 = "qqploti.nl.ps", file4 = "qqnormre.nl.ps", 
        file5 = "covre.nl.ps", file6 = "diagtrplti.nl.ps", file7 = "diagtrpltp.nl.ps") {
        cat("\\documentclass{article}\n", file = file)
        cat("\\usepackage{epsfig}\n", file = file, append = TRUE)
        cat("\\textwidth=6in\n", file = file, append = TRUE)
        cat("\\textheight=8.5in\n", file = file, append = TRUE)
        cat("\\oddsidemargin=.1in\n", file = file, append = TRUE)
        cat("\\evensidemargin=.1in\n", file = file, append = TRUE)
        cat("\\headheight=-.5in\n", file = file, append = TRUE)
        cat("\\begin{document}\n", file = file, append = TRUE)
        cat("\\title{Results from NLME}\n", file = file, append = TRUE)
        cat("\\maketitle\n", file = file, append = TRUE)
        print(xtable(r$tTable, "Summary Population mean Parameter estimate from NLME"), 
            file = file, append = TRUE)
        print(xtable(var.nlme, "Intersubject Variance-Covariance Matrix from NLME"), 
            file = file, append = TRUE)
        print(xtable(sigma.nlme, "Intrasubject Variance from NLME"), 
            file = file, append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
        print(xtable(AIC.table, "Model Selection Criteria"), 
            file = file, append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
        print(xtable(coef.table, "Descriptive statistics for individual Parameters"), 
            file = file, append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
 if (!is.null(x$cov.id)) print(xtable(cov.desc, "Covariates"), file = file, append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file1, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Trellis plot of Concentration vs Time Data}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file2, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Predicted vs Observed and Residual vs Predicted Plots}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file3, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{QQPlots: Stage I Model, QQplots provide a visual assessment for normality.}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file4, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{QQPlots: Stage II Model}\n", file = file, 
            append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
if (!is.null(x$cov.id)){
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file5, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Individual Random effects vs Covariates}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
}
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file6, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Individual Predicted and Observed Values vs Time}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file7, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Population Predicted and Observed Values vs Time}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\end{document}\n", file = file, append = TRUE)
        cat(paste("\ncompleted export to", file))
    }
    format.NLME(file = nameFile$file, file1 = "trplt.ps", file2 = nameFile$file2, 
        file3 = nameFile$file3, file4 = nameFile$file4, file5 = "re1.ps", 
        file6 = "diagtrplti.ps", file7 = "diagtrpltp.ps")
#   return(list(MM = MM))
}



tex.NONMEM<-function (x=x, nameDir=nameDir,
    nameFile = nameFile, descStructure = descStructure,...) 
{


cat ("Results written to", nameDir, "\n")

    
#AIC
    n <- nrow(x$pkdat$dat)
    pages <- round(length(unique(x$pkdat$dat$id))/16)
    K <- sum(x$param[, 1] != 0)
    loglike <- -0.5 * x$of
    AIC.table <- data.frame(PKtools.AIC(loglike = loglike, n = n, 
        K = K), row.names="")

#covariates
if (!is.null(x$cov)){
    cov.desc <- t(apply(x$cov[, ], 2, desc, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    cov.desc <- data.frame(cov.desc[-c(1),,drop=FALSE])
    names(cov.desc) <- gsub("X", " ", names(cov.desc))
#    row.names(cov.desc) <- x$nameData$covnames
#    cov.desc <- as.matrix(cov.desc)
#    cov.desc <- formatC(cov.desc, format = "fg", digits = 3)
}

#individual parameters
    ip.desc <- t(apply(x$ip[, ], 2, desc, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    ip.desc <- data.frame(ip.desc[-c(1),,drop=FALSE])
    names(ip.desc) <- gsub("X", " ", names(ip.desc))
#    row.names(ip.desc) <- x$nameData$reparams
#    ipdesc <- as.matrix(ip.desc)
#    ipdesc <- formatC(ipdesc, format = "fg", digits = 3)

    postscript(file = nameFile$file1, paper = "letter")
    library(lattice)
    print(trplot(x, xvarlab=x$nameData$xvarlab, yvarlab=x$nameData$yvarlab, pages=pages))
    dev.off()
    postscript(file = "trplt.ps", paper = "letter")
    print(trplot(x, xvarlab=x$nameData$xvarlab, yvarlab=x$nameData$yvarlab, pages=1))
    dev.off()
    postscript(file = nameFile$file2, paper = "letter")
    diagplot(x)
    dev.off()
    postscript(file = nameFile$file3, paper = "letter")
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 5))
    qqnorm(x$pred$WRES, main = "Weighted residuals from the population model")
    qqline(x$pred$WRES)
    qqnorm(x$pred$IWRE, main = "Weighted residuals from the individual model")
    qqline(x$pred$IWRE)
    dev.off()

    postscript(file = nameFile$file4, paper = "letter")
    re <- x$re
    nre <- ncol(re)
    par(mfrow = c(nre - 1, 1), mar = c(4, 5, 3, 5))
    for (i in 2:nre) {
        qqnorm(re[, i], main = paste("Random effects for log(", 
            x$nameData$reparams[i - 1], ")", sep = ""))
        qqline(re[, i])
    }
    dev.off()

if (!is.null(x$cov)){
    postscript(file = nameFile$file5, paper = "letter")
    ncov <- ncol(x$cov)
    re <- x$re
    nre <- ncol(re)
    par(mfrow = c(1, 1), mar = c(4, 5, 3, 5))
    for (j in 2:ncov) {
        for (i in 2:nre) {
            plot(x$cov[, j], re[, i], xlab = x$nameData$covnames[j - 
                1], ylab = paste("Random effects for log(", x$nameData$reparams[i - 
                1], ")", sep = ""))
            abline(h = 0)
        }
    }
    dev.off()

    postscript(file = "re1.ps", paper = "letter")
    par(mar=c(4, 5, 3, 5), bg="grey")
    if (ncov - 1 <= 16) 
        par(mfrow = c(ceiling((ncov - 1)/4), 4))
    if (ncov - 1 <= 12) 
        par(mfrow = c(ceiling((ncov - 1)/3), 3))
    if (ncov - 1 <= 8) 
        par(mfrow = c(ceiling((ncov - 1)/2), 2))
    if (ncov - 1 <= 4) 
        par(mfrow = c(ncov - 1, 1))
    for (j in 2:ncov) {
        plot(x$cov[, j], re[, 2], xlab = x$nameData$covnames[j - 
            1], ylab = paste("Random effects for log(", x$nameData$reparams[1], 
            ")", sep = ""))
        abline(h = 0)
    }
    dev.off()
}
    postscript(file = nameFile$file6, paper = "letter")
    diagtrplot(x, "i", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = pages) 
    dev.off()
    
    postscript(file = "diagtrplti.ps", paper = "letter")
    diagtrplot(x, "i", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = 1) 
    dev.off()
    
    postscript(file = nameFile$file7, paper = "letter")
    diagtrplot(x, "p", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = pages) 
    dev.off()
    postscript(file = "diagtrpltp.ps", paper = "letter")
    diagtrplot(x, "p", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = 1) 
    dev.off()
    NM.param <- data.frame(x$param)
    names(NM.param) <- c("Estimate", "Standard Error")
    rownames(NM.param) <- rbind(cbind(x$nameData$tparams), cbind(x$nameData$varnames), 
        "sigma2")
    NMparam <- as.matrix(NM.param)
    NMparam <- formatC(NMparam, format = "fg", digits = 3)

    library(xtable)
    format.NONMEM <- function(file = "mynm.tex", file1 = nameFile$file1, 
        file2 = nameFile$file2, file3 = nameFile$file3, file4 = nameFile$file4, 
        file5 = nameFile$file5, file6 = nameFile$file6, file7 = nameFile$file7) {
       cat("\\documentclass{article}\n", file = file)
        cat("\\usepackage{epsfig}\n", file = file, append = TRUE)
        cat("\\textwidth=6in\n", file = file, append = TRUE)
        cat("\\textheight=8.5in\n", file = file, append = TRUE)
        cat("\\oddsidemargin=.1in\n", file = file, append = TRUE)
        cat("\\evensidemargin=.1in\n", file = file, append = TRUE)
        cat("\\headheight=-.5in\n", file = file, append = TRUE)
        cat("\\begin{document}\n", file = file, append = TRUE)
        cat("\\title{Results from NONMEM}\n", file = file, append = TRUE)
        cat("\\maketitle\n", file = file, append = TRUE)
        print(xtable(NMparam, "Population parameter estimates"), 
            file = file, append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
        print(xtable(AIC.table, "Model Selection Criteria"), 
            file = file, append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
        print(xtable(ip.desc, "Individual Parameters"), file = file, 
            append = TRUE)
        if (!is.null(x$cov))  print(xtable(cov.desc, "Covariates"), file = file, append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file1, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Trellis plot of Concentration vs Time Data}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file2, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Predicted vs Observed and Residual vs Predicted Plots}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file3, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{QQPlots: Stage I Model, QQplots provide a visual assessment for normality.}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file4, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{QQPlots: Stage II Model}\n", file = file, 
            append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
if (!is.null(x$cov)){ 
       cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file5, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Individual Random effects vs Covariates}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
}
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file6, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Individual Predicted and Observed Values vs Time}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file7, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Population Predicted and Observed Values vs Time}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\end{document}\n", file = file, append = TRUE)
        cat(paste("\ncompleted export to", file))
    }
    format.NONMEM(file = nameFile$file, file1 = "trplt.ps", file2 = nameFile$file2, 
        file3 = nameFile$file3, file4 = nameFile$file4, file5 = "re1.ps", 
        file6 = "diagtrplti.ps", file7 = "diagtrpltp.ps")
#         return(list(NM = NM))

}


tex.WinBUGS<-function (x=x, nameDir=nameDir, nameFile = nameFile, descStructure = descStructure,...) 
{


cat ("Results written to", nameDir, "\n")

#Covariates
if (!is.null(x$cov)){
    ncov <- ncol(x$cov)
    cov.desc <- t(apply(x$cov[, ], 2, desc, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    cov.desc <- data.frame(cov.desc[-c(1),,drop=FALSE])
    names(cov.desc) <- gsub("X", " ", names(cov.desc))
    row.names(cov.desc) <- x$nameData$covnames
#    covdesc <- as.matrix(cov.desc)
#    covdesc <- formatC(covdesc, format = "fg", digits = 3)
}
    postscript(file = "qqploti.ps", paper = "letter")
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 5))
    qqnorm(x$pred$presid, main = "Residuals from population model")
    qqline(x$pred$presid)
    qqnorm(x$pred$iresid, main = "Residuals from individual model")
    qqline(x$pred$iresid)
    dev.off()
    postscript(file = nameFile$file0, paper = "letter")
    count <- ncol(x$sims.list$mu)
    par(mfrow = c(count, 1), mar = c(4, 5, 3, 5))
    for (i in 1:count) {
        plot(density(x$sims.list$mu[, i]), main = x$nameData$coef[i])
        abline(v = x$mean$mu[i])
    }
    dev.off()
    pages <- round(length(unique(x$pred$id))/16)
    postscript(file = nameFile$file1, paper = "letter")
    library(lattice)
    print(trplot(x,xvarlab=x$nameData$xvarlab, yvarlab=x$nameData$yvarlab, pages=pages))
    dev.off()
    postscript(file = "trplt.ps", paper = "letter")
    print(trplot(x,xvarlab=x$nameData$xvarlab, yvarlab=x$nameData$yvarlab, pages=1))
    dev.off()
    postscript(file = nameFile$file2, paper = "letter")
    diagplot(x=x)
    dev.off()
    postscript(file = nameFile$file3, paper = "letter")
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 5))
    qqnorm(x$pred$presid, main = "Residuals from population model")
    qqline(x$pred$presid)
    qqnorm(x$pred$iresid, main = "Residuals from individual model")
    qqline(x$pred$iresid)
    dev.off()
    postscript(file = nameFile$file4, paper = "letter")
    re <- x$mean$re
    nre <- ncol(re)
    par(mfrow = c(nre, 1), mar = c(4, 5, 3, 5))
    for (i in 1:nre) {
        qqnorm(re[, i], main = paste("log(", 
            x$nameData$reparams[i], ")", sep = ""))
        qqline(re[, i])
    }
    dev.off()
if (!is.null(x$cov)){
    postscript(file = nameFile$file5, paper = "letter")
    re <- x$mean$re
    ncov <- ncol(x$cov.id)
    nre <- ncol(re)
    par(mfrow = c(1, 1), mar = c(4, 5, 3, 5))
    for (i in 1:nre) {
        for (j in 2:ncov) {
            plot(x$cov.id[, j], re[, i], xlab = x$nameData$covnames[j - 
                1], ylab = paste("Random effects for log(", x$nameData$reparams[i], 
                ")", sep = ""))
            abline(h = 0)
        }
    }
    dev.off()
    postscript(file = "re1.ps", paper = "letter")
    par(mar=c(4, 5, 3, 5), bg="grey")
    if (ncov <= 16) 
        par(mfrow = c(ceiling(ncov/4), 4))
    if (ncov <= 12) 
        par(mfrow = c(ceiling(ncov/3), 3))
    if (ncov <= 8) 
        par(mfrow = c(ceiling(ncov/2), 2))
    if (ncov <= 4) 
        par(mfrow = c(ncov, 1))
    for (j in 2:ncov) {
        plot(x$cov.id[, j], re[, 1], xlab = x$nameData$covnames[j - 
            1], ylab = paste("Random effects for log(", x$nameData$reparams[1], 
            ")", sep = ""))
        abline(h = 0)
    }
    dev.off()
}
    postscript(file = nameFile$file6, paper = "letter")
    diagtrplot(x, "i", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = pages) 
    dev.off()
    postscript(file = "diagtrplti.ps", paper = "letter")
    diagtrplot(x, "i", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = 1) 
    dev.off()
    postscript(file = nameFile$file7, paper = "letter")
    diagtrplot(x, "p", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = pages) 
    dev.off()
    postscript(file = "diagtrpltp.ps", paper = "letter")
    diagtrplot(x, "p", xvarlab = x$nameData$xvarlab, yvarlab = x$nameData$yvarlab, pages = 1) 
    dev.off()
    cnt <- 6 + length(descStructure$pcts)
    lp <- length(x$nameData$params)
    i <- 1
    j <- 1
    var.table <- desc(x$sims.list$itau[, i, j], pcts = descStructure$pcts, 
        nsig = descStructure$nsig)
    for (i in 1:lp) {
        for (j in 1:lp) {
            if (i == 1 && j == 1) 
                var.table <- desc(x$sims.list$itau[, i, j], pcts = descStructure$pcts, 
                  nsig = descStructure$nsig)
            else var.table <- rbind(var.table, desc(x$sims.list$itau[, 
                i, j], pcts = descStructure$pcts))
        }
    }
    row.names(var.table) <- x$nameData$varnames
    var.table <- data.frame(var.table)
    names(var.table) <- gsub("X", " ", names(var.table))
    var.table <- as.matrix(var.table)
    var.table <- formatC(var.table, format = "fg", digits = 3)

    sigma.table <- t(desc(x$sims.list$sigma2, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    sigma.table <- data.frame(sigma.table[,, drop=FALSE])
    names(sigma.table) <- gsub("X", " ", names(sigma.table))
    row.names(sigma.table) <- c("Sigma2")
    sigma.table <- as.matrix(cbind(sigma.table))
    sigma.table <- formatC(sigma.table, format = "fg", digits = 3)

    mu.table <- t(apply(x$sims.list$mu[, ], 2, desc, pcts = descStructure$pcts, 
        nsig = descStructure$nsig))
    mu.table <- data.frame(mu.table)
    names(mu.table) <- gsub("X", " ", names(mu.table))
    row.names(mu.table) <- x$nameData$coef
    mu.table <- as.matrix(mu.table)
    mu.table <- formatC(mu.table, format = "fg", digits = 3)
    tableList <- list()
    for (k in 1:lp) {
        tableList[[x$nameData$tparams[k]]] <- t(apply(x$sims.list$beta[, 
            , k], 2, desc, pcts = descStructure$pcts, nsig = descStructure$nsig))
    }
    for (k in 1:lp) {
        tableList[[x$nameData$tparams[k]]] <- data.frame(tableList[[x$nameData$tparams[k]]])
        names(tableList[[x$nameData$tparams[k]]]) <- gsub("X", 
            " ", names(tableList[[x$nameData$tparams[k]]]))
    }
    library(xtable)
    format.WB <- function(file = "wb.tex", file0 = "hist.ps", 
        file1 = "trplt.wb.ps", file2 = "diagplt.wb.ps", file3 = "qqploti.wb.ps", 
        file4 = "qqnormre.wb.ps", file5 = "covre.wb.ps", file6 = "diagtrplti.wb.ps", 
        file7 = "diagtrpltp.wb.ps") {
        cat("\\documentclass{article}\n", file = file)
        cat("\\usepackage{epsfig}\n", file = file, append = TRUE)
        cat("\\textwidth=6in\n", file = file, append = TRUE)
        cat("\\textheight=8.5in\n", file = file, append = TRUE)
        cat("\\oddsidemargin=.1in\n", file = file, append = TRUE)
        cat("\\evensidemargin=.1in\n", file = file, append = TRUE)
        cat("\\headheight=-.5in\n", file = file, append = TRUE)
        cat("\\begin{document}\n", file = file, append = TRUE)
        cat("\\title{Results from WinBUGS}\n", file = file, append = TRUE)
        cat("\\maketitle\n", file = file, append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
        print(xtable(mu.table, "Population Mean Estimates"), 
            file = file, append = TRUE)
        print(xtable(var.table, "Population Variance Estimates"), 
            file = file, append = TRUE)
        print(xtable(t(sigma.table), "Intrasubject Variance"), 
            file = file, append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file0, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Densities of the population parameter estimates}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        for (k in 1:lp) {
            mn.ind <- min(50, data$data$n.ind)
            print(xtable(tableList[[x$nameData$tparams[k]]][1:mn.ind, 
                ], paste("Individual(", x$nameData$tparams[k], 
                ")", sep = "")), file = file, append = TRUE)
        }
        cat("\\clearpage\n", file = file, append = TRUE)
        if (!is.null(x$cov.id)) print(xtable(cov.desc, "Covariates"), file = file, 
            append = TRUE)
        cat("\\clearpage\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file1, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Trellis plot of Concentration vs Time Data}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file2, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Predicted vs Observed and Residual vs Predicted Plots}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file3, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{QQPlots: Stage I Model, QQplots provide a visual assessment for normality.}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file4, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{QQPlots: Stage II Model}\n", file = file, 
            append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
if (!is.null(x$cov)){
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file5, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Individual Random effects vs Covariates}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
}
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file6, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Individual Predicted and Observed Values vs Time}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\begin{figure}\n", file = file, append = TRUE)
        cat(paste("\\includegraphics[height=6in,angle=270]{", 
            file7, "}\n", sep = ""), file = file, append = TRUE)
        cat("\\caption{Population Predicted and Observed Values vs Time}\n", 
            file = file, append = TRUE)
        cat("\\end{figure}\n", file = file, append = TRUE)
        cat("\\end{document}\n", file = file, append = TRUE)
        cat(paste("\ncompleted export to", file))
    }
    format.WB(file = nameFile$file, file0 = nameFile$file0, file1 = "trplt.ps", 
        file2 = nameFile$file2, file3 = nameFile$file3, file4 = nameFile$file4, 
        file5 = "re1.ps", file6 = "diagtrplti.ps", file7 = "diagtrpltp.ps")
#    return(list(WB = WB))
}
 


tex <- function (x=x, nameDir=nameDir, nameFile = nameFile, descStructure = descStructure,...) 
UseMethod("tex")

tex.default <- function (x=x, nameDir=nameDir, nameFile = nameFile, descStructure = descStructure,...) 
stop(paste("tex method does not yet exist for instances of class",
class(x)))    

HTMLtools.PKNLME <-function (x=x, nameDir=nameDir, nameFile = nameFile, descStructure=list(pcts=c(0.025,0.05,0.95,0.975),nsig=4),...) 
{ 
library(R2HTML)

cat ("Results written to", nameDir, "\n")

#Covariates
if (!is.null(x$cov.id)){
    cov.desc <- t(apply(x$cov.id[, ], 2, desc, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    cov.desc <- data.frame(cov.desc[-c(1),,drop=FALSE])
    names(cov.desc) <- gsub("X", " ", names(cov.desc))
    row.names(cov.desc) <- x$nameData$covnames
#    covdesc <- as.matrix(cov.desc)
#    covdesc <- formatC(covdesc, format = "fg", digits = 3)
}

#NLME parameters
    r <- summary(x$mm)
    r.tTable <- as.matrix(r$tTable)
    r.tTable <- formatC(r.tTable, format = "fg", digits = 3)
    sigma.nlme <- data.frame(r$sigma^2)
    names(sigma.nlme) <- c("sigma2")
    sigma.nlme <- as.matrix(sigma.nlme)
    sigma.nlme <- data.frame(formatC(sigma.nlme, format = "fg", digits = 3), row.names="")
    var.nlme <- getPsi(r)
    var.nlme <- as.matrix(var.nlme)
    var.nlme <- formatC(var.nlme, format = "fg", digits = 3)
    n <- nrow(data$pkvar)

#AIC table
    AIC.table <- data.frame(PKtools.AIC(loglike = logLik(x$mm), n = n, 
        K = attr(logLik(x$mm), "df")), row.names="")
    colnames(AIC.table) <- c("AIC", "AICc", "OF", "K")

#individual PK parameters
    pages = round(length(unique(data$pkvar$id))/16)
    coef.table<- t(apply(coef(x$mm)[, ], 2, desc, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    coef.table <- data.frame(coef.table[,,drop=FALSE])
    names(coef.table) <- gsub("X", " ", names(coef.table))

plotFun<-function(x, f) f(x)

plot.TrPre<- function(x) function(){
    id <- x$pkdata$id
    CONC <- x$pkdata$conc
    TIME <- x$pkdata$time
    dat <- data.frame(cbind(id, CONC, TIME))
    library(lattice)
    ID <- as.factor(dat$id)
    print(xyplot(CONC ~ TIME | ID, data = dat, xlab = x$nameData$xvarlab, xlim = c(-1, 
        max(TIME)), ylab = x$nameData$yvarlab, ylim = c(min(CONC), max(CONC)), 
        span = 2, layout = c(4, 4, pages=1), aspect = 1, panel = function(x, 
            y, span) {
            panel.xyplot(x, y, col = "black", cex = 1.2, type = "b")
        }))
}
plot.Tr = plot.TrPre(x)

qplot.Pre<-function(x) function(){
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 5), bg="grey")
    qqnorm(x$mm$residuals[, 1], main = "Residuals from population model")
    qqline(x$mm$residuals[, 1])
    qqnorm(x$mm$residuals[, 2], main = "Residuals from individual model")
    qqline(x$mm$residuals[, 2])}
qplot<-qplot.Pre(x)
    
qplot.rePre<-function(x) function(){
    re <- x$mm$coef$random$id
    nre <- ncol(re)
    par(mfrow = c(nre, 1), mar = c(4, 5, 3, 5), bg="grey")
    for (i in 1:nre) {
        qqnorm(re[, i], main = paste("Random effects for log(", 
            x$nameData$reparams[i], ")", sep = ""))
        qqline(re[, i])
    }
}
qplot.re<-qplot.rePre(x)

if (!is.null(x$cov.id)){      
re.covPre<-function(x)function(){ 
    re <- x$mm$coef$random$id
    nre <- ncol(re)
    ncov <- ncol(x$cov.id)
    par(bg = "grey")
    if (ncov-1<=16)
       par(mfrow=c(ceiling((ncov-1)/4),4))
    if (ncov-1<=12)
       par(mfrow=c(ceiling((ncov-1)/3),4))
    if (ncov-1<=8)
       par(mfrow=c(ceiling((ncov-1)/2),4))
    if (ncov-1<=4)
       par(mfrow=c(ceiling((ncov-1)),1))
    for (i in 1:nre){
    for (j in 2:ncov) {
            plot(x$cov.id[, j], re[, i], xlab = x$nameData$covnames[j - 
                1], ylab = x$nameData$reparams[i])
            abline(h = 0)
        }
   }
   } 
re.cov<-re.covPre(x)   
   }

dplot.Pre<-function(x) function() diagplot(x)
dplot<-dplot.Pre(x)

diagtrploti.Pre<- function (x)  function()
{
    ID <- x$pkdata$id
    CONC <- x$pkdata$conc
    TIME <- x$pkdata$time
    PRED <- x$mm$fitted[, 1]
    IPRE <- x$mm$fitted[, 2]
    RES <- x$mm$resid[, 1]
    IRES <- x$mm$resid[, 2]
    dat <- data.frame(cbind(ID, CONC, TIME, PRED, IPRE, RES, 
        IRES))
    library(lattice)
    typedv <- rep(1, length(dat$ID))
    typepred <- rep(2, length(dat$ID))
    dvdata <- cbind(dat$ID, dat$TIME, dat$CONC, typedv)
    #    preddata <- cbind(dat$ID, dat$TIME, dat$PRED, typepred)
    preddata <- cbind(dat$ID, dat$TIME, dat$IPRE, typepred)
    plotdata <- as.data.frame(rbind(dvdata, preddata))
    names(plotdata) <- c("id", "time", "dv.pred", "typepred")
    newid <- factor(plotdata$id)
    sps <- trellis.par.get("superpose.symbol")
    sps$pch <- 1:7
    trellis.par.set("superpose.symbol", sps)
    print(xyplot(dv.pred ~ time | newid, data = plotdata, groups = typepred, 
        xlab = x$nameData$xvarlab, ylab = x$nameData$yvarlab, span = 2, layout = c(4, 
            4, pages=1), aspect = 1, type = "b", panel = "panel.superpose", 
        panel.group = "panel.xyplot", key = list(columns = 2, 
            text = list(paste(c("Values:  Observed", " Predicted"))), 
            points = Rows(sps, 1:2))))
}
diagtrploti<-diagtrploti.Pre(x)

diagtrplotp.Pre<- function (x)  function()
{
    ID <- x$pkdata$id
    CONC <- x$pkdata$conc
    TIME <- x$pkdata$time
    PRED <- x$mm$fitted[, 1]
    IPRE <- x$mm$fitted[, 2]
    RES <- x$mm$resid[, 1]
    IRES <- x$mm$resid[, 2]
    dat <- data.frame(cbind(ID, CONC, TIME, PRED, IPRE, RES, 
        IRES))
    library(lattice)
    typedv <- rep(1, length(dat$ID))
    typepred <- rep(2, length(dat$ID))
    dvdata <- cbind(dat$ID, dat$TIME, dat$CONC, typedv)
    preddata <- cbind(dat$ID, dat$TIME, dat$PRED, typepred)
    #preddata <- cbind(dat$ID, dat$TIME, dat$IPRE, typepred)
    plotdata <- as.data.frame(rbind(dvdata, preddata))
    names(plotdata) <- c("id", "time", "dv.pred", "typepred")
    newid <- factor(plotdata$id)
    sps <- trellis.par.get("superpose.symbol")
    sps$pch <- 1:7
    trellis.par.set("superpose.symbol", sps)
    print(xyplot(dv.pred ~ time | newid, data = plotdata, groups = typepred, 
        xlab = x$nameData$xvarlab, ylab = x$nameData$yvarlab, span = 2, layout = c(4, 
            4, pages=1), aspect = 1, type = "b", panel = "panel.superpose", 
        panel.group = "panel.xyplot", key = list(columns = 2, 
            text = list(paste(c("Values:  Observed", " Predicted"))), 
            points = Rows(sps, 1:2))))
}
diagtrplotp<-diagtrplotp.Pre(x)

HTMLInitFile(outdir=nameDir, file=nameFile$file, extension="html", BackGroundColor="#E9967A")
 file.app<-paste(nameFile$file,".html",sep="")
         HTML.title("Results from NLME", HR=3,file=file.app)
         HTML(r$tTable, caption="Mean Population Parameter Estimates",file = file.app)
         HTML(var.nlme,caption="Intersubject Variance-Covariance Matrix", file = file.app)
         HTML(sigma.nlme, caption="Intrasubject Variance", file = file.app)
         HTML(AIC.table, caption="Model Selection Criteria", file = file.app)
         HTML(coef.table, caption="Descriptive statistics for Individual Parameters", file=file.app)
             if (!is.null(x$cov.id)) HTML(cov.desc, caption="Covariates",file = file.app)
         HTMLplot(Caption="Trellis Plots", plotFunction=plot.Tr, GraphFileName="plot.Tr", GraphSaveAs="png", append=TRUE, file=file.app)    
         HTMLplot(Caption="Predicted vs Observed and Residual vs Predicted Plots", plotFunction=dplot, GraphFileName="dplot", GraphSaveAs="png", append=TRUE, file=file.app)
         HTMLplot(Caption="QQplots:Stage I Model, QQplots provide a visual assessment of normality", plotFunction=qplot, GraphFileName="qplot", GraphSaveAs="png", append=TRUE, file=file.app)
         HTMLplot(Caption="QQPlots: Stage II Model", plotFunction=qplot.re, GraphFileName="qplot.re", GraphSaveAs="png", append=TRUE, file=file.app)
 if (!is.null(x$cov.id)) HTMLplot(Caption="Individual Random effects vs Covariates", plotFunction=re.cov, GraphFileName="re.cov", GraphSaveAs="png", append=TRUE, file=file.app)        
         HTMLplot(Caption ="Individual Predicted and Observed Values vs Time",plotFunction=diagtrploti,GraphFileName="diagtrploti", GraphSaveAs="png", append=TRUE, file=file.app)    
         HTMLplot(Caption = "Population Predicted and Observed Values vs Time",plotFunction=diagtrplotp,GraphFileName="diagtrplotp", GraphSaveAs="png", append=TRUE, file=file.app)    
HTMLEndFile()
     }

HTMLtools.NONMEM <-function (x = x, nameDir = nameDir, nameFile = nameFile, descStructure = list(pcts = c(0.025, 0.05, 0.95, 0.975),
        nsig = 4),...){
    library(R2HTML)
    cat("Results written to", nameDir, "\n")
    n <- nrow(x$pkdat$dat)
    pages <- round(length(unique(x$pkdat$dat$id))/16)
    K <- sum(x$param[, 1] != 0)
    loglike <- -0.5 * x$of
    NM.param <- data.frame(x$param)
    names(NM.param) <- c("Estimate", "Standard Error")
    rownames(NM.param) <- rbind(cbind(x$nameData$tparams), cbind(x$nameData$varnames),
        "sigma2")
    NMparam <- as.matrix(NM.param)
    NMparam <- formatC(NMparam, format = "fg", digits = 3)
         
    AIC.table <- data.frame(PKtools.AIC(loglike = loglike, n = n,
        K = K), row.names = "")
    if (!is.null(x$cov)) {
        cov.desc <- t(apply(x$cov[, ], 2, desc, pcts = descStructure$pcts,
            nsig = descStructure$nsig))
        cov.desc <- data.frame(cov.desc[-c(1), , drop = FALSE])
        names(cov.desc) <- gsub("X", " ", names(cov.desc))
        row.names(cov.desc) <- x$nameData$covnames
    }
    ip.desc <- t(apply(x$ip[, ], 2, desc, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    ip.desc <- data.frame(ip.desc[-c(1), , drop = FALSE])
    names(ip.desc) <- gsub("X", " ", names(ip.desc))
    
    plotFun<-function(x, f) f(x)
        
    plot.TrPre<- function(x) function(){
    dat <- x$pred
    library(lattice)
    Id <- as.factor(dat$ID)
    print(xyplot(CONC ~ TIME | Id, data = dat, xlab =x$nameData$xvarlab, xlim = c(-1, 
        max(time)), ylab = x$nameData$yvarlab, ylim = c(min(dat$CONC), max(dat$CONC)), 
        span = 2, layout = c(4, 4, 1), aspect = 1, panel = function(x, 
            y, span) {
            panel.xyplot(x, y, col = "black", cex = 1.2, type = "b")
        }))
    }
    plot.Tr = plot.TrPre(x)
    
dplot.Pre<-function(x) function() diagplot(x)
dplot<-dplot.Pre(x)

   
qplot.Pre<-function(x) function(){
    par(mfrow = c(2, 1), mar = c(4, 5, 3, 5), bg="grey") 
    qqnorm(x$pred$WRES, main = "Weighted residuals from the population model")
    qqline(x$pred$WRES)
    qqnorm(x$pred$IWRE, main = "Weighted residuals from the individual model")
    qqline(x$pred$IWRE)
    }
qplot<-qplot.Pre(x)
    
    
qplot.rePre<-function(x) function(){
    re <- x$re
    nre <- ncol(re)
    par(mfrow = c(nre - 1, 1), mar = c(4, 2, 3, 2), bg = "grey")
    for (i in 2:nre) {
        qqnorm(re[, i], main = paste("Random effects for log(",
            x$nameData$reparams[i - 1], ")", sep = ""))
        qqline(re[, i])
    }
}    
qplot.re<-qplot.rePre(x)    


    if (!is.null(x$cov)) {
re.covPre<-function(x)function(){             
        par(bg = "grey")
        ncov <- ncol(x$cov)
        re <- x$re
        nre <- ncol(re)
        if (ncov-1<=16)
           par(mfrow=c(ceiling((ncov-1)/4),4))
        if (ncov-1<=12)
           par(mfrow=c(ceiling((ncov-1)/3),4))
        if (ncov-1<=8)
           par(mfrow=c(ceiling((ncov-1)/2),4))
        if (ncov-1<=4)
           par(mfrow=c(ceiling((ncov-1)),1))
        for (i in 2:nre) {
           for (j in 2:ncov) {         
                plot(x$cov[, j], re[, i], xlab = x$nameData$covnames[j -
                  1], ylab = x$nameData$reparams[i - 1])
                abline(h = 0)
            }
        }
      }
re.cov<-re.covPre(x)   
   }


diagtrplotp.Pre<- function (x)  function(){
    dat <- x$pred
    library(lattice)
    typedv <- rep(1, length(dat$ID))
    typepred <- rep(2, length(dat$ID))
    dvdata <- cbind(dat$ID, dat$TIME, dat$CONC, typedv) 
    preddata <- cbind(dat$ID, dat$TIME, dat$PRED, typepred)
    #preddata <- cbind(dat$ID, dat$TIME, dat$IPRE, typepred)
    plotdata <- as.data.frame(rbind(dvdata, preddata))
    names(plotdata) <- c("id", "time", "dv.pred", "typepred")
    newid <- factor(plotdata$id)
    sps <- trellis.par.get("superpose.symbol")
    sps$pch <- 1:7
    trellis.par.set("superpose.symbol", sps)
    print(xyplot(dv.pred ~ time | newid, data = plotdata, groups = typepred, 
        xlab = x$nameData$xvarlab, ylab = x$nameData$yvarlab, span = 2, layout = c(4, 
            4, pages), aspect = 1, type = "b", panel = "panel.superpose", 
        panel.group = "panel.xyplot", key = list(columns = 2, 
            text = list(paste(c("Values:  Observed", " Predicted"))), 
            points = Rows(sps, 1:2))))
}
diagtrplotp<-diagtrplotp.Pre(x)

diagtrploti.Pre<- function (x)  function(){
    dat <- x$pred
    library(lattice)
    typedv <- rep(1, length(dat$ID))
    typepred <- rep(2, length(dat$ID))
    dvdata <- cbind(dat$ID, dat$TIME, dat$CONC, typedv) 
    #preddata <- cbind(dat$ID, dat$TIME, dat$PRED, typepred)
    preddata <- cbind(dat$ID, dat$TIME, dat$IPRE, typepred)
    plotdata <- as.data.frame(rbind(dvdata, preddata))
    names(plotdata) <- c("id", "time", "dv.pred", "typepred")
    newid <- factor(plotdata$id)
    sps <- trellis.par.get("superpose.symbol")
    sps$pch <- 1:7
    trellis.par.set("superpose.symbol", sps)
    print(xyplot(dv.pred ~ time | newid, data = plotdata, groups = typepred, 
        xlab = x$nameData$xvarlab, ylab = x$nameData$yvarlab, span = 2, layout = c(4, 
            4, pages), aspect = 1, type = "b", panel = "panel.superpose", 
        panel.group = "panel.xyplot", key = list(columns = 2, 
            text = list(paste(c("Values:  Observed", " Predicted"))), 
            points = Rows(sps, 1:2))))
}
diagtrploti<-diagtrploti.Pre(x)

HTMLInitFile(outdir=nameDir, file=nameFile$file, extension="html", BackGroundColor="#E9967A")
 file.app<-paste(nameFile$file,".html",sep="")
         HTML.title("Results from NONMEM", HR=3,file=file.app)
         HTML(NMparam, caption="Population Parameter Estimates",file = file.app)
         HTML(AIC.table, caption="Model Selection Criteria", file = file.app)
         HTML(ip.desc, caption = "Descriptive statistics for Individual Parameters",file = file.app)
    if (!is.null(x$cov))
        HTML(cov.desc, caption = "Covariates", file = file.app)
         HTMLplot(Caption="Trellis Plots", plotFunction=plot.Tr, GraphFileName="plot.Tr", GraphSaveAs="png", append=TRUE, file=file.app)    
         HTMLplot(Caption="Predicted vs Observed and Residual vs Predicted Plots", plotFunction=dplot, GraphFileName="dplot", GraphSaveAs="png", append=TRUE, file=file.app)
         HTMLplot(Caption="QQplots:Stage I Model, QQplots provide a visual assessment of normality", plotFunction=qplot, GraphFileName="qplot", GraphSaveAs="png", append=TRUE, file=file.app)
         HTMLplot(Caption="QQPlots: Stage II Model", plotFunction=qplot.re, GraphFileName="qplot.re", GraphSaveAs="png", append=TRUE, file=file.app)
 if (!is.null(x$cov.id)) HTMLplot(Caption="Individual Random effects vs Covariates", plotFunction=re.cov, GraphFileName="re.cov", GraphSaveAs="png", append=TRUE, file=file.app)        
         HTMLplot(Caption ="Individual Predicted and Observed Values vs Time",plotFunction=diagtrploti,GraphFileName="diagtrploti", GraphSaveAs="png", append=TRUE, file=file.app)    
         HTMLplot(Caption = "Population Predicted and Observed Values vs Time",plotFunction=diagtrplotp,GraphFileName="diagtrplotp", GraphSaveAs="png", append=TRUE, file=file.app)    
HTMLEndFile()
     }


HTMLtools.WinBUGS<-function (x = x, nameDir = nameDir, nameFile = nameFile,
    descStructure = list(pcts = c(0.025, 0.05, 0.95, 0.975),
        nsig = 4), ...)
{
    library(R2HTML)
    cat("Results written to", nameDir, "\n")
    pages <- round(length(unique(x$pred$id))/16)
    if (!is.null(x$cov.id)) {
        cov.desc <- t(apply(x$cov.id[, ], 2, desc, pcts = descStructure$pcts,
            nsig = descStructure$nsig))
        cov.desc <- data.frame(cov.desc[-c(1), , drop = FALSE])
        names(cov.desc) <- gsub("X", " ", names(cov.desc))
        row.names(cov.desc) <- x$nameData$covnames
    }
    
    cnt <- 6 + length(descStructure$pcts)
    lp <- length(x$nameData$params)
    i <- 1
    j <- 1
    var.table <- desc(x$sims.list$itau[, i, j], pcts = descStructure$pcts,
        nsig = descStructure$nsig)
    for (i in 1:lp) {
        for (j in 1:lp) {
            if (i == 1 && j == 1)
                var.table <- desc(x$sims.list$itau[, i, j], pcts = descStructure$pcts,
                  nsig = descStructure$nsig)
            else var.table <- rbind(var.table, desc(x$sims.list$itau[,
                i, j], pcts = descStructure$pcts))
        }
    }
    row.names(var.table) <- x$nameData$varnames
    var.table <- data.frame(var.table)
    names(var.table) <- gsub("X", " ", names(var.table))
    var.table <- as.matrix(var.table)
    var.table <- formatC(var.table, format = "fg", digits = 3)
    sigma.table <- t(desc(x$sims.list$sigma2, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    sigma.table <- data.frame(sigma.table[, , drop = FALSE])
    names(sigma.table) <- gsub("X", " ", names(sigma.table))
    row.names(sigma.table) <- c("Sigma2")
    sigma.table <- as.matrix(cbind(sigma.table))
    sigma.table <- formatC(sigma.table, format = "fg", digits = 3)
    mu.table <- t(apply(x$sims.list$mu[, ], 2, desc, pcts = descStructure$pcts,
        nsig = descStructure$nsig))
    mu.table <- data.frame(mu.table)
    names(mu.table) <- gsub("X", " ", names(mu.table))
    row.names(mu.table) <- x$nameData$coef
    mu.table <- as.matrix(mu.table)
    mu.table <- formatC(mu.table, format = "fg", digits = 3)
    tableList <- list()
    for (k in 1:lp) {
        tableList[[x$nameData$tparams[k]]] <- t(apply(x$sims.list$beta[,
            , k], 2, desc, pcts = descStructure$pcts, nsig = descStructure$nsig))
    }
    for (k in 1:lp) {
        tableList[[x$nameData$tparams[k]]] <- data.frame(tableList[[x$nameData$tparams[k]]])
        names(tableList[[x$nameData$tparams[k]]]) <- gsub("X",
            " ", names(tableList[[x$nameData$tparams[k]]]))
    }
        
plotFun<-function(x, f) f(x)

densities.Pre<-function(x) function(){
count <- ncol(x$sims.list$mu)
    if (count <= 16)
        par(mfrow = c(ceiling(count/4), 4), mar = c(4, 2, 3,
            2), bg = "grey")
    if (count <= 12)
        par(mfrow = c(ceiling(count/3), 3), mar = c(4, 2, 3,
            2), bg = "grey")
    if (count <= 8)
        par(mfrow = c(ceiling(count/2), 2), mar = c(4, 2, 3,
            2), bg = "grey")
    if (count <= 4)
        par(mfrow = c(count, 1), mar = c(4, 2, 3, 2), bg = "grey")
    for (i in 1:count) {
        plot(density(x$sims.list$mu[, i]), main = x$nameData$coef[i])
        abline(v = x$mean$mu[i])
    }
}

plot.densities<-densities.Pre(x)

plot.TrPre<- function(x) function(){
    library(lattice)
    dat <- x$pred
    ID <- as.factor(dat$id)
    print(xyplot(conc ~ time | ID, data = dat, xlab = x$nameData$xvarlab, xlim = c(-1, 
        max(time)), ylab = x$nameData$yvarlab, ylim = c(min(dat$conc), max(dat$conc)), 
        span = 2, layout = c(4, 4, 1), aspect = 1, panel = function(x, 
            y, span) {
            panel.xyplot(x, y, col = "black", cex = 1.2, type = "b")
        }))
}
plot.Tr = plot.TrPre(x)

qplot.Pre<-function(x) function(){
  par(mfrow = c(2, 1), mar = c(4, 5, 3, 5), bg = "grey")
    qqnorm(x$pred$presid, main = "Residuals from population model")
    qqline(x$pred$presid)
    qqnorm(x$pred$iresid, main = "Residuals from individual model")
    qqline(x$pred$iresid)
}
qplot<-qplot.Pre(x)

qplot.rePre<-function(x) function(){
    re <- x$mean$re
    nre <- ncol(re)
    par(mfrow = c(nre, 1), mar = c(4, 5, 3, 5), bg="grey")
    for (i in 1:nre) {
        qqnorm(re[, i], main = paste("Random effects for log(",
            x$nameData$reparams[i], ")", sep = ""))
        qqline(re[, i])
    }
}
qplot.re<-qplot.rePre(x)

if (!is.null(x$cov.id)){
re.covPre<-function(x)function(){
    par(bg="grey")
    re <- x$mean$re
    nre <- ncol(re)
    ncov <- ncol(x$cov.id)
    if (ncov-1<=16)
       par(mfrow=c(ceiling((ncov-1)/4),4))
    if (ncov-1<=12)
       par(mfrow=c(ceiling((ncov-1)/3),4))
    if (ncov-1<=8)
       par(mfrow=c(ceiling((ncov-1)/2),4))
    if (ncov-1<=4)
       par(mfrow=c(ceiling((ncov-1)),1))
    for (i in 1:nre) {
       for (j in 2:ncov) {  
            plot(x$cov.id[, j], re[, i], xlab = x$nameData$covnames[j -
                1], ylab = x$nameData$reparams[i])
            abline(h = 0)
        }
   }
   }
re.cov<-re.covPre(x)
   }

dplot.Pre<-function(x) function() diagplot(x)
dplot<-dplot.Pre(x)

diagtrploti.Pre<- function (x)  function()
{
  trdata <- x$pred
    library(lattice)
    typedv <- rep(1, length(trdata$id))
    typepred <- rep(2, length(trdata$id))
    dvdata <- cbind(trdata$id, trdata$time, trdata$conc, typedv)
    preddata <- cbind(trdata$id, trdata$time, trdata$ipred, typepred)
    plotdata <- as.data.frame(rbind(dvdata, preddata))
    names(plotdata) <- c("id", "time", "dv.pred", "typepred")
    newid <- factor(plotdata$id)
    sps <- trellis.par.get("superpose.symbol")
    sps$pch <- 1:7
    trellis.par.set("superpose.symbol", sps)
    print(xyplot(dv.pred ~ time | newid, data = plotdata, groups = typepred, 
        xlab = x$nameData$xvarlab, ylab = x$nameData$yvarlab, span = 2, layout = c(4, 
            4, pages), aspect = 1, type = "b", panel = "panel.superpose", 
        panel.group = "panel.xyplot", key = list(columns = 2, 
            text = list(paste(c("Values:  Observed", " Predicted"))), 
            points = Rows(sps, 1:2))))
}
diagtrploti<-diagtrploti.Pre(x)

diagtrplotp.Pre<- function (x)  function()
{
  trdata <- x$pred
    library(lattice)
    typedv <- rep(1, length(trdata$id))
    typepred <- rep(2, length(trdata$id))
    dvdata <- cbind(trdata$id, trdata$time, trdata$conc, typedv)
    preddata <- cbind(trdata$id, trdata$time, trdata$ppred, typepred)
    plotdata <- as.data.frame(rbind(dvdata, preddata))
    names(plotdata) <- c("id", "time", "dv.pred", "typepred")
    newid <- factor(plotdata$id)
    sps <- trellis.par.get("superpose.symbol")
    sps$pch <- 1:7
    trellis.par.set("superpose.symbol", sps)
    print(xyplot(dv.pred ~ time | newid, data = plotdata, groups = typepred, 
        xlab = x$nameData$xvarlab, ylab = x$nameData$yvarlab, span = 2, layout = c(4, 
            4, pages), aspect = 1, type = "b", panel = "panel.superpose", 
        panel.group = "panel.xyplot", key = list(columns = 2, 
            text = list(paste(c("Values:  Observed", " Predicted"))), 
            points = Rows(sps, 1:2))))
}
diagtrplotp<-diagtrplotp.Pre(x)

HTMLInitFile(outdir=nameDir, file=nameFile&file, extension="html", BackGroundColor="#E9967A")
 file.app<-paste(nameFile$file,".html",sep="")
 HTML.title("Results from WinBUGS", HR = 3, file = file.app)
    HTML(mu.table, caption = "Population Mean Estimates", file = file.app)
    HTML(var.table, caption = "Population Variance Estimates",
        file = file.app)
    HTML(sigma.table, caption = "Intrasubject Variance Estimate",
        file = file.app)
    for (k in 1:lp) {
        mn.ind <- min(50, data$data$n.ind)
        HTML(tableList[[x$nameData$tparams[k]]][1:mn.ind, ], caption = paste("Individual(",
            x$nameData$tparams[k], ")", sep = ""), file = file.app)
    }
    if (!is.null(x$cov.id))
        HTML(cov.desc, caption = "Covariates", file = file.app)
HTMLplot(Caption = "Densities of the population parameter estimates",plotFunction=plot.densities, GraphFileName="plot.densities", GraphSaveAs="png", append=TRUE, 
         file=file.app)
HTMLplot(Caption="Trellis Plots", plotFunction=plot.Tr, GraphFileName="plot.Tr", GraphSaveAs="png", append=TRUE, 
         file=file.app)
HTMLplot(Caption="Predicted vs Observed and Residual vs Predicted Plots", plotFunction=dplot, GraphFileName="dplot", 
         GraphSaveAs="png", append=TRUE, file=file.app)
HTMLplot(Caption="QQplots:Stage I Model, QQplots provide a visual assessment of normality", plotFunction=qplot, 
         GraphFileName="qplot", GraphSaveAs="png", append=TRUE, file=file.app)
HTMLplot(Caption="QQPlots: Stage II Model", plotFunction=qplot.re, GraphFileName="qplot.re", GraphSaveAs="png", 
         append=TRUE, file=file.app)
if (!is.null(x$cov.id)) HTMLplot(Caption="Individual Random effects vs Covariates", plotFunction=re.cov, GraphFileName="re.cov", 
 GraphSaveAs="png", append=TRUE, file=file.app)
HTMLplot(Caption ="Individual Predicted and Observed Values vs Time",plotFunction=diagtrploti,GraphFileName="diagtrploti", 
         GraphSaveAs="png", append=TRUE, file=file.app)
HTMLplot(Caption = "Population Predicted and Observed Values vs Time",plotFunction=diagtrplotp,GraphFileName="diagtrplotp", GraphSaveAs="png", append=TRUE, file=file.app)
HTMLEndFile()    
}
    
HTMLtools<-function (x = x, nameDir= nameDir, nameFile = nameFile, descStructure = descStructure,...)
UseMethod("HTMLtools")

HTMLtools.default <- function (x = x, nameDir= nameDir, 
nameFile = nameFile, descStructure = descStructure,...)
stop(paste("HTMLtools method does not yet exist for instances of class",
class(x)))


#####################################################
# getPsi function provided by Jose Pinheiro
# permission to include in runPK (renamed PKtools)
# given 2004 no copyright claim by the present author
#####################################################                                                                                


getPsi <- function (object, level = Q)
{
    val <- as.matrix(object$modelStruct$reStruct)
    Q <- length(val)
    if (length(level) == 1) {
        val <- object$sigma^2 * val[[level]]
    }
    else {
        val <- lapply(val[level], function(x, sig2) sig2 * x,
            sig2 = object$sigma^2)
    }
    val
}

#####################################################
# interface to R by A.Gelman 
# permission to include in runPK given 2004
# no copyright claim by the present author
#####################################################                                                                                


attach.all <- function (a, overwrite = TRUE, name = "attach.all")
{
    if (overwrite) {
        for (j in 1:length(a)) {
            if (names(a)[j] %in% ls(.GlobalEnv))
                remove(list = names(a)[j], envir = .GlobalEnv)
        }
    }
    attach(a, name = name)
}
attach.bugs <- function (a)
{
    attach.all(a, name = "bugs.all")
    attach.all(a$sims.list, name = "bugs.sims")
}
bugs <- function (data, inits, parameters.to.save, model.file = "model.bug",
    n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter/2),
    n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/1000)),
    debug = FALSE, attach.sims = TRUE, print.summary = TRUE,
    plot.summary = TRUE, digits.summary = 1, display.parallel = FALSE,
    DIC = TRUE, bugs.directory = "c:/Program Files/WinBUGS14/",
    dos.location = "c:/progra~1/winbug~1/winbug~1")
{
    start.time <- Sys.time()
    bugs.data.inits(data, inits, n.chains, digits = 5)
    if (DIC)
        parameters.to.save <- c(parameters.to.save, "deviance")
    bugs.script(parameters.to.save, n.chains, n.iter, n.burnin,
        n.thin, bugs.directory, model.file, debug = debug)
    bugs.run(n.burnin, bugs.directory, dos.location)
    sims <- bugs.sims(parameters.to.save, n.chains, n.iter, n.burnin,
        n.thin, attach.sims, DIC)
    bugs.plot(sims, model.file, n.chains, n.iter, n.burnin, n.thin,
        print.summary, plot.summary, digits.summary, display.parallel,
        DIC, start.time)
    return(sims)
}
bugs.data.inits <- function (data, inits, n.chains, digits)
{
    if (is.numeric(unlist(data)))
        write.datafile(round.bugs.list(data, digits), "data.txt")
    else {
        data.list <- lapply(data, get, pos = 1)
        names(data.list) <- data
        write.datafile(round.bugs.list(data.list, digits), "data.txt")
    }
    for (i in 1:n.chains) {
        if (is.function(inits))
            write.datafile(round.bugs.list(inits(), digits),
                paste("inits", i, ".txt", sep = ""))
        else write.datafile(round.bugs.list(inits[[i]], digits),
            paste("inits", i, ".txt", sep = ""))
    }
}
bugs.plot <- function (sims, model.file, n.chains, n.iter, n.burnin, n.thin,
    print.summary = TRUE, plot.summary = TRUE, digits.summary,
    display.parallel = FALSE, DIC = TRUE, start.time)
{
    working.directory <- gsub("\\\\", "/", paste(getwd(), "/",
        sep = ""))
    model <- paste(working.directory, model.file, sep = "")
    end.time <- Sys.time()
    if (print.summary) {
        cat("Inference for Bugs model at \"", model, "\"\n ",
            n.chains, " chains, each with ", n.iter, " iterations (first ",
            n.burnin, " discarded)", sep = "")
        if (n.thin > 1)
            cat(", n.thin =", n.thin)
        cat("\n n.sims =",sims$n.sims, "iterations saved\n")
        print(difftime(end.time, start.time))
        print(round(sims$summary, digits.summary))
        if (!is.null(sims$DIC)) {
            pD <- sims$pD
            DIC <- sims$DIC
            cat(" pD =", round(pD, 1), "and DIC =", round(DIC,
                1), "(using the rule, pD = var(deviance)/2)\n")
            cat("\n For each parameter, n.eff is a crude measure of effective sample size,\n and Rhat is the potential scale reduction factor (at convergence, Rhat=1).\n DIC is an estimate of expected predictive error (lower deviance is better).\n")
        }
        flush.console()
    }
    if (plot.summary) {
        mar.old <- par("mar")
        par(pty = "m")
        layout(matrix(c(1, 2), 1, 2))
        bugs.plot.summary(sims, DIC)
        bugs.plot.inferences(sims, display.parallel)
        mtext(paste("Bugs model at \"", model, "\", ", n.chains,
            " chains, each with ", n.iter, " iterations", sep = ""),
            outer = TRUE, line = -1, cex = 0.7)
        if (.Device == "windows")
            bringToTop()
        par(mar = mar.old)
    }
}
bugs.plot.inferences <- function (sims, display.parallel)
{
    if (.Device == "windows" | (.Device == "null device" & options("device") ==
        "windows")) {
        cex.names <- 0.7
        cex.axis <- 0.6
        cex.tiny <- 0.4
        cex.points <- 0.7
        standard.width <- 30
        max.width <- 40
        min.width <- 0.02
    }
    else {
        cex.names <- 0.7
        cex.axis <- 0.6
        cex.tiny <- 0.4
        cex.points <- 0.3
        standard.width <- 30
        max.width <- 40
        min.width <- 0.01
    }
    rootnames <- sims$root.short
    n.roots <- length(rootnames)
    sims.array <- sims$sims.array
    n.chains <- sims$n.chains
    dimension.short <- sims$dimension.short
    indexes.short <- sims$indexes.short
    long.short <- sims$long.short
    height <- 0.6
    par(mar = c(0, 0, 1, 0))
    plot(c(0, 1), c(-n.roots - 0.5, -0.4), ann = FALSE, bty = "n",
        xaxt = "n", yaxt = "n", type = "n")
    W <- max(strwidth(rootnames, cex = cex.names))
    B <- (1 - W)/3.8
    A <- 1 - 3.5 * B
    if (display.parallel)
        text(A, -0.4, "80% interval for each chain", adj = 0,
            cex = cex.names)
    else text(A, -0.4, "medians and 80% intervals", adj = 0,
        cex = cex.names)
    num.height <- strheight(1:9, cex = cex.tiny)
    for (k in 1:n.roots) {
        text(0, -k, rootnames[k], adj = 0, cex = cex.names)
        J <- min(length(long.short[[k]]), max.width)
        if (k == 1)
            index <- 1:J
        else index <- sum(unlist(lapply(long.short, length))[1:(k -
            1)]) + 1:J
        spacing <- 3.5/max(J, standard.width)
        med <- rep(NA, J)
        i80 <- array(NA, c(J, 2))
        med.chains <- array(NA, c(J, sims$n.chains))
        i80.chains <- array(NA, c(J, sims$n.chains, 2))
        for (j in 1:J) {
            med[j] <- median(sims.array[, , index[j]])
            i80[j, ] <- quantile(sims.array[, , index[j]], c(0.1,
                0.9))
            for (m in 1:n.chains) {
                med.chains[j, m] <- quantile(sims.array[, m,
                  index[j]], 0.5)
                i80.chains[j, m, ] <- quantile(sims.array[, m,
                  index[j]], c(0.1, 0.9))
            }
        }
        rng <- range(i80, i80.chains)
        p.rng <- pretty(rng, n = 2)
        b <- height/(max(p.rng) - min(p.rng))
        a <- -(k + height/2) - b * p.rng[1]
        lines(A + c(0, 0), -k + c(-height/2, height/2))
        if (min(p.rng) < 0 & max(p.rng) > 0)
            lines(A + B * spacing * c(0, J + 1), rep(a, 2), lwd = 0.5,
                col = "gray")
        for (x in p.rng) {
            text(A - B * 0.2, a + b * x, x, cex = cex.axis)
            lines(A + B * c(-0.05, 0), rep(a + b * x, 2))
        }
        for (j in 1:J) {
            if (display.parallel) {
                for (m in 1:n.chains) {
                  interval <- a + b * i80.chains[j, m, ]
                  if (interval[2] - interval[1] < min.width)
                    interval <- mean(interval) + c(-1, 1) * min.width/2
                  lines(A + B * spacing * rep(j + 0.6 * (m -
                    (n.chains + 1)/2)/n.chains, 2), interval,
                    lwd = 0.5, col = m + 1)
                }
            }
            else {
                lines(A + B * spacing * rep(j, 2), a + b * i80[j,
                  ], lwd = 0.5)
                for (m in 1:n.chains) points(A + B * spacing *
                  j, a + b * med.chains[j, m], pch = 20, cex = cex.points,
                  col = m + 1)
            }
            dk <- dimension.short[k]
            if (dk > 0) {
                for (m in 1:dk) {
                  index0 <- indexes.short[[k]][[j]][m]
                  if (j == 1)
                    text(A + B * spacing * j, -k - height/2 -
                      0.05 - num.height * (m - 1), index0, cex = cex.tiny)
                  else if (index0 != indexes.short[[k]][[j -
                    1]][m] & (index0%%(floor(log10(index0) +
                    1)) == 0))
                    text(A + B * spacing * j, -k - height/2 -
                      0.05 - num.height * (m - 1), index0, cex = cex.tiny)
                }
            }
        }
        if (J < length(long.short[[k]]))
            text(-0.015, -k, "*", cex = cex.names, col = "red")
    }
}
bugs.plot.summary <- function (sims, DIC)
{
    if (.Device == "windows" | (.Device == "null device" & options("device") ==
        "windows")) {
        cex.names <- 0.7
        cex.top <- 0.7
        cex.points <- 0.7
        max.length <- 50
        min.width <- 0.01
    }
    else {
        cex.names <- 0.7
        cex.top <- 0.7
        cex.points <- 0.3
        max.length <- 80
        min.width <- 0.005
    }
    summ <- sims$summary
    sims.array <- sims$sims.array
    n.chains <- sims$n.chains
    n.parameters <- nrow(summ)
    J0 <- unlist(lapply(sims$long.short, length))
    if (DIC)
        J0 <- J0[1:(length(J0) - 1)]
    J <- J0
    total <- ceiling(sum(J + 0.5))
    while ((total > max.length) & max(J) > 1) {
        J[J == max(J)] <- max(J) - 1
        total <- ceiling(sum(J + 0.5))
    }
    pos <- -1
    ypos <- NULL
    id <- NULL
    ystart <- NULL
    jj <- 1:J[1]
    n.roots <- length(sims$root.short)
    if (DIC)
        n.roots <- n.roots - 1
    for (k in 1:n.roots) {
        ystart <- c(ystart, pos)
        for (j in 1:J[k]) {
            ypos <- c(ypos, pos)
            pos <- pos - 1
            id <- c(id, j)
        }
        pos <- pos - 0.5
        if (k > 1)
            jj <- c(jj, sum(J0[1:(k - 1)]) + (1:J[k]))
    }
    bottom <- min(ypos) - 1
    med <- rep(NA, sum(J))
    i80 <- array(NA, c(sum(J), 2))
    i80.chains <- array(NA, c(sum(J), n.chains, 2))
    for (j in 1:sum(J)) {
        med[j] <- median(sims.array[, , jj[j]])
        i80[j, ] <- quantile(sims.array[, , jj[j]], c(0.1, 0.9))
        for (m in 1:n.chains) i80.chains[j, m, ] <- quantile(sims.array[,
            m, jj[j]], c(0.1, 0.9))
    }
    rng <- range(i80, i80.chains)
    p.rng <- pretty(rng, n = 2)
    b <- 2/(max(p.rng) - min(p.rng))
    a <- -b * p.rng[1]
    par(mar = c(0, 0, 1, 3))
    plot(c(0, 1), c(min(bottom, -max.length) - 3, 2.5), ann = FALSE,
        bty = "n", xaxt = "n", yaxt = "n", type = "n")
    W <- max(strwidth(unlist(dimnames(summ)[[1]]), cex = cex.names))
    B <- (1 - W)/3.6
    A <- 1 - 3.5 * B
    B <- (1 - A)/3.5
    b <- B * b
    a <- A + B * a
    text(A + B * 1, 2.5, "80% interval for each chain", cex = cex.top)
    text(A + B * 3, 2.6, "R-hat", cex = cex.top)
    lines(A + B * c(0, 2), c(0, 0))
    lines(A + B * c(2.5, 3.5), c(0, 0))
    lines(A + B * c(0, 2), rep(bottom, 2))
    lines(A + B * c(2.5, 3.5), rep(bottom, 2))
    if (min(p.rng) < 0 & max(p.rng) > 0)
        lines(rep(a, 2), c(0, bottom), lwd = 0.5, col = "gray")
    for (x in p.rng) {
        text(a + b * x, 1, x, cex = cex.names)
        lines(rep(a + b * x, 2), c(0, -0.2))
        text(a + b * x, bottom - 1, x, cex = cex.names)
        lines(rep(a + b * x, 2), bottom + c(0, 0.2))
    }
    for (x in seq(1, 2, 0.5)) {
        text(A + B * (1.5 + seq(1, 2, 0.5)), rep(1, 3), c("1",
            "1.5", "2+"), cex = cex.names)
        lines(A + B * rep(1.5 + x, 2), c(0, -0.2))
        text(A + B * (1.5 + seq(1, 2, 0.5)), rep(bottom - 1,
            3), c("1", "1.5", "2+"), cex = cex.names)
        lines(A + B * rep(1.5 + x, 2), bottom + c(0, 0.2))
    }
    for (j in 1:sum(J)) {
        name <- dimnames(summ)[[1]][jj[j]]
        if (id[j] == 1)
            text(0, ypos[j], name, adj = 0, cex = cex.names)
        else {
            pos <- as.vector(regexpr("[[]", name))
            text(strwidth(substring(name, 1, pos - 1), cex = cex.names),
                ypos[j], substring(name, pos, nchar(name)), adj = 0,
                cex = cex.names)
        }
        for (m in 1:n.chains) {
            interval <- a + b * i80.chains[j, m, ]
            if (interval[2] - interval[1] < min.width)
                interval <- mean(interval) + c(-1, 1) * min.width/2
            lines(interval, rep(ypos[j] - 0.1 * (m - (n.chains +
                1)/2), 2), lwd = 1, col = m + 1)
            points(A + B * (1.5 + min(max(summ[jj[j], "Rhat"],
                1), 2)), ypos[j], pch = 20, cex = cex.points)
        }
    }
    for (k in 1:n.roots) {
        if (J[k] < J0[k])
            text(-0.015, ystart[k], "*", cex = cex.names, col = "red")
    }
    if (sum(J != J0) > 0)
        text(0, bottom - 3, "*  array truncated for lack of space",
            adj = 0, cex = cex.names, col = "red")
}
bugs.return.settings <- function (bugs.directory)
{
    registry.default <- readBin(paste(bugs.directory, "System/Rsrc/Registry_default.odc",
        sep = ""), "character", 1e+05, size = 1)
    writeBin(registry.default, paste(bugs.directory, "System/Rsrc/Registry.odc",
        sep = ""))
}
bugs.run <- function (n.burnin, bugs.directory, dos.location)
{
    bugs.update.settings(n.burnin, bugs.directory)
    inopt=options()
    on.exit(options(inopt))
#   options.save <<- options("warn")
#    options.save options("warn")
    options(warn = 2, error = bugs.run.error)
    system(paste(dos.location, "/par", "script.txt"))
    options(error = NULL)
#    options(options.save)
    bugs.return.settings(bugs.directory)
    if (length(grep("Bugs did not run correctly", scan("coda1.txt",
        character(), quiet = TRUE, sep = "\n"))) > 0)
        stop("Look at the log file and\ntry again with debug=T and figure out what went wrong within Bugs.")
}
bugs.run.error <- function ()
{
    cat("Error in the bugs.run() function.\n", "Check that Winbugs14 is in the directory c:\\Program Files\\WinBUGS14\\ and\n",
        "if you have a WinBUGS14Beta directory on your computer, delete it.\n",
        "If this error message still comes up, try editing the bugs.R file:\n",
        "In the definition of the bugs() function, change the assigment,\n",
        "  dos.location = \"c:/progra~1/winbug~1/winbug~1\"\n",
        "to\n", "  dos.location = \"c:/progra~1/winbug~2/winbug~1\"\n",
        "and save the change.\n", "Now exit R (saving your workspace), restart R, and try again.\n",
        "If that does not work, try\n", "  dos.location = \"c:/progra~1/winbug~3/winbug~1\"\n",
        "(There is some problem with accessing files in MS-DOS and Windows.)\n")
    options(error = NULL)
#    options(options.save)
    inopt=options()
    on.exit(options(inopt))
}
bugs.script <- function (parameters.to.save, n.chains, n.iter, n.burnin, n.thin,
    bugs.directory, model.file, debug = FALSE)
{
    nch <- nchar(model.file)
    suffix <- substr(model.file, nch - 3, nch)
    working.directory <- gsub("\\\\", "/", paste(getwd(), "/",
        sep = ""))
    if (suffix == ".bug") {
        model.file.bug <- model.file
        model.file <- paste(substr(model.file, 1, nch - 4), ".txt",
            sep = "")
        system(paste(Sys.getenv("COMSPEC"), "/c copy", model.file.bug,
            model.file))
    }
    else if (suffix != ".txt") {
        stop("model file must be a .bug or .txt file.")
    }
    if (n.chains < 2)
        stop("n.chains must be at least 2")
    n.keep <- ceiling(n.iter/n.thin) - ceiling(n.burnin/n.thin)
    if (n.keep < 2)
        stop("(n.iter-n.burnin)/n.thin must be at least 2")
    script <- paste(bugs.directory, "script.txt", sep = "")
    model <- paste(working.directory, model.file, sep = "")
    data <- paste(working.directory, "data.txt", sep = "")
    history <- paste(working.directory, "history.odc", sep = "")
    coda <- paste(working.directory, "coda", sep = "")
    logfile <- paste(working.directory, "log.odc", sep = "")
    inits <- rep(NA, n.chains)
    for (i in 1:n.chains) {
        inits[i] <- paste(working.directory, "inits", i, ".txt",
            sep = "")
    }
    initlist <- paste("inits (", 1:n.chains, ", '", inits, "')\n",
        sep = "")
    savelist <- paste("set (", parameters.to.save, ")\n", sep = "")
    cat("display ('log')\n", "check ('", model, "')\n", "data ('",
        data, "')\n", "compile (", n.chains, ")\n", initlist,
        "gen.inits()\n", "beg (", ceiling(n.burnin/n.thin) +
            1, ")\n", "thin.updater (", n.thin, ")\n", savelist,
        "update (", ceiling(n.iter/n.thin), ")\n", "stats (*)\n",
        "history (*, '", history, "')\n", "coda (*, '", coda,
        "')\n", "save ('", logfile, "')\n", file = script, sep = "",
        append = FALSE)
    if (!debug)
        cat("quit ()\n", file = script, append = TRUE)
    sims.files <- paste("coda", 1:n.chains, ".txt", sep = "")
    for (i in 1:n.chains) cat("Bugs did not run correctly.\n",
        file = sims.files[i], append = FALSE)
}
bugs.sims <- function (parameters.to.save, n.chains, n.iter, n.burnin, n.thin,
    attach.sims = TRUE, DIC = TRUE)
{
    sims.files <- paste("coda", 1:n.chains, ".txt", sep = "")
    index <- read.table("codaIndex.txt", header = FALSE, sep = "\t")
    parameter.names <- as.vector(index[, 1])
    n.keep <- index[1, 3] - index[1, 2] + 1
    n.parameters <- length(parameter.names)
    n.sims <- n.keep * n.chains
    sims <- array(NA, c(n.sims, n.parameters))
    sims.array <- array(NA, c(n.keep, n.chains, n.parameters))
    root.long <- rep(NA, n.parameters)
    indexes.long <- as.list(rep(NA, n.parameters))
    for (i in 1:n.parameters) {
        temp <- decode.parameter.name(parameter.names[i])
        root.long[i] <- temp$root
        indexes.long[[i]] <- temp$indexes
    }
    n.roots <- length(parameters.to.save)
    left.bracket.short <- as.vector(regexpr("[[]", parameters.to.save))
    right.bracket.short <- as.vector(regexpr("[]]", parameters.to.save))
    root.short <- ifelse(left.bracket.short == -1, parameters.to.save,
        substring(parameters.to.save, 1, left.bracket.short -
            1))
    dimension.short <- rep(0, n.roots)
    indexes.short <- as.list(rep(NA, n.roots))
    n.indexes.short <- as.list(rep(NA, n.roots))
    long.short <- as.list(rep(NA, n.roots))
    length.short <- rep(NA, n.roots)
    for (j in 1:n.roots) {
        long.short[[j]] <- (1:n.parameters)[root.long == root.short[j]]
        length.short[j] <- length(long.short[[j]])
        if (length.short[j] == 0)
            stop(paste("parameter", root.short[[j]], "is not in the model"))
        else if (length.short[j] > 1) {
            dimension.short[j] <- length(indexes.long[[long.short[[j]][1]]])
            n.indexes.short[[j]] <- rep(NA, dimension.short[j])
            for (k in 1:dimension.short[j]) n.indexes.short[[j]][k] <- length(unique(unlist(lapply(indexes.long[long.short[[j]]],
                .subset, k))))
            length.short[j] <- prod(n.indexes.short[[j]])
            if (length(long.short[[j]]) != length.short[j])
                stop(paste("error in parameter", root.short[[j]],
                  "in parameters.to.save"))
            indexes.short[[j]] <- as.list(rep(NA, length.short[j]))
            for (k in 1:length.short[j]) indexes.short[[j]][[k]] <- indexes.long[[long.short[[j]][k]]]
        }
    }
    rank.long <- unlist(long.short)
    for (i in 1:n.chains) {
        sims.i <- scan(sims.files[i], quiet = TRUE)[2 * (1:(n.keep *
            n.parameters))]
        sims[(n.keep * (i - 1) + 1):(n.keep * i), ] <- sims.i
        sims.array[, i, ] <- sims.i
    }
    dimnames(sims) <- list(NULL, parameter.names)
    dimnames(sims.array) <- list(NULL, NULL, parameter.names)
    summary <- monitor(sims.array, keep.all = TRUE)
    last.values <- as.list(rep(NA, n.chains))
    for (i in 1:n.chains) {
        n.roots.0 <- ifelse(DIC, n.roots - 1, n.roots)
        last.values[[i]] <- as.list(rep(NA, n.roots.0))
        names(last.values[[i]]) <- root.short[1:n.roots.0]
        for (j in 1:n.roots.0) {
            if (dimension.short[j] <= 1) {
                last.values[[i]][[j]] <- sims.array[n.keep, i,
                  long.short[[j]]]
                names(last.values[[i]][[j]]) <- NULL
            }
            else last.values[[i]][[j]] <- aperm(array(sims.array[n.keep,
                i, long.short[[j]]], rev(n.indexes.short[[j]])),
                dimension.short[j]:1)
        }
    }
    sims <- sims[sample(n.sims), ]
    sims.list <- as.list(rep(NA, n.roots))
    summary.mean <- as.list(rep(NA, n.roots))
    summary.sd <- as.list(rep(NA, n.roots))
    summary.median <- as.list(rep(NA, n.roots))
    names(sims.list) <- root.short
    names(summary.mean) <- root.short
    names(summary.sd) <- root.short
    names(summary.median) <- root.short
    for (j in 1:n.roots) {
        if (length.short[j] == 1) {
            sims.list[[j]] <- sims[, long.short[[j]]]
            summary.mean[[j]] <- summary[long.short[[j]], "mean"]
            summary.sd[[j]] <- summary[long.short[[j]], "sd"]
            summary.median[[j]] <- summary[long.short[[j]], "50%"]
        }
        else {
            sims.list[[j]] <- aperm(array(sims[, long.short[[j]]],
                c(n.sims, rev(n.indexes.short[[j]]))), c(1, (dimension.short[j] +
                1):2))
            summary.mean[[j]] <- aperm(array(summary[long.short[[j]],
                "mean"], rev(n.indexes.short[[j]])), dimension.short[j]:1)
            summary.sd[[j]] <- aperm(array(summary[long.short[[j]],
                "sd"], rev(n.indexes.short[[j]])), dimension.short[j]:1)
            summary.median[[j]] <- aperm(array(summary[long.short[[j]],
                "50%"], rev(n.indexes.short[[j]])), dimension.short[j]:1)
        }
    }
    summary <- summary[rank.long, ]
    all <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
        n.thin = n.thin, n.keep = n.keep, n.sims = n.sims, sims.array = sims.array[,
            , rank.long], sims.list = sims.list, sims.matrix = sims[,
            rank.long], summary = summary, mean = summary.mean,
        sd = summary.sd, median = summary.median, root.short = root.short,
        long.short = long.short, dimension.short = dimension.short,
        indexes.short = indexes.short, last.values = last.values)
    if (DIC) {
        deviance <- all$sims.array[, , dim(sims.array)[3]]
        pD <- rep(NA, n.chains)
        DIC <- rep(NA, n.chains)
        for (i in 1:n.chains) {
            pD[i] <- var(deviance[, i])/2
            DIC[i] <- mean(deviance[, i]) + pD[i]
        }
        all <- c(all, list(pD = mean(pD), DIC = mean(DIC)))
    }
    if (attach.sims) {
        attach.all(all, name = "bugs.all")
        attach.all(sims.list, name = "bugs.sims")
    }
    return(all)
}
bugs.update.settings <- function (n.burnin, bugs.directory)
{
    char.burnin <- as.character(n.burnin)
    registry <- readBin(paste(bugs.directory, "System/Rsrc/Registry_default.odc",
        sep = ""), "character", 1e+05, size = 1)
    info <- registry[regexpr("Int", registry) > 0]
    while (regexpr("\r", info) > 0) {
        newline <- regexpr("\r", info)
        info <- substring(info, newline + 1)
        line <- substring(info, 1, regexpr("\r", info) - 1)
        if (regexpr("AdaptivePhase", line) > 0) {
            numpos <- regexpr("Int", line) + 4
            num <- substring(line, numpos)
            if (as.numeric(num) > n.burnin) {
                num.new <- paste(paste(rep(" ", nchar(num) -
                  nchar(char.burnin)), sep = "", collapse = ""),
                  char.burnin, sep = "")
                line.new <- sub(num, num.new, line)
                registry <- sub(line, line.new, registry)
            }
        }
    }
    writeBin(registry, paste(bugs.directory, "System/Rsrc/Registry.odc",
        sep = ""))
}
chisqdf <- function (A, varA)
{
    2 * (A^2/varA)
}
conv.par <- function (x)
{
    alpha <- 0.05
    m <- ncol(x)
    n <- nrow(x)
    xdot <- apply(x, 2, mean)
    s2 <- apply(x, 2, var)
    W <- mean(s2)
    B <- n * var(xdot)
    muhat <- mean(xdot)
    varW <- var(s2)/m
    varB <- B^2 * 2/(m - 1)
    covWB <- (n/m) * (cov(s2, xdot^2) - 2 * muhat * cov(s2, xdot))
    sig2hat <- ((n - 1) * W + B)/n
    quantiles <- quantile(as.vector(x), probs = c(0.025, 0.25,
        0.5, 0.75, 0.975))
    if (W > 1e-08) {
        postvar <- sig2hat + B/(m * n)
        varpostvar <- max(0, (((n - 1)^2) * varW + (1 + 1/m)^2 *
            varB + 2 * (n - 1) * (1 + 1/m) * covWB)/n^2)
        post.df <- min(chisqdf(postvar, varpostvar), 1000)
        post.range <- muhat + sqrt(postvar) * qt(1 - alpha/2,
            post.df) * c(-1, 0, 1)
        varlo.df <- chisqdf(W, varW)
        confshrink.range <- sqrt(c(postvar/W, (n - 1)/n + (1 +
            1/m) * (1/n) * (B/W) * qf(0.975, m - 1, varlo.df)) *
            (post.df + 3)/(post.df + 1))
        n.eff <- m * n * min(sig2hat/B, 1)
        list(post = post.range, quantiles = quantiles, confshrink = confshrink.range,
            n.eff = n.eff)
    }
    else {
        list(post = muhat * c(1, 1, 1), quantiles = quantiles,
            confshrink = c(1, 1), n.eff = 1)
    }
}
conv.par.log <- function (r)
{
    conv.p <- conv.par(log(r))
    list(post = exp(conv.p$post), quantiles = exp(conv.p$quantiles),
        confshrink = conv.p$confshrink, n.eff = conv.p$n.eff)
}
conv.par.logit <- function (r)
{
    conv.p <- conv.par(logit(r))
    list(post = invlogit(conv.p$post), quantiles = invlogit(conv.p$quantiles),
        confshrink = conv.p$confshrink, n.eff = conv.p$n.eff)
}
cov <- function (a, b)
{
    m <- length(a)
    ((mean((a - mean(a)) * (b - mean(b)))) * m)/(m - 1)
}
decode.parameter.name <- function (a)
{
    left.bracket <- regexpr("[[]", a)
    if (left.bracket == -1) {
        root <- a
        dimension <- 0
        indexes <- NA
    }
    else {
        root <- substring(a, 1, left.bracket - 1)
        indexes <- NULL
        right.bracket <- regexpr("[]]", a)
        a <- substring(a, left.bracket + 1, right.bracket - 1)
        while (nchar(a) > 0) {
            comma <- regexpr(",", a)
            if (comma == -1) {
                indexes <- c(indexes, as.numeric(a))
                a <- ""
            }
            else {
                indexes <- c(indexes, as.numeric(substring(a,
                  1, comma - 1)))
                a <- substring(a, comma + 1, nchar(a))
            }
        }
        dimension <- length(indexes)
    }
    return(list(root = root, dimension = dimension, indexes = indexes))
}

format.data <- function (x,...) {
    datalist=x
    if (!is.list(datalist) || is.data.frame(datalist))
        stop("Argument to format.data must be a list.")
    n <- length(datalist)
    datalist.string <- as.list(rep(NA, n))
    for (i in 1:n) {
        if (length(datalist[[i]]) == 1)
            datalist.string[[i]] <- paste(names(datalist)[i],
                "=", as.character(datalist[[i]]), sep = "")
        if (is.vector(datalist[[i]]) & length(datalist[[i]]) >
            1)
            datalist.string[[i]] <- paste(names(datalist)[i],
                "=c(", paste(as.character(datalist[[i]]), collapse = ", "),
                ")", sep = "")
        if (is.array(datalist[[i]]))
            datalist.string[[i]] <- paste(names(datalist)[i],
                "= structure(.Data= c(", paste(as.character(as.vector(aperm(datalist[[i]]))),
                  collapse = ", "), "), .Dim=c(", paste(as.character(dim(datalist[[i]])),
                  collapse = ", "), "))", sep = "")
    }
    datalist.tofile <- paste("list(", paste(unlist(datalist.string),
        collapse = ", "), ")", sep = "")
    return(datalist.tofile)
}


invlogit <- function (x)
{
    1/(1 + exp(-x))
}
logit <- function (x)
{
    log(x/(1 - x))
}

monitor <- function (a, trans = NULL, keep.all = FALSE, Rupper.keep = FALSE)
{
    output <- NULL
    nparams <- ifelse(length(dim(a)) < 3, 1, dim(a)[length(dim(a))])
    if (length(dim(a)) == 2)
        a <- array(a, c(dim(a), 1))
    if (!keep.all) {
        n <- floor(dim(a)[1]/2)
        a <- a[(n + 1):(2 * n), , ]
    }
    if (is.null(trans))
        trans <- ifelse((apply(a <= 0, 3, sum)) == 0, "log",
            "")
    for (i in 1:nparams) {
        if (trans[i] == "log")
            conv.p <- conv.par.log(a[, , i])
        else if (trans[i] == "logit")
            conv.p <- conv.par.logit(a[, , i])
        else conv.p <- conv.par(a[, , i])
        output <- rbind(output, c(mean(a[, , i]), sqrt(var(as.vector(a[,
            , i]))), conv.p$quantiles, conv.p$confshrink, round.sci(conv.p$n.eff,
            2)))
    }
    dimnames(output) <- list(dimnames(a)[[3]], c("mean", "sd",
        "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "Rupper",
        "n.eff"))
    if (!Rupper.keep)
        return(output[, -(ncol(output) - 1)])
    else return(output)
}



round.bugs <- function (x, digits)
{
    a=x
    r <- round(a, digits)
    exponent <- floor(log10(abs(a)))
    number <- abs(a) * 10^-exponent
    int <- floor(number)
    int.frac <- formatC(number, format = "f", digits = digits)
    prefix <- ifelse(a > 0, " ", "-")
    zeroes <- rep("0", length(a))
    scientific <- paste(prefix, int.frac, "E", exponent, sep = "")
    return(ifelse(a == 0, zeroes, scientific))
}
round.bugs.list <- function (x, digits)
{
    inits=x
    lapply(inits, round.bugs, digits)
}
round.sci <- function (x, digits)
{
    a=x
    round(a, min(0, digits - 1 - floor(log10(a))))
}



write.datafile <- function (datalist, towhere, fill = TRUE) {
    if (!is.list(datalist) || is.data.frame(datalist))
        stop("First argument to write.datafile must be a list.")
    cat(format.data(datalist), file = towhere, fill = fill)
}
