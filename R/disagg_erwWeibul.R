library(cmaes)

erwWeibul <- function(t, par=c(10, 1, 1, 50, -0.01)) {
  par[1] + par[2] * par[3] * par[4] * ((t - par[5])*par[3])^(par[4]-1) * exp(-par[3]*(t-par[5])^par[4])
}

startWeibul <- list(latePeak = function() c(1, 1, 30, -0.01),
                     default = function() c(1,1,1.2, runif(1, -1, 0)))

optFun <- function(par, qStart, TM, qEnd, HQ) {
  eVal <- 1e6# +runif(1,-1e2,1e2)

  if(abs(par[2]) < 1e-3)
    return(eVal)
  
  p1 <- qStart - par[1] * par[2] * par[3] * ((0 - par[4])*par[2])^(par[3]-1) * exp(-par[2]*(0-par[4])^par[3])

  if(is.nan(p1) | is.na(p1))
    return(eVal)
  
  if(par[4] > 0)
    return(eVal)
  
  dayDis <- erwWeibul(0:23/24, c(p1,par))

  if(any(is.nan(dayDis)) | any(is.infinite(dayDis)) | any(is.na(dayDis)))
    return(eVal)
  if (any(dayDis < 0))
    return(eVal)
  # if (dayDis[24] <= qEnd)
  #   return(eVal)
  # if (dayDis[24] >= HQ)
  #   return(eVal)

  # if (((dayDis[2]-p1+1e-4)/(max(dayDis)-p1)) < 0.01)
  #   return(eVal)

  diffQEnd <- ifelse(dayDis[24] <= qEnd, qEnd - dayDis[24], 0)
  
  absDiffHQ <- abs(max(dayDis)-HQ)
  absDiffTM <- abs(sum(dayDis)/24 - TM)^2

  res <- absDiffHQ+absDiffTM+diffQEnd*10
  
  if(is.infinite(res))
    return(eVal)
  
  return(res)
}

lorentz <- function(t, par=c(10,3.079638, 0.8639837, 0.9580365)) {
  par[1] + 2 * par[2]/pi*par[3]/(4*(t-par[4])^2+par[3]^2)
}

startLorentz <- list(latePeak = function() c(1,.02,0.95),
                     default = function() c(0,1,runif(1)))

optFunLor <- function(par, qStart, TM, qEnd, HQ) {
  eVal <- 1e6# +runif(1,-1e2,1e2)
  
  if(abs(par[2]) < 1e-3 & abs(par[3]) < 1e-3)
    return(eVal)

  p1 <- qStart - 2 * par[1]/pi*par[2]/(4*(0-par[3])^2+par[2]^2)
  
  if(is.nan(p1))
    return(eVal)
  
  # ensure minimum width of Lorentz function
  # if(abs(par[2]) < 0.2)
  #   return(eVal)
  
  dayDis <- lorentz(0:23/24, c(p1,par))
  
  if(any(is.nan(dayDis)) | any(is.infinite(dayDis)))
    return(eVal)
  if (any(dayDis < 0))
    return(eVal)
  # if (dayDis[24] < qEnd)
  #   return(eVal)
  # if (dayDis[24] >= HQ)
  #   return(eVal)

  diffQEnd <- ifelse(dayDis[24] <= qEnd, qEnd - dayDis[24], 0)
  
  absDiffHQ <- abs(max(dayDis)-HQ) 
  absDiffTM <- abs(sum(lorentz(0:23/24, c(p1,par)))/24 - TM)^2

  # cat(absDiffHQ, absDiffTM, diffQEnd, "\n")
  return(absDiffHQ+absDiffTM+diffQEnd*10)
}

#############
lorentzErwWeibul <- function(qStart, TM, qEnd, HQ, qPrec, decision) {
  stopifnot(decision %in% c("wagner", "bestfit"))
  
  if(any(is.na(c(qStart, TM, qEnd, HQ, qPrec))) | TM > HQ ) {
    diagRow <- matrix(c(qStart, TM, qEnd, HQ, rep(NA, 9)), nrow = 1)
    colnames(diagRow) <- c("qStart", "TM", "qEnd", "HQ",
                           "wei1", "wei2", "wei3", "wei4", "wei5", 
                           "lor1", "lor2", "lor3", "lor4")
    res <- rep(NA,24)
    attr(res, "diagnose") <- diagRow
    return(res)
  }
  
  # if (HQ - TM > TM & qEnd > TM)
  #   peakType <- "latePeak"
  # else 
    peakType <- "default"
  
  # Lorentzfunktion
  parFitLor <- optim(par=startLorentz[[peakType]](), fn=optFunLor, 
                     qStart=qStart, TM=TM, qEnd=qEnd, HQ=HQ,
                    control = list(maxit=1e3), method="L-BFGS-B")
  count <- 0
  # cat(parFitLor$value, "\n")
  
  while (parFitLor$value > (HQ-TM)/5 & count < 50) {
    parFitLor <- optim(par=startLorentz[[peakType]](),
                        fn=optFunLor, 
                       qStart=qStart, TM=TM, qEnd=qEnd, HQ=HQ,
                       control = list(maxit=1e3), 
                       method="L-BFGS-B",
                       lower = c(-Inf, 0, 0),
                       upper = c(Inf, 1, 1))
    count <- count + 1
    # cat(parFitLor$value, "\n")
  }
  
  if (count > 0 & count < 50 )
    warning("Retried Lorentz optimisation ", count, " times.")
  
  if(count == 50) {
    decision <- "bestfit"
    warning("Lorentz optimisation aborted for: ", 
            paste(c("qStart", " TM", " qEnd", " HQ"), round(c(qStart, TM, qEnd, HQ), 2), sep=": "),
            " with error:", parFitLor$value)
  }
  
  par <- parFitLor$par
  p1 <- qStart - 2 * par[1]/pi*par[2]/(4*(0-par[3])^2+par[2]^2)
  dayDisLor <- lorentz(0:23/24, c(p1,par))
  
  p1Lor <- p1
    
  ## erweiterte Weibull
  parFit <- optim(par=startWeibul[[peakType]](), fn=optFun,
                  qStart=qStart, TM=TM, qEnd=qEnd, HQ=HQ,
                  method = "L-BFGS-B", control = list(maxit=1e4),
                  lower = c(0, 0, -Inf, -1),
                  upper = c(Inf, Inf, Inf, 0))
  
  count <- 0
  while (parFit$value > (HQ-TM)/5 & count < 50) {
    parFit <- optim(par=startWeibul[[peakType]](), fn=optFun, 
                    qStart=qStart, TM=TM, qEnd=qEnd, HQ=HQ,
                    method = "L-BFGS-B", control = list(maxit=1e4),
                    lower = c(0, -10, -Inf, -1),
                    upper = c(Inf, Inf, Inf, 0))
    count <- count + 1
  }
  
  if (count > 0 & count < 50)
    warning("Retried Weibull optimisation ", count, " times.")
  
  if(count == 50) {
    decision <- "bestfit"
    warning("Weibull optimisation aborted for: ",
            paste(c("qStart", " TM", " qEnd", " HQ"), round(c(qStart, TM, qEnd, HQ), 2), sep=": "),
            " with error:", parFit$value)
  }
  
  par <- parFit$par
  p1 <- qStart - par[1] * par[2] * par[3] * ((0 - par[4])*par[2])^(par[3]-1) * exp(-par[2]*(0-par[4])^par[3])
  dayDisWei <- erwWeibul(0:23/24, c(p1,par))
  
  diagRow <- matrix(c(qStart, TM, qEnd, HQ, p1, parFit$par, p1Lor, parFitLor$par), nrow = 1)
  colnames(diagRow) <- c("qStart", "TM", "qEnd", "HQ",
                         "wei1", "wei2", "wei3", "wei4", "wei5", 
                         "lor1", "lor2", "lor3", "lor4")
  
  attr(dayDisWei, "diagnose") <- diagRow
  attr(dayDisLor, "diagnose") <- diagRow

  if (decision == "wagner") {
    A <- TM > qPrec | TM > qEnd # wachsend - fallend
    B <- (TM > dayDisLor[24] & dayDisLor[24] > qEnd) | (TM < dayDisLor[24] & dayDisLor[24] < qEnd) # zwischen aktuell und Folgetag
    
    if (A) {
      if (B) {
        cat("[Selected Lorentz.]\n")
        return(dayDisLor)
      }
      else {
        cat("[Selected Weibull.]\n")
        return(dayDisWei)
      }
    }
    
    C <- qPrec <= TM & TM <= qEnd # monoton steigend
    
    if (C) {
      cat("[Selected Weibull.]\n")
      return(dayDisWei)
    }
  
    D <- qPrec >= TM & TM >= qEnd # monoton fallend
  
    if (D) {
      if (B) {
        cat("[Selected Lorentz.]\n")
        return(dayDisLor)
      } else {
        cat("[Selected Weibull.]\n")
        return(dayDisWei)
      }
    }
    
    E <- qPrec > TM & TM < qEnd # fallend - wachsend
    
    if (E) {
      cat("[Selected Lorentz.]\n")
      return(dayDisLor)
    }
  } 
  if (decision == "bestfit") {
    if (parFitLor$value < parFit$value) {
      cat("[Selected Lorentz.]\n")
      return(dayDisLor)
    } else {
      cat("[Selected Weibull.]\n")
      return(dayDisWei)
    }
  }
}


###

a3<-function(qs,q1,q2,q3){
  a3t<- -1/9*(6*qs-2*q3+7*q2-11*q1)
  return(a3t)
}

a2<-function(t,qs,q1,q2,q3){
  a2t<-1/6*((12*t+12)*qs-q3+t*(-4*q3+14*q2-22*q1)+8*q2-19*q1)
  return(a2t)
}

a1<-function(t,qs,q1,q2,q3){
  a1t<- -1/36*((72*t^2+144*t+42)*qs+4*q3+t*(-12*q3+96*q2-228*q1)+t^2*(-24*q3+84*q2-132*q1)-23*q2-23*q1)
  return(a1t)
}

a0<-function(t,qs,q1,q2,q3){
  a0t<- 1/72*((48*t^3+144*t^2+84*t-12)*qs+t*(8*q3-46*q2-46*q1)+q3+t^2*(-12*q3+96*q2-228*q1)+t^3*(-16*q3+56*q2-88*q1)-8*q2+91*q1)
  return(a0t)
}

disaggpoly<-function(t,a0,a1,a2,a3){
  return(a3*t^3+a2*t^2+a1*t+a0)
}

disagg<-function(q, scheitelTag=c(), decision="wagner", diagnose=FALSE){
  
  q[,1]<-as.Date(q[,1], "%d.%m.%Y")
  if(length(scheitelTag) > 0)
    scheitelTag[,1]<-as.Date(scheitelTag[,1], "%d.%m.%Y")
  
  n <- length(q[,1])
  outv<-numeric()
  
  ###Startwerte berechnen mit t=1###
  a30<-a3(q[1,2],q[1,2],q[2,2],q[3,2])
  a20<-a2(0.5,q[1,2],q[1,2],q[2,2],q[3,2])
  a10<-a1(0.5,q[1,2],q[1,2],q[2,2],q[3,2])
  a00<-a0(0.5,q[1,2],q[1,2],q[2,2],q[3,2])
  
  outv <- disaggpoly(1:24/24,a00,a10,a20,a30)
  
  #für jeden Ursprungszeitschritt neue Koeffizienten
  
  pb <- txtProgressBar(0, n-3, style = 3)
  
  diagList <- NULL
  
  for (j in 1:(n-3)) {
    setTxtProgressBar(pb, j)

    qS <- outv[24*j]
    if(is.na(qS))
      qS <- q[j+1,2]
    
    if (q[j+1,1] %in% scheitelTag[,1]) {
      sInd <- which(scheitelTag[,1] == q[j+1,1])
      resDisAg <- lorentzErwWeibul(qStart=qS, TM=q[j+1,2],
                                   qEnd=q[j+2,2], HQ=scheitelTag[sInd,2],
                                   qPrec=q[j,2], decision)
      outv <- c(outv, resDisAg)
      if(diagnose)
        diagList <- rbind(diagList, attr(resDisAg, "diagnose"))
    } else {
      a3j<-a3(qS,q[j+1,2],q[j+2,2],q[j+3,2])
      a2j<-a2(0.5,qS,q[j+1,2],q[j+2,2],q[j+3,2])
      a1j<-a1(0.5,qS,q[j+1,2],q[j+2,2],q[j+3,2])
      a0j<-a0(0.5,qS,q[j+1,2],q[j+2,2],q[j+3,2])
      
      ##nutze Startwerte als Startwertbedingung##
      outv <- c(outv, disaggpoly(1:24/24,a0j,a1j,a2j,a3j))
    }
    
    if (any(is.na(tail(outv, 24))) | any(is.infinite(tail(outv, 24)))) {
      # eretze mit Polynom
      a3j<-a3(qS,q[j+1,2],q[j+2,2],q[j+3,2])
      a2j<-a2(0.5,qS,q[j+1,2],q[j+2,2],q[j+3,2])
      a1j<-a1(0.5,qS,q[j+1,2],q[j+2,2],q[j+3,2])
      a0j<-a0(0.5,qS,q[j+1,2],q[j+2,2],q[j+3,2])
      
      
      outv <- c(outv[-c(length(outv)-0:23)], 
                disaggpoly(1:24/24,a0j,a1j,a2j,a3j))
      
      warning("Neither Weibull nor Lorentz could disaggregate the timeseries; using polynomial fallback.")
    }
  }
  
  close(pb)
  
  thour <- seq(from=ISOdate(format(q[1,1], "%Y"),format(q[1,1], "%m"),format(q[1,1], "%d"), hour=0), length.out=length(outv), by="hour")
  
  res <- data.frame(thour=thour,
                    dissAgg = outv)
  
  attr(res, "diagnose") <- diagList
  
  return(res)
}
