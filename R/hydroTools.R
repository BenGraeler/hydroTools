# aggregate time series for different duration levels
durationLevels <- function(series, timeRes = 5/60, durLevels = c(1,3,6,12,24),
                           align = "left", endrule = "NA", ..., NaN.rm=FALSE) {
  stopifnot(is.numeric(series))
  stopifnot(length(timeRes) == 1)
  
  if(any(durLevels/timeRes < ceiling(durLevels/timeRes)))
    warning("The durLevels are not a multiple of the temporal resolution timeRes. Multiplicies are rounded to the next larger integer.")
  
  warnOnce <- FALSE
  res <- matrix(NA, length(series), ncol = length(durLevels))
  for(i in 1:length(durLevels)) {
    dl <- durLevels[i]
    res[,i] <- runmean(series, dl/timeRes, align = align,
                       endrule = endrule, ...) * dl/timeRes
    if (sum(is.nan(res[,i])) > 0) {
      if(!NaN.rm)
        warnOnce <- TRUE
      else
        res[is.nan(res)] <- NA
    }
  }
  
  if(warnOnce)
    warning("NaNs have been produced, but not treated.")
  colnames(res) <- paste("dl", durLevels, sep="")
  
  return(res)
}

# select the hydrological year/season 
hydroYears <- function(start, end, breaks=c("11-01", "10-31")) {
  stopifnot(end - start > 0)
  n.breaks <- length(breaks)
  stopifnot(n.breaks %% 2 == 0)

  jmp.breaks <- which(diff(rank(sapply(strsplit(breaks, "-", fixed=T),
                                       function(x) x[1]))) != 1)
  
  if (end - start == 1)
      yrs <- c(rep(start, jmp.breaks),
               rep(end, n.breaks - jmp.breaks))
  else 
    yrs <- c(rep(start, jmp.breaks),
             rep((start+1):(end-1), each=n.breaks),
             rep(end, n.breaks - jmp.breaks))
  apply(matrix(paste(yrs, breaks, sep="-"),
               ncol = 2, byrow = T),
        1, function(x) paste(x[1], x[2],sep="/"))
}


.aggInt <- function(ts, int, FUN, ...) {
  t(sapply(int, function(x) apply(ts[x], 2, FUN, ...)))
}

## extract events from time series

# input: ts, lower threshold, time gap, tc/area
# output: data.frame with: time, endTime, Volume, peak intensity, duration

ts2events <- function(ts, var, AE = 0, tc = 0.5*AE^0.6, 
                      gap = max(4*60, 2*tc*60),
                      thd = 0.01, tunit = "mins") {
  stopifnot(inherits(ts, "xts"))
  
  if(ncol(ts) > 1 & missing(var)) {
    warning("The first column of xts time series is used. Specify 'var' to select a diffrent one.")
    ts <- ts[,1]
  }
  if (!missing(var))
    ts <- ts[,var]
  
  aboveThd <- which(ts > thd)
  eStart <- c(aboveThd[1], aboveThd[which(diff(aboveThd) > 1)+1])
  eEnd <- c(aboveThd[which(diff(aboveThd) > 1)], tail(aboveThd,1))
  
  timeTs <- index(ts)
  boolInd <- which(difftime(timeTs[eStart[-1]], timeTs[eEnd[-length(eEnd)]],
                            units = tunit) > gap)
  eStart <- eStart[c(1, boolInd+1)]
  eEnd <- eEnd[c(boolInd, length(eEnd))]
  
  aboveZero <- which(ts > 0)
  eStartZero <- c(aboveZero[1], aboveZero[which(diff(aboveZero) > 1)+1])
  eEndZero <- c(aboveZero[which(diff(aboveZero) > 1)],tail(aboveZero,1))
  
  for (i in 1:length(eStart)) {
    eStart[i] <- min(eStart[i], tail(eStartZero[eStartZero <= eStart[i] & eStartZero > c(0,eEnd)[i]],1))
    eEnd[i] <- max(eEnd[i], eEndZero[eEndZero >= eEnd[i] & eEndZero < c(eStart,length(ts))[i+1]][1], na.rm = T)
  }
  
  res <- data.frame(time=timeTs[eStart], endTime=timeTs[eEnd])
  stats <- do.call(rbind, 
                   .aggInt(ts, paste(res$time,res$endTime,sep="/"), 
                           function(x) data.frame(volume=sum(x),
                                                  peakTime=attr(x,"names")[which.max(x)], peak=max(x))))
  res$volume <- stats$volume
  res$peakTime <- as.POSIXct(as.character(stats$peakTime),
                             format = "%Y-%m-%d %H:%M:%S",
                             tz = attributes(index(ts))$tzone)
  res$peak <- stats$peak
  res$duration <- difftime(res$endTime, res$time, units = tunit)
  
  return(res)
}

###################
# # calculate mowing window quantiles
# ts <- cbind(ts, runquantile(ts, w.width, probs, endrule = "NA", align = w.align))


## count quantile exceedances for a given window
exQt <- function(ts, probs = 0.5, qt,
                 w.width = 10, w.align = "right", 
                 exceed = TRUE, relative = TRUE) {
  if(missing(qt))
    qt <- quantile(ts, probs)
  
  # check whether either qt or w.width has length 1
  stopifnot(length(qt) == 1 | length(w.width) == 1)
  
  exQtSingle <- function(.qt, .w.width) {
    # count exceedances
    if(exceed)
      cExc <- apply(embed(ts, .w.width), 1, function(x) sum(x > .qt))
    else
      cExc <- apply(embed(ts, .w.width), 1, function(x) sum(x < .qt))
    
    if(relative)
      cExc <- cExc/.w.width
    
    switch(w.align,
           left = c(cExc, rep(NA, .w.width-1)),
           center = c(rep(NA, floor((.w.width-1)/2)), cExc, rep(NA, ceiling((.w.width-1)/2))),
           right = c(rep(NA, .w.width-1), cExc),
           stop("Alignment method \"", w.align, "\" is not implemented!"))
  }
  
  if(length(w.width) == 1) {
    res <- do.call(cbind, lapply(qt, function(x) exQtSingle(x, w.width)))
    colnames(res) <- names(qt)
  } else {
    res <- do.call(cbind, lapply(w.width, function(x) exQtSingle(qt, x)))
    colnames(res) <- as.character(w.width)
  }
  
  return(xts(res, index(ts)))
}