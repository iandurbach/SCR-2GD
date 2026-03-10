##############################################################################
## package 'secrdesign'
## optimalSpacing.R
## 2017-07-15 moved from Lambda.R
## 2017-08-09 finished (?) fiddling with plot etc.
## 2018-11-28 distribution; openCR
## 2019-02-15 try-error catch bad uniroot in interpRSE
## 2020-01-20 2.6.0 removed openCR
## 2024-04-28 "bg" added to dotargs in plot method
## 2026-01-01 validate step in simRSE
## 2026-01-01 minsimRSE checks for NA
##############################################################################

interpRSE <- function (values) {
  ## find minimum at n = r
  nminr <- function(RR) {
    n <- approx(values$R, values$n, RR)$y
    r <- approx(values$R, values$r, RR)$y
    n-r
  }
  rangeR <- range(values$R)
  if (diff(rangeR)>0) {
    ur <- try(uniroot(nminr, interval = rangeR))
    if (inherits(ur, 'try-error'))
      list(minimum = NA, objective = NA)
    else
      list(minimum = ur$root, objective = approx(values$R, values$RSE, ur$root)$y)
  }
  else         list(minimum = NA, objective = NA)
  
}

oneRSE <- function (R, k, traps, xsigma, detectpar, noccasions, nrepeats, detectfn,
                    maskargs = list(), CF) {
  # hold traps constant, vary sigma, D, mask
  tmpsigma <- spacing(traps) / R
  tmpD <- (k / tmpsigma)^2
  detectpar$sigma <- tmpsigma
  maskargs$buffer <- xsigma * tmpsigma
  maskargs$spacing <- maskargs$spacing / R
  mask <- do.call(make.mask, maskargs)
  if (nrow(mask) < 600) warning ("mask has only ", nrow(mask), " points")
  minnrRSE(D = tmpD * nrepeats, traps, mask, detectpar, noccasions, detectfn, CF)
}
##############################################################################

getrotRSE <- function (R, k, traps, xsigma, detectpar, noccasions, nrepeats,
                       detectfn, maskargs = list(), CF, 
                       distribution = "poisson") {
  
  tmpsigma <- spacing(traps) / R
  tmpD <- (k / tmpsigma)^2
  maskargs$buffer <- tmpsigma * xsigma
  maskargs$spacing <- maskargs$spacing / R
  mask <- do.call(make.mask, maskargs)
  detectpar$sigma <- tmpsigma
  
  nrm <- Enrm(tmpD, traps = traps, mask = mask, detectpar = detectpar,
              noccasions = noccasions, detectfn = detectfn) * nrepeats
  
  rotRSE <- minnrRSE(tmpD * nrepeats, traps = traps, mask = mask, detectpar = detectpar,
                     noccasions = noccasions, detectfn = detectfn, CF = 1.0)
  
  if (distribution == "binomial") {
    # 2019-01-08 blocked as doesn't seem to do anything
    # Pxy <- pdot(mask, traps, detectfn, detectpar, noccasions)
    # esa <- sum(Pxy) * attr(mask, 'area')
    # 2019-01-08 inserted * nrepeats 
    rotRSE <- sqrt(rotRSE^2 - 1 / (tmpD * nrepeats * maskarea(mask)))
  }
  c(R, nrm, rotRSE*CF)
}
##############################################################################

simRSEfn <- function (R, k, traps, xsigma, detectpar, noccasions, nrepeats, detectfn,
                      nx, nrepl, fitfunction, allargs) {
  tmpsigma <- spacing(traps) / R  # vector same length as R
  tmpD <- (k / tmpsigma)^2        # vector same length as R
  
  ## scenarios
  scen1 <- make.scenarios(trapsindex = 1, noccasions = noccasions, nrepeats = nrepeats, D = tmpD,
                          sigma = 1, lambda0 = detectpar$lambda0, detectfn = detectfn)
  scen1$sigma <- tmpsigma  ## specify here to avoid crossing with D
  
  ## general arguments
  defaultargs <- list(nrepl = nrepl, trapset = traps, scenarios = scen1, fit = TRUE,
                      fit.function = fitfunction, byscenario = FALSE, ncores = NULL, 
                      xsigma = xsigma, nx = nx)
  dotsargs <- allargs[names(allargs) %in% c("seed", "ncores")]
  runargs <- replacedefaults(defaultargs, dotsargs)
  
  ## sim.popn arguments: only Ndist
  defaultpopargs <- list(Ndist = "poisson")
  popargs <- allargs[names(allargs) %in% c("Ndist")]
  runargs$pop.args <- replacedefaults(defaultpopargs, popargs)
  
  ## arguments for model fitting function
  defaultfitargs <- list(start = "true", 
                         detectfn = detectfn, 
                         details = list())
  fitargs <- allargs[!(names(allargs) %in% c("seed", "ncores", "Ndist"))]
  fitargs <- fitargs[names(fitargs) %in% names(formals(fitfunction))]   ## 2019-02-15
  runargs$fit.args <- replacedefaults(defaultfitargs, fitargs)    
  
  ## run
  sims1 <- do.call(run.scenarios, runargs)
  allstats <- select.stats(sims1, parameter = "D", c("estimate", "RSE", "RB","ERR"))
  
  # optional validate step 2026-01-01
  defaultvalidargs <- list(x = allstats, test = "RSE", validrange = c(0, Inf), quietly = TRUE)
  validargs <- allargs[names(allargs) %in% c("test", "validrange")]
  if (!is.null(validargs) && length(validargs)>0) {
    validargs <- replacedefaults(defaultvalidargs, validargs)
    allstats <- do.call(validate, validargs)
  }
  
  # save each replicate
  eachRSE <- do.call(rbind, allstats$output)[, c("estimate", "RSE")]
  eachRSE <- cbind(R = rep(R, each = nrepl), eachRSE)
  
  # summarise
  tmp <- summary(allstats, fields = c("n","mean","se","rms"))$OUTPUT
  tmp2 <- lapply(tmp, unlist)
  simout <- as.data.frame(t(sapply(tmp2, '[', c('n2','mean2','se2', 'mean3','se3','rms4'))))
  names(simout) <- c("n","RSE.mean","RSE.se","RB.mean", "RB.se","rRMSE")
  simout$rRMSE <- simout$rRMSE / scen1$D
  
  list(eachRSE = eachRSE, summary = cbind(data.frame(R = R), simout))
}
##############################################################################

optimalSpacing <- function (
    D, 
    traps, 
    detectpar, 
    noccasions,
    nrepeats = 1,
    detectfn = c('HHN', 'HHR', 'HEX', 'HAN', 'HCG', 'HN', 'HR', 'EX'),
    fittedmodel = NULL,
    xsigma = 4,
    R = seq(0.2, 4, 0.2),
    CF = 1.0,
    distribution = c("poisson", "binomial"),
    fit.function = c("none", "secr.fit"),
    simulationR = seq(0.4, 4, 0.4),
    nrepl = 10,
    plt = FALSE,
    ...) {
  
  ## if a fitted model is provided then all preceding arguments are overridden
  ## (D, traps, detectpar, noccasions, detectfn)
  
  ## This function seeks the detector spacing that minimises RSE for a given geometry.
  ## It uses the trick of
  ## (i)  expressing spacing as a multiple of sigma, and
  ## (ii) varying (relative) spacing by scaling sigma and other spatial parameters
  ##      while keeping the detector array fixed.
  ## For relative spacing = R, the scaling factors are
  ##    sigma = sigma / R
  ##    D = (k / sigma)^2
  ##    buffer = xsigma * sigma / R
  ##    maskspacing = maskspacing / R
  
  ## note that if trap spacing is doubled,
  ## - the number of animals per trap square is quadrupled
  ## - sigma / spacing is halved
  
  ## suppress simulations with nrepl = 0
  ## suppress optimisation with CF = NA
  ## control values with R
  fit.function <- match.arg(fit.function)
  distribution <- match.arg(distribution)
  if (!is.null(fittedmodel)) {
    if (ms(fittedmodel$capthist))
      stop ("optimalSpacing requires single-session data")
    pred <- predict(fittedmodel)
    detectfn <- fittedmodel$detectfn
    detectpar = list(g0 = pred['g0','estimate'],
                     lambda0 = pred['lambda0','estimate'],
                     sigma = pred['sigma','estimate'],
                     z = pred['z','estimate'],
                     w = pred['w','estimate'])
    
    traps <- attr(fittedmodel$capthist, 'traps')
    if (any(detector(traps)=='single'))
      warning ("results are unreliable for single-catch traps when lambda0 ",
               "inferred from multi-catch model")
    noccasions <- ncol(fittedmodel$capthist)
    
    if (fittedmodel$CL) {
      D <- derived(fittedmodel)['D','estimate']
    }
    else {
      D <- pred['D','estimate']
    }
  }
  else {
    if (is.character(detectfn)) detectfn <- match.arg(detectfn)
    detectfn <- secr:::secr_valid.detectfn(detectfn, valid = c(0,1,2,14:19))
    
  }
  dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 2 to 16
  detectfn <- dfc$detectfn
  detectpar <- dfc$detectpar
  
  args <- list(...)
  k <- detectpar$sigma * D^0.5
  
  # prepare make.mask arguments
  defaultmaskargs <- list (traps = traps, buffer = xsigma * spacing(traps),
                           nx = 64, type = "trapbuffer")
  dotsargs <- args[names(args) %in% c('nx', 'type', 'poly','poly.habitat')]
  maskargs <- replacedefaults(defaultmaskargs, dotsargs)
  tmpmask <- do.call(make.mask, maskargs) # baseline mask constructed at R = 1
  maskargs$spacing <- spacing(tmpmask)
  
  if (any(detector(traps) == 'single')) {
    warning ("treating single-catch traps as multi-catch", call. = FALSE)
    detector(traps) <- 'multi'
  }
  
  #################
  values <- sapply(R, getrotRSE, k, traps, xsigma, detectpar,
                   noccasions, nrepeats, detectfn, maskargs, CF, distribution)
  values <- as.data.frame(t(values))
  names(values) <- c("R", "n", "r", "m", "RSE")
  #################
  
  if (!is.na(CF)) {
    if (!is.null(values))
      opt <- try(interpRSE(values))
    else
      ## this path is never taken in current version 2017-09-27
      opt <- optimize(oneRSE,  interval = range(R), k = k, traps = traps,
                      xsigma = xsigma, detectpar = detectpar, noccasions = noccasions,
                      nrepeats = nrepeats, detectfn = detectfn, maskargs = maskargs,
                      CF = CF)
  }
  #################
  if (fit.function %in% c("secr.fit")) {
    args$details <- list(distribution = distribution)   ## blocks other supplied details
    simRSE <- simRSEfn (simulationR, k, traps, xsigma, detectpar, noccasions, nrepeats,
                        detectfn, nx = maskargs$nx, nrepl, fit.function, args)
  }
  else simRSE <- NULL
  #################
  
  rotRSE <- list(values = values)
  if (!is.na(CF)) {
    rotRSE$optimum.spacing <- opt$minimum * detectpar$sigma # spacing(traps)
    rotRSE$optimum.R <- opt$minimum
    rotRSE$minimum.RSE <- opt$objective
  }
  else {
    rotRSE$optimum.spacing <- rotRSE$optimum.R <- rotRSE$minimum.RSE <- NA
  }
  
  out <- list(rotRSE = rotRSE,
              simRSE = simRSE)
  attr(out, "noccasions") <- noccasions
  attr(out, "traps") <- traps
  attr(out, "nrepeats") <- nrepeats
  attr(out, "detectpar") <- detectpar
  attr(out, "detectfn") <- detectfn
  
  class(out) <- c("optimalSpacing", "list")
  
  if (plt) {
    ## fixed 2017-08-24
    args$x <- out
    if (is.null(args$add)) args$add <- FALSE
    if (is.null(args$plottype)) { 
      if (is.null(out$simRSE)) 
        args$plottype <- "nrm"
      else 
        args$plottype <- "RSE"
    }
    do.call(plot, args)
    invisible(out)
  }
  else {
    out
  }
}
##############################################################################

plot.optimalSpacing <- function (x, add = FALSE, plottype = c("both", "RSE", "nrm", "RB", "RMSE"), 
                                 xtype = c('relative','absolute'), xoffset = 0, ...) {
  ## need to define missing cases
  args <- list(...)
  plottype <- match.arg(plottype)
  xtype <- match.arg(xtype)
  if (plottype == "nrm") {
    y <- x$rotRSE$values$n + x$rotRSE$values$r
  }
  else if (plottype == "RMSE") {
    if(is.null(x$simRSE)) stop ("simulations required for RMSE")
    y <- x$simRSE$summary$rRMSE
  }
  else if (plottype == "RB") {
    if(is.null(x$simRSE)) stop ("simulations required for RB")
    y <- x$simRSE$summary$RB.mean
  }
  else {
    y <- x$rotRSE$values$RSE
    if(is.null(x$simRSE) && all(is.na(y))) {
      warning ("RSE all NA")
    }
  }
  
  R <- x$rotRSE$values$R
  if (xtype == 'absolute') {
    sigma <- attr(x, 'detectpar')$sigma
  }
  else {
    sigma <- 1
  }
  if (!add) {
    miny <- 0
    maxy <- 0.5
    if (!all(is.na(y))) {
      miny <- min(c(0,y), na.rm = TRUE)*1.3
      maxy <- max(y, na.rm = TRUE)*1.3
    }
    maxx <- max(R, na.rm = TRUE)
    minx <- min(R, na.rm = TRUE)
    if (minx < 0.2 * (maxx-minx)) minx <- 0
    defaultargs <- list(x = 0, y = 0, type = "n", las = 1,
                        xlab = expression(paste("Spacing -  ", sigma, "  units")),
                        ylab = expression(paste("RSE ", hat(italic(D)))),
                        ylim = c(miny, maxy),
                        xlim = c(minx, maxx))
    if (plottype == 'nrm') defaultargs$ylab <- "Number"
    if (plottype == 'RB') defaultargs$ylab <- expression(paste("RB ", hat(italic(D))))
    if (plottype == 'RMSE') defaultargs$ylab <- expression(paste("rRMSE ", hat(italic(D))))
    if (xtype == 'absolute') {
      defaultargs$xlab <- "Spacing -  m"
      minx <- minx * sigma
      maxx <- maxx * sigma
      defaultargs$ylim <- c(0, maxy)
      defaultargs$xlim <- c(minx, maxx)
    }
    dotsargs <- args[names(args) %in% c("xlab", "ylab", "xlim", "ylim", "las",
                                        "xaxs", "yaxs")]
    plotargs <- replacedefaults(defaultargs, dotsargs)
    do.call(plot, plotargs)
  }
  
  if (plottype %in% c("both","RSE","RB", "RMSE")) {
    defaultargs <- list(col = "black", lwd = 1, cex = 1, pch = 21)
    dotsargs <- args[names(args) %in% c("col", "lwd", "lty", "cex", "pch", "bg", "type")]
    plotargs <- replacedefaults(defaultargs, dotsargs)
    
    # approximate RSE
    if (plottype == "both") {
      plotargs$x <- R * sigma + xoffset
      plotargs$y <- y
      do.call(lines, plotargs)
    }
    
    # simulated RSE, RB, RMSE
    if (!is.null(x$simRSE)) {
      plotargs$x <- x$simRSE$summary$R * sigma  + xoffset
      if (plottype %in% c("both","RSE")) {
        plotargs$y <- x$simRSE$summary$RSE.mean
        segments(plotargs$x, plotargs$y - 2 * x$simRSE$summary$RSE.se,
                 plotargs$x, plotargs$y + 2 * x$simRSE$summary$RSE.se)
        do.call(points, plotargs)
      }
      if (plottype == "RB") {
        plotargs$y <- x$simRSE$summary$RB.mean
        segments(plotargs$x, plotargs$y - 2 * x$simRSE$summary$RB.se,
                 plotargs$x, plotargs$y + 2 * x$simRSE$summary$RB.se)
        do.call(points, plotargs)
      }
      if (plottype == "RMSE") {
        plotargs$y <- x$simRSE$summary$rRMSE
        do.call(points, plotargs)
      }
    }
  }
  if (plottype == "nrm") {
    defaultargs <- list(col = "blue", lwd = 1, cex = 1, pch = 16, type = 'l')
    dotsargs <- args[names(args) %in% c("col", "lwd", "lty", "cex", "pch", "type")]
    plotargs <- replacedefaults(defaultargs, dotsargs)
    
    plotargs$x <- R * sigma
    plotargs$y <- x$rotRSE$values$n
    do.call(lines, plotargs)
    
    plotargs$col <- "red"
    plotargs$y <- x$rotRSE$values$r
    do.call(lines, plotargs)
    
    plotargs$lty <- 2
    plotargs$y <- x$rotRSE$values$m
    do.call(lines, plotargs)
    
    legend (x = "top", legend = c("n","r","m"), lty=c(1,1,2), col=c("blue","red","red"), horiz = TRUE,
            cex = par()$cex * 1.2)
  }
  
}

print.optimalSpacing <- function (x, ...) {
  attr(x,"class") <- NULL
  print(x)
}

##############################################################################

minsimRSE <- function (object, ...) UseMethod("minsimRSE")
minsimRSE.optimalSpacing <- function (object, cut = 0.2, plt = FALSE, 
                                      verbose = FALSE, incr = 0.1, ...) {
  if (is.null(object$simRSE)) stop ("requires optimalSpacing object with simulations")
  object$simRSE$summary <- object$simRSE$summary[
    !(is.na(object$simRSE$summary$RSE.mean) | object$simRSE$summary$RSE.mean>=10),]
  sumy <- object$simRSE$summary   
  ok <- sumy$RSE.mean <= min(sumy$RSE.mean)*(1+cut)
  sumy <- sumy[ok,]
  lm1 <- lm(RSE.mean ~ R + I(R^2), data = sumy)
  cba <- coef(lm1)
  R <- c(-cba[2]/2/cba[3])  # -b/2a
  rse <- predict(lm1, newdata=data.frame(R=R))
  newR <- seq(min(sumy$R), max(sumy$R), incr)
  newRSE <- predict(lm1, newdata = data.frame(R = newR))
  if (plt) {
    object$rotRSE$values$RSE[] <- NA  # suppress ROT values for plot
    plot(object, ...)
    lines (newR, newRSE)
  }
  if (verbose) {
    out <- list(
      model  = lm1, 
      fitted = data.frame(R = newR, RSE = newRSE),
      R      = unname(R), 
      RSE    = unname(rse))
  }
  else {
    out <- c(R = unname(R), RSE = unname(rse))
  }
  if (plt) invisible(out) else out
}
##############################################################################