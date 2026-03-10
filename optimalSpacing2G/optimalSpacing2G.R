##############################################################################
## Two-group optimal spacing (absolute spacing search)
## Poisson only, no simulation block
## Objective: maximise either min(En1 + En2, Er1 + Er2) or min(En1, En2, Er1, Er2)
##############################################################################

optimalSpacing2G <- function(
    D1, D2,
    traps0,
    detectpar1, detectpar2,
    noccasions,
    nrepeats = 1,
    detectfn1 = "HHN",
    detectfn2 = "HHN",
    xsigma = 4,
    spacing_m = seq(200, 4000, 200),
    criterion = c("sum_min", "all_min"),
    CF = 1.0,
    ...
) {

  criterion <- match.arg(criterion)

  detectfn1 <- secr:::secr_valid.detectfn(detectfn1, valid = c(0, 1, 2, 14:19))
  detectfn2 <- secr:::secr_valid.detectfn(detectfn2, valid = c(0, 1, 2, 14:19))
  dfc1 <- secrdesign:::dfcast(detectfn1, detectpar1)
  dfc2 <- secrdesign:::dfcast(detectfn2, detectpar2)
  detectfn1 <- dfc1$detectfn
  detectfn2 <- dfc2$detectfn
  detectpar1 <- dfc1$detectpar
  detectpar2 <- dfc2$detectpar

  base_spacing <- spacing(traps0)

  args <- list(...)
  defaultmaskargs <- list(nx = 64, type = "trapbuffer")
  dotsmask <- args[names(args) %in% c("nx", "type", "poly", "poly.habitat")]
  maskargs_base <- secrdesign:::replacedefaults(defaultmaskargs, dotsmask)

  if (any(detector(traps0) == "single")) {
    warning("treating single-catch traps as multi-catch", call. = FALSE)
    detector(traps0) <- "multi"
  }

  values <- lapply(
    spacing_m,
    getCrit2G_abs,
    traps0,
    base_spacing,
    xsigma,
    detectpar1,
    detectpar2,
    D1,
    D2,
    noccasions,
    nrepeats,
    detectfn1,
    detectfn2,
    maskargs_base,
    CF,
    criterion
  )
  values <- do.call(rbind, values)
  values <- as.data.frame(values)
  names(values) <- c("spacing", "En1", "En2", "Er1", "Er2", "crit")

  opt <- if (!is.na(CF) && nrow(values) > 0 && diff(range(values$spacing)) > 0) {
    interpCritMax(values)
  } else {
    list(minimum = NA, objective = NA)
  }

  out <- list(
    values = values,
    optimum.spacing = opt$minimum,
    maximum.crit = opt$objective,
    criterion = criterion,
    traps.base = traps0,
    detectpar1 = detectpar1,
    detectpar2 = detectpar2
  )
  class(out) <- c("optimalSpacing2G", "list")
  out
}

##############################################################################

getCrit2G_abs <- function(
    S,
    traps0,
    base_spacing,
    xsigma,
    detectpar1,
    detectpar2,
    D1,
    D2,
    noccasions,
    nrepeats,
    detectfn1,
    detectfn2,
    maskargs_base,
    CF,
    criterion
) {

  scalefac <- S / base_spacing
  trapS <- scale_traps(traps0, scalefac)

  maskargs <- maskargs_base
  maskargs$traps <- trapS
  maskargs$buffer <- xsigma * max(detectpar1$sigma, detectpar2$sigma) + S
  mask <- do.call(make.mask, maskargs)

  En1 <- Enrm(D1, trapS, mask, detectpar1, noccasions, detectfn1) * nrepeats
  En2 <- Enrm(D2, trapS, mask, detectpar2, noccasions, detectfn2) * nrepeats
  Er1 <- En1[2]
  Er2 <- En2[2]
  En1 <- En1[1]
  En2 <- En2[1]

  critval <- switch(
    criterion,
    sum_min = min(En1 + En2, Er1 + Er2),
    all_min = min(En1, En2, Er1, Er2)
  )

  c(S, En1, En2, Er1, Er2, critval * CF)
}

##############################################################################

interpCritMax <- function(values) {
  f <- approxfun(values$spacing, values$crit, rule = 2)
  opt <- optimize(f, interval = range(values$spacing), maximum = TRUE)
  list(minimum = opt$maximum, objective = opt$objective)
}

##############################################################################

scale_traps <- function(traps, scalefac) {
  tr <- traps
  tr$x <- tr$x * scalefac
  tr$y <- tr$y * scalefac
  for (nm in c("spacing", "spacex", "spacey")) {
    if (!is.null(attr(tr, nm))) {
      attr(tr, nm) <- attr(tr, nm) * scalefac
    }
  }
  tr
}

##############################################################################
