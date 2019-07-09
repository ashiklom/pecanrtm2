#' @title Calculate burnin value
#' 
#' @description Automatically detect burnin based on one of several methods.
#' @param jags_out List of MCMC sample matrices or `mcmc.list` object
#' @param threshold Maximum value of Gelman diagnostic
#' @param method Character string indicating method. Options are 
#' "moving.window" (default) or "gelman.plot".
#' @param use.confidence Logical. If TRUE (default), use 95% confidence 
#' interval for Gelman Diagnostic. If FALSE, use the point estimate.
#' @param ... Other parameters to methods
#' 
#' @details 
#' See "gelman_diag_mw" and "gelman_diag_gelmanPlot"
#' 
#' @examples
#'      z1 <- coda::mcmc(c(rnorm(2500, 5), rnorm(2500, 0)))
#'      z2 <- coda::mcmc(c(rnorm(2500, -5), rnorm(2500, 0)))
#'      z <- coda::mcmc.list(z1, z2)
#'      burnin <- getBurnin(z, threshold = 1.05)
#' @author Alexey Shiklomanov, Michael Dietze
#' @export

getBurnin <- function(jags_out,
                      threshold = 1.1, 
                      use.confidence = TRUE,
                      method = "moving.window",
                      plotfile = "/dev/null",
                      ...) {
  if (method == "moving.window") {
    GBR <- try(gelman_diag_mw(jags_out, ...))
  } else if (method == "gelman.plot") {
    GBR <- try(gelman_diag_gelmanPlot(jags_out, ...))
  } else {
    stop("Unknown method: ", method)
  }
  if (class(GBR) == "try-error") {
    message("Unable to calculate Gelman diagnostic. Assuming no convergence.")
    return(1)
  }
  column <- ifelse(use.confidence, 2, 1)
  gbr_values <- GBR[, -(1:2), column, drop = FALSE]
  gbr_exceed <- gbr_values > threshold
  if (all(!gbr_exceed)) {
    # Chains converged instantly -- no burnin required
    burnin <- 2     # This isn't 1 to allow testing for convergence with `burnin == 1`
  } else {
    index <- utils::tail(which(rowSums(gbr_exceed) > 0), 1) + 1
    stopifnot(length(index) == 1,
              class(index) %in% c("numeric", "integer"))
    if (index > dim(GBR)[1]) {
      burnin <- NA
    } else {
      burnin <- GBR[index, "Start", column]
    }
  }
  if (is.na(burnin)) {
    message("*** Chains have not converged yet ***")
    mvals <- as.data.frame(matrix(gbr_values, nrow(gbr_values), ncol(gbr_values)))
    colnames(mvals) <- colnames(gbr_values)
    mex <- as.data.frame(matrix(gbr_exceed, nrow(gbr_exceed), ncol(gbr_exceed)))
    colnames(mex) <- sprintf("PSRF %s > %.2f", colnames(gbr_exceed), threshold)
    print(cbind(utils::tail(mvals), utils::tail(mex)))
    burnin <- 1
  }
  return(burnin)
} # getBurnin

#' @title Automatically calculate and apply burnin value
#'
#' @author Michael Dietze, Alexey Shiklomanov
#' @param return.burnin Logical. If `TRUE`, return burnin value in addition to 
#' samples (as list). Default = FALSE.
#' @param ... Additional arguments for \code{getBurnin}, \code{gelman_diag_mw}, 
#' and \code{gelman.diag}.
#' @inheritParams getBurnin
#' @examples
#'      z1 <- coda::mcmc(c(rnorm(2500, 5), rnorm(2500, 0)))
#'      z2 <- coda::mcmc(c(rnorm(2500, -5), rnorm(2500, 0)))
#'      z <- coda::mcmc.list(z1, z2)
#'      z_burned <- autoburnin(z)
#' @export
autoburnin <- function(jags_out, return.burnin = FALSE, ...) {
  burnin <- getBurnin(jags_out, ...)
  if (burnin == 1) {
    samples <- jags_out
  } else if (burnin > 1) {
    samples <- stats::window(jags_out, start = burnin)
  } else {
    stop("Bad return value for burnin: \n",
         burnin)
  }
  if (return.burnin) {
    out <- list(samples = samples, burnin = burnin)
  } else {
    out <- samples
  }
  return(out)
} # autoburnin

#' @title Make MCMC list from samples list
#'
#' @param samps samples list (output from invert.custom)
#' @export
makeMCMCList <- function(samps) {
  samps.mcmc <- lapply(samps, coda::mcmc)
  stopifnot(all(sapply(samps.mcmc, coda::is.mcmc)))
  samps.mcmc.list <- coda::mcmc.list(samps.mcmc)
  stopifnot(coda::is.mcmc.list(samps.mcmc.list))
  return(samps.mcmc.list)
} # makeMCMCList

##' Create priors for BayesianTools
##' 
##' Helper function for creating log-priors compatible with BayesianTools package
##'
##' @param prior.sel `data.frame` containing prior distributions of the selected parameters
##'
##' @return out Prior class object for BayesianTools package
##' @details `prior.sel` must contain the following columns:
##'   * `distn` -- String describing a distribution; e.g. `norm` for `dnorm`, `rnorm`, etc.
##'   * `parama`, `paramb` -- First and second parameters, respectively, of the corresponding distribution
##'
##' Optionally, `prior.sel` may also contain the following columns:
##'   * `param_name` -- Parameter name, which will be carried through to the prior object and sampler
##'   * `lower`, `upper` -- Lower and upper bounds, respectively. These can be leveraged by the BayesianTools samplers.
##'   * `best` -- Best guess for a parameter estimate. BayesianTools can also use this, though I'm not sure how...
##'
##' @author Istem Fer, Alexey Shiklomanov
##' @export
pda.create.btprior <- function(prior.sel) {
  
  # TODO: test exponential -- it only has one argument, so this won't work

  # Returns a function that calculates the density of the specified 
  # distribution given the parameters
  ddist_generator <- function(distn, a, b) {
    fun_string <- paste0('d', distn)
    f <- match.fun(fun_string)
    out <- function(x) f(x, a, b, log = TRUE)
    return(out)
  }

  # Returns a function that draws from the specified distribution with the 
  # specified parameters
  rdist_generator <- function(distn, a, b) {
    fun_string <- paste0('r', distn)
    f <- match.fun(fun_string)
    out <- function(n = 1) f(n, a, b)
    return(out)
  }

  # Create a list of density and random draw functions
  ddist_funs <- with(prior.sel, mapply(ddist_generator, distn, parama, paramb))
  rdist_funs <- with(prior.sel, mapply(rdist_generator, distn, parama, paramb))
  if ('param_name' %in% names(prior.sel)) {
    names(ddist_funs) <- names(rdist_funs) <- prior.sel[['param_name']]
  }

  # `mapply` statement returns
  density <- function(params) {
    dens_vec <- mapply(function(f, x) f(x), ddist_funs, params)   # Returns vector of log densities
    out <- sum(dens_vec)
    return(out)
  }

  # Returns vector of random draws
  sampler <- function() {
    out <- vapply(rdist_funs, function(f) f(), numeric(1))
    return(out)
  }

  # BayesianTools lower and upper bounds and best guess, if specified in data.frame
  lower <- NULL
  if ('lower' %in% names(prior.sel)) {
      lower <- prior.sel[['lower']]
  }
  upper <- NULL
  if ('upper' %in% names(prior.sel)) {
      upper <- prior.sel[['upper']]
  }
  best <- NULL
  if ('best' %in% names(prior.sel)) {
      best <- prior.sel[['best']]
  }
  
  # Use createPrior{BayesianTools} function to create prior class object compatible
  # with rest of the functions
  out <- BayesianTools::createPrior(density = density, sampler = sampler, 
                                    lower = lower, upper = upper, best = best)
  return(out)
} # pda.create.btprior

#' @title Calculate Gelman diagnostic on moving window
#'
#' @author Alexey Shiklomanov
#' @param x MCMC samples, of class \code{mcmc} or \code{mcmc.list}
#' @param width_fraction Fractional width of moving window. Default=0.1.
#' @param width Width of moving window. Default is niter(x)*width_fraction
#' @param njump Number of windows to calculate over
#' @param include.mpsrf Whether to calculate multivariate PSRF and include in output (default = FALSE).
#' @return Gelman Diagnostic 3D array. First dim -- mean (1) and 95% confidence (2). Second dim -- iteration
#' @export
gelman_diag_mw <- function(x,
                           width_fraction = 0.1,
                           width = ceiling(coda::niter(x)*width_fraction),
                           njump = 50,
                           include.mpsrf = TRUE,
                           ...) {

  stopifnot(class(x) %in% c("mcmc", "mcmc.list"))
  stopifnot(width %% 1 == 0)
  stopifnot(njump %% 1 == 0)
  startx <- stats::start(x)
  endx <- stats::end(x)
  a <- floor(seq(startx, endx - width + 1, length.out = njump))
  b <- ceiling(seq(startx + width - 1, endx, length.out = njump))
  if (length(a) < 1) {
    stop("Start index vector has length 0")
  }
  if (length(b) < 1) {
    stop("End index vector has length 0")
  }
  if (length(a) != length(b)) {
    stop("Start and end index vector length mismatch.\n",
         "Start length = ", length(a), "\n",
         "End length = ", length(b))
  }
  n_row <- length(a)
  n_col <- coda::nvar(x) + 2
  vnames <- coda::varnames(x)
  if (is.null(vnames)) {
    vnames <- paste0("V", seq_len(coda::nvar(x)))
  }
  col_names <- c("Start", "End", vnames)
  if (include.mpsrf) {
    n_col <- n_col + 1
    col_names <- c(col_names, "mpsrf")
  }
  gdmat <- array(numeric(), c(n_row, n_col, 2))
  dimnames(gdmat)[[2]] <- col_names
  gdmat[,1,] <- a
  gdmat[,2,] <- b
  for (i in seq_len(n_row)) {
    xsub <- stats::window(x, start=a[i], end=b[i])
    gd_raw <- coda::gelman.diag(xsub, 
                                autoburnin=FALSE,
                                multivariate = include.mpsrf)
    gd <- gd_raw$psrf
    if (include.mpsrf) {
      gd <- rbind(gd, "mpsrf" = rep(gd_raw$mpsrf, 2))
    }
    gdmat[i, -(1:2), ] <- gd
  }
  return (gdmat)
} # gelman_diag_mw

#' @title Calculate Gelman Diagnostic using coda::gelman.plot
#' 
#' @author Alexey Shiklomanov
#' @inheritParams gelman_diag_mw
#' @description Calculates Gelman diagnostic cumulatively. This is a much 
#' more conservative approach than the moving-window method.
#' @export
gelman_diag_gelmanPlot <- function(x, ...) {
  grDevices::pdf(file = NULL)
  GBR_raw <- coda::gelman.plot(x)
  grDevices::dev.off()
  GBR <- array(numeric(), dim(GBR_raw$shrink) + c(0, 2, 0))
  dimnames(GBR)[[2]] <- c("Start", "End", dimnames(GBR_raw$shrink)[[2]])
  GBR[,-(1:2),] <- GBR_raw$shrink
  GBR[, 2, ] <- GBR_raw$last.iter
  GBR[, 1, 1] <- GBR[, 1, 2] <- c(1, GBR[-nrow(GBR), 2, 1] + 1)
  return(GBR)
}
