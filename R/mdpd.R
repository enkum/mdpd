#' Robust Estimation of Mean Vector and Covariance Matrix for Contaminated Gaussian Datasets via Density Power Divergence Minimization
#'
#' Computes MDPD estimator via custom M-iteration
#'
#' @usage covMdpd(x)
## raw.only = FALSE, nsamp = "usual",
#' covMdpd(x, alpha = 2.0/ncol(x), h.min = h.bdp.n.p(0.5, nrow(x), ncol(x)), rcond.min = 1E-6, tolIter = 1E-9, maxiter = 500L)
#'
#' ## S3 method for class 'mdpd'
#' print.mdpd(x, digits = getOption("digits"), prefix = "\t",...)
#'
#' ## S3 method for class 'mdpd'
#' plot.mdpd(x, which = c("all", "distance", "qqchi2"), classic = FALSE)
#'
#' @param x multivariate data set; matrix.
## @param raw.only logical; if TRUE,
##@param nsamp A logical value
#' @param alpha numeric value between 0 and 1; represents DPD power.
#' @param h.min numeric; represents the smallest desired best subset size.
#' @param rcond.min numeric; condition number tolerance.
#' @param tolIter numeric; iteration error tolerance.
#' @param maxiter numeric; maximum iteration number. Default value set to 500.
#'
#'
#' @details Consider the \emph{p}-variate DPD objective function
#'
#' \deqn{H_n(\boldsymbol{\theta}) = \int_{\mathbb{R}^{p}} f^{1+\alpha}(\boldsymbol{x}|\boldsymbol{\theta}) \mathrm{d}\boldsymbol{x} -
#' \left(1 + \frac{1}{\alpha} \right) n^{-1} \sum_{i=1}^{n} f^{\alpha}(\boldsymbol{x}_{i}|\boldsymbol{\theta}),}
#' where \eqn{f(\boldsymbol{x}|\boldsymbol{\theta})} is a \emph{p}-variate probability density.
#' Here, we choose the Gaussian density
#' \eqn{f(\boldsymbol{x}|\boldsymbol{\mu}, \boldsymbol{\Sigma}) =  (2\pi)^{-\frac{p}{2}} |\boldsymbol{\Sigma}|^{-\frac{1}{2}}
#' \exp(-\frac{1}{2} \big\|\boldsymbol{\Sigma}^{-1/2}(\boldsymbol{x} - \boldsymbol{\mu})\big\|^{2})} with
#' \eqn{\boldsymbol{\theta} = (\boldsymbol{\mu}, \boldsymbol{\Sigma})}, where
#' \eqn{\boldsymbol{\mu} \in \mathbb{R}^{p}} and \eqn{\boldsymbol{\Sigma} \in \mathrm{SPD}(p \times p)} are the location vector and
#' scatter matrix respectively.
#' For a given sample \eqn{ \boldsymbol{X} = \{x_1, \ldots, x_n\} } from \eqn{ (1- \epsilon)f(\boldsymbol{x}|\boldsymbol{\theta}) +
#' \epsilon \mathrm{d}G(\boldsymbol{x}) } with \eqn{\epsilon \in [0, 1/2) } and arbitrary unknown cdf \eqn{G(\boldsymbol{x})},
#' the minimum density power divergence estimator is given as:
#' \deqn{\boldsymbol{\hat\theta}^{\mathrm{DPD}}_{n} := \underset{\boldsymbol{\theta} \in \Theta}{\arg \min} \quad H_{n}(\boldsymbol{\theta})}
#'
#'
#' @return A list with components
#' \item{raw.center }{ Location vector}
#' \item{raw.cov }{ Scatter matrix}
#' \item{raw.weights }{ Weights}
#' \item{raw.mah }{ Squared Mahalanobis distances}
#' \item{crit }{ DPD objective value}
#' \item{iter }{ Number of iterations}
#'
#' @author Michael Pokojovy (\email{michael.pokojovy@@gmail.com}), Andrews T. Anum, Ebenezer Nkum, Abhijit Mandal
#'
#' @keywords Provide keywords here
#'
#' @seealso \code{\link{covMcd}}, \code{\link{bdp.max}}, \code{\link{bdp.p.alpha}}, \code{\link{alpha.p.bdp}}, \code{\link{h.bdp.n.p}}
#'
#' @references Basu, A., Harris, I. R., Hjort, N. L., and Jones, M. (1998). Robust and efficient estimation by minimising a density power divergence. \emph{Biometrika}, \strong{85}(3):549–559.
#'
#' Hubert M., Debruyne M., Rousseeuw P.J. (2018). Minimum covariance determinant and extensions. \emph{WIREs Computational Statistics}, 10:e1421.
#'
#' Rousseeuw, P. J., and van Driessen, K. (1999). A fast algorithm for the minimum covariance determinant estimator. \emph{Technometrics}, \strong{41}, 212–223.
#'
#' @name covMdpd
#'
#' @examples
#' n = 100
#' p = 5
#' x = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p))
#' covMdpd(x) #prints out a list of outputs including the location vector and scatter matrix
#'
#' n = 250
#' p = 10
#' x = MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = diag(p))
#' res = covMdpd(x, alpha = 0.5)
#' print(res)
#' plot(res, which = "distance", classic = TRUE)
#'
#' ## example with outliers
#' set.seed(1)
#' eps = 0.1
#' n = 250
#' p = 5
#' ncp = 150
#' n.out = eps*n
#' out.index = sample(x = n, n.out)
#' x.cont = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p))
#' x.cont[out.index, ] = x.cont[out.index, ] + sqrt(ncp/p)
#' res = covMdpd(x.cont)
#' print(res)
#' plot(res, which = "all", classic = TRUE)

#' @export
covMdpd <- function(x, raw.only = FALSE, nsamp = "usual",
                    alpha = 2.0/ncol(x), h.min = h.bdp.n.p(0.5, nrow(x), ncol(x)),
                    rcond.min = 1E-6, tolIter = 1E-9, maxiter = 500L) {
  p = ncol(x)
  n = nrow(x)

  C1 = (1.0 + alpha)^(-0.5*p)
  C2 = (1.0 + 1.0/alpha)

  mahal <- function(x, mu, iroot.Sigma) {
    rowSums((sweep(x, 2L, mu, check.margin = TRUE) %*% iroot.Sigma)^2) # instead of mahalanobis(x, center = mu, cov = Sigma)
  }

  Hn.crit <- function(log.det, md, alpha) {
    # Computes the DPD objective function
    p = ncol(x)

    wt = exp(-0.5*alpha*md)

    beta = -0.5*alpha*(p*log(2*pi) + log.det)
    C    = C1 - C2*mean(wt)

    return(C*exp(beta))
  }

  mean.shift.MDPD <- function(x, mu0 = colMeans(x), Sigma0 = cov(x)) {
    # Check and amend the initial value (if necessary)
    spec   = svd(Sigma0)
    spec$d = pmax(spec$d, max(spec$d)*rcond.min)

    iroot.Sigma = spec$u %*% diag(1.0/sqrt(spec$d)) %*% t(spec$v)
    md = mahal(x, mu0, iroot.Sigma) # instead of mahalanobis(x, center = mu, cov = Sigma)

    wt0 = rep(0.0, n)

    C1 = (1.0 + alpha)^(-0.5*p)
    C2 = (1.0 + 1.0/alpha)

    ## Iterate
    for (iter in 1:maxiter) {
      wt = exp(-0.5*alpha*md)

      # adjust wt if necessary
      if (sum(wt) < h.min) {
        mult = uniroot(function(mult) sum(exp(-0.5*alpha*md*mult)) - h.min,
                       lower = 0.0, upper = 2.0)$root # upper = 1.0 is theoretically enough, but just to be numerically sure
        #wt = exp(-0.5*alpha*md*mult)
      }

      if (mean((wt - wt0)^2) < tolIter*tolIter) {
        break
      } else {
        wt0 = wt
      }

      C3 = C2*mean(wt) - C1 # Note C > 0 since mean(wt) >= 0.5 and 0 <= alpha <= 2/p
      C  = (C2/n)/C3

      mu    = colMeans(wt*x)/mean(wt)
      Sigma = C*crossprod(sqrt(wt)*sweep(x, 2L, mu, check.margin = TRUE))

      spec   = svd(Sigma)
      spec$d = pmax(spec$d, max(spec$d)*rcond.min)
      iroot.Sigma = spec$u %*% diag(1.0/sqrt(spec$d)) %*% t(spec$v)

      md = mahal(x, mu, iroot.Sigma)
    }

    crit = Hn.crit(sum(log(spec$d)), md, alpha)

    return(list(raw.center = mu, raw.cov = Sigma, raw.weights = wt/sum(wt), raw.mah = md,
                crit = crit, iter = iter))
  }

  #nsamp = 500L

  if (alpha == 0.0) {
    raw.center = colMeans(x)
    raw.cov    = cov(x)

    spec   = Sigma0
    spec$d = pmax(spec$d, max(spec$d)*rcond.min)
    iroot.Sigma = spec$u %*% diag(1.0/sqrt(spec$d)) %*% t(spec$v)
    md = mahal(x, mu0, iroot.Sigma)

    ans = list(raw.center = raw.center, raw.cov = raw.cov, raw.weights = rep(1.0, n)/n, raw.mah = md,
               crit = Hn.crit(sum(log(spec$d)), md, alpha), iter = 0L)
  } else {
    if (nsamp == "usual") {
      ans = mean.shift.MDPD(x, colMeans(x), cov(x))
    } else {
      mcd = robustbase::covMcd(x, alpha = 1.0 - bdp.p.alpha(p, alpha), nsamp = nsamp)
      ans = mean.shift.MDPD(x, mcd$raw.center, mcd$raw.cov)
    }
  }

  METHOD = "Computing MDPD estimator via custom M-iteration"
  DNAME = paste(deparse(substitute(x)))
  names(ans$crit) <- "DPD objective value"
  names(ans$iter) <- "number of iterations"

  rval = list(raw.center = ans$raw.center,
              raw.cov = ans$raw.cov,
              raw.weights	= ans$raw.weights,
              raw.mah = ans$raw.mah,
              method = METHOD,
              crit = ans$crit,
              iter = ans$iter,
              data.name = DNAME,
              Xdata = x
              )
  class(rval) <- "mdpd"
  rval
}

#' @export
#'
print <- function(x, ...){
  UseMethod("print", x)
  NextMethod()
}

#' @export
#'

print.mdpd <- function(x, digits = getOption("digits"), prefix = "\t",...){
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")

  cat("data:  ", x$data.name, "\n", sep = "")
  cat("\n")
  out <- character()
  if(!is.null(x$crit)){
    out <- c(out, paste(names(x$crit), "=",
                        format(signif(x$crit, max(1L, digits - 2L)))))
  }

  if(!is.null(x$iter)){
    out <- c(out, paste(names(x$iter), "=",
                        format(signif(x$iter, max(1L, digits - 2L)))))
  }

  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")

  if(!is.null(x$raw.center)) {
    cat("\n")
    cat("Location vector:\n")
    print(x$raw.center, digits=digits, ...)
  }

  if(!is.null(x$raw.cov)) {
    cat("\n")
    cat("Scatter matrix:\n")
    print(x$raw.cov, digits=digits, ...)
  }

  if(!is.null(x$raw.weights)) {
    cat("\n")
    cat("Weights:\n")
    print(x$raw.weights, digits=digits, ...)
  }

  if(!is.null(x$raw.mah)) {
    cat("\n")
    cat("Mahalanobis distance:\n")
    print(x$raw.mah, digits=digits, ...)
  }

  cat("\n")
  invisible(x)

}

#' @export
#'
plot.mdpd <- function(x, which = c("all", "distance", "qqchi2"),
                      classic = FALSE, cutoff = NULL, id.n, labels.id = 1:nrow(x$Xdata),
                      tol = 1e-7, cex.id = 0.75, label.pos = c(4,2),...){

  mydistplot <- function(distvec, cutoff, classic = FALSE, id.n, ...) {

    n <- length(distvec)
    if (missing(id.n)){
      id.n <- length(which(distvec > cutoff))
    }
    ylab <- paste("Square Root of", if (classic)
      "Mahalanobis"
      else "Robust", "distance")
    plot(distvec, type = "p", ylab = ylab, xlab = "Index", main = "Distance Plot")
    label(1:n, distvec, id.n)
    abline(h = cutoff)
  }

  qqplot <- function(x, p, cutoff = sqrt(qchisq(0.975, p)),
                     classic = FALSE, id.n) {
    n <- length(x)
    if (missing(id.n))
      id.n <- length(which(x > cutoff))
    qq <- sqrt(qchisq(((1:n) - 1/3)/(n + 1/3), p))
    x <- sort(x, index.return = TRUE)
    ix <- x$ix
    x <- x$x
    ylab <- paste(if (classic)
      "Mahalanobis"
      else "Robust", "distance")
    xlab <- "Square root of the quantiles of the Chi-Squared distribution"
    plot(qq, x, xlab = xlab, ylab = ylab, main = "Chi-Squared QQ-Plot")
    label(qq, x, id.n, ind = (n - id.n + 1):n, labs = ix)
    abline(0, 1, lty = 2)
  }

  label <- function(x, y, id.n, ind = sort.list(y, decreasing = TRUE)[1:id.n],
                    labs = labels.id, adj.x = TRUE) {
    if (id.n > 0) {
      labpos <- if (adj.x)
        label.pos[1 + as.numeric(x > mean(range(x)))]
      else 3
      text(x[ind], y[ind], labs[ind], cex = cex.id, xpd = TRUE,
           pos = labpos, offset = 0.25)
    }
  }

  if (!inherits(x, "mdpd"))
    stop("Object must be of class \"mdpd\"")
  if(!is.character(which))
    stop("'which' must be set to either 'all', 'distance' or 'qqchi2' ")

  if (!is.matrix(x$Xdata) || !is.numeric(x$Xdata))
    stop("x is not a numeric dataframe or matrix.")
  n <- dim(x$Xdata)[1]
  p <- dim(x$Xdata)[2]

  if (is.null(cutoff)){
    cutoff <- sqrt(qchisq(0.975, p))
  }

  which <- match.arg(which)
  md <- sqrt(mahalanobis(x$Xdata, colMeans(x$Xdata), var(x$Xdata), tol = tol))
  rres = covMdpd(x$Xdata)
  rd <- sqrt(mahalanobis(x$Xdata, rres$raw.center, rres$raw.cov, tol = tol))

  if (which == "all" || which == "distance") {
    if (classic) {
      opr <- if (prod(par("mfrow")) == 1)
        par(mfrow = c(1, 2), pty = "m")
      else list()
    }
    mydistplot(rd, cutoff, id.n = id.n)
    if (classic) {
      mydistplot(md, cutoff, classic = TRUE, id.n = id.n)
      par(opr)
    }
    on.exit(par(opr))
  }

  if (which == "all" || which == "qqchi2") {
    if (classic) {
      opr <- if (prod(par("mfrow")) == 1)
        par(mfrow = c(1, 2), pty = "m")
      else list()
    }
    qqplot(rd, p, cutoff = cutoff, id.n = id.n)
    if (classic) {
      qqplot(md, p, cutoff = cutoff, classic = TRUE, id.n = id.n)
      par(opr)
    }
    on.exit(par(opr))
  }
}
