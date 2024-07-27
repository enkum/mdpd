#' Computes \eqn{\alpha}, power of DPD
#'
#' Computes \eqn{\alpha} for DPD given dimension \emph{p} and nominal breakdown point \emph{bdp}
#'
#' @usage alpha.p.bdp(p, bdp)
#'
#' @param p numeric; number of variables.
#' @param bdp numeric; nominal breakdown point.
#'
#' @return numeric scalar
#'
#' @seealso \code{\link{covMdpd}}, \code{\link{bdp.max}}, \code{\link{bdp.p.alpha}}, \code{\link{h.bdp.n.p}}
#'
#' @examples
#' alpha = alpha.p.bdp(p = 5, bdp = 0.12)
#' alpha
#' alpha = alpha.p.bdp(p = 10, bdp = 0.03)
#' alpha

alpha.p.bdp <- function(p, bdp) {
  if ((bdp < 0.0) || (bdp > bdp.max(2.0))) {
    stop(paste("bdp must be between 0 and", bdp.max(p), "for p =", p))
  }

  alpha = uniroot(function(alpha) bdp.p.alpha(p, alpha) - bdp, lower = 0.0, upper = 2.0/p)$root
}



