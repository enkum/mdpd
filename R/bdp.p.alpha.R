#' Computes nominal breakdown point for DPD
#'
#' Computes the nominal breakdown point for DPD given \emph{p} and \eqn{\alpha}
#'
#' @usage bdp.p.alpha(p, alpha)
#'
#' @param p numeric; number of variables.
#' @param alpha numeric value between 0 and 1; represents DPD power.
#'
#' @return numeric scalar
#'
#' @seealso \code{\link{covMdpd}}, \code{\link{bdp.max}}, \code{\link{alpha.p.bdp}}, \code{\link{h.bdp.n.p}}
#'
#' @examples
#' bdp.p.alpha(p = 5, alpha = 0.5)
#' bdp.p.alpha(p = 8, alpha = 0.75)

bdp.p.alpha <- function(p, alpha) {
  alpha/((1.0 + alpha)^(p/2.0 + 1.0))
}



