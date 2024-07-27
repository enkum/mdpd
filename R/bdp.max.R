#' Computes maximum nominal breakdown point for DPD
#'
#' Computes the maximum nominal breakdown point for DPD given \emph{p}.
#'
#' @usage bdp.max(p)
#'
#' @param p numeric; number of variables.
#'
#' @return numeric scalar
#'
#' @seealso \code{\link{covMdpd}}, \code{\link{bdp.p.alpha}}, \code{\link{alpha.p.bdp}}, \code{\link{h.bdp.n.p}}
#'
#' @examples
#' bdp.max(p = 5)
#' bdp.max(p = 15)

bdp.max <- function(p) {
  return(2.0*(p/(p + 2.0))^(p/2.0)/(p + 2.0))
}



