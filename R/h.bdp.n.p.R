#' Computes smallest necessary subsample size \emph{h} for multivariate robust estimators
#'
#' Computes the smallest necessary uncontaminated subsample size to ensure desired bdp
#' for given sample size \emph{n} and number of variable \emph{p}
#'
#' @usage h.bdp.n.p(bdp, n, p)
#'
#' @param p numeric; number of columns.
#' @param n numeric; sample size.
#' @param bdp numeric value between 0 and 0.5; nominal breakdown point.
#'
#' @return numeric scalar
#'
#' @seealso \code{\link{covMdpd}}, \code{\link{h.alpha.n}}, \code{\link{bdp.max}}, \code{\link{bdp.p.alpha}}, \code{\link{alpha.p.bdp}}
#'
#' @examples
#' h.bdp.n.p(bdp = 0.01, n = 100, p = 20)
#' h.bdp.n.p(bdp = 0.18, n = 50, p = 5)

h.bdp.n.p <- function(bdp, n, p) {
  return(robustbase::h.alpha.n(1 - bdp, n, p))
}
