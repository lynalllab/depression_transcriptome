#no longer available on cran or updated on github, so downloading relevant function from:
#https://github.com/unistbig/metapro/blob/master/R/metapro.R
#' @title wZ
#' @description P-value combination based on weighted Z-method
#' @param p A numeric vector of p-values
#' @param weight A numeric vector of weights (e.g., sample sizes)
#' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effect, and vice versa. Default: TRUE.
#' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
#' @return p : Combined p-value
#' @return overall.eff.direction : The direction of combined effects.
#' @return sumz : Sum of transformed z-score
#' @examples wZ(p=c(0.01,0.2,0.8), weight = c(20,10,40), is.onetail=FALSE, eff.sign=c(1,-1,1))
#' @import "metap"
#' @references Becker BJ (1994). "Combining significance levels." In Cooper H, Hedges LV (eds.), A handbook of research synthesis, 215-230. Russell Sage, New York.
#' @references Stouffer SA, Suchman EA, DeVinney LC, Star SA, Williams RMJ (1949). The American soldier, vol 1: Adjustment during army life. Princeton University Press, Princeton.
#' @references Mosteller, F. & Bush, R.R. (1954). Selected quantitative techniques. In: Handbook of Social Psychology, Vol. 1 (G. Lindzey, ed.), pp. 289-334. Addison-Wesley, Cambridge, Mass.
#' @export
wZ = function(p, weight = NULL, is.onetail = TRUE, eff.sign)
{
  if(length(p) == 0){stop("No input p-values exist.")}
  if(is.null(weight)){weight = rep(1, length(p))}
  idx.na = which(is.na(p))
  if(length(idx.na)>0){
    p = p[-idx.na];
    weight = weight[-idx.na]
    if(!is.onetail)
    {
      eff.sign = eff.sign[-idx.na]
    }
  }
  if(is.onetail){res = sumz(p = p, weights = weight)
  }else{
    p1 = p2 = p
    idx_pos = which(eff.sign >= 0)
    idx_neg = which(eff.sign < 0)
    # positive direction
    p1[idx_pos] = p[idx_pos]/2
    p1[idx_neg] = 1 - p[idx_neg]/2
    res1 = sumz(p = p1, weights = weight)
    # negative direction
    p2[idx_pos] = 1 - p[idx_pos]/2
    p2[idx_neg] = p[idx_neg]/2
    res2 = sumz(p = p2, weights = weight)
    
    if(res1$p <= res2$p){res = res1;overall.eff.direction = "+"}else{res = res2; overall.eff.direction = "-"}
  }
  if(is.onetail)
  {
    RES = list(p = res$p[,1], sumz = res$z[,1])
  }else{
    RES = list(p = res$p[,1] * 2 , overall.eff.direction = overall.eff.direction, sumz = res$z[,1])
  }
  
  return(RES)
}


sumz <-
  function(p, weights = NULL, data = NULL, subset = NULL,
           na.action = na.fail, log.p = FALSE, log.input = FALSE)  {
    if(is.null(data)) data <- sys.frame(sys.parent())
    mf <- match.call()
    mf$data <- NULL
    mf$subset <- NULL
    mf$na.action <- NULL
    mf[[1]] <- as.name("data.frame")
    mf <- eval(mf, data)
    if(!is.null(subset)) mf <- mf[subset,]
    mf <- na.action(mf)
    p <- as.numeric(mf$p)
    weights <- mf$weights
    noweights <- is.null(weights)
    if(noweights) weights <- rep(1, length(p))
    if(length(p) != length(weights)) warning("Length of p and weights differ")
    if(log.input) {
      keep <- p < 0
    } else {
      keep <- (p > 0) & (p < 1)
    }
    invalid <- sum(1L * keep) < 2
    if(invalid) {
      warning("Must have at least two valid p values")
      res <- list(z = NA_real_, p = NA_real_,
                  validp = p[keep], weights = weights)
    } else {
      if(sum(1L * keep) != length(p)) {
        warning("Some studies omitted")
        omitw <- weights[!keep]
        if((sum(1L * omitw) > 0) & !noweights)
          warning("Weights omitted too")
      }
      zp <- (qnorm(p[keep], lower.tail = FALSE, log.p = log.input) %*%
               weights[keep]) / sqrt(sum(weights[keep]^2))
      res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE,
                                    log.p = log.p),
                  validp = p[keep], weights = weights)
    }
    class(res) <- c("sumz", "metap")
    res
  }
print.sumz <- function(x, ...) {
  cat("sumz = ", x$z, "p = ", x$p, "\n")
  invisible(x)
}