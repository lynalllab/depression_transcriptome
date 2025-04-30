
  
# source = https://slowkow.com/notes/ggplot2-qqplot/

#' @param ps Vector of p-values.
#' @param ci Size of the confidence interval, 95% by default.
#' @return A ggplot2 plot.
#' @examples
#' library(ggplot2)
#' gg_qqplot(runif(1e2)) + theme_grey(base_size = 24)
gg_qqplot <- function(ps, ci = 0.95, name) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) + 
    labs(title = name)
}

pdf("qq_plots.pdf")
gg_qqplot(ps = as.vector(bacon_meta_dge_pooled$pvalue.pval_BacWeightedZ), name = "Pooled")
gg_qqplot(ps = as.vector(bacon_meta_dge_male$pvalue.pval_BacWeightedZ), name = "Male")
gg_qqplot(ps = as.vector(bacon_meta_dge_female$pvalue.pval_BacWeightedZ), name = "Female")

gg_qqplot(ps = as.vector(bacon_meta_ncc_dge_pooled$pvalue.pval_BacWeightedZ), name = "Pooled ncc")
gg_qqplot(ps = as.vector(bacon_meta_ncc_dge_male$pvalue.pval_BacWeightedZ), name = "Male ncc")
gg_qqplot(ps = as.vector(bacon_meta_ncc_dge_female$pvalue.pval_BacWeightedZ), name = "Female ncc")

gg_qqplot(ps = as.vector(bacon_meta_sens_dge_pooledl$pvalue.pval_BacWeightedZ), name = "Pooled sensitivity")
gg_qqplot(ps = as.vector(bacon_meta_sens_dge_male$pvalue.pval_BacWeightedZ), name = "Male sensitivity")
gg_qqplot(ps = as.vector(bacon_meta_sens_dge_female$pvalue.pval_BacWeightedZ), name = "Female sensitivity")

gg_qqplot(ps = as.vector(bacon_meta_dte_pooledl$pvalue.pval_BacWeightedZ), name = "Pooled DTE")
gg_qqplot(ps = as.vector(bacon_meta_dte_male$pvalue.pval_BacWeightedZ), name = "Male DTE")
gg_qqplot(ps = as.vector(bacon_meta_dte_fem$pvalue.pval_BacWeightedZ), name = "Female DTE")

dev.off()


