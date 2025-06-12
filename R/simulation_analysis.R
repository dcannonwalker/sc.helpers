#' Get group labels
#' @param y A `SingleCellExperiment` or `DGEGLM` object
#' @export
get_coefs <- function(y, ...) {
    UseMethod("get_coefs")
}

#' Get simulated coefficient values for cluster and condition
#' @inheritParams get_groups.SingleCellExperiment
#' @export
get_coefs.SingleCellExperiment <- function(y) {
    x <- rowData(y)
    cl <- log(x$GroupDE.Group2) - log(x$GroupDE.Group1)
    cn <- log(x$ConditionDE.Condition2) - log(x$ConditionDE.Condition1)
    cbind(clusterGroup2 = cl, conditionCondition2 = cn)
}

#' Get estimated coefficient values for cluster and condition
#' @param y A `DGEGLM` object
#' @export
get_coefs.DGEGLM <- function(y) {

}

#' Compare estimated to simulated parameters
#' @param sim A `splatter` pop sim
#' @param fit An `edgeR` fit object
#' @param coef.name Name of the parameter to compare
#' @param plot Logical; make a plot?
#' @export
compare_coef <- function(sim, fit, coef.name, plot = TRUE) {
    est <- fit$coefficients[, coef.name]
    act <- get_coefs(sim)[, coef.name]
    if (plot) plot(act, est,
                   xlab = paste0("Simulated ", coef.name),
                   ylab = paste0("Estimated ", coef.name))
    out <- cbind(est, act)
    colnames(out) <- paste0(coef.name, c(".est", ".act"))
}
