#' Get the treatment condition for each sample
#' @param sim A `splatter` pop sim
#' @export
get_conditions <- function(sim) {
    conditions.df <- data.frame(sample = names(sim$Condition),
                                condition = sim$Condition)
    tapply(conditions.df, conditions.df$sample,
           function(s) unique(s$condition))
}

#' Get group labels
#' @param y A `SingleCellExperiment` or `Seurat` object
#' @export
get_groups <- function(y, ...) {
    UseMethod("get_groups")
}

#' Get group labels from a `SingleCellExperiment`
#' @param y A `SingleCellExperiment`
get_groups.SingleCellExperiment <- function(y) {
    x <- colData(y)
    sg <- paste0(x[["Sample"]], x[["Group"]])
    factor(sg, labels = seq(1, length(unique(sg))))
}
