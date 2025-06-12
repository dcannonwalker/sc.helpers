#' Get the treatment condition for each sample
#' @param sim A `splatter` pop sim
#' @export
get_conditions <- function(sim) {
    conditions.df <- data.frame(sample = names(sim$Condition),
                                condition = sim$Condition)
    tapply(conditions.df, conditions.df$sample,
           function(s) unique(s$condition))
}

#' Get the design variables for pseudobulk samples
#' @param sim A `splatter` pop sim
#' @export
get_pseudobulk_variables <- function(sim) {
    x <- colData(sim)
    x <- data.frame(Sample = x$Sample,
                     Group = x$Group,
                     Condition = x$Condition)
    x <- tapply(x, ~Sample + Group, function(grp) unique(grp$Condition))
    x <- data.frame(sample = rownames(x), x)
    reshape(x, direction = "long",
            idvar = "sample",
            varying = list(2:3), times = c("Group1", "Group2"),
            v.names = "condition", timevar = "group")
}

#' Get group labels
#' @param y A `SingleCellExperiment` or `Seurat` object
#' @export
get_groups <- function(y, ...) {
    UseMethod("get_groups")
}

#' Get group labels from a `SingleCellExperiment`
#' @param y A `SingleCellExperiment`
#' @export
get_groups.SingleCellExperiment <- function(y) {
    x <- colData(y)
    factor(paste0(x[["Sample"]], ":", x[["Group"]]))
}
