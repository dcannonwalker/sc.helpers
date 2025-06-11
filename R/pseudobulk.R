#' Aggregate single cell counts to create pseudobulk samples
#'
#' @param y A `DGEList` object
#' @export
create_pseudobulk_samples <- function(y) {
    # consider rewriting this to use a new S4 class that inherits
    # from DGEList
    grp <- y$samples$group
    group_mat <- Matrix::sparse.model.matrix(~0 + grp)
    counts.pb <- y$pseudodata$pseudo.counts %*% group_mat
    s.pb <- levels(grp)
    DGEList(counts = as.matrix(counts.pb), samples = s.pb)
}
