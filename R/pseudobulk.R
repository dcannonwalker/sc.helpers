#' Aggregate single cell counts to create pseudobulk samples
#' @param y A `DGEList` object of single cell counts with pseudodata attached
#' @export
create_pb_samples <- function(y) {
    # consider rewriting this to use a new S4 class that inherits
    # from DGEList
    grp <- y$samples$group
    group_mat <- Matrix::sparse.model.matrix(~0 + grp)
    colnames(group_mat) <- gsub("^grp", "", colnames(group_mat))
    counts.pb <- y$pseudodata$pseudo.counts %*% group_mat
    s.pb <- levels(grp)
    DGEList(counts = as.matrix(counts.pb), samples = s.pb)
}

#' Create pseudobulk dispersion matrix
#' @param y A `DGEList` object of single cell counts with fitted dispersions
#' @export
create_pb_dispersion <- function(y) {
    ncells <- table(y$samples$group)
    t(sapply(y$tagwise.dispersion, function(d) d / ncells))
}

#' Get pseudobulk normalization factors
#' @inheritParams create_pb_dispersion
#' @export
get_pb_nf <- function(y) {
    ncells <- table(y$samples$group)
    ncells / median(ncells)
}
