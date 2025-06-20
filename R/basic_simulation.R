estimate_parameters <- function(counts) {

}

#' Simulate counts from the negative binomial
#' @param mu A matrix of means - matches the dimension of the count matrix
#' @param phi A vector of tagwise dispersions - length matches `nrow(mu)`
#' @param S A vector of samplewise normalization factors -
#' length matches `ncol(mu)`
simulate_counts <- function(mu, phi, S) {
    N <- ncol(mu)
    G <- nrow(mu)
    if (!length(phi) %in% c(1, G)) stop("phi must be length 1 or nrow(mu)")
    if (!length(S) %in% c(1, N)) stop("S must be length 1 or ncol(mu)")
    counts <- matrix(nrow = G, ncol = N)
    # generate samplewise counts
    for (n in seq(1, N)) {
        counts[, n] <- rnbinom(n = G, size = 1 / phi, mu = mu[, n] * S)
    }
    return(counts)
}

#' Create matrix of means
#' @param X_g Tagwise design matrix
#' @param beta Matrix of mean parameters -
#' `nrow(beta) == G` and `ncol(beta) == ncol(X_g)`
create_means <- function(X_g, beta) {
    if (ncol(beta) != ncol(X_g)) stop("ncol(beta) != ncol(X_g)")
    mu <- matrix(nrow = nrow(beta), ncol = ncol(beta))
    for (g in seq(1, nrow(beta))) {
        mu[g, ] <- X_g %*% beta[g, ]
    }
    return(mu)
}
