#' Generate pseudodata from single cell counts
#' @param y A matrix or `DGEList` object
generate_sc_pseudodata <- function(y, ...) {
    UseMethod("generate_sc_pseudodata")
}

#' Generate pseudodata from single cell counts stored in a `DGEList`
#' @param y A `DGEList` object
generate_sc_pseudodata.DGEList <- function(y, prior.df=10,
                                           niter = 4,
                                           tol=1e-06, verbose=FALSE, ...) {
    ##### Basically a complete rip of `edgeR::estimateTagwiseDisp.DGEList()` #####
    y <- validDGEList(y)
    group <- y$samples$group
    lib.size <- y$samples$lib.size * y$samples$norm.factors
    dispersion <- y$common.dispersion
    if(is.null(dispersion)) stop("No common.dispersion found in the DGEList object. Run estimateCommonDisp first.")

    y$prior.df <- prior.df
    nlibs <- ncol(y$counts)
    ngroups <- length(unique(group))
    y$prior.n <- prior.df/(nlibs - ngroups)
    ##### End Copy #####

    # edit line
    out <- generate_sc_pseudodata(y$counts, group=group,
                                  lib.size=lib.size,
                                  dispersion=dispersion,
                                  prior.df=prior.df,
                                  tol=tol, niter = niter,
                                  verbose=verbose)
}

#' Generate pseudodata from a matrix of single cell counts
#' @param y A matrix of single cell counts
generate_sc_pseudodata.default <- function(y, group=NULL, lib.size=NULL,
                                           dispersion,
                                           prior.df=10,
                                           tol=1e-06, niter = 4,
                                           verbose=FALSE, ...) {
    ##### Basically a complete rip of `edgeR::estimateTagwiseDisp.default()` #####
    #	Check y
    y <- as.matrix(y)
    ntags <- nrow(y)
    nlibs <- ncol(y)

    #	Check group
    if(is.null(group)) group <- rep(1, nlibs)
    if(length(group)!=nlibs) stop("Incorrect length of group.")
    group <- dropEmptyLevels(group)

    #	Check lib.size
    if(is.null(lib.size)) lib.size <- colSums(y)
    if(length(lib.size)!=nlibs) stop("Incorrect length of lib.size.")

    for (i in seq(1, niter)) {
        disp.out <- do_iter(
            y = y, group = group, dispersion = dispersion,
            lib.size = lib.size, ntags = ntags, prior.df = prior.df,
            nlibs = nlibs,
            tol = tol, verbose = verbose
        )
        if (i > 1) if (verbose) message(
            "iter: ", i,
            "\n",
            "norm of difference in dispersion: ",
            sqrt(sum((dispersion - disp.out)^2))
        )
        dispersion <- disp.out
    }

    list(
        tagwise.dispersion = disp.out,
        pseudodata = equalizeLibSizes(y, group=group,
                                      dispersion=disp.out,
                                      lib.size=lib.size)
    )
}

do_iter <- function(y, group, dispersion, lib.size,
                    ntags, prior.df, nlibs,
                    tol,
                    verbose) {
    eq <- equalizeLibSizes(y, group=group, dispersion=dispersion, lib.size=lib.size)
    u <- splitIntoGroups(eq$pseudo.counts, group=group)
    delta <- rep(0, ntags)
    ngroups <- length(unique(group))
    prior.n <- prior.df/(nlibs - ngroups)

    if(verbose) message("Tagwise dispersion optimization begun, may be slow, progress reported every 100 tags")
    for(tag in seq_len(ntags)) {
        delta.this <- optimize(weightedCondLogLikDerDelta,
                               interval=c(1e-4,100/(100+1)),
                               tol=tol, maximum=TRUE, y=u,
                               tag=tag, ntags=ntags, prior.n=prior.n, der=0L)
        delta[tag] <- delta.this$maximum
        if(verbose) if(tag%%100==0) message("tag ",tag)
    }
    tagwise.dispersion <- delta/(1-delta)
    if(verbose) cat("\n")

    tagwise.dispersion
}
