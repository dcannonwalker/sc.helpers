#' Get the treatment condition for each sample
#' @param sim A `splatter` pop sim
get_conditions <- function(sim) {
    conditions.df <- data.frame(sample = names(sim$Condition),
                                condition = sim$Condition)
    tapply(conditions.df, conditions.df$sample,
           function(s) unique(s$condition))
}
