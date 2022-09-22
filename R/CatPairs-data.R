#' Categorical Cause-Effect Pairs
#'
#' Cause-effect pairs extracted from R packages MASS and datasets for which the pairwise causal relationships are clear from the context, and at least one of the variables in each pair is categorical. For non-categorical variable, we discretized it at 5 evenly spaced quantiles.The current version contains 33 categorical cause-effect pairs.
#'
#' @docType data
#'
#' @usage data(CatPairs)
#'
#' @format A list of length 2. The first element is a list of 33 cause-effect pairs as data frames with the first column being the cause and the second column being the effect. The second element is a list of sources of each pair.
#'
#' @keywords datasets
"CatPairs"
