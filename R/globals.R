utils::globalVariables(c(
  "maindat",
  "interp_n",
  "iterations",
  "ref_u",
  "ref_Ru",
  "xunik",
  "bstAUC",
  "cumdrop",
  "surv"
))

#' @useDynLib TwoSPC, .registration = TRUE
NULL

