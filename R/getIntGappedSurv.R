#' INTERNAL: data management not for user
#'
#' @description
#' This function gets survival function for each group for interval censored data.
#'
#' @param npmle_df pass from intROCsurv.
#' @param stepsize pass from intROCsurv.
#'
#'@importFrom data.table between
#'@keywords internal
#'@noRd
#'
getIGsurv <- function(npmle_df, stepsize){

  jump = min(stepsize)/2
  t <- seq(0, npmle_df$L[1], jump)
  surv <- rep(1, length(t))
  skm <- cbind(t=t, surv=surv)

  for (j in 2:nrow(npmle_df)) {
    if (npmle_df$R[j-1] != npmle_df$L[j]){
      t_temp <- seq(npmle_df$R[j-1], npmle_df$L[j], jump)
      s_temp <- rep(1-npmle_df$cumdrop[j-1], length(t_temp))
      skm <- rbind(skm, cbind(t_temp, s_temp))
    }
  }

  return(skm)
}

