#' Estimate SE and CI for AUC based on 1000 sample bootstrap
#'
#' @param time passed from ROCsurv
#' @param event passed from ROCsurv
#' @param group passed from ROCsurv
#' @param method passed from ROCsurv
#' @param level passed from ROCsurv
#' @param B number of bootstrap resamples
#'
#' @return SE and CI for indicated level
#'
#' @importFrom stats sd
#' @keywords internal
#' @noRd

btsp <- function(time, event, group, method, B, level) {

    id <- 1:length(time)
    maindat <- cbind(id, time, event, group)
    BTest <- rep(0, B)

    for (b in 1:B) {
      set.seed(b+7679)
      idsamps <- sample(id, length(id), replace = TRUE)
      btdat <- data.frame(maindat[idsamps,])
      KMests <- with(btdat, getKMtab(time, event, group))
      if (KMests[[2]]==0) {BTest[b] <- completeROC(KMests[[1]], silent = TRUE)}
      if (KMests[[2]]!=0 & method=="restrict") {BTest[b] <- restrictROC(KMests[[1]], silent = TRUE)}
    }

    bootstrapCIs <- c(sd(BTest), sort(BTest)[(1-level)*B], sort(BTest)[level*B])
    lowstring <- paste(level*100, "%", " Bootstrap Lower",sep="")
    upstring <- paste(level*100, "%", " Bootstrap Upper",sep="")
    names(bootstrapCIs) <- c("Bootstrap SE", lowstring, upstring)
    return(bootstrapCIs)
}
