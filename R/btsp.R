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
    rownames(maindat) <- id
    g1dat <- maindat[maindat$group==1,]
    g2dat <- maindat[maindat$group==0,]

    BTest <- rep(0, B)

    for (b in 1:B) {
      set.seed(b+7679)
      g1idsamps <- sample(as.integer(rownames(g1dat)), nrow(g1dat), replace = TRUE)
      g21idsamps <- sample(as.integer(rownames(g2dat)), nrow(g2dat), replace = TRUE)
      btdat <- data.frame(rbind(g1dat[g1idsamps, ], g2dat[g2idsamps, ]))
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
