

#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom data.table data.table

ROCsurvRec <- function(id, episode_j, time, group, mi,
                       checkPH = FALSE, compare=FALSE, area=TRUE, silent=FALSE, abtwc=FALSE,
                       xlab=NULL, ylab=NULL, main=NULL, MWCR, cex.axis = 1.5, cex.lab = 1.5,
                       legend.inset=0.02, legend.cex=1.5, lty = c(2,1,3), lwd = 1.5){

    dat <- data.frame(id, episode_j, time, group, mi)
    dat$delta <- ifelse(mi>1, 1, 0)
    dat <- dat[order(id,j),]

    MWCR <- getKMtabIntCens(dat)

    res<-onlyROC(MWCR[[1]], xlab, ylab, main, cex.axis = cex.axis,
            cex.lab = cex.lab, lty = lty, label.inset = label.inset,
            label.cex = label.cex, lwd = lwd)

}
