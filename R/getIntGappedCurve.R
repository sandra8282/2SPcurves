
#' INTERNAL: data management not for user
#'
#' @description
#' This function gets the roc for each group for interval censored data.
#'
#' @param npmle_0 passed from intROCsurv.
#' @param npmle_1 passed from intROCsurv.
#' @param xlab passed from intROCsurv.
#' @param ylab passed from intROCsurv.
#' @param main passed from intROCsurv.
#' @param cex.axis passed from intROCsurv.
#' @param cex.lab passed from intROCsurv.
#' @param lwd passed from intROCsurv.
#'
#'@importFrom data.table setDT
#'@importFrom data.table setkey
#'@importFrom data.table foverlaps
#'@importFrom DescTools AUC
#'
#'@keywords internal
#'@noRd
#'
getIGroc <- function(npmle_0, npmle_1, xlab, ylab, main, cex.axis,
                     cex.lab, lwd){

  dt0.L = dt1.R = dt0.R = dt1.L = L = NULL

  #search for overlapping intervals
  setDT(npmle_0); setDT(npmle_1);
  n0 = nrow(npmle_0); n1 = nrow(npmle_1);
  setkey(npmle_1)
  overlaps_ind <- foverlaps(npmle_0[,1:2], npmle_1[,1:2], type="any",
                            nomatch=0L, which = TRUE)
  combined_dt <- cbind(dt0 = npmle_0[overlaps_ind$xid,],
                       dt1 = npmle_1[overlaps_ind$yid,])
  notoverlap <- c(which(combined_dt$dt0.L == combined_dt$dt1.R),
                    which(combined_dt$dt0.R == combined_dt$dt1.L))
  temp <- combined_dt[notoverlap,]
  combined_dt <- combined_dt[dt0.L != dt1.R]
  combined_dt <- combined_dt[dt0.R != dt1.L]
  combined_dt$move <- "overlap"

  ind0 <- which(temp$dt0.L %in% combined_dt$dt0.L)
  ind1 <- which(temp$dt1.L %in% combined_dt$dt1.L)

  temp0o <- temp[-ind0, 1:4]
  temp1o <- temp[-ind1, 5:8]

  #find non-overlapping intervals for controls
  temp0 <- npmle_0[!(L %in% combined_dt$dt0.L),]
  if (nrow(temp0o)>0){temp0 <- rbind(temp0, unname(temp0o))}
  temp0 <- cbind(temp0, matrix(ncol = ncol(temp0), nrow = nrow(temp0)),
                 rep("control", nrow(temp0)))

  #find non-overlapping intervals for treated
  temp1 <- npmle_1[!(L %in% combined_dt$dt1.L),]
  if (nrow(temp1o)>0){temp1 <- rbind(temp1, unname(temp1o))}
  temp1 <- cbind(matrix(ncol = ncol(temp1), nrow = nrow(temp1)),
                 temp1, rep("trt", nrow(temp1)))

  #combine overlapping and non-overlapping intervals (with respective key)
  colnames(temp0) = colnames(temp1) = colnames(combined_dt)
  combined_dt <- rbind(combined_dt, temp0, temp1)
  combined_dt <- combined_dt[!duplicated(combined_dt),]
  combined_dt$dt0.L <- ifelse(is.na(combined_dt$dt0.L),
                              combined_dt$dt1.L, combined_dt$dt0.L)
  combined_dt$dt1.L <- ifelse(is.na(combined_dt$dt1.L),
                              combined_dt$dt0.L, combined_dt$dt1.L)
  combined_dt <- combined_dt[order(dt1.L, dt0.L),]

  #calculate survival and keep only needed columns
  combined_dt$S_0 <- 1 - combined_dt$dt0.cumdrop
  combined_dt$S_1 <- 1 - combined_dt$dt1.cumdrop
  combined_dt = combined_dt[,c(1:2, 5:6, 10:11, 9)]

  #find intervals that have uncertainty at end point and skip plotting those
  skip_ind <- c(which(duplicated(combined_dt[,1:2])),
                which(duplicated(combined_dt[,3:4])))
  combined_dt$skip <- rep("no", nrow(combined_dt))

  if (length(skip_ind)>0) {combined_dt$skip[skip_ind-1] = "yes"}

  #plot
  plot(NULL, type="n", las=1,
       xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
       xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)

  forplot = coord_new = c(1,1); type = "";
  for (i in 1:nrow(combined_dt)){
    if (combined_dt$skip[i]=="no"){
      if (combined_dt$move[i] == "overlap"){
        coord_new2 = c(combined_dt$S_0[i], combined_dt$S_1[i])
        rect(xright = coord_new[1], ytop = coord_new[2],
             xleft = coord_new2[1], ybottom = coord_new2[2],
             col = "grey", border = "grey")
        type = c(type, "overlap")
        colt = "grey"
      } else {
        colt = "black"
        if (combined_dt$move[i] == "trt"){
          coord_new2 = c(coord_new[1], combined_dt$S_1[i])
          type = c(type, "trt")
        } else {
          coord_new2 = c(combined_dt$S_0[i], coord_new[2])
          type = c(type, "control")
        }
      }
      segments(x0=coord_new[1], y0=coord_new[2],
               x1=coord_new2[1], y1=coord_new2[2], lwd = lwd, lty = 1, col = colt)
      points(rbind(coord_new, coord_new2), pch=20, cex=0.75)
      forplot = rbind(forplot, coord_new2)
      coord_new = coord_new2
    }
  }
  abline(c(0,1), col = "grey", lty=1, lwd = lwd-0.25)
  auc = AUC(forplot[,1], forplot[,2])
  rownames(forplot) = NULL

  colnames(combined_dt) = c("L0", "R0", "L1", "R1",
                            "S0", "S1", "Move", "Skip")
  combined_dt$L0 <- ifelse(is.na(combined_dt$R0), NA, combined_dt$L0)
  combined_dt$L1 <- ifelse(is.na(combined_dt$R1), NA, combined_dt$L1)
  combined_dt$order <- 1:nrow(combined_dt)

  forplot = data.frame(forplot, type)
  colnames(forplot) = c("u", "R_u", "type")
  return(list(curve = forplot, auc = auc))
}
