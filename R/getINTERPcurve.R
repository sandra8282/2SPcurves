#' Create an ROC curve for survival data from clinical trials.
#'
#' @description
#' This function interpolates missing portions of ROC curve
#'
#' @param control_pf passed from ROCsurv
#' @param trt_pf passed from ROCsurv
#' @param res passed from ROCsurv
#' @return Interpolated portions of ROC
#' \itemize{
#'  \item Data frame with the estimated curve.
#' }
#'
#' @importFrom Icens EMICM
#' @importFrom data.table data.table
#' @importFrom dplyr %>% filter select
#' @importFrom stats complete.cases
#' @importFrom DescTools AUC
#' @importFrom stats approx
#'
#' @keywords internal
#'
#' @noRd

getINTERPcurve <- function(control_pf, trt_pf, res, interp_n=1000) {

    ### survival estimations for missing portions through interpolation
    control_surv <- control_pf %>% mutate(surv = 1-cumdrop, surv0 = lag(surv, default = 1),
                                          group = 0, interpolated = "No")
    control_surv$interpolated <- ifelse(control_surv$R==Inf, "Yes", control_surv$interpolated)
    control_surv$R <- ifelse(control_surv$R==Inf, control_surv$R[nrow(control_surv)-1]*5, control_surv$R)
    control_surv_new <- data.frame(t=0, s=1)
    for (i in 1:nrow(control_surv)) {
      if (control_surv[i,1]!=control_surv[i, 2]){
        ti = approx(control_surv[i,1:2], c(control_surv$surv0[i],control_surv$surv[i]), n = interp_n)
        sfct = data.frame(t = ti$x, s = ti$y)
        control_surv_new = rbind(control_surv_new, sfct)
      }
    }

    trt_surv <- trt_pf %>% mutate(surv = 1- cumdrop, surv0 = lag(surv, default = 1),
                                  group = 0, interpolated = "No")
    trt_surv$interpolated <- ifelse(trt_surv$R==Inf, "Yes", trt_surv$interpolated)
    trt_surv$R <- ifelse(trt_surv$R==Inf, trt_surv$R[nrow(trt_surv)-1]+trt_surv$R[nrow(trt_surv)-1]*2, trt_surv$R)
    trt_surv_new <- data.frame(t=0, s=1)
    for (i in 1:nrow(trt_surv)) {
      if (trt_surv[i,1]!=trt_surv[i, 2]){
        ti = approx(trt_surv[i,1:2], c(trt_surv$surv0[i],trt_surv$surv[i]), n = interp_n)
        sfct = data.frame(t = ti$x, s = ti$y)
        trt_surv_new = rbind(trt_surv_new, sfct)
      }
    }


    # add missing portions to two sample curve ######################################################################################
    control_surv_new$group = 0
    trt_surv_new$group = 1
    surv_new = distinct(rbind(control_surv_new, trt_surv_new))
    surv_new = surv_new[order(surv_new$t, -surv_new$s),]
    KMests <- data.frame(getIntSKM(surv_new))
    res_temp <- data.frame(get4plot(skm = KMests))
    if (res_temp$x[nrow(res_temp)]==0|res_temp$y[nrow(res_temp)]==0){
      res_temp[nrow(res_temp)+1,]=c(0, 0, 1)
    }
    res_temp = distinct(res_temp)
    res_temp <- dplyr::left_join(res_temp, res, by = c("x"="u", "y"="R_u"))
    #res_temp <- res_temp %>% filter(!(nrow(res_temp)!=1&res_temp$y==1&is.na(res_temp$type)))
    if (length(table(duplicated(res_temp[,(1:2)])))>1){
        res_temp <- res_temp[-duplicated(res_temp[,(1:2)]),]
    }

    return(list(res_temp=res_temp, control_surv_new=control_surv_new,
                trt_surv_new=trt_surv_new))
}
