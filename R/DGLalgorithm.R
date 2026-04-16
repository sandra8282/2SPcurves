####### Based on the implementation in doubcens by James McVittie
####### made some small corrections to accommodate right censored onset times,
####### inputs for censoring to be something besides Inf and improve speed.
####### Added option to get cdf and survival of T.
####### Sandra Castro-Pearson - 8/8/21

#' De Grottule and Legakos algorithm
#'
#' @param o_l passed from other functions.
#' @param o_r passed from other functions.
#' @param e_l passed from other functions.
#' @param e_r passed from other functions.
#' @param maxtime passed from other functions.
#' @param ctime passed from other functions.
#'
#' @return onset and gap estimates
#'
#' @keywords internal
#' @noRd
#'
#'

DGLalg = function(o_l, o_r, e_l, e_r, maxtime, ctime) {

  if (any(e_l - o_l < 0)) stop("at least one o_l > e_l")
  if (any(e_r - o_l < 0)) stop("at least one o_l > e_r")
  if (any(e_r - o_r < 0)) stop("at least one o_r > e_r")
  if (any(e_l == -Inf))   stop("at least one e_l = -Inf")

  X_L = o_l; X_R = o_r; Z_L = e_l; Z_R = e_r
  N   = length(X_L)

  # build x_val and t_val
  x_val = NULL; t_val = NULL
  for (i in 1:N) {
    x_val = c(x_val, X_L[i], X_R[i])
    t_val = c(t_val,
              Z_L[i] - X_L[i],
              ifelse(Z_L[i] - X_R[i] < 0, -Inf, Z_L[i] - X_R[i]),
              Z_R[i] - X_L[i],
              ifelse(abs(Z_R[i]) == Inf & abs(X_R[i]) == Inf,
                     Inf, Z_R[i] - X_R[i]))
  }
  x_val = sort(unique(x_val[x_val > 0]))
  t_val = sort(unique(t_val[t_val > -Inf]))
  x_val = x_val[-which(x_val == ctime)]

  tvind  = which(t_val > maxtime)
  ltvind = length(tvind)
  if (ltvind > 4) {
    newltvind = floor(length(tvind) / 2)
    tvind = tvind[newltvind:ltvind]
  } else {
    tvind = max(tvind)
  }
  t_val = t_val[-tvind]

  nx = as.integer(length(x_val))
  nt = as.integer(length(t_val))

  res = .Fortran("dglalg_f",
                 x_l     = as.double(X_L),
                 x_r     = as.double(X_R),
                 z_l     = as.double(Z_L),
                 z_r     = as.double(Z_R),
                 n       = as.integer(N),
                 maxtime = as.double(maxtime),
                 ctime   = as.double(ctime),
                 x_val   = as.double(x_val),
                 nx      = nx,
                 t_val   = as.double(t_val),
                 nt      = nt,
                 w_out   = double(nx),
                 f_out   = double(nt),
                 counter = integer(1),
                 info    = integer(1))

  if (res$info ==  1) stop("algorithm failed to converge in 1000 iterations")
  if (res$info == -1) stop("zero denominator in EM update")

  w_new = res$w_out
  f_new = res$f_out

  onset = cbind(t = x_val, f_t = w_new)
  gap   = cbind(t = t_val, f_t = f_new,
                F_t = cumsum(f_new),
                S_t = 1 - cumsum(f_new))
  onset = onset[onset[, 1] <= maxtime, ]
  gap   = gap[gap[, 1] <= maxtime, ]

  return(list(onset = onset, gap = gap, iterations = res$counter))
}
