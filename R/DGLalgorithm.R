####### Based on the implementation in doubcens by James McVittie
####### made some small corrections to accomodate right censored onset times,
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

  #data set-up
  if (length(which(e_l - o_l<0))>0) {stop("at least one o_l > e_l")}
  if (length(which(e_r - o_l<0))>0) {stop("at least one o_l > e_r")}
  if (length(which(e_r - o_r<0))>0) {stop("at least one o_r > e_r")}
  if(length(which(e_l == -Inf))>0) {stop("at least one e_l = -Inf")}

  #rename variable to match paper
  X_L = o_l; X_R = o_r; Z_L = e_l; Z_R = e_r;
  N = length(X_L); x_val = t_val = NULL;

  for (i in 1:N) {
    x_val = c(x_val, X_L[i], X_R[i])
    t_val = c(t_val,
              Z_L[i] - X_L[i],
              ifelse(Z_L[i] - X_R[i]<0, -Inf, Z_L[i] - X_R[i]),
              Z_R[i] - X_L[i],
              ifelse(abs(Z_R[i])==Inf&abs(X_R[i])==Inf, Inf, Z_R[i] - X_R[i])
              )
  }

  x_val = sort(unique(x_val[x_val > 0]))
  t_val = sort(unique(t_val[t_val > -Inf]))

  x_val = x_val[-which(x_val==ctime)]
  tvind <- which(t_val>maxtime)
  ltvind <- length(tvind)
  if (ltvind>4){
    newltvind <- floor(length(tvind)/2)
    tvind <- tvind[newltvind:ltvind]
  } else {tvind <- max(tvind)}
  t_val = t_val[-tvind]
  alpha = matrix(0, nrow=length(x_val)*length(t_val), ncol=N)
  combos = expand.grid(1:length(x_val), 1:length(t_val))
  jth_index = combos[,1]; kth_index = combos[,2]
  for (i in 1:N) {
    for (j in 1:length(combos$Var1)) {
      if (X_L[i] <= x_val[combos[j,1]] & x_val[combos[j,1]] <= X_R[i] &
          Z_L[i] <= x_val[combos[j,1]] + t_val[combos[j,2]] &
          x_val[combos[j,1]] + t_val[combos[j,2]] <= Z_R[i])
      {
        alpha[j,i]=1
      }
    }
  }
  w_new = rep(1/length(x_val), length(x_val))
  f_new = rep(1/length(t_val), length(t_val))
  mu = matrix(0, nrow=length(alpha[,1]), ncol=length(alpha[1,]))
  for (i in 1:N) {
    w_group = w_new[jth_index]
    f_group = f_new[kth_index]
    for (j in 1:length(combos$Var1)) {
      mu[j,i] = (alpha[j,i]*w_new[combos[j,1]]*f_new[combos[j,2]]) /
        (sum(alpha[,i]*w_group*f_group))
    }
  }

  eps = 0.0005
  diff_onset = w_new
  diff_failure = f_new
  tol = sum(abs(diff_onset)) + sum(abs(diff_failure))
  counter = 0

  while (tol > eps) {
    if (counter == 1000) {stop("algorithm failed to converge in 1000 iterations")}
    w_old = w_new
    f_old = f_new
    for (j in 1:length(x_val)) {
      w_new[j] = sum(mu[combos$Var1==j,]) / N
    }
    for (k in 1:length(t_val)) {
      f_new[k] = sum(mu[combos$Var2==k,]) / N
    }
    diff_onset = w_old - w_new
    diff_failure = f_old - f_new
    tol = sum(abs(diff_onset)) + sum(abs(diff_failure))
    for (i in 1:N) {
      w_group = w_new[jth_index]
      f_group = f_new[kth_index]
      for (j in 1:length(combos$Var1)) {
        mu[j,i] = (alpha[j,i]*w_new[combos[j,1]]*f_new[combos[j,2]]) /
          (sum(alpha[,i]*w_group*f_group))
      }
    }
    counter = counter + 1
  }

  # if (counter < 1000) {
  #   print(paste("Algorithm converged in", counter, "iterations", sep = " "))
  # }

  onset <- cbind(t = x_val, f_t = w_new)
  gap <- cbind(t = t_val, f_t = f_new, F_t =  cumsum(f_new), S_t = 1-cumsum(f_new))

  onset <- onset[which(onset[,1]<=maxtime),]
  gap <- gap[which(gap[,1]<=maxtime),]

  return(list(onset = onset, gap=gap, iterations = counter))

}
