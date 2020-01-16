#' Simulate a clinical trial with survival outcomes for treatment and control groups based on various proportional hazards models.
#'
#' @description
#' This function uses the simsurv and coxed packages to simulate event and censoring times for a proportional hazards
#' model with either: parametric distribution assumptions, or a user defined hazard function, or in using a flexible-hazard
#' method that needs no distributions to be specified. The user can provide treatment assignments for the subjects or have
#' them randomly assigned to each subject.
#'
#' @param size Numeric value specifying the number of subjects or sample size.
#' @param followup Numeric value specifying the simulated study's follow-up time in years.
#' @param beta Numeric value fot the treatment effect.
#' @param X Optional data frame with two columns for subject id's and their assigned treatment.
#' The column for treatment assignment must be a logical or numeric vector where values of TRUE or 1 will indicate individual was assigned to the treatment
#' arm and FALSE or 0 to control arms. Any other values are invalid.
#' If $X$ is not provided the treatment will be assigned randomly with probability 0.5 and subject id's will be a sequence from 1 to the specified sample size.
#' @param dist Optional character string specifying the parametric survival distribution to be used. Options include "weibull" (the default), "exponential", or "gompertz"
#'  with parameters \code{alpha = 0.5, gamma = 1.5}. This is ignored if the \code{hazard} argument is provided. See \strong{Details} and \strong{Examples}.
#' @param hazard Optional user-defined hazard function, with arguments \code{t}, \code{X}, and \code{beta}. This function should return
#' the hazard at time \code{t} for an individual with treatment assigment supplied via \code{X} and treatment effect supplied via \code{beta}.
#' See \strong{Details} and \strong{Examples}.
#' @param censprop Optional the proportion of observations to designate as right-censored when the flexible-hazards model is used. This is ignored if the \code{hazard}
#' or \code{dist} arguments are provided. See \strong{Details} and \strong{Examples}.
#'
#' @return A data frame containing:
#' \itemize{
#'  \item subject id's
#'  \item Treatment assignment 0-control, 1-treatment
#'  \item Time to event (time)
#'  \item Event indicator (event: 1-event occurred, 0-no event/censored).
#' }
#'
#' @import survival
#' @import simsurv
#' @importFrom coxed sim.survdata
#' @importFrom stats na.omit
#' @importFrom stats rbinom
#'
#' @examples
#'
#' #Simulate a clinical trial with 500 subjects that were followed for 2 years.
#' #where the true hazard ratio of 0.2 indicates a large treatment effect.
#'
#' n = 500
#' maxt = 2
#' beta = 0.2
#'
#' #1) Assuming data follows Wiebull(0.5, 1.5).
#'     set.seed(28)
#'     simdata <- simulate(size = n, followup = maxt, beta = beta, dist = "weibull")
#'     # Check the maximum survival time in the data the censoring rate.
#'      maxt <- max(simdata$time)
#'      maxt
#'      cens <- table(simdata$event)[1]/500 #estimate proportion censored
#'      cens
#'
#' #2) Assuming no distribution and using flexible hazards model with 25% censoring rate.
#'     set.seed(28)
#'     simdata2 <- simulate(size = n, followup = maxt, beta = beta, censprop = 0.25)
#'     # Check the maximum survival time in the data the censoring rate.
#'      maxt <- max(simdata2$time)
#'      maxt
#'      cens <- table(simdata2$event)[1]/500 #estimate proportion censored
#'      cens
#'
#' @export

simSurvTrial <- function(size, followup, beta, X=NULL, dist=NULL, hazard=NULL, censprop=NULL) {

  if (is.null(X)) {X <- data.frame(id = 1:size, trt = stats::rbinom(size, 1, 0.5))}

  check1 <- c(dist, hazard, censprop)
  if (is.null(check1)) {stop("Error: Must provide one of the following arguments - dist, hazard or censprop.")}

  if (!is.null(dist)) {
    s1 <- simsurv(dist = dist, lambdas = 0.5, gammas = 1.5, betas = c(trt = beta), x = X, maxt = followup)
    simdata <- data.frame(cbind(id = s1$id, treatment = X$trt, time = s1$eventtime, event = s1$status))
    return(simdata)
  }

  if (!is.null(censprop)) {
    oldw <- getOption("warn")
    options(warn = -1)
    res <- sim.survdata(N=size, T=followup*365, X=data.frame(treatment = X[,2]), beta = log(beta),
                        num.data.frames = 1, censor = censprop)
    simdata <- data.frame(id = seq(1, size, 1))
    simdata <- cbind(simdata, res$data)
    simdata[,4] <- as.numeric(simdata[,4])
    colnames(simdata) <- c("id","treatment", "time", "event")
    options(warn = oldw)
    return(simdata)
  }

}
