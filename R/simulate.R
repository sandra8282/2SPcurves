#' Simulate a clinical trial with survival outcomes for treatment and control groups based on proportional hazards.
#'
#' @description
#' This function randomly assigns individuals to treatment or control groups then uses
#' the coxed package to simulate event and censoring times based on a proportional hazards
#' model without any other parametric asumptions based on the flexible-hazard method.
#' See Harden and Kropko (2018) for further details.
#'
#' @param time The simulated study's follow-up time (maximum time a subject could survive without censoring).
#' @param size The total number of subjects or sample size.
#' @param hazardratio The treatment to control hazard ratio.
#' @param censprop Proportion of observations to designate as right-censored.
#'
#' @return A data frame containing:
#' \itemize{
#'  \item Treatment assignment 0-control, 1-treatment
#'  \item Time to event (time)
#'  \item Event indicator (event: 1-event occurred, 0-no event/censored).
#' }
#'
#' @import survival
#' @importFrom coxed sim.survdata
#' @importFrom stats na.omit
#' @importFrom stats rbinom
#'
#' @examples
#'
#' #Simulate a study of 500 subjects followed for 365 days where the hazard ratio of 0.2 shows treatment prolongs life.
#' set.seed(28)
#' simdata <- simulate(500, 365, 0.2, 0.15)
#'
#' #Check proportion censored and hazard ratio estimates and compare to true.
#' cens <- unname(table(simdata$event)[1]/500) #estimate proportion censored
#' treatment_hr <- unname(exp(coxph(Surv(time, event) ~ treatment, data=simdata, ties="breslow")$coefficients))
#' cens
#' treatment_hr
#'
#' @export

simulate <-function(size, time, hazardratio, censprop) {
  treatment = rbinom(size, 1, 0.5)
  oldw <- getOption("warn")
  options(warn = -1)
  res <- sim.survdata(T=time, X=data.frame(treatment), beta = log(hazardratio),
                      num.data.frames = 1, censor = censprop)
  simdata <- data.frame(id = seq(1, size, 1))
  simdata <- cbind(simdata, res$data)
  simdata[,4] <- as.numeric(simdata[,4])
  colnames(simdata) <- c("id","treatment", "time", "event")
  options(warn = oldw)
  return(simdata)
}
