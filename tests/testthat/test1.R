#This will be used for a test

# n = 500
# maxt = 2
# beta = 0.2
#
# #1) Assuming data follows Wiebull(0.5, 1.5).
#     set.seed(28)
#     simdata <- simulate(size = n, followup = maxt, beta = beta, dist = "weibull")
#     # Check the maximum survival time in the data the censoring rate.
#      maxt <- max(simdata$time)
#      maxt
#      cens <- table(simdata$event)[1]/500 #estimate proportion censored
#      cens
