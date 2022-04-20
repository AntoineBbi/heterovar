#' longitudinal_data is a R function to generate data from the following linear mixed model
#' Y_ij|b_i ~ N(intercept + b_0i + (slope + b_1i) x t_ij , sigma_e_0 + sigma_e_1 x t_ij)
#' with b_i ~ N(0, Sigma2_b)
#'
#'
#' @param n number of subjects
#' @param tmin numerical positive value indicating the minimum of time windows (baseline, by default 0)
#' @param tmax numerical positive value indicating the maximum of the follow-up (by default 1)
#' @param sd_tmin positive value of error standard deviation at tmin
#' @param sd_tmax positive value of error standard deviation at tmax
#' @param pas time between two successive measurements
#' @param intercept value of intercept fixed effect at tmin for the linear trajectory
#' @param slope value of slope fixed effect for the linear trajectory
#' @param sigma2_b covariance matrix of random effects b_i
#' @param group name of the subject group (type: characters)
#'
#' @return A dataframe with generated data
#'
#' @export
#'
#' @examples
#'
#' n <- 200
#' data_constant <- longitudinal_data(n = n,
#'                                    tmin = 0,
#'                                    tmax = 5,
#'                                    pas = 0.5,
#'                                    sd_tmin = 10*2/3,
#'                                    sd_tmax = 10*2/3,
#'                                    intercept = 0,
#'                                    slope = 1,
#'                                    sigma2_b = matrix(c(sigma2_b_00, sigma2_b_01, sigma2_b_01, sigma2_b_11), ncol = 2),
#'                                    group = "Constant variability")
#'
#' data_decrease <- longitudinal_data(n = n,
#'                                    tmin = 0,
#'                                    tmax = 5,
#'                                    pas = 0.5,
#'                                    sd_tmin = 10,
#'                                    sd_tmax = 10/3,
#'                                    intercept = 0,
#'                                    slope = 1,
#'                                    sigma2_b = matrix(c(sigma2_b_00, sigma2_b_01, sigma2_b_01, sigma2_b_11), ncol = 2),
#'                                    group = "Decreasing variability")
#'
#' data_increase <- longitudinal_data(n = n,
#'                                    tmin = 0,
#'                                    tmax = 5,
#'                                    pas = 0.5,
#'                                    sd_tmin = 10/3,
#'                                    sd_tmax = 10,
#'                                    intercept = 0,
#'                                    slope = 1,
#'                                    sigma2_b = matrix(c(sigma2_b_00, sigma2_b_01, sigma2_b_01, sigma2_b_11), ncol = 2),
#'                                    group = "Increasing variability")
#'
#' # wide data
#' data_trans <- data_increase$data_trans
#' data_trans$id <- data_trans$id + n
#' data_trans <- rbind(data_trans, data_constant$data_trans)
#' data_trans$id <- data_trans$id + n
#' data_trans <- rbind(data_trans, data_decrease$data_trans)
#' data_trans <- data_trans[order(data_trans$id), ]
#' data_trans$id
#' table(data_trans$variability)
#'
#'
#' # longitudinal data
#' data_long <- data_increase$data_long
#' data_long$id <- data_long$id + n
#' data_long <- rbind(data_long, data_constant$data_long)
#' data_long$id <- data_long$id + n
#' data_long <- rbind(data_long, data_decrease$data_long)
#' data_long <- data_long[order(data_long$id), ]
#' unique(data_long$id)
#' table(data_long$variability)
#'
#'
#' # Generation of Survival data
#' surv.data <- survival_data(data_trans,
#'                            censoring_param = c(5,2),
#'                            tmax = 5,
#'                            lambda = -6,
#'                            covariate = matrix(rbinom(nrow(data_trans), 1, prob = 0.4), ncol = 1),
#'                            alpha_cov = c(0.3),
#'                            alpha_sigma = 0.5,
#'                            alpha_mu = 0.1)

longitudinal_data <- function(n,
                              tmin = 0,
                              tmax = 1,
                              sd_tmin,
                              sd_tmax,
                              pas,
                              intercept,
                              slope,
                              sigma2_b,
                              group){

  # Simulated model with variance time varying
  # Y_ij|b_i ~ N(intercept + b_0i + (slope + b_1i) * t_ij , sigma_e_0 + sigma_e_1 * t_ij)
  # b_i ~ N(0, Sigma2_b)
  ## n          : number of subjects
  ## tmin       : minimum of time windows (baseline, by default 0)
  ## tmax       : maximum of the follow-up (by default 1)
  ## sd_tmin    : value of error standard deviation at tmin
  ## sd_tmax    : value of error standard deviation at tmax
  ## pas        : time between two successive measurements
  ## intercept  : value of intercept fixed effect at tmin for the linear trajectory
  ## slope      : value of slope fixed effect for the linear trajectory
  ## sigma2_b   : covariance matrix of random effects b_i
  ## group      : name of the subject group (type: characters)

  ### define the visit vector
  visit <- seq(tmin, tmax, by=pas)

  ### variance of errors with linear tends: sigma_e(t) = sigma_e_0 + sigma_e_1 * t
  sigma_e_0 <- sd_tmin
  sigma_e_1 <- (sd_tmax - sigma_e_0)/tmax
  sigma2_b_00 <- 210
  sigma2_b_01 <- -18
  sigma2_b_11 <- 10
  sigma2_b <- matrix(c(sigma2_b_00, sigma2_b_01, sigma2_b_01, sigma2_b_11), ncol = 2)

  ## simulation of individual trajectories
  ### transversal data
  data_trans <- data.frame(cbind(id = seq(1,n),
                                 intercept = intercept,
                                 slope = slope,
                                 sigma_e_0 = sigma_e_0,
                                 sigma_e_1 = sigma_e_1,
                                 sigma2_b_00 = sigma2_b_00,
                                 sigma2_b_01 = sigma2_b_01,
                                 sigma2_b_11 = sigma2_b_11))
  ### generate random effects
  mat_b <- mvtnorm::rmvnorm(n, c(0,0), sigma2_b)
  colnames(mat_b) <- c("b0", "b1")
  data_trans <- cbind(data_trans, mat_b)
  ### group of variability
  data_trans$variability <- group

  ### longitudinal data
  data_long <- as.data.frame(cbind(rep(1:n, length(visit)), # ID patient
                                   sort(rep(visit, n)), # value of visit
                                   matrix(rep(t(mat_b), length(visit)), ncol=ncol(mat_b), byrow=TRUE) # specific random effect
  )
  )

  colnames(data_long) <- c("id","visit", paste("b", 0:(ncol(mat_b)-1), sep = ""))
  data_long$sd_error <- sigma_e_0 + sigma_e_1 * data_long$visit
  data_long <- data_long[order(data_long$id, data_long$visit), ]
  data_long$error <- rnorm(nrow(data_long), 0, sd = data_long$sd_error)
  ## latent process
  data_long$y <- intercept + data_long$b0 + data_long$visit * ( slope + data_long$b1 ) + data_long$error

  data_long$variability <- group

  # output of function
  out <- list(data_trans = data_trans,
              data_long = data_long)
  return(out)

}
