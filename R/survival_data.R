#' survival_data is a R function used to generate survival datasets from exponential model
#' with standard deviation and current value from longitudinal_data function included as covariates
#'
#' @param data dataframe (one raw by id) return by longitudinal_data function
#' @param censoring_param 2-length vector including parameters of censoring beta distribution (by default, uniform distribution between 0 and tmax given 'censoring_param = c(1,1)')
#' @param tmax maximum of the follow-up (by default 1)
#' @param lambda parameter of baseline risk function as h_0(t)=exp(lambda)
#' @param covariate matrix (nxp) with n is nraw(data) and p the number of time-independent covariates
#' @param alpha_cov p-length vector of parameters associated with p covariates
#' @param alpha_sigma parameter associated with the subject-specific standard deviation
#' @param alpha_mu parameter associated with the subject-specific current value
#'
#' @return a dataframe with generated data
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
#'
survival_data <- function(data,
                          censoring_param = c(1,1),
                          tmax = 1,
                          lambda = 1,
                          covariate = NULL,
                          alpha_cov = NULL,
                          alpha_sigma = 0.03,
                          alpha_mu = 0){

  # Simulated model with variance time varying
  ## data          : data.frame (one raw by id) return by longitudinal_data function
  ## censure_param : 2-length vector including parameters of censoring beta distribution (by default, uniform distribution between 0 and tmax given 'censoring_param = c(1,1)')
  ## tmax          : maximum of the follow-up (by default 1)
  ## lambda        : parameter of baseline risk function as h_0(t)=exp(lambda)
  ## covariate     : matrix (nxp) with n is nraw(data) and p the number of time-independent covariates
  ## alpha_cov     : p-length vector of parameters associated with p covariates
  ## alpha_sigma   : parameter associated with the subject-specific standard deviation
  ## alpha_mu      : parameter associated with the subject-specific current value


  out <- data

  if(!("matrix" %in% class(covariate)))
    stop("'covariate' is not a matrix object.\n")

  if(nrow(covariate)!=nrow(data))
    stop("'data' and 'covariate' have different number of raws.\n")

  if(ncol(covariate)!=length(alpha_cov))
    stop("number of colums of 'covariate' and length of alpha_cov are different.\n")

  if(is.null(covariate) || is.null(alpha_cov)){
    pred_baseline <- 0
  }else{
    pred_baseline <- as.numeric(covariate%*%alpha_cov)
    covariate <- as.data.frame(covariate)
    colnames(covariate) <- paste("cov",1:ncol(covariate), sep = "")
  }



  # initialization of constants
  c1 <- alpha_sigma*data$sigma_e_0 + alpha_mu*(data$intercept+data$b0) + lambda + pred_baseline
  c2 <- alpha_sigma*data$sigma_e_1 + alpha_mu*(data$slope+data$b1)

  # Generation of time-to-event
  u <- runif(nrow(data))
  out$time_event <- log(-log(u) * c2 * exp(-c1) + 1) / c2

  # Generation of censoring time
  out$time_censoring <- tmax*rbeta(n = nrow(data),
                                   shape1 = censoring_param[1],
                                   shape2 = censoring_param[2])

  # Observed time
  out$time <- apply(out[, c("time_event", "time_censoring")], 1, min)
  out$time[which(is.na(out$time_event))] <- out$time_censoring[which(is.na(out$time_event))]
  out$event <- as.numeric(out$time_event<=out$time_censoring)
  out$event[which(is.na(out$event))] <- 0

  # Add covariate in data
  out <- cbind(out, covariate)

  # output of function
  out <- list(data = out,
              control = list(draws.id = data.frame(cbind(out$id, c1 = c1, c2 = c2, u = u)),
                             censoring_param = censoring_param,
                             tmax = tmax,
                             lambda = lambda,
                             covariate = covariate,
                             alpha_cov = alpha_cov,
                             alpha_sigma = alpha_sigma,
                             alpha_mu = alpha_mu)
  )
  return(out)

}
