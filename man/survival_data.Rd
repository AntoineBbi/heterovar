% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_data.R
\name{survival_data}
\alias{survival_data}
\title{survival_data is a R function used to generate survival datasets from exponential model
with standard deviation and current value from longitudinal_data function included as covariates}
\usage{
survival_data(
  data,
  censoring_param = c(1, 1),
  tmax = 1,
  lambda = 1,
  covariate = NULL,
  alpha_cov = NULL,
  alpha_sigma = 0.03,
  alpha_mu = 0
)
}
\arguments{
\item{data}{dataframe (one raw by id) return by longitudinal_data function}

\item{censoring_param}{2-length vector including parameters of censoring beta distribution (by default, uniform distribution between 0 and tmax given 'censoring_param = c(1,1)')}

\item{tmax}{maximum of the follow-up (by default 1)}

\item{lambda}{parameter of baseline risk function as h_0(t)=exp(lambda)}

\item{covariate}{matrix (nxp) with n is nraw(data) and p the number of time-independent covariates}

\item{alpha_cov}{p-length vector of parameters associated with p covariates}

\item{alpha_sigma}{parameter associated with the subject-specific standard deviation}

\item{alpha_mu}{parameter associated with the subject-specific current value}
}
\value{
a dataframe with generated data
}
\description{
survival_data is a R function used to generate survival datasets from exponential model
with standard deviation and current value from longitudinal_data function included as covariates
}
\examples{

n <- 200
data_constant <- longitudinal_data(n = n,
                                   tmin = 0,
                                   tmax = 5,
                                   pas = 0.5,
                                   sd_tmin = 10*2/3,
                                   sd_tmax = 10*2/3,
                                   intercept = 0,
                                   slope = 1,
                                   sigma2_b = matrix(c(sigma2_b_00, sigma2_b_01, sigma2_b_01, sigma2_b_11), ncol = 2),
                                   group = "Constant variability")

data_decrease <- longitudinal_data(n = n,
                                   tmin = 0,
                                   tmax = 5,
                                   pas = 0.5,
                                   sd_tmin = 10,
                                   sd_tmax = 10/3,
                                   intercept = 0,
                                   slope = 1,
                                   sigma2_b = matrix(c(sigma2_b_00, sigma2_b_01, sigma2_b_01, sigma2_b_11), ncol = 2),
                                   group = "Decreasing variability")

data_increase <- longitudinal_data(n = n,
                                   tmin = 0,
                                   tmax = 5,
                                   pas = 0.5,
                                   sd_tmin = 10/3,
                                   sd_tmax = 10,
                                   intercept = 0,
                                   slope = 1,
                                   sigma2_b = matrix(c(sigma2_b_00, sigma2_b_01, sigma2_b_01, sigma2_b_11), ncol = 2),
                                   group = "Increasing variability")

# wide data
data_trans <- data_increase$data_trans
data_trans$id <- data_trans$id + n
data_trans <- rbind(data_trans, data_constant$data_trans)
data_trans$id <- data_trans$id + n
data_trans <- rbind(data_trans, data_decrease$data_trans)
data_trans <- data_trans[order(data_trans$id), ]
data_trans$id
table(data_trans$variability)


# longitudinal data
data_long <- data_increase$data_long
data_long$id <- data_long$id + n
data_long <- rbind(data_long, data_constant$data_long)
data_long$id <- data_long$id + n
data_long <- rbind(data_long, data_decrease$data_long)
data_long <- data_long[order(data_long$id), ]
unique(data_long$id)
table(data_long$variability)


# Generation of Survival data
surv.data <- survival_data(data_trans,
                           censoring_param = c(5,2),
                           tmax = 5,
                           lambda = -6,
                           covariate = matrix(rbinom(nrow(data_trans), 1, prob = 0.4), ncol = 1),
                           alpha_cov = c(0.3),
                           alpha_sigma = 0.5,
                           alpha_mu = 0.1)

}
