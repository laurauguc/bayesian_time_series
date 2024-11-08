# bayesian
bayesian_arima <- function(x, order, h) { # seperate out obtaining predictions (no "h")
  d = order[2]
  if (d) {
    x_diff = difference(x, d)
    stan_data <- list(num_obs = length(x_diff), y = x_diff, p = order[1], q = order[3], h = h)
  } else {
    stan_data <- list(num_obs = length(x), y = x, p = order[1], q = order[3], h = h)
  }
  stan_fit <- rstan::stan(file = 'stan_models/ARMA_simple.stan', data = stan_data)
  return(stan_fit)
}

comparison_plot <- function(frequentist_fit, bayesian_fit, x_observed, x_future) {
  h = length(x_future)
  size = length(x_observed) + h
  d = length(frequentist_fit[["model"]][["Delta"]])
  x <- c(x_observed, x_future)

  # inputs for plot
  frequentist_y_fit <- c(fitted(frequentist_fit), rep_len(NA, h))
  frequentist_y_forecast <- c(rep_len(NA, size - h - 1), frequentist_y_fit[size-h],
                              as.numeric(predict(frequentist_fit, n.ahead = h)$pred))

  # note on command above: in order to connect the fit line with the forecast line, included one less NA and added x_observed[size-h]
  frequentist_pred_interval <- forecast(frequentist_fit, h=h, level=90)
  frequentist_ub <- c(fitted(frequentist_fit) + 1.64*sqrt(frequentist_fit$sigma2),
             frequentist_pred_interval[["upper"]])
  frequentist_lb <- c(fitted(frequentist_fit) - 1.64*sqrt(frequentist_fit$sigma2),
             frequentist_pred_interval[["lower"]])

  y_pred_mean <- obtain_mean_pred(bayesian_fit, x_observed, d, h)

  y_fit <- rstan::extract(bayesian_fit, 'y_fit')$y_fit
  y_pred <- rstan::extract(bayesian_fit, 'y_pred')$y_pred

  # reverse differencing: to do: d in funciton not used. removed. only performing one stsep.
  if (d) {
    for (d_i in d:1) {
      x_observed_diff_above <- difference(x_observed, d_i-1)
      y_fit <- reverse_differencing_y_fit(y_fit, x_observed_diff_above)
      y_pred <- reverse_differencing_y_forecast(y_pred, x_observed_diff_above)
      }
    }

  y_pred_sims <- cbind(y_fit, y_pred)
  row.names(y_pred_sims)  <- 1:dim(y_pred_sims)[1]

  y_pred_lb <- apply(y_pred_sims, MARGIN = 2, FUN = function(x) {quantile(x, probs = c(.05))})
  y_pred_ub <- apply(y_pred_sims, MARGIN = 2, FUN = function(x) {quantile(x, probs = c(.95))})
  y_fit_mean <- replace(y_pred_mean, (size - h+1):length(y_pred_mean), NA)
  y_forecast_mean <- replace(y_pred_mean, 1:(size - h - 1), NA)
  # note on line above: in order to connect the y_forecast_mean to the y_fit_mean line, the y_forecast_mean contains h + 1 non-NA values, the first one of which corresponds to the latest actual.

  # plot
  p <- ggplot(mapping = aes(x = 1:size)) +
    geom_ribbon(
      aes(ymin = y_pred_lb, ymax = y_pred_ub, fill = "Bayesian credible interval"),
      alpha = .2
    ) +
    geom_ribbon(
      aes(ymin = frequentist_lb, ymax = frequentist_ub, fill = "Prediction interval for frequentist model"),
      alpha = .2
    ) +
    geom_line(mapping = aes(y = y_fit_mean,
                            linetype = "In-sample (Fitted values)", color = "Bayesian")) +
    geom_line(mapping = aes(y = y_forecast_mean,
                            linetype = "Out-of-sample (Forecasted values)",
                            color = "Bayesian")) +
    geom_line(mapping = aes(y = frequentist_y_fit,
                            linetype = "In-sample (Fitted values)", color = "Frequentist")) +
    geom_line(mapping = aes(y = frequentist_y_forecast,
                            linetype = "Out-of-sample (Forecasted values)",
                            color = "Frequentist")) +
    geom_point(
      aes(y = x_observed, x = 1:(size-h), shape = 'Training')
    ) +
    geom_point(
      aes(y = x_future, x = (size-h+1):size, shape = 'Testing')
    ) +
    xlab("Simulation #") + ylab("Value") +
    scale_fill_manual(name = "90% Interval",
                      values = c("Bayesian credible interval" = "blue",
                      "Prediction interval for frequentist model" = "orange")) +
    scale_shape_manual(name = "Simulated datasets", values = c("Training" = 20, "Testing" = 21)) +
    scale_linetype_manual(name = "Prediction type", values = c("solid", "dashed")) +
    scale_colour_manual(name = "Statistical paradigm", values =c('blue', "orange"))

  return(p)
}



obtain_mean_pred <- function(bayesian_fit, x, d, h) {
  p = bayesian_fit@par_dims[["phi"]]
  q = bayesian_fit@par_dims[["theta"]]

  sigma_sims <- rstan::extract(bayesian_fit, "sigma")$sigma
  if (p){
    phi_sims <- rstan::extract(bayesian_fit, "phi")$phi
    phi_sims_mean <- apply(phi_sims, 2, mean)
  }
  if (q) {
    theta_sims = rstan::extract(bayesian_fit, "theta")$theta
    theta_sims_mean <- apply(theta_sims, 2, mean)
  }
  if (d) {
    x_without_diff <- x
    x <- difference(x, d)
  }
  num_obs = length(x)

  # fitted values
  err_mean <- c(0)
  x_fit_mean <- c(0)
  for (t in 2:num_obs) {
    eta <- 0
    if (p) for (i in 1:min(p, t-1)) {eta =+ x[t-i] * phi_sims_mean[i]}
    if (q) for (i in 1:min(q, t-1)) {eta =+ err_mean[t-i] * theta_sims_mean[i]}
    eta_show <- eta
    x_fit_mean[t] <- eta
    err_mean[t] <- x[t] - x_fit_mean[t]
  }
  if (d) for (i in d:1) {
    x_actual_above <- difference(x_without_diff, (i-1))
    x_fit_mean <- reverse_diff_other_vector(x_fit_mean, x_actual_above)
  }

  # forecasted values
  x_forecast_mean <- c()
  # No errors for mean values.

  for (t in 1:h) {
    p_new = min(t-1, p) # new data
    p_rem = p - p_new 
    q_new = min(t-1, q)
    q_rem = q - q_new
    d_new <- min(t-1, d)
    d_rem <- d - d_new
    eta <- 0
    if (p_new) for (i in 1:p_new) {eta =+ x_forecast_mean[t-i] * phi_sims_mean[i]}
    if (p_rem) for (i in 1:p_rem) {eta =+ x[num_obs+1-i] * phi_sims_mean[i + p_new]}
    # for mean values: error = 0
    if (q_rem) for (i in 1:min(q_rem, t)) {eta =+ err_mean[num_obs + 1 - i] * theta_sims_mean[i+ q_new]}
    x_forecast_mean[t] <- eta
    }

  if (d) for (i in d:1) {
    x_actual_above <- difference(x_without_diff, i-1)
    x_forecast_mean <- reverse_diff_forecast_mean(x_forecast_mean, x_actual_above[length(x_actual_above)])
  }
  y_pred_mean <- c(x_fit_mean, x_forecast_mean)
  return(y_pred_mean)
}


######### Differencing and Reverse differencing functions ############

reverse_diff_forecast_mean <- function(x_diff, init_value) {
  x_diff[1] =+ init_value
  x <- cumsum(x_diff)
  return(x)
}

difference <- function(x, d) {
  if (d > 0) for (i in 1:d) {x <- diff(x)}
  return(x)
}

reverse_diff_other_vector <- function(x_diff_fit, x_actual) {

  x_fit <- c(0, x_diff_fit) + c(x_actual[1], x_actual[1:(length(x_actual)-1)])

  return(x_fit)
}

reverse_differencing_y_fit <- function(y_fit, y) { # y_above

  y_fit <- sweep(cbind(rep_len(0,dim(y_fit)[1]),y_fit),
                     2,
                     c(y[1],y[1:(length(y)-1)]),
                     "+")
  return(y_fit)
  }

reverse_differencing_y_forecast <- function(y_pred, y) { # init value is: y[length(y)]
    y_pred[,1] = y_pred[,1] + y[length(y)]

    for (i in 2:dim(y_pred)[2]) {
      y_pred[,i] <- y_pred[,i] + y_pred[,(i-1)]
    }
    return(y_pred)
  }
