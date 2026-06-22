
process_fit <- function(fit, var = "log_annual_precipitation") {

  if(var == "log_annual_precipitation") {
    x_transform <- exp
    global_c_center <- 6.7
  } else if (var == "average_temperature") {
    x_transform <- identity
    global_c_center <- 19.5
  } else {
    stop("Unsupported variable: ", var)
  }

  var_between <- paste0(var, "_mn")
  var_within <- paste0(var, "_c")

  message("\t\033[32mComputing global mean predictions by ", var_between, "...\033[0m")

  data_fit_global_mn <-
    estimate_means(fit, by = var_between) |>
    as_tibble() |>
    rename(estimate = Mean, conf.low = CI_low, conf.high = CI_high) |>
    mutate(
      taxon_name = "average_mn",
      annual_precipitation = x_transform(.data[[var_between]]),
      across(c("estimate", "conf.low", "conf.high"), ~ exp(.x)),
    )

  message("\t\033[32mComputing global centered predictions by ", var_within, "...\033[0m")
  data_fit_global_c <-
    estimate_means(fit, by = var_within) |>
    as_tibble() |>
    rename(estimate = Mean, conf.low = CI_low, conf.high = CI_high) |>
    mutate(
      taxon_name = "average_c",
      annual_precipitation = x_transform(global_c_center + .data[[var_within]]),
      across(c("estimate", "conf.low", "conf.high"), ~ exp(.x))
    )

  message("\t\033[32mComputing individual predictions for all observations...\033[0m")
  data_fit_ind <-
    marginaleffects::predictions(fit) |>
    as_tibble() |>
    select(
      any_of(c("identifier", "taxon_name", "decimal_latitude", "decimal_longitude",
               "annual_precipitation", var)),
      all_of(c(var_between, var_within)),
      any_of(c("value", "log_value")),
      estimate, conf.low, conf.high
    ) |>
    mutate(
      log_residual = estimate - log_value,
      residual = exp(estimate) - exp(log_value),
      across(c("estimate", "conf.low", "conf.high"), ~ exp(.x))
    )

  data_fit_ind_range <-
    # use the original data but add predictions from the fitted model
    # =only keep first and last record for each species to show the range of rainfall
    data_fit_ind |>
    group_by(taxon_name) |>
    arrange(.data[[var]]) |>
    slice(c(1, n())) |>
    ungroup()

  message("\t\033[32mComputing predictions at species mean...\033[0m")
  data_fit_mn <-
    marginaleffects::predictions(fit,
      newdata = data_fit_ind |>
        mutate("{var_within}" := 0) |>
        select(all_of(c(var_between, var_within, "taxon_name"))) |>
        distinct()
    ) |>
    as_tibble() |>
     ## Keep the log estimates but also convert to original scale for plotting
     mutate(
       # Convert back to original scale
       across(c("estimate", "conf.low", "conf.high"), ~ exp(.x))
     )

    if(var ==  "log_annual_precipitation") {
    data_fit_mn <-
      data_fit_mn |>
      mutate(
        annual_precipitation = x_transform(.data[[var_between]])
      )
    }

  message("\t\033[32mDone processing fit!\033[0m")

  list(
    global_mn = data_fit_global_mn,
    global_c = data_fit_global_c,
    ind = data_fit_ind,
    ind_range = data_fit_ind_range,
    mn = data_fit_mn
  )
}
