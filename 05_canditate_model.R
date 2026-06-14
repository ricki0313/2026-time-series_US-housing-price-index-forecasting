# ============================================================
# Function
# ============================================================

check_required_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required. Please install it first."))
    }
  }
}

load_housing_data <- function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  required_cols <- c("date", "hpi")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  df$date <- as.Date(df$date)
  df <- df[order(df$date), ]
  
  if (!"log_hpi" %in% names(df)) {
    df$log_hpi <- log(df$hpi)
  }
  
  if (!"diff_log_hpi" %in% names(df)) {
    df$diff_log_hpi <- c(NA, diff(df$log_hpi))
  }
  
  return(df)
}

make_monthly_ts <- function(values, dates) {
  start_year <- as.integer(format(min(dates), "%Y"))
  start_month <- as.integer(format(min(dates), "%m"))
  
  ts(
    values,
    start = c(start_year, start_month),
    frequency = 12
  )
}

add_dummy_variables <- function(df) {
  df$post2002 <- ifelse(df$date >= as.Date("2002-01-01"), 1, 0)
  
  df$crisis <- ifelse(
    df$date >= as.Date("2007-01-01") &
      df$date <= as.Date("2010-12-01"),
    1,
    0
  )
  
  return(df)
}

prepare_model_inputs <- function(df) {
  df <- add_dummy_variables(df)
  
  log_hpi_ts <- make_monthly_ts(
    values = df$log_hpi,
    dates = df$date
  )
  
  df_diff <- df[!is.na(df$diff_log_hpi), ]
  
  diff_log_hpi_ts <- make_monthly_ts(
    values = df_diff$diff_log_hpi,
    dates = df_diff$date
  )
  
  xreg_crisis <- as.matrix(df_diff[, c("crisis"), drop = FALSE])
  xreg_both <- as.matrix(df_diff[, c("post2002", "crisis"), drop = FALSE])
  
  list(
    df = df,
    df_diff = df_diff,
    log_hpi_ts = log_hpi_ts,
    diff_log_hpi_ts = diff_log_hpi_ts,
    xreg_crisis = xreg_crisis,
    xreg_both = xreg_both
  )
}

fit_candidate_models <- function(inputs) {
  models <- list(
    # Baseline model
    "Random Walk with Drift" = forecast::Arima(
      inputs$log_hpi_ts,
      order = c(0, 1, 0),
      include.drift = TRUE,
      method = "ML"
    ),
    
    # ARIMA models with drift
    "ARIMA(1,1,0) with drift" = forecast::Arima(
      inputs$log_hpi_ts,
      order = c(1, 1, 0),
      include.drift = TRUE,
      method = "ML"
    ),
    
    "ARIMA(1,1,1) with drift" = forecast::Arima(
      inputs$log_hpi_ts,
      order = c(1, 1, 1),
      include.drift = TRUE,
      method = "ML"
    ),
    
    "ARIMA(2,1,1) with drift" = forecast::Arima(
      inputs$log_hpi_ts,
      order = c(2, 1, 1),
      include.drift = TRUE,
      method = "ML"
    ),
    
    # SARIMA models with drift
    "SARIMA(1,1,1)(1,0,0)[12] with drift" = forecast::Arima(
      inputs$log_hpi_ts,
      order = c(1, 1, 1),
      seasonal = list(order = c(1, 0, 0), period = 12),
      include.drift = TRUE,
      method = "ML"
    ),
    
    "SARIMA(1,1,1)(0,0,1)[12] with drift" = forecast::Arima(
      inputs$log_hpi_ts,
      order = c(1, 1, 1),
      seasonal = list(order = c(0, 0, 1), period = 12),
      include.drift = TRUE,
      method = "ML"
    ),
    
    # ARIMAX models with dummy variables
    # These models use dummy variables to capture changes in average growth.
    "ARIMAX(1,1,1) with crisis dummy" = forecast::Arima(
      inputs$diff_log_hpi_ts,
      order = c(1, 0, 1),
      xreg = inputs$xreg_crisis,
      include.mean = TRUE,
      method = "ML"
    ),
    
    "ARIMAX(1,1,1) with post-2002 dummy and crisis dummy" = forecast::Arima(
      inputs$diff_log_hpi_ts,
      order = c(1, 0, 1),
      xreg = inputs$xreg_both,
      include.mean = TRUE,
      method = "ML"
    ),
    
    # SARIMAX models with dummy variables
    # These models add the 12-month dependence structure and use dummy variables
    # to capture changes in average growth.
    "SARIMAX(1,1,1)(1,0,0)[12] with crisis dummy" = forecast::Arima(
      inputs$diff_log_hpi_ts,
      order = c(1, 0, 1),
      seasonal = list(order = c(1, 0, 0), period = 12),
      xreg = inputs$xreg_crisis,
      include.mean = TRUE,
      method = "ML"
    ),
    
    "SARIMAX(1,1,1)(1,0,0)[12] with post-2002 dummy and crisis dummy" = forecast::Arima(
      inputs$diff_log_hpi_ts,
      order = c(1, 0, 1),
      seasonal = list(order = c(1, 0, 0), period = 12),
      xreg = inputs$xreg_both,
      include.mean = TRUE,
      method = "ML"
    )
  )
  
  return(models)
}

extract_ic_table <- function(models) {
  table <- do.call(
    rbind,
    lapply(names(models), function(model_name) {
      fit <- models[[model_name]]
      
      data.frame(
        Model = model_name,
        AIC = as.numeric(AIC(fit)),
        AICc = as.numeric(fit$aicc),
        BIC = as.numeric(BIC(fit)),
        row.names = NULL
      )
    })
  )
  
  table <- table[order(table$AICc), ]
  
  return(table)
}

extract_parameter_table <- function(models) {
  table <- do.call(
    rbind,
    lapply(names(models), function(model_name) {
      fit <- models[[model_name]]
      coefs <- stats::coef(fit)
      
      if (length(coefs) == 0) {
        return(NULL)
      }
      
      se <- rep(NA_real_, length(coefs))
      names(se) <- names(coefs)
      
      if (!is.null(fit$var.coef)) {
        available_se <- sqrt(diag(fit$var.coef))
        se[names(available_se)] <- available_se
      }
      
      z_value <- coefs / se
      p_value <- 2 * (1 - stats::pnorm(abs(z_value)))
      
      data.frame(
        Model = model_name,
        Parameter = names(coefs),
        Estimate = as.numeric(coefs),
        Std_Error = as.numeric(se),
        z_value = as.numeric(z_value),
        p_value = as.numeric(p_value),
        row.names = NULL
      )
    })
  )
  
  return(table)
}

print_model_tables <- function(ic_table, parameter_table) {
  cat("\n================ Information Criteria Table ================\n")
  print(ic_table, row.names = FALSE, digits = 5)
  
  cat("\n================ Parameter Estimation Table ================\n")
  print(parameter_table, row.names = FALSE, digits = 5)
}

print_key_summaries <- function(models) {
  cat("\n================ ARIMA(1,1,1) with drift Summary ================\n")
  print(summary(models[["ARIMA(1,1,1) with drift"]]))
  
  cat("\n===== ARIMAX(1,1,1) with post-2002 dummy and crisis dummy Summary =====\n")
  print(summary(models[[
    "ARIMAX(1,1,1) with post-2002 dummy and crisis dummy"
  ]]))
  
  cat("\n===== SARIMAX(1,1,1)(1,0,0)[12] with post-2002 dummy and crisis dummy Summary =====\n")
  print(summary(models[[
    "SARIMAX(1,1,1)(1,0,0)[12] with post-2002 dummy and crisis dummy"
  ]]))
}

run_model_estimation <- function(file_path) {
  set.seed(123)
  
  check_required_packages(c("forecast"))
  
  df <- load_housing_data(file_path)
  inputs <- prepare_model_inputs(df)
  models <- fit_candidate_models(inputs)
  
  ic_table <- extract_ic_table(models)
  parameter_table <- extract_parameter_table(models)
  
  print_model_tables(
    ic_table = ic_table,
    parameter_table = parameter_table
  )
  
  print_key_summaries(models)
  
  results <- list(
    data = df,
    inputs = inputs,
    models = models,
    ic_table = ic_table,
    parameter_table = parameter_table
  )
  
  return(invisible(results))
}

# ============================================================
# Main
# ============================================================

# Run the analysis and print results in the console.
invisible(run_model_estimation(
  file_path = "C:/Users/use/Desktop/NCCU_STAT/Time Series/final/data/processed/monthly_housing_train.csv"
))