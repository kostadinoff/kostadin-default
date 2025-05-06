###### ---Libraries---######

if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, haven, modelsummary, MKinfer, rstatix, finalfit, tinytable, monochromeR, ggstats, epitools, ggsurvfit, broom, rstan, brms, gtsummary, ggplot2, quantreg, patchwork, tidymodels, gt, epiR, readxl, scales, marginaleffects, moments, ggthemes, entropy, emmeans, janitor, easystats, showtext, tseries, rugarch, forecast, fractal, brglm2, sysfonts, pracma, MASS, zoo, detectseparation)

###### ---Options---######

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
options(ggplot2.messages = FALSE)
options(dplyr.width = Inf)
###### ---Functions---######

# Function to round numeric columns in a tibble

mutate_round <- function(data, digits = 2) {
  if (!is.data.frame(data)) {
    stop("Input must be a data frame or tibble.")
  }
  data |>
    dplyr::mutate(across(where(is.numeric), ~ janitor::round_half_up(., digits)))
}

# Function to calculate confidence intervals for proportions

pcit <- function(data, conf.level = 0.95) {
  # Identify numeric columns
  numeric_cols <- sapply(data, is.numeric)

  # Ensure there are at least two numeric columns
  if (sum(numeric_cols) < 2) {
    stop("The dataset must contain at least two numeric columns.")
  }

  # Get indices of the first and second numeric columns
  numeric_col_indices <- which(numeric_cols)

  if (length(numeric_col_indices) < 2) {
    stop("Not enough numeric columns found.")
  }

  # Extract the first and second numeric columns
  x <- data[[numeric_col_indices[1]]]
  n <- data[[numeric_col_indices[2]]]

  # Ensure that n >= x for all cases
  if (any(n < x)) {
    stop("Each trial count must be greater than or equal to the success count.")
  }

  # Initialize vectors for confidence intervals
  lower <- numeric(nrow(data))
  upper <- numeric(nrow(data))

  # Calculate confidence intervals for each row
  for (i in seq_len(nrow(data))) {
    ci <- binom.test(
      x = x[i],
      n = n[i],
      conf.level = conf.level
    )$conf.int
    lower[i] <- ci[1]
    upper[i] <- ci[2]
  }

  # Identify non-numeric columns for grouping
  non_numeric_cols <- names(data)[!numeric_cols]

  # Create a tibble including all original columns and calculated results
  result <- tibble(
    # Include non-numeric columns
    data[, non_numeric_cols, drop = FALSE],

    # Calculated columns with distinct names
    successes = x,
    trials = n,
    proportion = x / n,
    lower = lower,
    upper = upper,
    conf.level = conf.level
  )

  return(result)
}

# Function to compare proportions between groups

compare_proportions <- function(data,
                                conf.level = 0.95,
                                method = "holm") {
  # Validate the input data
  required_cols <- c("proportion", "trials")
  if (!all(required_cols %in% names(data))) {
    stop("The input data must contain 'proportion' and 'trials' columns.")
  }
  # Identify group columns dynamically (assuming non-numeric columns are groups)
  non_numeric_cols <- names(data)[!sapply(data, is.numeric)]
  if (length(non_numeric_cols) == 0) {
    stop("No non-numeric columns found for grouping.")
  }
  # Get all combinations of non-numeric columns for group identification
  group_cols <- non_numeric_cols
  # Ensure there are at least two groups to compare
  if (nrow(data) < 2) {
    stop("The dataset must contain at least two groups for comparison.")
  }
  # Create a grid of pairwise combinations of the groups
  combinations <- combn(seq_len(nrow(data)), 2, simplify = FALSE)
  # Function to calculate differences, p-values, and confidence intervals for each pair
  calculate_difference <- function(index_pair) {
    i <- index_pair[1]
    j <- index_pair[2]

    group1 <- data[i, ]
    group2 <- data[j, ]

    # Proportion difference
    prop_diff <- group1$proportion - group2$proportion
    # Standard error for the difference in proportions
    pooled_se <- sqrt(
      group1$proportion * (1 - group1$proportion) / group1$trials +
        group2$proportion * (1 - group2$proportion) / group2$trials
    )

    small_constant <- 1e-10
    pooled_se <- max(pooled_se, small_constant)

    # Prevent division by zero
    if (pooled_se == 0) {
      stop("Standard error is zero; cannot compute z-score and p-value.")
    }
    z_score <- prop_diff / pooled_se
    # Calculate p-value for the proportion difference
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    # Calculate confidence interval for the proportion difference
    z_critical <- qnorm(1 - (1 - conf.level) / 2)
    ci_lower <- prop_diff - z_critical * pooled_se
    ci_upper <- prop_diff + z_critical * pooled_se

    result <- tibble(
      !!!setNames(
        lapply(group_cols, function(col) {
          group1[[col]]
        }),
        paste0("gr1_", group_cols)
      ),
      !!!setNames(
        lapply(group_cols, function(col) {
          group2[[col]]
        }),
        paste0("gr2_", group_cols)
      ),
      prop_diff = prop_diff,
      z_score = z_score,
      p_value = p_value,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )

    return(result)
  }

  # Apply the calculate_difference function to each pair
  results <- map_dfr(combinations, calculate_difference)

  # Adjust p-values for multiple comparisons
  results <- results %>%
    mutate(adj_p_value = p.adjust(p_value, method = method))

  return(results)
}


# Function to calculate the proportion stratified by a grouping variable

compare_proportions_by <- function(data,
                                   conf.level = 0.95,
                                   method = "holm") {
  # Validate the input data
  required_cols <- c("proportion", "trials")
  if (!all(required_cols %in% names(data))) {
    stop("The input data must contain 'proportion' and 'trials' columns.")
  }

  # Automatically determine group and subgroup variables
  group_var <- names(data)[1]
  subgroup_var <- names(data)[2]

  # Split the data by group
  grouped_data <- split(data, data[[group_var]])

  # Function to calculate differences, p-values, and confidence intervals for each subgroup
  calculate_difference_within_group <- function(df) {
    # Generate pairwise comparisons for subgroups
    pairwise_comparisons <- combn(unique(df[[subgroup_var]]), 2, simplify = FALSE)

    # Initialize an empty result dataframe
    result_list <- list()

    for (pair in pairwise_comparisons) {
      subgroup1 <- df %>% filter(!!sym(subgroup_var) == pair[1])
      subgroup2 <- df %>% filter(!!sym(subgroup_var) == pair[2])

      # Ensure there is exactly one row in each subgroup
      if (nrow(subgroup1) != 1 || nrow(subgroup2) != 1) {
        warning("Each subgroup should have exactly one entry for comparison.")
        next
      }

      # Calculate proportion difference
      prop_diff <- subgroup1$proportion - subgroup2$proportion

      # Standard error for the difference in proportions
      pooled_se <- sqrt(
        subgroup1$proportion * (1 - subgroup1$proportion) / subgroup1$trials +
          subgroup2$proportion * (1 - subgroup2$proportion) / subgroup2$trials
      )
      small_constant <- 1e-10
      pooled_se <- max(pooled_se, small_constant)
      # Prevent division by zero
      if (pooled_se == 0) {
        warning("Standard error is zero; cannot compute z-score and p-value.")
        next
      }

      z_score <- prop_diff / pooled_se

      # Calculate p-value for the proportion difference
      p_value <- 2 * (1 - pnorm(abs(z_score)))

      # Calculate confidence interval for the proportion difference
      z_critical <- qnorm(1 - (1 - conf.level) / 2)
      ci_lower <- prop_diff - z_critical * pooled_se
      ci_upper <- prop_diff + z_critical * pooled_se

      # Append the results for this pairwise comparison
      result_list[[length(result_list) + 1]] <- tibble(
        group = unique(df[[group_var]]),
        subgroup1 = pair[1],
        subgroup2 = pair[2],
        prop_diff = prop_diff,
        z_score = z_score,
        p_value = as.numeric(p_value),
        # Ensure p_value is numeric
        ci_lower = ci_lower,
        ci_upper = ci_upper
      )
    }

    # Combine all results into one dataframe
    bind_rows(result_list)
  }

  # Apply the calculation function to each group
  results <- map_dfr(grouped_data, calculate_difference_within_group)

  # Adjust p-values for multiple comparisons if p_value column is correctly numeric
  if (nrow(results) > 0) {
    results <- results %>%
      mutate(adj_p_value = p.adjust(p_value, method = method)) %>%
      dplyr::select(
        group,
        subgroup1,
        subgroup2,
        prop_diff,
        z_score,
        p_value,
        ci_lower,
        ci_upper,
        adj_p_value
      )
  } else {
    results <- results %>%
      mutate(adj_p_value = NA) %>%
      dplyr::select(
        group,
        subgroup1,
        subgroup2,
        prop_diff,
        z_score,
        p_value,
        ci_lower,
        ci_upper,
        adj_p_value
      )
  }

  return(results)
}

# Univariate categorical plot

univariate_cat_plot <- function(data, variable) {
  variable <- rlang::ensym(variable)
  title <- paste("Univariate Categorical Plot of", rlang::as_name(variable))

  # Filter out non-finite values (NA, NaN)
  variable_data <- dplyr::pull(data, !!variable)
  variable_data <- variable_data[!is.na(variable_data) & !is.nan(variable_data)]

  data %>%
    count({{ variable }}) %>%
    filter(!is.na({{ variable }})) %>%
    mutate(prop = n / sum(n)) %>%
    ggplot(aes(y = fct_reorder({{ variable }}, prop), x = prop)) +
    geom_col(
      alpha = 0.6,
      fill = "gray60",
      color = "black"
    ) +
    geom_label(
      aes(label = paste0(n, " (", scales::percent(prop), ")")),
      color = "black",
      size = 5,
      family = "Roboto Condensed",
      hjust = -0.1
    ) +
    scale_x_continuous(labels = scales::percent) +
    expand_limits(x = 1) +
    guides(
      x = guide_axis(cap = "both"),
      y = guide_axis(cap = "both")
    ) +
    labs(
      x = "Proportion",
      y = "Category",
      title = title
    )
}

# Univariate continuous plot
univariate_cont_plot <- function(data, variable) {
  variable <- rlang::ensym(variable)
  title <- paste("Univariate Continuous Plot of", rlang::as_name(variable))

  data %>%
    ggplot(aes(x = !!variable)) +
    geom_density(
      adjust = 1 / 2,
      fill = "gray90",
      color = "black",
      alpha = 0.6
    ) +
    geom_vline(
      aes(xintercept = mean({{ variable }}, na.rm = TRUE)),
      color = "red",
      linetype = 1,
      linewidth = 1
    ) +
    geom_vline(
      aes(xintercept = median({{ variable }}, na.rm = TRUE)),
      color = "blue",
      linetype = "dashed",
      linewidth = 1,
    ) +
    annotate(
      "label",
      x = mean(data[[rlang::as_name(variable)]], na.rm = TRUE),
      y = 0.9 * max(density(data[[rlang::as_name(variable)]], na.rm = TRUE)$y),
      label = paste0("Mean: ", round(mean(data[[rlang::as_name(variable)]], na.rm = TRUE), 2)),
      color = "red",
      size = 5,
      family = "Roboto Condensed",
      hjust = -0.1
    ) +
    annotate(
      "label",
      x = median(data[[rlang::as_name(variable)]], na.rm = TRUE),
      y = 0.8 * max(density(data[[rlang::as_name(variable)]], na.rm = TRUE)$y),
      label = paste0("Median: ", round(median(data[[rlang::as_name(variable)]], na.rm = TRUE), 2)),
      color = "blue",
      size = 5,
      family = "Roboto Condensed",
      hjust = -0.1
    ) +
    expand_limits(x = c(min(data[[rlang::as_name(variable)]], na.rm = TRUE), max(data[[rlang::as_name(variable)]], na.rm = TRUE))) +
    guides(
      x = guide_axis(cap = "both"),
      y = guide_axis(cap = "both")
    ) +
    labs(
      x = "Value",
      y = "Density",
      title = title
    )
}

# Univariate plot: Check variable type and apply the appropriate plot
univariate_plot <- function(data, variable) {
  variable <- rlang::ensym(variable)

  # Check if the variable is numeric or categorical
  if (is.numeric(data[[rlang::as_name(variable)]])) {
    # Call the continuous plot function
    univariate_cont_plot(data, !!variable)
  } else if (is.factor(data[[rlang::as_name(variable)]]) || is.character(data[[rlang::as_name(variable)]])) {
    # Call the categorical plot function
    univariate_cat_plot(data, !!variable)
  } else {
    stop("The variable must be either numeric or categorical (factor/character).")
  }
}

# Function for regression

kk_reg <- function(data, outcome, predictors, log_outcome = FALSE, custom_formula = NULL, ...) {
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }
  if (!rlang::is_string(outcome)) {
    stop("`outcome` must be a character string.")
  }
  if (!is.character(predictors)) {
    stop("`predictors` must be a character vector.")
  }
  if (!is.logical(log_outcome)) {
    stop("`log_outcome` must be a logical value (TRUE or FALSE).")
  }

  # Capture outcome as symbol
  outcome <- rlang::ensym(outcome)

  # Ensure outcome exists in the data
  if (!rlang::as_name(outcome) %in% colnames(data)) {
    stop("Outcome variable '", rlang::as_name(outcome), "' not found in the data.")
  }

  # Ensure predictors exist in the data
  missing_predictors <- setdiff(predictors, colnames(data))
  if (length(missing_predictors) > 0) {
    stop("Predictor variable(s) not found in the data: ", paste(missing_predictors, collapse = ", "))
  }

  # If log_outcome is TRUE, transform the outcome variable
  if (log_outcome) {
    data <- data %>% mutate(!!outcome := log(!!outcome))
  }

  # Determine the type of outcome variable
  outcome_type <- case_when(
    is.factor(data[[rlang::as_name(outcome)]]) && nlevels(data[[rlang::as_name(outcome)]]) == 2 ~ "binary", # Binary outcome
    is.ordered(data[[rlang::as_name(outcome)]]) ~ "ordinal", # Ordinal outcome
    is.numeric(data[[rlang::as_name(outcome)]]) ~ "continuous", # Continuous outcome
    TRUE ~ NA_character_
  )

  if (is.na(outcome_type)) {
    stop("Outcome must be numeric, binary factor, or ordinal factor.")
  }

  # Function to fit a model
  fit_model <- function(formula) {
    if (outcome_type == "binary") {
      return(glm(formula, data = data, family = binomial(), ...)) # Logistic regression for binary outcome
    } else if (outcome_type == "ordinal") {
      return(MASS::polr(formula, data = data, Hess = TRUE, ...)) # Ordinal regression
    } else if (outcome_type == "continuous") {
      return(lm(formula, data = data, ...)) # Linear regression for continuous outcome
    }
  }

  # Function to process model output
  process_results <- function(model, model_type) {
    results <- broom::tidy(model, conf.int = TRUE) %>%
      mutate(
        model_type = model_type,
        outcome_type = outcome_type,
        AIC = AIC(model),
        BIC = BIC(model)
      )

    if (outcome_type == "binary") {
      # Convert to odds ratios (OR) for logistic regression
      results <- results %>%
        mutate(
          estimate = exp(estimate), # Odds ratio for logistic regression
          conf.low = exp(conf.low), # Confidence interval for OR
          conf.high = exp(conf.high),
          percent_change = (estimate - 1) * 100, # Percentage change for OR
          percent_change_low = (conf.low - 1) * 100, # Lower CI for percentage change
          percent_change_high = (conf.high - 1) * 100 # Upper CI for percentage change
        )
    } else if (outcome_type == "continuous" && log_outcome) {
      # Convert to percentage change if linear regression with log outcome
      results <- results %>%
        mutate(
          percent_change = (exp(estimate) - 1) * 100,
          percent_change_low = (exp(conf.low) - 1) * 100,
          percent_change_high = (exp(conf.high) - 1) * 100
        )
    } else if (outcome_type == "continuous") {
      # Calculate percentage change for continuous (no log transformation)
      results <- results %>%
        mutate(
          percent_change = estimate * 100, # Simple percentage change for continuous variable
          percent_change_low = conf.low * 100,
          percent_change_high = conf.high * 100
        )
    } else {
      # No transformation for normal linear regression
      results <- results %>%
        mutate(percent_change = NA, percent_change_low = NA, percent_change_high = NA)
    }

    return(results)
  }

  # Fit univariate models
  univariate_results <- map_dfr(predictors, function(predictor) {
    formula <- as.formula(paste(rlang::as_name(outcome), "~", predictor))
    model <- fit_model(formula)
    process_results(model, "univariate") %>%
      mutate(predictor = predictor)
  })

  # Fit multivariate model
  if (is.null(custom_formula)) {
    multivariate_formula <- as.formula(paste(rlang::as_name(outcome), "~", paste(predictors, collapse = " + ")))
  } else {
    multivariate_formula <- custom_formula
  }
  multivariate_model <- fit_model(multivariate_formula)
  multivariate_results <- process_results(multivariate_model, "multivariate")

  # Combine results
  results <- bind_rows(univariate_results, multivariate_results)

  return(results)
}
# Additional function for regression

krk_reg <- function(data, outcome, predictors, log_outcome = FALSE, custom_formula = NULL, include_diagnostics = TRUE, ...) {
  # Load required packages
  library(dplyr)
  library(rlang)
  library(broom)
  library(purrr)

  # Validate inputs
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!is_string(outcome)) stop("`outcome` must be a character string.")
  if (!is.character(predictors)) stop("`predictors` must be a character vector.")
  if (!is.logical(log_outcome)) stop("`log_outcome` must be a logical value (TRUE or FALSE).")
  if (!is.logical(include_diagnostics)) stop("`include_diagnostics` must be a logical value.")

  # Capture outcome as symbol
  outcome <- ensym(outcome)
  outcome_name <- as_name(outcome)

  # Check if outcome exists
  if (!outcome_name %in% colnames(data)) {
    stop("Outcome variable '", outcome_name, "' not found in the data.")
  }

  # Check if predictors exist
  missing_predictors <- setdiff(predictors, colnames(data))
  if (length(missing_predictors) > 0) {
    stop("Predictor variable(s) not found in the data: ", paste(missing_predictors, collapse = ", "))
  }

  # Transform outcome if log_outcome is TRUE
  if (log_outcome) {
    if (!is.numeric(data[[outcome_name]]) || any(data[[outcome_name]] <= 0, na.rm = TRUE)) {
      stop("`log_outcome = TRUE` requires a positive numeric outcome.")
    }
    data <- data %>% mutate(!!outcome := log(!!outcome))
  }

  # Determine outcome type
  outcome_type <- case_when(
    is.factor(data[[outcome_name]]) && nlevels(data[[outcome_name]]) == 2 ~ "binary",
    is.ordered(data[[outcome_name]]) ~ "ordinal",
    is.numeric(data[[outcome_name]]) ~ "continuous",
    TRUE ~ NA_character_
  )
  if (is.na(outcome_type)) stop("Outcome must be numeric, binary factor, or ordered factor.")

  # Function to fit a model
  fit_model <- function(formula) {
    tryCatch(
      {
        if (outcome_type == "binary") {
          glm(formula, data = data, family = binomial(), ...)
        } else if (outcome_type == "ordinal") {
          if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required for ordinal regression.")
          MASS::polr(formula, data = data, Hess = TRUE, ...)
        } else if (outcome_type == "continuous") {
          lm(formula, data = data, ...)
        }
      },
      error = function(e) {
        stop("Model fitting failed: ", e$message)
      }
    )
  }

  # Function to calculate diagnostics (including model significance)
  calculate_diagnostics <- function(model) {
    if (!include_diagnostics) {
      return(tibble())
    }

    diagnostics <- list()
    summary_model <- summary(model)

    # Add model significance test
    if (outcome_type == "continuous") {
      diagnostics$r_squared <- summary_model$r.squared
      diagnostics$adj_r_squared <- summary_model$adj.r.squared
      diagnostics$residual_std_error <- summary_model$sigma
      diagnostics$model_p_value <- pf(summary_model$fstatistic[1],
        summary_model$fstatistic[2],
        summary_model$fstatistic[3],
        lower.tail = FALSE
      )
    } else if (outcome_type == "binary") {
      null_model <- update(model, ~1)
      loglik_model <- as.numeric(logLik(model))
      loglik_null <- as.numeric(logLik(null_model))
      diagnostics$pseudo_r_squared <- as.numeric(1 - (loglik_model / loglik_null))
      diagnostics$nagelkerke_r_squared <- as.numeric((1 - exp(-2 * (loglik_model - loglik_null))) /
        (1 - exp(2 * loglik_null / nrow(data))))
      diagnostics$model_p_value <- pchisq(2 * (loglik_model - loglik_null),
        df = length(coef(model)) - 1,
        lower.tail = FALSE
      )
      if (requireNamespace("pROC", quietly = TRUE)) {
        roc_curve <- pROC::roc(model$y, fitted(model), quiet = TRUE)
        diagnostics$auc_roc <- as.numeric(pROC::auc(roc_curve))
      }
    } else if (outcome_type == "ordinal") {
      null_model <- update(model, ~1)
      loglik_model <- as.numeric(logLik(model))
      loglik_null <- as.numeric(logLik(null_model))
      diagnostics$pseudo_r_squared <- as.numeric(1 - (loglik_model / loglik_null))
      diagnostics$nagelkerke_r_squared <- as.numeric((1 - exp(-2 * (loglik_model - loglik_null))) /
        (1 - exp(2 * loglik_null / nrow(data))))
      diagnostics$model_p_value <- pchisq(2 * (loglik_model - loglik_null),
        df = length(coef(model)),
        lower.tail = FALSE
      )
    }
    as_tibble(diagnostics)
  }

  # Function to process model output
  process_results <- function(model, model_type) {
    estimate_label <- if (outcome_type == "binary") {
      "odds_ratio"
    } else if (outcome_type == "continuous" && log_outcome) {
      "exp_coef"
    } else if (outcome_type == "ordinal") {
      "odds_ratio"
    } else {
      "coef"
    }

    results <- broom::tidy(model, conf.int = TRUE) %>%
      mutate(
        model_type = model_type,
        outcome_type = outcome_type,
        AIC = AIC(model),
        BIC = BIC(model),
        estimate_label = estimate_label,
        coef.type = if_else(grepl("\\|", term), "scale", "coefficient")
      )

    if (outcome_type %in% c("binary", "ordinal") || (outcome_type == "continuous" && log_outcome)) {
      results <- results %>%
        mutate(
          estimate = if_else(coef.type == "coefficient", exp(estimate), estimate),
          conf.low = if_else(coef.type == "coefficient", exp(conf.low), conf.low),
          conf.high = if_else(coef.type == "coefficient", exp(conf.high), conf.high),
          percent_change = if_else(coef.type == "coefficient", (estimate - 1) * 100, NA_real_),
          percent_change_low = if_else(coef.type == "coefficient", (conf.low - 1) * 100, NA_real_),
          percent_change_high = if_else(coef.type == "coefficient", (conf.high - 1) * 100, NA_real_)
        )
    } else if (outcome_type == "continuous") {
      results <- results %>%
        mutate(
          percent_change = estimate * 100,
          percent_change_low = conf.low * 100,
          percent_change_high = conf.high * 100
        )
    }

    diagnostics <- calculate_diagnostics(model)
    results %>% bind_cols(diagnostics)
  }

  # Fit univariate models
  univariate_results <- map_dfr(predictors, function(predictor) {
    formula <- as.formula(paste(outcome_name, "~", predictor))
    model <- fit_model(formula)
    process_results(model, "univariate") %>% mutate(predictor = predictor)
  })

  # Fit multivariate model
  if (is.null(custom_formula)) {
    multivariate_formula <- as.formula(paste(outcome_name, "~", paste(predictors, collapse = " + ")))
  } else {
    multivariate_formula <- custom_formula
  }
  multivariate_model <- fit_model(multivariate_formula)
  multivariate_results <- process_results(multivariate_model, "multivariate")

  # Combine results
  bind_rows(univariate_results, multivariate_results)
}

kkonehot <- function(data, column) {
  # Check if input is a tibble
  if (!inherits(data, "tbl_df")) {
    stop("Input 'data' must be a tibble")
  }

  # Check if column exists in data
  if (!column %in% names(data)) {
    stop("Specified column not found in the dataset")
  }

  # Get the variable
  var <- data[[column]]

  # Check if variable is character or factor
  if (!(is.character(var) || is.factor(var))) {
    stop("Column must be character or factor type for one-hot encoding")
  }

  # Convert character to factor if needed
  if (is.character(var)) {
    var <- as.factor(var)
  }

  # Get levels of the factor
  levels <- levels(var)

  # Create one-hot encoded columns
  for (level in levels) {
    new_col_name <- paste0(column, "_", level)
    data[[new_col_name]] <- ifelse(var == level, 1, 0)
  }

  # Return the modified dataset
  return(data)
}

###### ---Define theme---######

set_plot_font <- function(font = "Roboto Condensed", size = 18) {
  showtext_auto()

  # Try to add the Google font dynamically
  tryCatch(
    {
      sysfonts::font_add_google(name = font, family = font, db_cache = FALSE)
      message("Successfully loaded font: ", font)
    },
    error = function(e) {
      message("Font '", font, "' not found on Google Fonts. Falling back to 'Arial'.")
      # Use a pre-installed system font as fallback (no file path needed)
      font <<- "Arial" # Update font variable to use Arial
    }
  )

  # Define relative font sizes based on the `size` parameter
  title_size <- size + 4
  subtitle_size <- size + 2
  caption_size <- size - 2
  axis_title_size <- size
  axis_text_size <- size
  strip_text_size <- size

  theme_nice <- ggthemes::theme_tufte() +
    theme(
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.ticks.length = unit(4, "mm"),
      plot.title = element_text(family = font, size = title_size, hjust = 0, vjust = 2),
      plot.subtitle = element_text(family = font, size = subtitle_size),
      plot.caption = element_text(family = font, hjust = 0.5, vjust = 1, size = caption_size),
      plot.caption.position = "plot",
      axis.title = element_text(family = font, size = axis_title_size),
      axis.text = element_text(family = font, size = axis_text_size),
      axis.text.x = element_text(margin = margin(5, b = 10)),
      strip.text = element_text(family = font, size = strip_text_size),
      axis.line = element_line()
    )

  theme_set(theme_nice)
}


kkplot <- function(...) {
  ggplot(...) +
    guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both"))
}


# Function to compute full summary statistics with grouped data and pipe support
# Parameters:
# - data: Data frame containing the column of interest
# - col: Unquoted column name to analyze
# - var_name: Optional name for the variable (defaults to column name)
# - verbose: Output level ("full", "basic", or "custom")
# - stats: For verbose="custom", list of stats to include (e.g., c("mean", "median"))
# - pairwise: Whether to compute pairwise comparisons for categorical variables
# - chi_probs: Optional expected probabilities for chi-square test

kk_summary <- function(data, col, var_name = NULL,
                       verbose = c("full", "basic", "custom"),
                       stats = NULL, pairwise = TRUE, chi_probs = NULL) {
  # Validate inputs
  verbose <- match.arg(verbose)
  if (!is.data.frame(data)) stop("Input 'data' must be a data frame")

  # Get column name and set var_name
  col_quo <- enquo(col)
  col_name <- quo_name(col_quo)
  if (is.null(var_name)) var_name <- col_name

  # Extract the column as a vector
  x <- pull(data, {{ col }})

  # Core function to compute stats for a single vector
  compute_stats <- function(x, var_name) {
    # Validate input
    if (is.null(x)) stop("Input cannot be NULL")
    if (is.data.frame(x) || is.matrix(x)) stop("Input must be a vector, not a data frame or matrix")
    if (!is.vector(x) && !is.factor(x)) stop("Input must be a vector or factor")

    # Convert factor to character
    if (is.factor(x)) x <- as.character(x)

    # Initialize result
    result <- list()

    # Handle missing values
    n_miss <- sum(is.na(x))
    result$var_name <- var_name
    result$type <- ifelse(is.numeric(x), "numeric", "categorical")
    result$n_total <- length(x)
    result$n_miss <- n_miss
    result$n_valid <- length(x) - n_miss
    result$miss_pct <- (n_miss / length(x)) * 100

    # Return basic info if all values are NA
    if (result$n_valid == 0) {
      result$note <- "All values are missing"
      return(result)
    }

    # Define stats to include based on verbose
    if (verbose == "basic") {
      stats <- if (result$type == "numeric") {
        c("mean", "median", "sd", "min", "max")
      } else {
        c("level_summary", "mode")
      }
    } else if (verbose == "custom" && !is.null(stats)) {
      stats <- stats
    } else {
      stats <- if (result$type == "numeric") {
        c(
          "mean", "huber_mean", "trim_mean", "median", "min", "max", "range",
          "variance", "sd", "se", "cv_pct", "mad", "iqr", "q1", "q3", "skewness",
          "kurtosis", "pct_5_95", "ci_mean_low", "ci_mean_up", "shapiro_p",
          "shapiro_int", "n_outliers", "outlier_values"
        )
      } else {
        c(
          "levels", "n_unique", "mode", "level_summary", "entropy", "evenness",
          "chi_p", "chi_int", "chi_note", "cramer_v", "pairwise_p", "pairwise_note"
        )
      }
    }

    # Numeric variable statistics
    if (result$type == "numeric") {
      x <- na.omit(x) # Remove NAs
      if ("mean" %in% stats) result$mean <- mean(x)
      if ("huber_mean" %in% stats) {
        # Huber M-estimator with tuning constant c=1.345
        huber_m <- function(x, c = 1.345, max_iter = 100, tol = 1e-6) {
          mu <- median(x)
          sigma <- mad(x, constant = 1.4826)
          for (i in 1:max_iter) {
            residuals <- x - mu
            weights <- pmin(c / abs(residuals / sigma), 1)
            weights[is.na(weights)] <- 1
            mu_new <- sum(weights * x) / sum(weights)
            if (abs(mu_new - mu) < tol) break
            mu <- mu_new
          }
          return(mu)
        }
        result$huber_mean <- huber_m(x)
      }
      if ("trim_mean" %in% stats) result$trim_mean <- mean(x, trim = 0.1)
      if ("median" %in% stats) result$median <- median(x)
      if ("min" %in% stats) result$min <- min(x)
      if ("max" %in% stats) result$max <- max(x)
      if ("range" %in% stats) result$range <- max(x) - min(x)
      if ("variance" %in% stats) result$variance <- var(x)
      if ("sd" %in% stats) result$sd <- sd(x)
      if ("se" %in% stats) result$se <- sd(x) / sqrt(length(x))
      if ("cv_pct" %in% stats) result$cv_pct <- sd(x) / mean(x) * 100
      if ("mad" %in% stats) result$mad <- mad(x, constant = 1.4826)
      if ("iqr" %in% stats) result$iqr <- IQR(x)
      if ("q1" %in% stats) result$q1 <- quantile(x, 0.25)
      if ("q3" %in% stats) result$q3 <- quantile(x, 0.75)
      if ("skewness" %in% stats) result$skewness <- skewness(x)
      if ("kurtosis" %in% stats) result$kurtosis <- kurtosis(x)
      if ("pct_5_95" %in% stats) result$pct_5_95 <- list(quantile(x, probs = c(0.05, 0.95)))
      if ("ci_mean_low" %in% stats || "ci_mean_up" %in% stats) {
        t_test <- t.test(x, conf.level = 0.95)
        result$ci_mean_low <- t_test$conf.int[1]
        result$ci_mean_up <- t_test$conf.int[2]
      }
      if ("shapiro_p" %in% stats || "shapiro_int" %in% stats) {
        if (length(x) >= 3 && length(x) <= 5000) {
          shapiro_test <- shapiro.test(x)
          result$shapiro_p <- shapiro_test$p.value
          result$shapiro_int <- ifelse(shapiro_test$p.value > 0.05,
            "normal (p > 0.05)",
            "non-normal (p <= 0.05)"
          )
        } else {
          result$shapiro_p <- NA
          result$shapiro_int <- "sample size out of Shapiro-Wilk range (3-5000)"
        }
      }
      if ("n_outliers" %in% stats || "outlier_values" %in% stats) {
        lower_fence <- result$q1 - 1.5 * result$iqr
        upper_fence <- result$q3 + 1.5 * result$iqr
        outliers <- x[x < lower_fence | x > upper_fence]
        result$n_outliers <- length(outliers)
        result$outlier_values <- list(outliers)
      }
    } else {
      # Categorical variable statistics
      x <- as.factor(x) # Convert to factor
      x <- x[!is.na(x)] # Remove NAs
      freq_table <- table(x)
      n_valid <- sum(freq_table)

      # Check for low counts
      low_counts <- all(freq_table < 5) && n_valid < 10
      if (low_counts && any(c("chi_p", "pairwise_p") %in% stats)) {
        result$note <- "Warning: Very low counts; chi-square and pairwise tests may be unreliable"
      }

      # Compute n_unique early for evenness
      if (any(c("n_unique", "evenness", "chi_p", "chi_int", "chi_note", "cramer_v", "pairwise_p") %in% stats)) {
        result$n_unique <- length(unique(x))
      }

      if ("levels" %in% stats) result$levels <- list(levels(x))
      if ("mode" %in% stats) result$mode <- names(sort(table(x), decreasing = TRUE))[1]

      if ("level_summary" %in% stats) {
        proportions <- prop.table(freq_table) * 100
        ci_list <- lapply(names(freq_table), function(level) {
          binom_result <- binom.test(freq_table[level], n_valid, conf.level = 0.95)
          tibble(
            level = level,
            count = as.numeric(freq_table[level]),
            prop_pct = proportions[level],
            ci_low_pct = binom_result$conf.int[1] * 100,
            ci_up_pct = binom_result$conf.int[2] * 100
          )
        })
        result$level_summary <- list(bind_rows(ci_list))
      }

      if ("entropy" %in% stats || "evenness" %in% stats) {
        freqs <- prop.table(freq_table)
        result$entropy <- entropy(freqs, unit = "log2")
      }

      if ("evenness" %in% stats) {
        result$evenness <- if (!is.null(result$n_unique) && result$n_unique > 1) {
          result$entropy / log2(result$n_unique)
        } else {
          0
        }
      }

      if (any(c("chi_p", "chi_int", "chi_note", "cramer_v") %in% stats) && result$n_unique > 1) {
        expected_prob <- if (!is.null(chi_probs)) {
          if (length(chi_probs) != result$n_unique || sum(chi_probs) != 1) {
            stop("chi_probs must match number of levels and sum to 1")
          }
          chi_probs
        } else {
          rep(1 / result$n_unique, result$n_unique)
        }
        if (low_counts) {
          result$chi_p <- NA
          result$chi_int <- "chi-square test skipped due to very low counts"
          result$chi_note <- "no test performed"
          result$cramer_v <- NA
        } else if (any(freq_table < 5) || n_valid < 10) {
          chi_test <- chisq.test(freq_table, p = expected_prob, simulate.p.value = TRUE, B = 2000)
          result$chi_note <- "simulated p-value used due to low counts or small sample size"
          result$chi_p <- chi_test$p.value
          result$chi_int <- ifelse(chi_test$p.value > 0.05,
            "uniform distribution (p > 0.05)",
            "non-uniform distribution (p <= 0.05)"
          )
          result$cramer_v <- sqrt(chi_test$statistic / (n_valid * (result$n_unique - 1)))
        } else {
          chi_test <- chisq.test(freq_table, p = expected_prob)
          result$chi_note <- "standard chi-square test"
          result$chi_p <- chi_test$p.value
          result$chi_int <- ifelse(chi_test$p.value > 0.05,
            "uniform distribution (p > 0.05)",
            "non-uniform distribution (p <= 0.05)"
          )
          result$cramer_v <- sqrt(chi_test$statistic / (n_valid * (result$n_unique - 1)))
        }
      } else if (any(c("chi_p", "chi_int", "chi_note", "cramer_v") %in% stats)) {
        result$chi_p <- NA
        result$chi_int <- "chi-square test not applicable (single level)"
        result$chi_note <- "no test performed"
        result$cramer_v <- NA
      }

      if ("pairwise_p" %in% stats && pairwise && result$n_unique > 2 && all(freq_table >= 1) && result$n_unique <= 10) {
        level_names <- names(freq_table)
        pairwise_results <- tibble(level1 = character(), level2 = character(), p_value = numeric())
        for (i in 1:(length(level_names) - 1)) {
          for (j in (i + 1):length(level_names)) {
            test_result <- binom.test(c(freq_table[i], freq_table[j]),
              n = c(freq_table[i] + freq_table[j], freq_table[i] + freq_table[j]),
              p = 0.5
            )
            pairwise_results <- add_row(pairwise_results,
              level1 = level_names[i],
              level2 = level_names[j],
              p_value = test_result$p.value
            )
          }
        }
        # Apply Holm-Bonferroni correction
        pairwise_results$p_value <- p.adjust(pairwise_results$p_value, method = "holm")
        result$pairwise_p <- list(pairwise_results)
        result$pairwise_note <- "Holm-adjusted p-values from exact binomial tests"
      } else if ("pairwise_p" %in% stats) {
        reason <- if (!pairwise) {
          "pairwise comparisons disabled"
        } else if (result$n_unique <= 2) {
          "2 or fewer levels"
        } else if (any(freq_table < 1)) {
          "zero counts in some levels"
        } else {
          "too many levels (>10)"
        }
        result$pairwise_p <- list("not performed")
        result$pairwise_note <- paste("pairwise comparisons skipped:", reason)
      }
    }

    return(result)
  }

  # Handle grouped data
  if (dplyr::is_grouped_df(data)) {
    result <- data %>%
      group_by(across(all_of(group_vars(.)))) %>%
      summarise(stats = list(compute_stats({{ col }}, var_name)), .groups = "keep") %>%
      unnest_wider(stats)
  } else {
    result <- as_tibble(compute_stats(x, var_name))
  }

  return(result)
}

# Function for time series

kk_time_series <- function(data, value_col = NULL, date_col = NULL, group_cols = NULL, date_format = "auto", date_pattern = NULL) {
  # Load required packages
  required_pkgs <- c("dplyr", "moments", "tseries", "pracma", "zoo", "rugarch", "forecast", "entropy", "fractal", "lubridate")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    warning(sprintf("Missing packages: %s", paste(missing_pkgs, collapse = ", ")))
  }

  # Validate data input
  if (missing(data)) {
    stop("Argument 'data' is missing, with no default")
  }

  # Validate value column
  if (is.null(value_col) && !"value" %in% names(data)) {
    stop("Please specify a value column using value_col parameter or ensure 'value' column exists")
  }

  # Convert value_col to symbol and extract
  if (!is.null(value_col)) {
    if (is.character(value_col)) value_col <- rlang::sym(value_col)
    data <- dplyr::mutate(data, value = !!value_col)
  }

  # Auto-detect date column if not specified
  if (is.null(date_col) && date_format == "auto") {
    potential_date_cols <- c("date", "year_week", "time", "year_month", "year_quarter")
    date_col <- potential_date_cols[potential_date_cols %in% names(data)][1]
    if (!is.na(date_col)) {
      date_col <- rlang::sym(date_col)
      # Infer date_format based on column name or content
      if (date_col == "date") {
        date_format <- "date"
      } else if (date_col == "year_week") {
        date_format <- "year_week_string"
      } else if (date_col == "year_month") {
        date_format <- "year_month_string"
      } else if (date_col == "year_quarter") {
        date_format <- "year_quarter_string"
      } else {
        date_format <- "date"
      }
    }
  }

  # Handle date formats
  if (date_format == "date" && !is.null(date_col)) {
    if (is.character(date_col)) date_col <- rlang::sym(date_col)
    data <- dplyr::mutate(data, date = !!date_col)
  } else if (date_format == "year_month" && all(c("year", "month") %in% names(data))) {
    data <- dplyr::mutate(data, date = lubridate::ymd(paste(year, month, "01", sep = "-")))
  } else if (date_format == "year_week" && all(c("year", "week") %in% names(data))) {
    data <- dplyr::mutate(data, date = lubridate::make_date(year) + lubridate::weeks(week - 1))
  } else if (date_format == "year_quarter" && all(c("year", "quarter") %in% names(data))) {
    data <- dplyr::mutate(data, date = lubridate::yq(paste(year, quarter, sep = "-Q")))
  } else if (date_format == "year_week_string" && !is.null(date_col)) {
    # Handle year_week as a single string column (e.g., "2020-1", "2020-01", "2020-W1")
    if (is.character(date_col)) date_col <- rlang::sym(date_col)
    data <- dplyr::mutate(data, date = {
      # Normalize and split year_week string
      year_week <- gsub("W|-", "-", as.character(!!date_col)) # Replace 'W' or '-' with '-'
      parts <- strsplit(year_week, "-")
      year <- as.integer(sapply(parts, `[`, 1))
      week <- as.integer(gsub("^0+", "", sapply(parts, `[`, 2))) # Remove leading zeros
      lubridate::make_date(year) + lubridate::weeks(week - 1)
    })
  } else if (date_format == "year_month_string" && !is.null(date_col)) {
    # Handle year_month as a single string column (e.g., "2020-01")
    if (is.character(date_col)) date_col <- rlang::sym(date_col)
    data <- dplyr::mutate(data, date = lubridate::ym(as.character(!!date_col)))
  } else if (date_format == "year_quarter_string" && !is.null(date_col)) {
    # Handle year_quarter as a single string column (e.g., "2020-Q1")
    if (is.character(date_col)) date_col <- rlang::sym(date_col)
    data <- dplyr::mutate(data, date = lubridate::yq(as.character(!!date_col)))
  } else if (date_format == "custom" && !is.null(date_col) && !is.null(date_pattern)) {
    # Handle custom date formats with user-specified pattern
    if (is.character(date_col)) date_col <- rlang::sym(date_col)
    data <- dplyr::mutate(data, date = lubridate::parse_date_time(!!date_col, orders = date_pattern))
  } else if ("date" %in% names(data)) {
    # Fallback to existing date column
  } else {
    stop("Unable to identify a valid date column. Specify date_col and date_format (auto, date, year_month, year_week, year_quarter, year_week_string, year_month_string, year_quarter_string, custom) or ensure 'date' column exists. Available columns: ", paste(names(data), collapse = ", "))
  }

  # Ensure date is proper format
  if (!inherits(data$date, c("Date", "POSIXct"))) {
    data$date <- try(as.Date(data$date), silent = TRUE)
    if (inherits(data$date, "try-error") || !inherits(data$date, c("Date", "POSIXct"))) {
      stop("Unable to convert 'date' to Date or POSIXct. Check date_format, date_col, or date_pattern. Sample values: ", paste(head(data[[as.character(date_col)]], 3), collapse = ", "))
    }
  }

  # Handle grouping
  if (!is.null(group_cols)) {
    group_cols <- rlang::syms(group_cols)
    data <- dplyr::group_by(data, !!!group_cols)
  }

  # Process time series for each group
  result <- data %>%
    dplyr::arrange(date) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::group_map(~ {
      # Count missing values
      missing_count <- sum(is.na(.x$value))

      # Create zoo object
      ts_zoo <- try(zoo::zoo(.x$value, order.by = .x$date), silent = TRUE)
      if (inherits(ts_zoo, "try-error")) stop("Failed to create zoo object")

      # Calculate statistics
      n <- length(.x$value)
      if (n == 0) stop("No valid data points after filtering")

      # Rate of change
      roc <- try({
        vals <- .x$value
        (vals[-1] - vals[-n]) / vals[-n] * 100
      }, silent = TRUE)
      roc_stats <- if (inherits(roc, "try-error")) rep(NA, 4) else c(
        mean(roc, na.rm = TRUE),
        stats::sd(roc, na.rm = TRUE),
        min(roc, na.rm = TRUE),
        max(roc, na.rm = TRUE)
      )

      # Time index for regression
      time_index <- as.numeric(difftime(.x$date, min(.x$date), units = "days"))

      # Trend strength
      trend_model <- try(stats::lm(value ~ time_index, data = .x), silent = TRUE)
      trend_slope <- trend_strength <- NA
      if (!inherits(trend_model, "try-error")) {
        trend_slope <- try(stats::coef(trend_model)[2], silent = TRUE)
        if (!inherits(trend_slope, "try-error")) {
          trend_fit <- try(stats::predict(trend_model), silent = TRUE)
          if (!inherits(trend_fit, "try-error")) {
            detrended <- .x$value - trend_fit
            trend_strength <- try(1 - stats::var(detrended) / stats::var(.x$value), silent = TRUE)
            if (inherits(trend_strength, "try-error")) trend_strength <- NA
          }
        }
      }

      # GARCH volatility
      garch_vol <- try({
        if (requireNamespace("rugarch", quietly = TRUE) && n >= 30) {
          suppressWarnings({
            spec <- rugarch::ugarchspec(
              mean.model = list(armaOrder = c(1, 0)),
              variance.model = list(model = "sGARCH", garchOrder = c(1, 1))
            )
            garch_fit <- rugarch::ugarchfit(spec, .x$value, solver = "hybrid")
            mean(rugarch::sigma(garch_fit), na.rm = TRUE)
          })
        } else {
          NA
        }
      }, silent = TRUE)
      if (inherits(garch_vol, "try-error")) garch_vol <- NA

      # FFT for dominant frequency
      dom_freq <- try({
        freq <- stats::spectrum(.x$value, plot = FALSE)
        freq$freq[which.max(freq$spec)]
      }, silent = TRUE)
      if (inherits(dom_freq, "try-error")) dom_freq <- NA

      # Entropy
      shannon_entropy <- try({
        if (requireNamespace("entropy", quietly = TRUE)) {
          breaks <- pretty(range(.x$value), n = min(10, n / 5))
          binned <- cut(.x$value, breaks, include.lowest = TRUE)
          counts <- table(binned)
          entropy::entropy(counts / sum(counts))
        } else {
          NA
        }
      }, silent = TRUE)
      if (inherits(shannon_entropy, "try-error")) shannon_entropy <- NA

      # ACF/PACF
      acf_lag1 <- try({
        acf_result <- stats::acf(.x$value, lag.max = 1, plot = FALSE)
        acf_result$acf[2]
      }, silent = TRUE)
      if (inherits(acf_lag1, "try-error")) acf_lag1 <- NA

      pacf_lag1 <- try({
        pacf_result <- stats::pacf(.x$value, lag.max = 1, plot = FALSE)
        pacf_result$acf[1]
      }, silent = TRUE)
      if (inherits(pacf_lag1, "try-error")) pacf_lag1 <- NA

      # Ljung-Box
      ljung_pval <- try({
        test <- stats::Box.test(.x$value, lag = min(10, n - 1), type = "Ljung-Box")
        test$p.value
      }, silent = TRUE)
      if (inherits(ljung_pval, "try-error")) ljung_pval <- NA

      # Stationarity
      adf_result <- try({
        if (requireNamespace("tseries", quietly = TRUE) && n >= 10) {
          suppressWarnings({
            test <- tseries::adf.test(.x$value)
            c(test$p.value, test$statistic)
          })
        } else {
          c(NA, NA)
        }
      }, silent = TRUE)
      if (inherits(adf_result, "try-error")) adf_result <- c(NA, NA)

      kpss_result <- try({
        if (requireNamespace("tseries", quietly = TRUE) && n >= 10) {
          suppressWarnings({
            test <- tseries::kpss.test(.x$value)
            c(test$p.value, test$statistic)
          })
        } else {
          c(NA, NA)
        }
      }, silent = TRUE)
      if (inherits(kpss_result, "try-error")) kpss_result <- c(NA, NA)

      # Hurst
      hurst <- try({
        if (requireNamespace("pracma", quietly = TRUE) && n >= 20) {
          temp <- capture.output({
            h_result <- pracma::hurstexp(.x$value, display = FALSE)
          })
          h_result$Hs
        } else {
          NA
        }
      }, silent = TRUE)
      if (inherits(hurst, "try-error")) hurst <- NA

      # Lyapunov
      lyap <- try({
        if (requireNamespace("fractal", quietly = TRUE) && n >= 30) {
          fractal::lyapunov(.x$value)
        } else {
          NA
        }
      }, silent = TRUE)
      if (inherits(lyap, "try-error")) lyap <- NA

      # Outliers
      outlier_count <- try({
        q1 <- stats::quantile(.x$value, 0.25)
        q3 <- stats::quantile(.x$value, 0.75)
        iqr <- q3 - q1
        sum(.x$value < (q1 - 1.5 * iqr) | .x$value > (q3 + 1.5 * iqr))
      }, silent = TRUE)
      if (inherits(outlier_count, "try-error")) outlier_count <- NA

      # Compile results
      dplyr::tibble(
        Metric = c(
          "Length of Series", "Mean", "Median", "Standard Deviation", "Variance", "Min", "Max", "Range",
          "Q1 (25th Percentile)", "Q3 (75th Percentile)", "Skewness", "Kurtosis", "CV (SD/Mean)",
          "Missing Count", "Outlier Count", "Mean Rate of Change (%)", "SD Rate of Change (%)",
          "Min ROC (%)", "Max ROC (%)", "Trend Slope", "Trend Strength", "ADF p-value", "ADF statistic",
          "KPSS p-value", "KPSS statistic", "ACF Lag 1", "PACF Lag 1", "Ljung-Box p-value",
          "Hurst Exponent", "GARCH Volatility", "Shannon Entropy", "Dominant Frequency", "Lyapunov Exponent"
        ),
        Value = tryCatch({
          c(
            n,
            mean(.x$value, na.rm = TRUE),
            stats::median(.x$value, na.rm = TRUE),
            stats::sd(.x$value, na.rm = TRUE),
            stats::var(.x$value, na.rm = TRUE),
            min(.x$value, na.rm = TRUE),
            max(.x$value, na.rm = TRUE),
            diff(range(.x$value, na.rm = TRUE)),
            stats::quantile(.x$value, 0.25, na.rm = TRUE),
            stats::quantile(.x$value, 0.75, na.rm = TRUE),
            if (requireNamespace("moments", quietly = TRUE)) moments::skewness(.x$value, na.rm = TRUE) else NA,
            if (requireNamespace("moments", quietly = TRUE)) moments::kurtosis(.x$value, na.rm = TRUE) else NA,
            if (mean(.x$value, na.rm = TRUE) != 0) stats::sd(.x$value, na.rm = TRUE) / abs(mean(.x$value, na.rm = TRUE)) else NA,
            missing_count,
            outlier_count,
            roc_stats[1], roc_stats[2], roc_stats[3], roc_stats[4],
            trend_slope,
            trend_strength,
            adf_result[1], adf_result[2],
            kpss_result[1], kpss_result[2],
            acf_lag1,
            pacf_lag1,
            ljung_pval,
            hurst,
            garch_vol,
            shannon_entropy,
            dom_freq,
            lyap
          )
        }, error = function(e) {
          warning("Error in calculating statistics: ", e$message)
          rep(NA, 33)
        })
      )
    }, .keep = TRUE) %>%
    dplyr::bind_rows(.id = "group") %>%
    dplyr::mutate(Value = round(as.numeric(Value), 4))

  return(result)
}