###### ---Libraries---######

if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, haven, modelsummary, MKinfer, rstatix, finalfit, tinytable, monochromeR, ggstats, epitools, ggsurvfit, broom, rstan, brms, gtsummary, quantreg, patchwork, tidymodels, gt, epiR, readxl, scales, marginaleffects, ggthemes, emmeans, janitor, easystats, showtext, brglm2, sysfonts, MASS, detectseparation)

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


###### ---Define theme---######
set_plot_font <- function(font = "Roboto Condensed", size = 18) {
  showtext::showtext_auto()

  # Try to add the Google font dynamically
  sysfonts::font_add_google(font, font)

  # Define relative font sizes based on the `size` parameter
  title_size <- size + 4
  subtitle_size <- size + 2
  caption_size <- size + 2
  axis_title_size <- size
  axis_text_size <- size
  strip_text_size <- size

  theme_nice <- ggthemes::theme_tufte() +
    theme(
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.ticks.length = unit(4, "mm"),
      plot.title = element_text(family = font, size = title_size, hjust = 0, vjust = 2),
      plot.subtitle = element_text(family = font, size = subtitle_size),
      plot.caption = element_text(family = font, size = caption_size, hjust = 1),
      axis.title = element_text(family = font, size = axis_title_size),
      axis.text = element_text(family = font, size = axis_text_size),
      axis.text.x = element_text(margin = margin(5, b = 10)),
      strip.text = element_text(family = font, size = strip_text_size),
      axis.line = element_line()
    )

  theme_set(theme_nice)
}
