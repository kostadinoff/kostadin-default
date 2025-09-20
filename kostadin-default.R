###### ---Libraries---######

if (!require(pacman)) install.packages("pacman")

# Core packages that are used throughout
pacman::p_load(tidyverse, rlang, broom)

# Helper function for conditional package loading
require_packages <- function(packages, context = "this analysis") {
  missing_pkgs <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(paste("The following packages are required for", context, "but not installed:", 
               paste(missing_pkgs, collapse = ", "), 
               "\nPlease install with: install.packages(c('", 
               paste(missing_pkgs, collapse = "', '"), "'))"))
  }
  invisible(TRUE)
}

# Load additional packages as needed within functions rather than globally

###### ---Options---######

options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())
options(ggplot2.messages = FALSE)
options(dplyr.width = Inf)
###### ---Functions---######

# Function to format a tibble with numeric values for display
format_tibble <- function(data, digits = 2) {
  data %>%
    mutate(
      Value_display = sapply(Value, function(x) {
        if (is.na(x)) "NA" else format(round(x, digits), nsmall = digits, scientific = FALSE, big.mark = ",")
      })
    )
}

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

# Needs: emmeans, sandwich, stats, dplyr, tibble, rlang
compare_proportions_glm <- function(
  data,
  group,                # grouping factor to compare
  x, n,                 # counts of "successes" and totals
  by = NULL,            # optional stratifier (within-stratum pairwise)
  covariates = NULL,    # optional vector of covariate column names for adjustment
  adjust = "holm",      # p-adjust method for multiple pairwise tests
  conf.level = 0.95,
  vcov_type = "HC3",    # robust SE type for sandwich::vcovHC
  drop_empty = TRUE     # drop strata with <2 groups
) {
  # Dependency checks
  for (pkg in c("emmeans", "sandwich")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required. Install it first.", pkg))
    }
  }

  group <- dplyr::ensym(group)
  x     <- dplyr::ensym(x)
  n     <- dplyr::ensym(n)
  by_   <- if (is.null(by)) NULL else dplyr::ensym(by)

  df <- data |>
    dplyr::mutate(
      .x = !!x,
      .n = !!n,
      .fail = .n - .x
    )

  if (any(df$.n < df$.x, na.rm = TRUE)) {
    stop("Each 'n' must be >= 'x'.")
  }

  # Build formula: cbind(x, n-x) ~ group (+ by) (+ covariates) (+ interactions if by)
  rhs_terms <- rlang::as_string(group)
  if (!is.null(by_)) rhs_terms <- paste(rhs_terms, rlang::as_string(by_), sep = " * ") # interaction allows within-stratum contrasts
  if (!is.null(covariates) && length(covariates)) {
    rhs_terms <- paste(rhs_terms, paste(covariates, collapse = " + "), sep = " + ")
  }
  fml <- stats::as.formula(paste0("cbind(.x, .fail) ~ ", rhs_terms))

  fit <- stats::glm(
    formula = fml,
    family  = stats::binomial(link = "logit"),
    data    = df
  )

  # Robust VCOV
  V <- sandwich::vcovHC(fit, type = vcov_type)

  # emmeans grid
  if (is.null(by_)) {
    emm <- emmeans::emmeans(fit, specs = ~ !!group, vcov. = V)
    # Regrid to response to get probabilities; then pairwise diffs on prob scale
    emm_resp <- emmeans::regrid(emm, transform = "response")
    cmp <- emmeans::contrast(emm_resp, method = "pairwise", adjust = adjust)
    out <- as.data.frame(emmeans::summary(cmp, infer = TRUE, level = conf.level))
    # Standardize column names
    out <- tibble::as_tibble(out) |>
      dplyr::transmute(
        group1   = contrast |> sub(" - .*", "", .),
        group2   = contrast |> sub(".* - ", "", .),
        estimate = estimate,      # difference in probabilities p1 - p2
        conf_low = lower.CL,
        conf_high= upper.CL,
        p_value  = p.value,
        adjust   = adjust,
        conf_level = conf.level
      )
  } else {
    by_name <- rlang::as_string(by_)
    emm <- emmeans::emmeans(fit, specs = list(group = ~ !!group, by = ~ !!by_), vcov. = V)
    emm_resp <- emmeans::regrid(emm, transform = "response")
    cmp <- emmeans::contrast(emm_resp, method = "pairwise", by = by_name, adjust = adjust)
    tmp <- as.data.frame(emmeans::summary(cmp, infer = TRUE, level = conf.level))
    out <- tibble::as_tibble(tmp) |>
      dplyr::transmute(
        !!by_name := .data[[by_name]],
        group1   = contrast |> sub(" - .*", "", .),
        group2   = contrast |> sub(".* - ", "", .),
        estimate = estimate,
        conf_low = lower.CL,
        conf_high= upper.CL,
        p_value  = p.value,
        adjust   = adjust,
        conf_level = conf.level
      )
    if (drop_empty) {
      out <- out |> dplyr::group_by(.data[[by_name]]) |>
        dplyr::filter(dplyr::n() > 0) |>
        dplyr::ungroup()
    }
  }

  out
}


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

    # Handle edge cases more appropriately
    if (pooled_se < .Machine$double.eps) {
      warning("Standard error is effectively zero; p-value and confidence intervals may be unreliable.")
      return(tibble(
        message = "Cannot compute reliable statistics with zero variance",
        prop_diff = prop_diff,
        z_score = NA_real_,
        p_value = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_
      ))
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

# Function for EGN
extract_age_from_egn <- function(egn, admission_date = Sys.Date()) {
  region_codes <- list(
    "000-043" = "Blagoevgrad", "044-093" = "Burgas", "094-139" = "Varna",
    "140-169" = "Veliko Tarnovo", "170-183" = "Vidin", "184-217" = "Vratsa",
    "218-233" = "Gabrovo", "234-281" = "Kardzhali", "282-301" = "Kyustendil",
    "302-319" = "Lovech", "320-341" = "Montana", "342-377" = "Pazardzhik",
    "378-395" = "Pernik", "396-435" = "Pleven", "436-501" = "Plovdiv",
    "502-527" = "Razgrad", "528-555" = "Ruse", "556-575" = "Silistra",
    "576-601" = "Sliven", "602-623" = "Smolyan", "624-721" = "Sofia City",
    "722-751" = "Sofia Province", "752-789" = "Stara Zagora", "790-821" = "Dobrich",
    "822-843" = "Targovishte", "844-871" = "Haskovo", "872-903" = "Shumen",
    "904-925" = "Yambol", "926-999" = "Other/Unknown"
  )
  weights <- c(2, 4, 8, 5, 10, 9, 7, 3, 6)

  egn <- as.character(egn)
  admission_date <- as.Date(admission_date)

  result <- data.frame(
    age = numeric(length(egn)),
    birth_date = as.Date(character(length(egn))),
    is_valid = logical(length(egn)),
    gender = character(length(egn)),
    region = character(length(egn)),
    birth_order = numeric(length(egn)),
    invalid_egn = character(length(egn)),
    invalid_reason = character(length(egn)), # New column for reason
    stringsAsFactors = FALSE
  )

  for (i in seq_along(egn)) {
    output <- list(
      age = NA_real_,
      birth_date = as.Date(NA_character_),
      is_valid = FALSE,
      gender = NA_character_,
      region = NA_character_,
      birth_order = NA_real_,
      invalid_egn = NA_character_,
      invalid_reason = NA_character_
    )

    if (is.na(egn[i]) || nchar(trimws(egn[i])) == 0) {
      output$invalid_egn <- egn[i]
      output$invalid_reason <- "Missing or empty EGN"
      result[i, ] <- output
      next
    }

    egn_length <- nchar(egn[i])

    # Helper function for date parsing (defined once)
    parse_egn_date <- function(egn_str, admission_date_val) {
      year <- as.numeric(substr(egn_str, 1, 2))
      month <- as.numeric(substr(egn_str, 3, 4))
      day <- as.numeric(substr(egn_str, 5, 6))
      
      # Determine century based on month indicator
      if (month >= 41 && month <= 52) {
        year <- year + 2000
        month <- month - 40
      } else if (month >= 21 && month <= 32) {
        year <- year + 1800
        month <- month - 20
      } else if (month >= 1 && month <= 12) {
        year <- year + 1900
      } else {
        return(list(valid = FALSE, reason = sprintf("Invalid month: %d", month)))
      }
      
      # Try to create valid date
      birth_date <- tryCatch(
        as.Date(paste(year, month, day, sep = "-")),
        error = function(e) {
          tryCatch(
            as.Date(paste(year, month, "01", sep = "-")),
            error = function(e2) as.Date(NA_character_)
          )
        }
      )
      
      if (is.na(birth_date)) {
        return(list(valid = FALSE, reason = "Invalid birth date"))
      }
      
      age <- floor(as.numeric(difftime(admission_date_val, birth_date, units = "days")) / 365.25)
      return(list(valid = TRUE, birth_date = birth_date, age = age))
    }
    
    if (egn_length == 6) {
      current_admission_date <- if (length(admission_date) == 1) admission_date else admission_date[i]
      date_result <- parse_egn_date(egn[i], current_admission_date)
      
      if (!date_result$valid) {
        output$invalid_egn <- egn[i]
        output$invalid_reason <- date_result$reason
      } else {
        output$birth_date <- date_result$birth_date
        output$age <- date_result$age
        output$is_valid <- TRUE
      }
      
      result[i, ] <- output
      next
    }

    if (egn_length == 10 && grepl("^\\d+$", egn[i])) {
      year <- as.numeric(substr(egn[i], 1, 2))
      month <- as.numeric(substr(egn[i], 3, 4))
      day <- as.numeric(substr(egn[i], 5, 6))

      month_indicator <- month
      if (month_indicator >= 41 && month_indicator <= 52) {
        year <- year + 2000
        month <- month - 40
      } else if (month_indicator >= 21 && month_indicator <= 32) {
        year <- year + 1800
        month <- month - 20
      } else if (month_indicator >= 1 && month_indicator <= 12) {
        year <- year + 1900
      } else {
        warning(sprintf("Invalid month in EGN at index %d: %s", i, egn[i]))
        output$invalid_egn <- egn[i]
        output$invalid_reason <- sprintf("Invalid month: %d", month_indicator)
        result[i, ] <- output
        next
      }

      output$birth_date <- tryCatch(
        as.Date(paste(year, month, day, sep = "-"), format = "%Y-%m-%d"),
        error = function(e) as.Date(NA_character_)
      )

      if (is.na(output$birth_date)) {
        output$invalid_egn <- egn[i]
        output$invalid_reason <- "Invalid birth date"
        result[i, ] <- output
        next
      }

      current_admission_date <- if (length(admission_date) == 1) admission_date else admission_date[i]
      output$age <- floor(as.numeric(difftime(current_admission_date, output$birth_date, units = "days")) / 365.25)

      digits <- as.numeric(strsplit(egn[i], "")[[1]])
      weighted_sum <- sum(digits[1:9] * weights)
      control_digit <- weighted_sum %% 11
      control_digit <- if (control_digit == 10) 0 else control_digit
      output$is_valid <- control_digit == digits[10]

      if (!output$is_valid) {
        output$invalid_reason <- "Invalid control digit"
      }

      ninth_digit <- as.numeric(substr(egn[i], 9, 9))
      output$gender <- if (ninth_digit %% 2 == 0) "Male" else "Female"

      three_digit_code <- as.numeric(substr(egn[i], 7, 9))
      if (ninth_digit %% 2 == 0) {
        output$birth_order <- (ninth_digit / 2) + 1
      } else {
        output$birth_order <- ((ninth_digit + 1) / 2)
      }

      output$region <- "Other/Unknown"
      for (range in names(region_codes)) {
        range_bounds <- as.numeric(unlist(strsplit(range, "-")))
        if (three_digit_code >= range_bounds[1] && three_digit_code <= range_bounds[2]) {
          output$region <- region_codes[[range]]
          break
        }
      }
    } else {
      output$invalid_egn <- egn[i]
      output$invalid_reason <- sprintf("Invalid length (%d) or non-numeric", egn_length)
    }

    result[i, ] <- output
  }

  return(result)
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


# Plot

###### ---Define theme---######

set_plot_font <- function(font = "Roboto Condensed", size = 18) {
  # Require necessary packages
  require_packages(c("sysfonts", "ggthemes"), "font management")
  
  # Determine which font to use
  font_to_use <- tryCatch({
    sysfonts::font_add_google(name = font, family = font, db_cache = FALSE)
    message("Successfully loaded font: ", font)
    font
  }, error = function(e) {
    message("Font '", font, "' not found on Google Fonts. Falling back to 'Arial'.")
    "Arial"  # Local assignment, not global
  })

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
      plot.title = element_text(family = font_to_use, size = title_size, hjust = 0, vjust = 2, margin = margin(t = 10, b = 10)),
      plot.subtitle = element_text(family = font_to_use, size = subtitle_size),
      plot.caption = element_text(family = font_to_use, hjust = 0.5, vjust = 1, size = caption_size),
      plot.caption.position = "plot",
      axis.title = element_text(family = font_to_use, size = axis_title_size),
      axis.text = element_text(family = font_to_use, size = axis_text_size),
      axis.text.x = element_text(margin = margin(5, b = 10)),
      strip.text = element_text(family = font_to_use, size = strip_text_size),
      axis.line = element_line()
    )

  theme_set(theme_nice)
}


kkplot <- function(...) {
  ggplot(...) +
    guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both"))
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
    kkplot(aes(y = fct_reorder({{ variable }}, prop), x = prop)) +
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
    kkplot(aes(x = !!variable)) +
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
      y = 0.8 * max(density(data[[rlang::as_name(variable)]], na.rm = TRUE)$y),
      label = paste0("Median: ", round(median(data[[rlang::as_name(variable)]], na.rm = TRUE), 2)),
      color = "blue",
      size = 5,
      family = "Roboto Condensed",
      hjust = -0.1a[[rlang::as_name(variable)]], na.rm = TRUE),
    ) + = 0.9 * max(density(data[[rlang::as_name(variable)]], na.rm = TRUE)$y),
    expand_limits(x = c(min(data[[rlang::as_name(variable)]], na.rm = TRUE), max(data[[rlang::as_name(variable)]], na.rm = TRUE))) +
    labs(or = "red",
      x = "Value",
      y = "Density",to Condensed",
      title = title
    ) +
}   annotate(
      "label",
# Univariate plot: Check variable type and apply the appropriate plot
univariate_plot <- function(data, variable) {ame(variable)]], na.rm = TRUE)$y),
  variable <- rlang::ensym(variable)nd(median(data[[rlang::as_name(variable)]], na.rm = TRUE), 2)),
      color = "blue",
  # Check if the variable is numeric or categorical
  if (is.numeric(data[[rlang::as_name(variable)]])) {
    # Call the continuous plot function
    univariate_cont_plot(data, !!variable)
  } else if (is.factor(data[[rlang::as_name(variable)]]) || is.character(data[[rlang::as_name(variable)]])) {e)]], na.rm = TRUE))) +
    # Call the categorical plot function
    univariate_cat_plot(data, !!variable)
  } else {"Density",
    stop("The variable must be either numeric or categorical (factor/character).")
  } )
}

# Function for regressionvariable type and apply the appropriate plot
univariate_plot <- function(data, variable) {
kk_reg <- function(data, outcome, predictors, log_outcome = FALSE, custom_formula = NULL, ...) {
  # Validate inputs
  if (!is.data.frame(data)) {numeric or categorical
    stop("`data` must be a data frame.")riable)]])) {
  } # Call the continuous plot function
  if (!rlang::is_string(outcome)) {riable)
    stop("`outcome` must be a character string.")ble)]]) || is.character(data[[rlang::as_name(variable)]])) {
  } # Call the categorical plot function
  if (!is.character(predictors)) {riable)
    stop("`predictors` must be a character vector.")
  } stop("The variable must be either numeric or categorical (factor/character).")
  if (!is.logical(log_outcome)) {
    stop("`log_outcome` must be a logical value (TRUE or FALSE).")
  }
# Function for regression
  # Capture outcome as symbol
  outcome <- rlang::ensym(outcome)predictors, log_outcome = FALSE, custom_formula = NULL, ...) {
  # Validate inputs
  # Ensure outcome exists in the data
  if (!rlang::as_name(outcome) %in% colnames(data)) {
    stop("Outcome variable '", rlang::as_name(outcome), "' not found in the data.")
  }f (!rlang::is_string(outcome)) {
    stop("`outcome` must be a character string.")
  # Ensure predictors exist in the data
  missing_predictors <- setdiff(predictors, colnames(data))
  if (length(missing_predictors) > 0) {ter vector.")
    stop("Predictor variable(s) not found in the data: ", paste(missing_predictors, collapse = ", "))
  }f (!is.logical(log_outcome)) {
    stop("`log_outcome` must be a logical value (TRUE or FALSE).")
  # If log_outcome is TRUE, transform the outcome variable
  if (log_outcome) {
    data <- data %>% mutate(!!outcome := log(!!outcome))
  }utcome <- rlang::ensym(outcome)

  # Determine the type of outcome variable
  outcome_type <- case_when(e) %in% colnames(data)) {
    is.factor(data[[rlang::as_name(outcome)]]) && nlevels(data[[rlang::as_name(outcome)]]) == 2 ~ "binary", # Binary outcome
    is.ordered(data[[rlang::as_name(outcome)]]) ~ "ordinal", # Ordinal outcome
    is.numeric(data[[rlang::as_name(outcome)]]) ~ "continuous", # Continuous outcome
    TRUE ~ NA_character_ist in the data
  )issing_predictors <- setdiff(predictors, colnames(data))
  if (length(missing_predictors) > 0) {
  if (is.na(outcome_type)) {(s) not found in the data: ", paste(missing_predictors, collapse = ", "))
    stop("Outcome must be numeric, binary factor, or ordinal factor.")
  }
  # If log_outcome is TRUE, transform the outcome variable
  # Function to fit a model
  fit_model <- function(formula) {ome := log(!!outcome))
    if (outcome_type == "binary") {
      return(glm(formula, data = data, family = binomial(), ...)) # Logistic regression for binary outcome
    } else if (outcome_type == "ordinal") {
      return(MASS::polr(formula, data = data, Hess = TRUE, ...)) # Ordinal regression
    } else if (outcome_type == "continuous") { && nlevels(data[[rlang::as_name(outcome)]]) == 2 ~ "binary", # Binary outcome
      return(lm(formula, data = data, ...)) # Linear regression for continuous outcome
    }s.numeric(data[[rlang::as_name(outcome)]]) ~ "continuous", # Continuous outcome
  } TRUE ~ NA_character_
  )
  # Function to process model output
  process_results <- function(model, model_type) {
    results <- broom::tidy(model, conf.int = TRUE) %>%rdinal factor.")
      mutate(
        model_type = model_type,
        outcome_type = outcome_type,
        AIC = AIC(model),ormula) {
        BIC = BIC(model)"binary") {
      )eturn(glm(formula, data = data, family = binomial(), ...)) # Logistic regression for binary outcome
    } else if (outcome_type == "ordinal") {
    if (outcome_type == "binary") {ta = data, Hess = TRUE, ...)) # Ordinal regression
      # Convert to odds ratios (OR) for logistic regression
      results <- results %>%a = data, ...)) # Linear regression for continuous outcome
        mutate(
          estimate = exp(estimate), # Odds ratio for logistic regression
          conf.low = exp(conf.low), # Confidence interval for OR
          conf.high = exp(conf.high),
          percent_change = (estimate - 1) * 100, # Percentage change for OR
          percent_change_low = (conf.low - 1) * 100, # Lower CI for percentage change
          percent_change_high = (conf.high - 1) * 100 # Upper CI for percentage change
        )odel_type = model_type,
    } else if (outcome_type == "continuous" && log_outcome) {
      # Convert to percentage change if linear regression with log outcome
      results <- results %>%
        mutate(
          percent_change = (exp(estimate) - 1) * 100,
          percent_change_low = (exp(conf.low) - 1) * 100,
          percent_change_high = (exp(conf.high) - 1) * 100n
        )ults <- results %>%
    } else if (outcome_type == "continuous") {
      # Calculate percentage change for continuous (no log transformation)
      results <- results %>%f.low), # Confidence interval for OR
        mutate(high = exp(conf.high),
          percent_change = estimate * 100, # Simple percentage change for continuous variable
          percent_change_low = conf.low * 100,* 100, # Lower CI for percentage change
          percent_change_high = conf.high * 100 * 100 # Upper CI for percentage change
        )
    } else {f (outcome_type == "continuous" && log_outcome) {
      # No transformation for normal linear regressionion with log outcome
      results <- results %>%
        mutate(percent_change = NA, percent_change_low = NA, percent_change_high = NA)
    }     percent_change = (exp(estimate) - 1) * 100,
          percent_change_low = (exp(conf.low) - 1) * 100,
    return(results)hange_high = (exp(conf.high) - 1) * 100
  }     )
    } else if (outcome_type == "continuous") {
  # Fit univariate modelsage change for continuous (no log transformation)
  univariate_results <- map_dfr(predictors, function(predictor) {
    formula <- as.formula(paste(rlang::as_name(outcome), "~", predictor))
    model <- fit_model(formula)mate * 100, # Simple percentage change for continuous variable
    process_results(model, "univariate") %>%0,
      mutate(predictor = predictor)f.high * 100
  })    )
    } else {
  # Fit multivariate modelfor normal linear regression
  if (is.null(custom_formula)) {
    multivariate_formula <- as.formula(paste(rlang::as_name(outcome), "~", paste(predictors, collapse = " + ")))
  } else {
    multivariate_formula <- custom_formula
  } return(results)
  multivariate_model <- fit_model(multivariate_formula)
  multivariate_results <- process_results(multivariate_model, "multivariate")
  # Fit univariate models
  # Combine resultss <- map_dfr(predictors, function(predictor) {
  results <- bind_rows(univariate_results, multivariate_results)edictor))
    model <- fit_model(formula)
  return(results)ts(model, "univariate") %>%
}     mutate(predictor = predictor)
# Additional function for regression

krk_reg <- function(data, outcome, predictors, log_outcome = FALSE, custom_formula = NULL, include_diagnostics = TRUE, ...) {
  # Load required packagesla)) {
  library(dplyr)_formula <- as.formula(paste(rlang::as_name(outcome), "~", paste(predictors, collapse = " + ")))
  library(rlang)
  library(broom)_formula <- custom_formula
  library(purrr)
  multivariate_model <- fit_model(multivariate_formula)
  # Validate inputslts <- process_results(multivariate_model, "multivariate")
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!is_string(outcome)) stop("`outcome` must be a character string.")
  if (!is.character(predictors)) stop("`predictors` must be a character vector.")
  if (!is.logical(log_outcome)) stop("`log_outcome` must be a logical value (TRUE or FALSE).")
  if (!is.logical(include_diagnostics)) stop("`include_diagnostics` must be a logical value.")
}
  # Capture outcome as symbolression
  outcome <- ensym(outcome)
  outcome_name <- as_name(outcome) predictors, log_outcome = FALSE, custom_formula = NULL, include_diagnostics = TRUE, ...) {
  # Load required packages
  # Check if outcome exists
  if (!outcome_name %in% colnames(data)) {
    stop("Outcome variable '", outcome_name, "' not found in the data.")
  }ibrary(purrr)

  # Check if predictors exist
  missing_predictors <- setdiff(predictors, colnames(data))me.")
  if (length(missing_predictors) > 0) {me` must be a character string.")
    stop("Predictor variable(s) not found in the data: ", paste(missing_predictors, collapse = ", "))
  }f (!is.logical(log_outcome)) stop("`log_outcome` must be a logical value (TRUE or FALSE).")
  if (!is.logical(include_diagnostics)) stop("`include_diagnostics` must be a logical value.")
  # Transform outcome if log_outcome is TRUE
  if (log_outcome) {as symbol
    if (!is.numeric(data[[outcome_name]]) || any(data[[outcome_name]] <= 0, na.rm = TRUE)) {
      stop("`log_outcome = TRUE` requires a positive numeric outcome.")
    }
    data <- data %>% mutate(!!outcome := log(!!outcome))
  }f (!outcome_name %in% colnames(data)) {
    stop("Outcome variable '", outcome_name, "' not found in the data.")
  # Determine outcome type
  outcome_type <- case_when(
    is.factor(data[[outcome_name]]) && nlevels(data[[outcome_name]]) == 2 ~ "binary",
    is.ordered(data[[outcome_name]]) ~ "ordinal",mes(data))
    is.numeric(data[[outcome_name]]) ~ "continuous",
    TRUE ~ NA_character_able(s) not found in the data: ", paste(missing_predictors, collapse = ", "))
  )
  if (is.na(outcome_type)) stop("Outcome must be numeric, binary factor, or ordered factor.")
  # Transform outcome if log_outcome is TRUE
  # Function to fit a model
  fit_model <- function(formula) {name]]) || any(data[[outcome_name]] <= 0, na.rm = TRUE)) {
    tryCatch(log_outcome = TRUE` requires a positive numeric outcome.")
      {
        if (outcome_type == "binary") {= log(!!outcome))
          glm(formula, data = data, family = binomial(), ...)
        } else if (outcome_type == "ordinal") {
          if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required for ordinal regression.")
          MASS::polr(formula, data = data, Hess = TRUE, ...)
        } else if (outcome_type == "continuous") {a[[outcome_name]]) == 2 ~ "binary",
          lm(formula, data = data, ...)"ordinal",
        }meric(data[[outcome_name]]) ~ "continuous",
      }, ~ NA_character_
      error = function(e) {
        stop("Model fitting failed: ", e$message)numeric, binary factor, or ordered factor.")
      }
    )unction to fit a model
  }it_model <- function(formula) {
    tryCatch(
  # Function to calculate diagnostics (including model significance)
  calculate_diagnostics <- function(model) {
    if (!include_diagnostics) {ata, family = binomial(), ...)
      return(tibble())come_type == "ordinal") {
    }     if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' required for ordinal regression.")
          MASS::polr(formula, data = data, Hess = TRUE, ...)
    diagnostics <- list()e_type == "continuous") {
    summary_model <- summary(model)...)
        }
    # Add model significance test
    if (outcome_type == "continuous") {
      diagnostics$r_squared <- summary_model$r.squared
      diagnostics$adj_r_squared <- summary_model$adj.r.squared
      diagnostics$residual_std_error <- summary_model$sigma
      diagnostics$model_p_value <- pf(summary_model$fstatistic[1],
        summary_model$fstatistic[2],
        summary_model$fstatistic[3],s (including model significance)
        lower.tail = FALSE function(model) {
      )(!include_diagnostics) {
    } else if (outcome_type == "binary") {
      null_model <- update(model, ~1)
      loglik_model <- as.numeric(logLik(model))
      loglik_null <- as.numeric(logLik(null_model))
      diagnostics$pseudo_r_squared <- as.numeric(1 - (loglik_model / loglik_null))
      diagnostics$nagelkerke_r_squared <- as.numeric((1 - exp(-2 * (loglik_model - loglik_null))) /
        (1 - exp(2 * loglik_null / nrow(data))))
      diagnostics$model_p_value <- pchisq(2 * (loglik_model - loglik_null),
        df = length(coef(model)) - 1,y_model$r.squared
        lower.tail = FALSEuared <- summary_model$adj.r.squared
      )iagnostics$residual_std_error <- summary_model$sigma
      if (requireNamespace("pROC", quietly = TRUE)) {statistic[1],
        roc_curve <- pROC::roc(model$y, fitted(model), quiet = TRUE)
        diagnostics$auc_roc <- as.numeric(pROC::auc(roc_curve))
      } lower.tail = FALSE
    } else if (outcome_type == "ordinal") {
      null_model <- update(model, ~1)y") {
      loglik_model <- as.numeric(logLik(model))
      loglik_null <- as.numeric(logLik(null_model))
      diagnostics$pseudo_r_squared <- as.numeric(1 - (loglik_model / loglik_null))
      diagnostics$nagelkerke_r_squared <- as.numeric((1 - exp(-2 * (loglik_model - loglik_null))) /
        (1 - exp(2 * loglik_null / nrow(data))))eric((1 - exp(-2 * (loglik_model - loglik_null))) /
      diagnostics$model_p_value <- pchisq(2 * (loglik_model - loglik_null),
        df = length(coef(model)),- pchisq(2 * (loglik_model - loglik_null),
        lower.tail = FALSEodel)) - 1,
      ) lower.tail = FALSE
    } )
    as_tibble(diagnostics)("pROC", quietly = TRUE)) {
  }     roc_curve <- pROC::roc(model$y, fitted(model), quiet = TRUE)
        diagnostics$auc_roc <- as.numeric(pROC::auc(roc_curve))
  # Function to process model output
  process_results <- function(model, model_type) {
    estimate_label <- if (outcome_type == "binary") {
      "odds_ratio" <- as.numeric(logLik(model))
    } else if (outcome_type == "continuous" && log_outcome) {
      "exp_coef"s$pseudo_r_squared <- as.numeric(1 - (loglik_model / loglik_null))
    } else if (outcome_type == "ordinal") {s.numeric((1 - exp(-2 * (loglik_model - loglik_null))) /
      "odds_ratio" * loglik_null / nrow(data))))
    } else {stics$model_p_value <- pchisq(2 * (loglik_model - loglik_null),
      "coef" length(coef(model)),
    }   lower.tail = FALSE
      )
    results <- broom::tidy(model, conf.int = TRUE) %>%
      mutate((diagnostics)
        model_type = model_type,
        outcome_type = outcome_type,
        AIC = AIC(model),odel output
        BIC = BIC(model),tion(model, model_type) {
        estimate_label = estimate_label,= "binary") {
        coef.type = if_else(grepl("\\|", term), "scale", "coefficient")
      )lse if (outcome_type == "continuous" && log_outcome) {
      "exp_coef"
    if (outcome_type %in% c("binary", "ordinal") || (outcome_type == "continuous" && log_outcome)) {
      results <- results %>%
        mutate(
          estimate = if_else(coef.type == "coefficient", exp(estimate), estimate),
          conf.low = if_else(coef.type == "coefficient", exp(conf.low), conf.low),
          conf.high = if_else(coef.type == "coefficient", exp(conf.high), conf.high),
          percent_change = if_else(coef.type == "coefficient", (estimate - 1) * 100, NA_real_),
          percent_change_low = if_else(coef.type == "coefficient", (conf.low - 1) * 100, NA_real_),
          percent_change_high = if_else(coef.type == "coefficient", (conf.high - 1) * 100, NA_real_)
        )utcome_type = outcome_type,
    } else if (outcome_type == "continuous") {
      results <- results %>%
        mutate(e_label = estimate_label,
          percent_change = estimate * 100,erm), "scale", "coefficient")
          percent_change_low = conf.low * 100,
          percent_change_high = conf.high * 100
        )utcome_type %in% c("binary", "ordinal") || (outcome_type == "continuous" && log_outcome)) {
    } results <- results %>%
        mutate(
    diagnostics <- calculate_diagnostics(model)ficient", exp(estimate), estimate),
    results %>% bind_cols(diagnostics) == "coefficient", exp(conf.low), conf.low),
  }       conf.high = if_else(coef.type == "coefficient", exp(conf.high), conf.high),
          percent_change = if_else(coef.type == "coefficient", (estimate - 1) * 100, NA_real_),
  # Fit univariate modelslow = if_else(coef.type == "coefficient", (conf.low - 1) * 100, NA_real_),
  univariate_results <- map_dfr(predictors, function(predictor) {", (conf.high - 1) * 100, NA_real_)
    formula <- as.formula(paste(outcome_name, "~", predictor))
    model <- fit_model(formula)"continuous") {
    process_results(model, "univariate") %>% mutate(predictor = predictor)
  })    mutate(
          percent_change = estimate * 100,
  # Fit multivariate modelow = conf.low * 100,
  if (is.null(custom_formula)) {conf.high * 100
    multivariate_formula <- as.formula(paste(outcome_name, "~", paste(predictors, collapse = " + ")))
  } else {
    multivariate_formula <- custom_formula
  } diagnostics <- calculate_diagnostics(model)
  multivariate_model <- fit_model(multivariate_formula)
  multivariate_results <- process_results(multivariate_model, "multivariate")

  # Combine resultsmodels
  bind_rows(univariate_results, multivariate_results)predictor) {
}   formula <- as.formula(paste(outcome_name, "~", predictor))
    model <- fit_model(formula)
kkonehot <- function(data, column) {te") %>% mutate(predictor = predictor)
  # Check if input is a tibble
  if (!inherits(data, "tbl_df")) {
    stop("Input 'data' must be a tibble")
  }f (is.null(custom_formula)) {
    multivariate_formula <- as.formula(paste(outcome_name, "~", paste(predictors, collapse = " + ")))
  # Check if column exists in data
  if (!column %in% names(data)) {m_formula
    stop("Specified column not found in the dataset")
  }ultivariate_model <- fit_model(multivariate_formula)
  multivariate_results <- process_results(multivariate_model, "multivariate")
  # Get the variable
  var <- data[[column]]
  bind_rows(univariate_results, multivariate_results)
  # Check if variable is character or factor
  if (!(is.character(var) || is.factor(var))) {
    stop("Column must be character or factor type for one-hot encoding")
  } Check if input is a tibble
  if (!inherits(data, "tbl_df")) {
  # Convert character to factor if needed
  if (is.character(var)) {
    var <- as.factor(var)
  } Check if column exists in data
  if (!column %in% names(data)) {
  # Get levels of the factorot found in the dataset")
  levels <- levels(var)

  # Create one-hot encoded columns
  for (level in levels) {
    new_col_name <- paste0(column, "_", level)
    data[[new_col_name]] <- ifelse(var == level, 1, 0)
  }f (!(is.character(var) || is.factor(var))) {
    stop("Column must be character or factor type for one-hot encoding")
  # Return the modified dataset
  return(data)
} # Convert character to factor if needed
  if (is.character(var)) {
    var <- as.factor(var)
# Function to compute full summary statistics with grouped data and pipe support
# Parameters:
# - data: Data frame containing the column of interest
# - col: Unquoted column name to analyze
# - var_name: Optional name for the variable (defaults to column name)
# - verbose: Output level ("full", "basic", or "custom")
# - stats: For verbose="custom", list of stats to include (e.g., c("mean", "median"))
# - pairwise: Whether to compute pairwise comparisons for categorical variables
# - chi_probs: Optional expected probabilities for chi-square test
  }
kk_summary <- function(data, col, var_name = NULL,
                       verbose = c("full", "basic", "custom"),
                       stats = NULL, pairwise = TRUE, chi_probs = NULL) {
  verbose <- match.arg(verbose)
  if (!is.data.frame(data)) stop("Input 'data' must be a data frame")

  col_quo <- enquo(col)ull summary statistics with grouped data and pipe support
  col_name <- quo_name(col_quo)
  if (is.null(var_name)) var_name <- col_name interest
# - col: Unquoted column name to analyze
  x <- pull(data, {{ col }})for the variable (defaults to column name)
# - verbose: Output level ("full", "basic", or "custom")
  compute_stats <- function(x, var_name) {tats to include (e.g., c("mean", "median"))
    if (is.null(x)) stop("Input cannot be NULL")isons for categorical variables
    if (is.data.frame(x) || is.matrix(x)) stop("Input must be a vector, not a data frame or matrix")
    if (!is.vector(x) && !is.factor(x)) stop("Input must be a vector or factor")
kk_summary <- function(data, col, var_name = NULL,
    if (is.factor(x)) x <- as.character(x) "basic", "custom"),
                       stats = NULL, pairwise = TRUE, chi_probs = NULL) {
    result <- list()rg(verbose)
  if (!is.data.frame(data)) stop("Input 'data' must be a data frame")
    n_miss <- sum(is.na(x))
    result$var_name <- var_name
    result$type <- ifelse(is.numeric(x), "numeric", "categorical")
    result$n_total <- length(x)me <- col_name
    result$n_miss <- n_miss
    result$n_valid <- length(x) - n_miss
    result$miss_pct <- (n_miss / length(x)) * 100
  compute_stats <- function(x, var_name) {
    if (result$n_valid == 0) {t cannot be NULL")
      result$note <- "All values are missing"p("Input must be a vector, not a data frame or matrix")
      return(result)) && !is.factor(x)) stop("Input must be a vector or factor")
    }
    if (is.factor(x)) x <- as.character(x)
    if (verbose == "basic") {
      stats <- if (result$type == "numeric") {
        c("mean", "median", "sd", "min", "max")
      } else {sum(is.na(x))
        c("level_summary", "mode")
      }ult$type <- ifelse(is.numeric(x), "numeric", "categorical")
    } else if (verbose == "custom" && !is.null(stats)) {
      stats <- stats n_miss
    } else {_valid <- length(x) - n_miss
      stats <- if (result$type == "numeric") {100
        c(
          "mean", "huber_mean", "trim_mean", "geometric_mean", "median", "min", "max", "range",
          "variance", "sd", "se", "cv_pct", "mad", "iqr", "q1", "q3", "skewness",
          "kurtosis", "pct_5_95", "ci_mean_low", "ci_mean_up", "shapiro_p",
          "shapiro_int", "ks_p", "ks_int", "n_outliers", "outlier_values"
        )
      } else {e == "basic") {
        c(s <- if (result$type == "numeric") {
          "levels", "n_unique", "mode", "level_summary", "entropy", "evenness",
          "gini_index", "chi_p", "chi_int", "chi_note", "cramer_v", "pairwise_p", "pairwise_note"
        )("level_summary", "mode")
      }
    } else if (verbose == "custom" && !is.null(stats)) {
      stats <- stats
    if (result$type == "numeric") {
      x <- na.omit(x)sult$type == "numeric") {
      if ("mean" %in% stats) result$mean <- mean(x)
      if ("huber_mean" %in% stats) {m_mean", "geometric_mean", "median", "min", "max", "range",
        huber_m <- function(x, c = 1.345, max_iter = 100, tol = 1e-6) {skewness",
          mu <- median(x)t_5_95", "ci_mean_low", "ci_mean_up", "shapiro_p",
          sigma <- mad(x, constant = 1.4826)n_outliers", "outlier_values"
          for (i in 1:max_iter) {
            residuals <- x - mu
            weights <- pmin(c / abs(residuals / sigma), 1)
            weights[is.na(weights)] <- 1"level_summary", "entropy", "evenness",
            mu_new <- sum(weights * x) / sum(weights)", "cramer_v", "pairwise_p", "pairwise_note"
            if (abs(mu_new - mu) < tol) break
            mu <- mu_new
          }
          return(mu)
        }esult$type == "numeric") {
        result$huber_mean <- huber_m(x)
      }f ("mean" %in% stats) result$mean <- mean(x)
      if ("trim_mean" %in% stats) result$trim_mean <- mean(x, trim = 0.1)
      if ("geometric_mean" %in% stats) result$geometric_mean <- if (all(x > 0)) exp(mean(log(x))) else NA
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
      if ("skewness" %in% stats) result$skewness <- moments::skewness(x)
      if ("kurtosis" %in% stats) result$kurtosis <- moments::kurtosis(x)
      if ("pct_5_95" %in% stats) result$pct_5_95 <- list(quantile(x, probs = c(0.05, 0.95)))
      if ("ci_mean_low" %in% stats || "ci_mean_up" %in% stats) {if (all(x > 0)) exp(mean(log(x))) else NA
        t_test <- t.test(x, conf.level = 0.95)- median(x)
        result$ci_mean_low <- t_test$conf.int[1]
        result$ci_mean_up <- t_test$conf.int[2])
      }f ("range" %in% stats) result$range <- max(x) - min(x)
      if ("shapiro_p" %in% stats || "shapiro_int" %in% stats) {
        if (length(x) >= 3 && length(x) <= 5000) {
          shapiro_test <- shapiro.test(x)d(x) / sqrt(length(x))
          result$shapiro_p <- shapiro_test$p.value(x) / mean(x) * 100
          result$shapiro_int <- ifelse(shapiro_test$p.value > 0.05,
            "normal (p > 0.05)",lt$iqr <- IQR(x)
            "non-normal (p <= 0.05)" <- quantile(x, 0.25)
          )q3" %in% stats) result$q3 <- quantile(x, 0.75)
        } else {ess" %in% stats) result$skewness <- moments::skewness(x)
          result$shapiro_p <- NA result$kurtosis <- moments::kurtosis(x)
          result$shapiro_int <- "sample size out of Shapiro-Wilk range (3-5000)".05, 0.95)))
        }("ci_mean_low" %in% stats || "ci_mean_up" %in% stats) {
      } t_test <- t.test(x, conf.level = 0.95)
      if ("ks_p" %in% stats || "ks_int" %in% stats) {
        if (length(x) > 5000) {test$conf.int[2]
          ks_test <- tryCatch(
            {apiro_p" %in% stats || "shapiro_int" %in% stats) {
              suppressWarnings(ks.test(x, "pnorm", mean(x), sd(x)))
            },iro_test <- shapiro.test(x)
            error = function(e) NULLo_test$p.value
          )esult$shapiro_int <- ifelse(shapiro_test$p.value > 0.05,
          if (!is.null(ks_test)) {
            result$ks_p <- ks_test$p.value
            result$ks_int <- ifelse(ks_test$p.value > 0.05,
              "normal (KS p > 0.05)",
              "non-normal (KS p <= 0.05)"
            )ult$shapiro_int <- "sample size out of Shapiro-Wilk range (3-5000)"
          } else {
            result$ks_p <- NA
            result$ks_int <- "KS test failed"stats) {
          }(length(x) > 5000) {
        } else {t <- tryCatch(
          result$ks_p <- NA
          result$ks_int <- "KS test not performed (n <= 5000)"(x)))
        }   },
      }     error = function(e) NULL
      if ("n_outliers" %in% stats || "outlier_values" %in% stats) {
        q1 <- quantile(x, 0.25)) {
        q3 <- quantile(x, 0.75)est$p.value
        iqr <- q3 - q1int <- ifelse(ks_test$p.value > 0.05,
        lower_fence <- q1 - 1.5 * iqr
        upper_fence <- q3 + 1.5 * iqr05)"
        outliers <- x[x < lower_fence | x > upper_fence]
        result$n_outliers <- length(outliers)
        result$outlier_values <- list(outliers)
      }     result$ks_int <- "KS test failed"
    } else {
      x <- as.factor(x)
      x <- x[!is.na(x)]- NA
      freq_table <- table(x)KS test not performed (n <= 5000)"
      n_valid <- sum(freq_table)
      low_counts <- all(freq_table < 5) && n_valid < 10
      if (low_counts && any(c("chi_p", "pairwise_p") %in% stats)) {
        result$note <- "Warning: Very low counts; chi-square and pairwise tests may be unreliable"
      } q3 <- quantile(x, 0.75)
      if (any(c("n_unique", "evenness", "chi_p", "chi_int", "chi_note", "cramer_v", "pairwise_p") %in% stats) && result$n_valid > 0) {
        result$n_unique <- length(unique(x))
      } upper_fence <- q3 + 1.5 * iqr
      if ("levels" %in% stats) result$levels <- list(levels(x))
      if ("mode" %in% stats) result$mode <- names(sort(table(x), decreasing = TRUE))[1]
      if ("level_summary" %in% stats) {utliers)
        proportions <- prop.table(freq_table) * 100
        ci_list <- lapply(names(freq_table), function(level) {
          binom_result <- binom.test(freq_table[level], n_valid, conf.level = 0.95)
          tibble::tibble(
            level = level,x)
            count = as.numeric(freq_table[level]),
            prop_pct = proportions[level], n_valid < 10
            ci_low_pct = binom_result$conf.int[1] * 100,% stats)) {
            ci_up_pct = binom_result$conf.int[2] * 100square and pairwise tests may be unreliable"
          )
        })any(c("n_unique", "evenness", "chi_p", "chi_int", "chi_note", "cramer_v", "pairwise_p") %in% stats)) {
        result$level_summary <- list(dplyr::bind_rows(ci_list))
      }
      if ("entropy" %in% stats || "evenness" %in% stats) {s(x))
        freqs <- prop.table(freq_table)e <- names(sort(table(x), decreasing = TRUE))[1]
        result$entropy <- entropy::entropy(freqs, unit = "log2")
      } proportions <- prop.table(freq_table) * 100
      if ("evenness" %in% stats) {eq_table), function(level) {
        result$evenness <- if (!is.null(result$n_unique) && result$n_unique > 1) {)
          result$entropy / log2(result$n_unique)
        } else {l = level,
          0 count = as.numeric(freq_table[level]),
        }   prop_pct = proportions[level],
      }     ci_low_pct = binom_result$conf.int[1] * 100,
      if ("gini_index" %in% stats) {$conf.int[2] * 100
        p <- prop.table(freq_table)
        result$gini_index <- 1 - sum(p^2)
      } result$level_summary <- list(dplyr::bind_rows(ci_list))
      if (any(c("chi_p", "chi_int", "chi_note", "cramer_v") %in% stats) && result$n_unique > 1) {
        expected_prob <- if (!is.null(chi_probs)) {tats) {
          if (length(chi_probs) != result$n_unique || sum(chi_probs) != 1) {
            stop("chi_probs must match number of levels and sum to 1")
          }
          chi_probs" %in% stats) {
        } else {venness <- if (!is.null(result$n_unique) && result$n_unique > 1) {
          rep(1 / result$n_unique, result$n_unique)
        } else {
        if (low_counts) {
          result$chi_p <- NA
          result$chi_int <- "chi-square test skipped due to very low counts"
          result$chi_note <- "no test performed"
          result$cramer_v <- NAble)
        } else if (any(freq_table < 5) || n_valid < 10) {
          chi_test <- chisq.test(freq_table, p = expected_prob, simulate.p.value = TRUE, B = 2000)
          result$chi_note <- "simulated p-value used due to low counts or small sample size" 1) {
          result$chi_p <- chi_test$p.valueprobs)) {
          result$chi_int <- ifelse(chi_test$p.value > 0.05,hi_probs) != 1) {
            "uniform distribution (p > 0.05)",of levels and sum to 1")
            "non-uniform distribution (p <= 0.05)"
          )hi_probs
          result$cramer_v <- sqrt(chi_test$statistic / (n_valid * (result$n_unique - 1)))
        } else {/ result$n_unique, result$n_unique)
          chi_test <- chisq.test(freq_table, p = expected_prob)
          result$chi_note <- "standard chi-square test"
          result$chi_p <- chi_test$p.value
          result$chi_int <- ifelse(chi_test$p.value > 0.05, very low counts"
            "uniform distribution (p > 0.05)",d"
            "non-uniform distribution (p <= 0.05)"
          )lse if (any(freq_table < 5) || n_valid < 10) {
          result$cramer_v <- sqrt(chi_test$statistic / (n_valid * (result$n_unique - 1)))B = 2000)
        } result$chi_note <- "simulated p-value used due to low counts or small sample size"
      } else if (any(c("chi_p", "chi_int", "chi_note", "cramer_v") %in% stats)) {
        result$chi_p <- NA- ifelse(chi_test$p.value > 0.05,
        result$chi_int <- "chi-square test not applicable (single level)"
        result$chi_note <- "no test performed"05)"
        result$cramer_v <- NA
      }   result$cramer_v <- sqrt(chi_test$statistic / (n_valid * (result$n_unique - 1)))
      if ("pairwise_p" %in% stats && pairwise && result$n_unique > 2 && all(freq_table >= 1) && result$n_unique <= 10) {
        level_names <- names(freq_table)ble, p = expected_prob)
        pairwise_results <- tibble::tibble(level1 = character(), level2 = character(), p_value = numeric())
        for (i in 1:(length(level_names) - 1)) {
          for (j in (i + 1):length(level_names)) {e > 0.05,
            test_result <- binom.test(c(freq_table[i], freq_table[j]),
              n = c(freq_table[i] + freq_table[j], freq_table[i] + freq_table[j]),
              p = 0.5
            )ult$cramer_v <- sqrt(chi_test$statistic / (n_valid * (result$n_unique - 1)))
            pairwise_results <- dplyr::add_row(pairwise_results,
              level1 = level_names[i],nt", "chi_note", "cramer_v") %in% stats)) {
              level2 = level_names[j],
              p_value = test_result$p.valuenot applicable (single level)"
            )t$chi_note <- "no test performed"
          }ult$cramer_v <- NA
        }
        pairwise_results$p_value <- p.adjust(pairwise_results$p_value, method = "holm")>= 1) && result$n_unique <= 10) {
        result$pairwise_p <- list(pairwise_results)
        result$pairwise_note <- "Holm-adjusted p-values from exact binomial tests"r(), p_value = numeric())
      } else if ("pairwise_p" %in% stats) {1)) {
        reason <- if (!pairwise) {(level_names)) {
          "pairwise comparisons disabled"req_table[i], freq_table[j]),
        } else if (result$n_unique <= 2) {able[j], freq_table[i] + freq_table[j]),
          "2 or fewer levels"
        } else if (any(freq_table < 1)) {
          "zero counts in some levels":add_row(pairwise_results,
        } else {vel1 = level_names[i],
          "too many levels (>10)"s[j],
        }     p_value = test_result$p.value
        result$pairwise_p <- list("not performed")
        result$pairwise_note <- paste("pairwise comparisons skipped:", reason)
      } }
    }   pairwise_results$p_value <- p.adjust(pairwise_results$p_value, method = "holm")
        result$pairwise_p <- list(pairwise_results)
    return(result)rwise_note <- "Holm-adjusted p-values from exact binomial tests"
  }   } else if ("pairwise_p" %in% stats) {
        reason <- if (!pairwise) {
  if (dplyr::is_grouped_df(data)) {abled"
    result <- data %>%ult$n_unique <= 2) {
      dplyr::group_by(dplyr::across(dplyr::all_of(dplyr::group_vars(.)))) %>%
      dplyr::summarise(stats = list(compute_stats({{ col }}, var_name)), .groups = "keep") %>%
      tidyr::unnest_wider(stats)evels"
  } else {else {
    result <- tibble::as_tibble(compute_stats(x, var_name))
  }     }
  return(result)airwise_p <- list("not performed")
}       result$pairwise_note <- paste("pairwise comparisons skipped:", reason)
      }
# Function for time series

# Load required packages
required_pkgs <- c("dplyr", "moments", "tseries", "pracma", "zoo", "rugarch", "forecast", "entropy")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {)) {
  warning(sprintf("Missing packages: %s. Some metrics may return NA.", paste(missing_pkgs, collapse = ", ")))
}     dplyr::group_by(dplyr::across(dplyr::all_of(dplyr::group_vars(.)))) %>%
      dplyr::summarise(stats = list(compute_stats({{ col }}, var_name)), .groups = "keep") %>%
# Define the functionider(stats)
kk_time_series <- function(data, value_col = NULL, date_col = "date", group_cols = NULL, wide_format = FALSE, skip_advanced = FALSE, round_digits = 4) {
  # Validate data inputs_tibble(compute_stats(x, var_name))
  if (missing(data)) {
    stop("Argument 'data' is missing, with no default")
  }

  # Validate value columns
  if (is.null(value_col) && !"value" %in% names(data)) {
    stop("Please specify a value column using value_col parameter or ensure 'value' column exists")
  }uired_pkgs <- c("dplyr", "moments", "tseries", "pracma", "zoo", "rugarch", "forecast", "entropy")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  # Convert value_col to symbol and extract
  if (!is.null(value_col)) {ackages: %s. Some metrics may return NA.", paste(missing_pkgs, collapse = ", ")))
    if (is.character(value_col)) value_col <- rlang::sym(value_col)
    data <- dplyr::mutate(data, value = !!value_col)
  }efine the function
kk_time_series <- function(data, value_col = NULL, date_col = "date", group_cols = NULL, wide_format = FALSE, skip_advanced = FALSE, round_digits = 4) {
  # Validate date column
  if (!date_col %in% names(data)) {
    stop(sprintf("Date column '%s' not found in data. Available columns: %s", date_col, paste(names(data), collapse = ", ")))
  }
  if (is.character(date_col)) date_col <- rlang::sym(date_col)
  # Validate value column
  # Ensure date is proper formatlue" %in% names(data)) {
  if (!inherits(data[[date_col]], c("Date", "POSIXct"))) {rameter or ensure 'value' column exists")
    stop(sprintf("Column '%s' must be of class Date or POSIXct. Current class: %s", date_col, class(data[[date_col]])[1]))
  }
  # Convert value_col to symbol and extract
  # Handle groupinge_col)) {
  if (!is.null(group_cols)) {l)) value_col <- rlang::sym(value_col)
    group_cols <- rlang::syms(group_cols)!value_col)
    data <- dplyr::group_by(data, !!!group_cols)
    # Check for duplicate dates within groups
    data <- data %>%lumn
      dplyr::group_by(!!!group_cols) %>%
      dplyr::mutate(.dup_dates = any(duplicated(!!date_col))) %>%olumns: %s", date_col, paste(names(data), collapse = ", ")))
      dplyr::ungroup()
    if (any(data$.dup_dates)) {ate_col <- rlang::sym(date_col)
      stop("Duplicate dates found within groups. Ensure unique dates per group or preprocess data to aggregate duplicates.")
    }nsure date is proper format
    data <- dplyr::select(data, -.dup_dates)"POSIXct"))) {
    data <- dplyr::group_by(data, !!!group_cols)ate or POSIXct. Current class: %s", date_col, class(data[[date_col]])[1]))
  } else {
    # Check for duplicate dates in ungrouped data
    if (any(duplicated(data[[date_col]]))) {
      stop("Duplicate dates found in data. Ensure unique dates or preprocess data to aggregate duplicates.")
    }roup_cols <- rlang::syms(group_cols)
  } data <- dplyr::group_by(data, !!!group_cols)
    # Check for duplicate dates within groups
  # Extract group keys for proper labeling
  group_keys <- if (!is.null(group_cols)) dplyr::group_keys(data) else NULL
  group_labels <- if (!is.null(group_keys)) {ed(!!date_col))) %>%
    apply(group_keys, 1, paste, collapse = "_")
  } else {y(data$.dup_dates)) {
    NULLop("Duplicate dates found within groups. Ensure unique dates per group or preprocess data to aggregate duplicates.")
  } }
    data <- dplyr::select(data, -.dup_dates)
  # Process time series for each grouproup_cols)
  result <- data %>%
    dplyr::arrange(!!date_col) %>% ungrouped data
    dplyr::filter(!is.na(value)) %>%l]]))) {
    dplyr::group_map(~ {tes found in data. Ensure unique dates or preprocess data to aggregate duplicates.")
      # Extract values and count observations
      y <- .x$value
      n <- length(y)
      if (n < 2) stop("Insufficient data points for time series analysis (n < 2)")
      if (n < 5) warning("Sample size is small (n = ", n, "). Some metrics may be unreliable or return NA.")
  group_labels <- if (!is.null(group_keys)) {
      # Check for negative valuesollapse = "_")
      if (any(y <= 0)) {
        warning("Non-positive values detected in group. Geometric mean growth rate may be unreliable.")
      }

      # Count missing values (before filtering)
      missing_count <- sum(is.na(.x$value))
    dplyr::arrange(!!date_col) %>%
      # Create zoo object for time series
      ts_zoo <- try(zoo::zoo(y, order.by = .x[[as.character(date_col)]]), silent = TRUE)
      if (inherits(ts_zoo, "try-error")) stop("Failed to create zoo object")
      y <- .x$value
      # Basic statistics
      mean_val <- mean(y, na.rm = TRUE)a points for time series analysis (n < 2)")
      median_val <- stats::median(y, na.rm = TRUE)= ", n, "). Some metrics may be unreliable or return NA.")
      sd_val <- stats::sd(y, na.rm = TRUE)
      var_val <- stats::var(y, na.rm = TRUE)
      min_val <- min(y, na.rm = TRUE)
      max_val <- max(y, na.rm = TRUE)detected in group. Geometric mean growth rate may be unreliable.")
      range_val <- diff(range(y, na.rm = TRUE))
      q1_val <- stats::quantile(y, 0.25, na.rm = TRUE)
      q3_val <- stats::quantile(y, 0.75, na.rm = TRUE)
      skewness_val <- if (requireNamespace("moments", quietly = TRUE) && n >= 3) {
        moments::skewness(y, na.rm = TRUE)
      } else { zoo object for time series
        warning("Sample size < 3 or 'moments' package missing for skewness. Returning NA.")
        NAinherits(ts_zoo, "try-error")) stop("Failed to create zoo object")
      }
      kurtosis_val <- if (requireNamespace("moments", quietly = TRUE) && n >= 4) {
        moments::kurtosis(y, na.rm = TRUE)
      } else {al <- stats::median(y, na.rm = TRUE)
        warning("Sample size < 4 or 'moments' package missing for kurtosis. Returning NA.")
        NAval <- stats::var(y, na.rm = TRUE)
      }in_val <- min(y, na.rm = TRUE)
      cv_val <- if (mean_val != 0 && n >= 2) {
        sd_val / abs(mean_val)y, na.rm = TRUE))
      } else {- stats::quantile(y, 0.25, na.rm = TRUE)
        warning("Mean is zero or sample size < 2 for CV. Returning NA.")
        NAness_val <- if (requireNamespace("moments", quietly = TRUE) && n >= 3) {
      } moments::skewness(y, na.rm = TRUE)
      } else {
      # Outlier count (IQR method)r 'moments' package missing for skewness. Returning NA.")
      outlier_count <- try(
        {
          iqr <- q3_val - q1_valeNamespace("moments", quietly = TRUE) && n >= 4) {
          sum(y < (q1_val - 1.5 * iqr) | y > (q3_val + 1.5 * iqr))
        },se {
        silent = TRUEle size < 4 or 'moments' package missing for kurtosis. Returning NA.")
      ) NA
      if (inherits(outlier_count, "try-error")) outlier_count <- NA
      cv_val <- if (mean_val != 0 && n >= 2) {
      # Chain-base metricsval)
      abs_increase_chain <- if (n > 1) diff(y) else NA
      t_y_chain <- if (n > 1 && all(y[1:(n - 1)] != 0)) {Returning NA.")
        (y[2:n] / y[1:(n - 1)])
      } else {
        warning("Zero values in denominator for chain-base growth rate. Returning NA.")
        rep(NA, n - 1)(IQR method)
      }utlier_count <- try(
      t_y_star_chain <- if (n > 1 && all(y[1:(n - 1)] != 0)) {
        ((y[2:n] / y[1:(n - 1)]) - 1) * 100
      } else {y < (q1_val - 1.5 * iqr) | y > (q3_val + 1.5 * iqr))
        warning("Zero values in denominator for chain-base rate of increase. Returning NA.")
        rep(NA, n - 1)
      }
      if (inherits(outlier_count, "try-error")) outlier_count <- NA
      # Fixed-base metrics (using first value as base)
      base_value <- y[1]cs
      abs_increase_fixed <- if (n > 0 && base_value != 0) {
        y - base_value(n > 1 && all(y[1:(n - 1)] != 0)) {
      } else {] / y[1:(n - 1)])
        warning("Base value is zero for fixed-base metrics. Returning NA.")
        rep(NA, n)ero values in denominator for chain-base growth rate. Returning NA.")
      } rep(NA, n - 1)
      t_y_fixed <- if (n > 0 && base_value != 0) (y / base_value) else rep(NA, n)
      t_y_star_fixed <- if (n > 0 && base_value != 0) ((y / base_value) - 1) * 100 else rep(NA, n)
        ((y[2:n] / y[1:(n - 1)]) - 1) * 100
      # Geometric mean growth rate (chain-base, coefficient form)
      t_y_chain_coeff <- t_y_chain[!is.na(t_y_chain) & t_y_chain > 0]crease. Returning NA.")
      geom_mean_growth <- if (length(t_y_chain_coeff) > 0) {
        exp(mean(log(t_y_chain_coeff), na.rm = TRUE))
      } else {
        NA
      }_value <- y[1]
      geom_mean_growth_pct <- (geom_mean_growth - 1) * 100 # Percentage formbs_increase_fixed <- if (n > 0 && base_value != 0) {

      # Mean rate of increase (T̅′) based on geometric mean growth rate      } else {
      mean_incrate <- geom_mean_growth - 1 # Coefficient form for fixed-base metrics. Returning NA.")
      mean_incrate_pct <- geom_mean_growth_pct # Percentage form

      # Summary statistics for chain-base and fixed-base metrics      t_y_fixed <- if (n > 0 && base_value != 0) (y / base_value) else rep(NA, n)
      mean_absinc_chain <- if (!all(is.na(abs_increase_chain))) mean(abs_increase_chain, na.rm = TRUE) else NA_value) - 1) * 100 else rep(NA, n)
      mean_devrate_chain <- if (!all(is.na(t_y_chain))) mean(t_y_chain, na.rm = TRUE) else NA
      mean_incrate_chain <- if (!all(is.na(t_y_star_chain))) mean(t_y_star_chain, na.rm = TRUE) else NA
      mean_absinc_fixed <- if (!all(is.na(abs_increase_fixed))) mean(abs_increase_fixed, na.rm = TRUE) else NA
      mean_devrate_fixed <- if (!all(is.na(t_y_fixed))) mean(t_y_fixed, na.rm = TRUE) else NA
      mean_incrate_fixed <- if (!all(is.na(t_y_star_fixed))) mean(t_y_star_fixed, na.rm = TRUE) else NA

      sd_absinc_chain <- if (length(abs_increase_chain[!is.na(abs_increase_chain)]) > 1) {        warning("No valid positive chain-base growth rates for geometric mean. Returning NA.")
        sd(abs_increase_chain, na.rm = TRUE)
      } else {
        NAn_growth_pct <- if (!is.na(geom_mean_growth)) (geom_mean_growth - 1) * 100 else NA
      }
      sd_devrate_chain <- if (length(t_y_chain[!is.na(t_y_chain)]) > 1) { Mean rate of increase (T̅′)
        sd(t_y_chain, na.rm = TRUE)else NA
      } else {a(geom_mean_growth)) geom_mean_growth_pct else NA
        NA
      }mmary statistics for chain-base and fixed-base metrics
      sd_incrate_chain <- if (length(t_y_star_chain[!is.na(t_y_star_chain)]) > 1) {ean_absinc_chain <- if (!all(is.na(abs_increase_chain))) mean(abs_increase_chain, na.rm = TRUE) else NA
        sd(t_y_star_chain, na.rm = TRUE)E) else NA
      } else {na(t_y_star_chain))) mean(t_y_star_chain, na.rm = TRUE) else NA
        NAinc_fixed <- if (!all(is.na(abs_increase_fixed))) mean(abs_increase_fixed, na.rm = TRUE) else NA
      }_devrate_fixed <- if (!all(is.na(t_y_fixed))) mean(t_y_fixed, na.rm = TRUE) else NA
      sd_absinc_fixed <- if (length(abs_increase_fixed[!is.na(abs_increase_fixed)]) > 1) {ean_incrate_fixed <- if (!all(is.na(t_y_star_fixed))) mean(t_y_star_fixed, na.rm = TRUE) else NA
        sd(abs_increase_fixed, na.rm = TRUE)
      } else {ease_chain[!is.na(abs_increase_chain)]) > 1) {
        NA_increase_chain, na.rm = TRUE)
      }se {
      sd_devrate_fixed <- if (length(t_y_fixed[!is.na(t_y_fixed)]) > 1) { NA
        sd(t_y_fixed, na.rm = TRUE)
      } else {h(t_y_chain[!is.na(t_y_chain)]) > 1) {
        NA_chain, na.rm = TRUE)
      }se {
      sd_incrate_fixed <- if (length(t_y_star_fixed[!is.na(t_y_star_fixed)]) > 1) { NA
        sd(t_y_star_fixed, na.rm = TRUE)
      } else {_star_chain[!is.na(t_y_star_chain)]) > 1) {
        NA_star_chain, na.rm = TRUE)
      }se {
 NA
      # Rate of change (existing ROC, for consistency with original function)      }
      roc <- if (n > 1 && all(y[1:(n - 1)] != 0)) {xed)]) > 1) {
        (y[2:n] - y[1:(n - 1)]) / y[1:(n - 1)] * 100
      } else {
        warning("Zero values in denominator for ROC. Returning NA.")
        rep(NA, n - 1)
      } <- if (length(t_y_fixed[!is.na(t_y_fixed)]) > 1) {
      roc_stats <- if (!all(is.na(roc))) { sd(t_y_fixed, na.rm = TRUE)
        c(
          mean(roc, na.rm = TRUE),
          sd(roc, na.rm = TRUE),
          min(roc, na.rm = TRUE),ngth(t_y_star_fixed[!is.na(t_y_star_fixed)]) > 1) {
          max(roc, na.rm = TRUE)= TRUE)
        )
      } else {A
        rep(NA, 4)
      }
 Rate of change (existing ROC, for consistency with original function)
      # Advanced metrics (skipped if skip_advanced = TRUE and n is small)      roc <- if (n > 1 && all(y[1:(n - 1)] != 0)) {
      trend_slope <- trend_strength <- NA
      if (!skip_advanced || n >= 10) {
        # Time index for regressionnator for ROC. Returning NA.")
        time_index <- as.numeric(difftime(.x[[as.character(date_col)]], min(.x[[as.character(date_col)]]), units = "days"))

        # Trend strength      roc_stats <- if (!all(is.na(roc))) {
        trend_model <- try(stats::lm(value ~ time_index, data = .x), silent = TRUE)
        if (!inherits(trend_model, "try-error")) {
          trend_slope <- try(stats::coef(trend_model)[2], silent = TRUE)
          if (!inherits(trend_slope, "try-error")) {
            trend_fit <- try(stats::predict(trend_model), silent = TRUE)
            if (!inherits(trend_fit, "try-error")) {
              detrended <- y - trend_fit
              trend_strength <- try(1 - stats::var(detrended) / stats::var(y), silent = TRUE)
              if (inherits(trend_strength, "try-error")) trend_strength <- NA
            }
          }ced metrics (skipped if skip_advanced = TRUE and n is small)
        }_slope <- trend_strength <- NA
      }(!skip_advanced || n >= 10) {
 # Time index for regression
      garch_vol <- NA        time_index <- as.numeric(difftime(.x[[as.character(date_col)]], min(.x[[as.character(date_col)]]), units = "days"))
      if (!skip_advanced || n >= 30) {
        garch_vol <- try(
          {y(stats::lm(value ~ time_index, data = .x), silent = TRUE)
            if (requireNamespace("rugarch", quietly = TRUE) && n >= 30) {(!inherits(trend_model, "try-error")) {
              suppressWarnings({
                # Difference data if non-stationaryope, "try-error")) {
                y_diff <- if (n > 1) diff(y) else yodel), silent = TRUE)
                spec <- rugarch::ugarchspec({
                  mean.model = list(armaOrder = c(1, 0)),
                  variance.model = list(model = "sGARCH", garchOrder = c(1, 1))ded) / stats::var(y), silent = TRUE)
                )
                garch_fit <- rugarch::ugarchfit(spec, y_diff, solver = "hybrid")
                if (is.null(garch_fit@fit$convergence) || garch_fit@fit$convergence != 0) {
                  # Try simpler model
                  spec <- rugarch::ugarchspec(
                    mean.model = list(armaOrder = c(0, 0)),
                    variance.model = list(model = "sGARCH", garchOrder = c(1, 1))
                  )
                  garch_fit <- rugarch::ugarchfit(spec, y_diff, solver = "hybrid")- try(
                }
                mean(rugarch::sigma(garch_fit), na.rm = TRUE)equireNamespace("rugarch", quietly = TRUE) && n >= 30) {
              })
            } else {# Difference data if non-stationary
              warning("Sample size < 30 or 'rugarch' package missing for GARCH. Returning NA.")ff <- if (n > 1) diff(y) else y
              NA
            }  mean.model = list(armaOrder = c(1, 0)),
          },     variance.model = list(model = "sGARCH", garchOrder = c(1, 1))
          silent = TRUE    )
        )it <- rugarch::ugarchfit(spec, y_diff, solver = "hybrid")
        if (inherits(garch_vol, "try-error") || is.null(garch_vol)) {       if (is.null(garch_fit@fit$convergence) || garch_fit@fit$convergence != 0) {
          warning("GARCH fitting failed. Returning NA.")
          garch_vol <- NA
        }model = list(armaOrder = c(0, 0)),
      }           variance.model = list(model = "sGARCH", garchOrder = c(1, 1))
           )
      dom_freq <- NA                  garch_fit <- rugarch::ugarchfit(spec, y_diff, solver = "hybrid")
      if (!skip_advanced || n >= 10) {
        dom_freq <- try(rch_fit), na.rm = TRUE)
          {
            freq <- stats::spectrum(y, plot = FALSE) } else {
            freq$freq[which.max(freq$spec)] package missing for GARCH. Returning NA.")
          },
          silent = TRUE}
        )
        if (inherits(dom_freq, "try-error")) dom_freq <- NA silent = TRUE
      }
 if (inherits(garch_vol, "try-error") || is.null(garch_vol)) {
      shannon_entropy <- NA          warning("GARCH fitting failed. Returning NA.")
      if (!skip_advanced || n >= 5) {
        shannon_entropy <- try(
          {
            if (requireNamespace("entropy", quietly = TRUE) && n >= 5) {
              breaks <- pretty(range(y), n = min(10, n / 5))
              binned <- cut(y, breaks, include.lowest = TRUE)
              counts <- table(binned)
              entropy::entropy(counts / sum(counts))
            } else {
              warning("Sample size < 5 or 'entropy' package missing for Shannon entropy. Returning NA.")q[which.max(freq$spec)]
              NA
            } = TRUE
          },
          silent = TRUEinherits(dom_freq, "try-error")) dom_freq <- NA
        )
        if (inherits(shannon_entropy, "try-error")) shannon_entropy <- NA
      }
f (!skip_advanced || n >= 5) {
      acf_lag1 <- NA        shannon_entropy <- try(
      pacf_lag1 <- NA
      if (!skip_advanced || n >= 3) {reNamespace("entropy", quietly = TRUE) && n >= 5) {
        acf_lag1 <- try(y), n = min(10, n / 5))
          {cut(y, breaks, include.lowest = TRUE)
            acf_result <- stats::acf(y, lag.max = 1, plot = FALSE, na.action = stats::na.pass)   counts <- table(binned)
            acf_result$acf[2]
          },
          silent = TRUE  warning("Sample size < 5 or 'entropy' package missing for Shannon entropy. Returning NA.")
        )
        if (inherits(acf_lag1, "try-error")) acf_lag1 <- NA   }

        pacf_lag1 <- try(          silent = TRUE
          {
            pacf_result <- stats::pacf(y, lag.max = 1, plot = FALSE, na.action = stats::na.pass)(inherits(shannon_entropy, "try-error")) shannon_entropy <- NA
            pacf_result$acf[1]
          },
          silent = TRUEg1 <- NA
        )
        if (inherits(pacf_lag1, "try-error")) pacf_lag1 <- NA(!skip_advanced || n >= 3) {
      }
   {
      ljung_pval <- NA            acf_result <- stats::acf(y, lag.max = 1, plot = FALSE, na.action = stats::na.pass)
      if (!skip_advanced || n >= 3) {$acf[2]
        ljung_pval <- try(
          {
            if (n >= 3) {
              test <- stats::Box.test(y, lag = min(10, n - 1), type = "Ljung-Box")lag1, "try-error")) acf_lag1 <- NA
              test$p.value
            } else {
              warning("Sample size < 3 for Ljung-Box test. Returning NA.")
              NAtion = stats::na.pass)
            }_result$acf[1]
          },
          silent = TRUElent = TRUE
        )
        if (inherits(ljung_pval, "try-error")) ljung_pval <- NAf (inherits(pacf_lag1, "try-error")) pacf_lag1 <- NA
      }

      adf_result <- c(NA, NA)      ljung_pval <- NA
      if (!skip_advanced || n >= 10) { >= 3) {
        adf_result <- try(
          {
            if (requireNamespace("tseries", quietly = TRUE) && n >= 10) { if (n >= 3) {
              suppressWarnings({ung-Box")
                test <- tseries::adf.test(y)
                c(test$p.value, test$statistic)
              })g-Box test. Returning NA.")
            } else {
              warning("Sample size < 10 or 'tseries' package missing for ADF test. Returning NA.")
              c(NA, NA)
            }
          },
          silent = TRUEinherits(ljung_pval, "try-error")) ljung_pval <- NA
        )
        if (inherits(adf_result, "try-error")) adf_result <- c(NA, NA)
      }
f (!skip_advanced || n >= 10) {
      kpss_result <- c(NA, NA)        adf_result <- try(
      if (!skip_advanced || n >= 10) {
        kpss_result <- try(ies", quietly = TRUE) && n >= 10) {
          {ngs({
            if (requireNamespace("tseries", quietly = TRUE) && n >= 10) {     test <- tseries::adf.test(y)
              suppressWarnings({
                test <- tseries::kpss.test(y)
                c(test$p.value, test$statistic)
              })ries' package missing for ADF test. Returning NA.")
            } else {NA, NA)
              warning("Sample size < 10 or 'tseries' package missing for KPSS test. Returning NA.")
              c(NA, NA)
            }
          },
          silent = TRUEinherits(adf_result, "try-error")) adf_result <- c(NA, NA)
        )
        if (inherits(kpss_result, "try-error")) kpss_result <- c(NA, NA)
      }
f (!skip_advanced || n >= 10) {
      hurst <- NA        kpss_result <- try(
      if (!skip_advanced || n >= 20) {
        hurst <- try(ies", quietly = TRUE) && n >= 10) {
          {sWarnings({
            if (requireNamespace("pracma", quietly = TRUE) && n >= 20) {     test <- tseries::kpss.test(y)
              temp <- capture.output({
                h_result <- pracma::hurstexp(y, display = FALSE)
              })
              h_result$Hsrning("Sample size < 10 or 'tseries' package missing for KPSS test. Returning NA.")
            } else {
              warning("Sample size < 20 or 'pracma' package missing for Hurst exponent. Returning NA.")
              NA
            } = TRUE
          },
          silent = TRUEinherits(kpss_result, "try-error")) kpss_result <- c(NA, NA)
        )
        if (inherits(hurst, "try-error")) hurst <- NA
      }
f (!skip_advanced || n >= 20) {
      # Compile results        hurst <- try(
      result_tibble <- dplyr::tibble(
        Metric = c(cma", quietly = TRUE) && n >= 20) {
          "Length of Series", "Mean", "Median", "Standard Deviation", "Variance", "Min", "Max", "Range",<- capture.output({
          "Q1 (25th Percentile)", "Q3 (75th Percentile)", "Skewness", "Kurtosis", "CV (SD/Mean)",
          "Missing Count", "Outlier Count",
          "Mean Abs Increase (Chain)", "Mean Growth Rate (Chain, Coeff)", "Mean Rate of Increase (Chain, %)",
          "Mean Abs Increase (Fixed)", "Mean Growth Rate (Fixed, Coeff)", "Mean Rate of Increase (Fixed, %)",
          "SD Abs Increase (Chain)", "SD Growth Rate (Chain, Coeff)", "SD Rate of Increase (Chain, %)",
          "SD Abs Increase (Fixed)", "SD Growth Rate (Fixed, Coeff)", "SD Rate of Increase (Fixed, %)",
          "Geometric Mean Growth Rate (Coeff)", "Geometric Mean Growth Rate (%)", "Mean Rate of Increase (T̅′, Coeff)", "Mean Rate of Increase (T̅′, %)",
          "Mean Rate of Change (%)", "SD Rate of Change (%)", "Min ROC (%)", "Max ROC (%)",
          "Trend Slope", "Trend Strength", "ADF p-value", "ADF statistic",
          "KPSS p-value", "KPSS statistic", "ACF Lag 1", "PACF Lag 1", "Ljung-Box p-value",
          "Hurst Exponent", "GARCH Volatility", "Shannon Entropy", "Dominant Frequency"
        ),
        Value = tryCatch(
          {
            c(t_tibble <- dplyr::tibble(
              n, mean_val, median_val, sd_val, var_val, min_val, max_val, range_val, = c(
              q1_val, q3_val, skewness_val, kurtosis_val, cv_val,in", "Max", "Range",
              missing_count, outlier_count,ss", "Kurtosis", "CV (SD/Mean)",
              mean_absinc_chain, mean_devrate_chain, mean_incrate_chain,
              mean_absinc_fixed, mean_devrate_fixed, mean_incrate_fixed,, "Mean Rate of Increase (Chain, %)",
              sd_absinc_chain, sd_devrate_chain, sd_incrate_chain,, "Mean Rate of Increase (Fixed, %)",
              sd_absinc_fixed, sd_devrate_fixed, sd_incrate_fixed,)", "SD Rate of Increase (Chain, %)",
              geom_mean_growth, geom_mean_growth_pct, mean_incrate, mean_incrate_pct,)", "SD Rate of Increase (Fixed, %)",
              roc_stats[1], roc_stats[2], roc_stats[3], roc_stats[4],an Rate of Increase (T̅′, Coeff)", "Mean Rate of Increase (T̅′, %)",
              trend_slope, trend_strength,C (%)", "Max ROC (%)",
              adf_result[1], adf_result[2], "ADF p-value", "ADF statistic",
              kpss_result[1], kpss_result[2], "ACF Lag 1", "PACF Lag 1", "Ljung-Box p-value",
              acf_lag1, pacf_lag1, ljung_pval,", "Shannon Entropy", "Dominant Frequency"
              hurst, garch_vol, shannon_entropy, dom_freq
            )
          },
          error = function(e) {c(
            warning("Error in calculating statistics: ", e$message)an_val, sd_val, var_val, min_val, max_val, range_val,
            rep(NA, 48)
          }ount, outlier_count,
        )   mean_absinc_chain, mean_devrate_chain, mean_incrate_chain,
      )     mean_absinc_fixed, mean_devrate_fixed, mean_incrate_fixed,
       sd_absinc_chain, sd_devrate_chain, sd_incrate_chain,
      return(result_tibble)              sd_absinc_fixed, sd_devrate_fixed, sd_incrate_fixed,
    }, .keep = TRUE)wth, geom_mean_growth_pct, mean_incrate, mean_incrate_pct,
ats[1], roc_stats[2], roc_stats[3], roc_stats[4],
  # Bind results with proper group labels              trend_slope, trend_strength,
  if (!is.null(group_cols)) {],
    result <- dplyr::bind_rows(result, .id = "group_idx") %>% kpss_result[2],
      dplyr::mutate(group = group_labels[as.integer(group_idx)]) %>%
      dplyr::select(-group_idx)
  } else {
    result <- dplyr::bind_rows(result)},
  }
         warning("Error in calculating statistics: ", e$message)
  # Convert to wide format if requested            rep(NA, 48)
  if (wide_format) {
    result <- result %>%
      tidyr::pivot_wider(names_from = Metric, values_from = Value, names_repair = "minimal")
  }
   return(result_tibble)
  return(result)    }, .keep = TRUE)
}
 # Bind results with proper group labels
# Function for time series metrics  if (!is.null(group_cols)) {
kk_time_metrics <- function(data, value_col = NULL, date_col = "date", group_cols = NULL) {ult, .id = "group_idx") %>%
  # Validate data input
  if (missing(data)) {oup_idx)
    stop("Argument 'data' is missing, with no default")
  }

  # Validate value column
  if (is.null(value_col) && !"value" %in% names(data)) {t if requested
    stop("Please specify a value column using value_col parameter or ensure 'value' column exists")
  }
   tidyr::pivot_wider(names_from = Metric, values_from = Value, names_repair = "minimal")
  # Convert value_col to symbol and extract  }
  if (!is.null(value_col)) {
    if (is.character(value_col)) value_col <- rlang::sym(value_col)
    data <- dplyr::mutate(data, value = !!value_col)
  }
unction for time series metrics
  # Validate date columnkk_time_metrics <- function(data, value_col = NULL, date_col = "date", group_cols = NULL) {
  if (!date_col %in% names(data)) {
    stop(sprintf("Date column '%s' not found in data. Available columns: %s", date_col, paste(names(data), collapse = ", ")))
  }
  if (is.character(date_col)) date_col <- rlang::sym(date_col)

  # Ensure date is proper format  # Validate value column
  if (!inherits(data[[date_col]], c("Date", "POSIXct"))) {lue" %in% names(data)) {
    stop(sprintf("Column '%s' must be of class Date or POSIXct. Current class: %s", date_col, class(data[[date_col]])[1]))rameter or ensure 'value' column exists")
  }

  # Handle grouping  # Convert value_col to symbol and extract
  if (!is.null(group_cols)) {e_col)) {
    group_cols <- rlang::syms(group_cols)l)) value_col <- rlang::sym(value_col)
    data <- dplyr::group_by(data, !!!group_cols)!value_col)
  }

  # Process time series for each group  # Validate date column
  result <- data %>%
    dplyr::arrange(!!date_col) %>%te column '%s' not found in data. Available columns: %s", date_col, paste(names(data), collapse = ", ")))
    dplyr::filter(!is.na(value)) %>%
    dplyr::group_map(~ {ol <- rlang::sym(date_col)
      # Extract values, time index, and period
      y <- .x$value
      n <- length(y)a[[date_col]], c("Date", "POSIXct"))) {
      if (n < 2) stop("Insufficient data points for time series analysis")lumn '%s' must be of class Date or POSIXct. Current class: %s", date_col, class(data[[date_col]])[1]))
      time_index <- 1:n
      period <- .x[[as.character(date_col)]][order(.x[[as.character(date_col)]])]

      # Autocorrelation coefficient r1  if (!is.null(group_cols)) {
      r1 <- if (n > 1) {ls)
        cor(y[1:(n - 1)], y[2:n], use = "pairwise.complete.obs", method = "pearson")_by(data, !!!group_cols)
      } else {
        NA
      }s time series for each group
t <- data %>%
      # Durbin-Watson Q_BP    dplyr::arrange(!!date_col) %>%
      q_bp <- if (n > 1) {alue)) %>%
        sum((y[2:n] - y[1:(n - 1)])^2) / sum(y^2)
      } else {
        NAvalue
      } length(y)
f (n < 2) stop("Insufficient data points for time series analysis")
      # Ljung-Box test with p-value      time_index <- 1:n
      ljung_box_p <- if (n > 1) {te_col)]][order(.x[[as.character(date_col)]])]
        test <- try(stats::Box.test(y, lag = 1, type = "Ljung-Box"), silent = TRUE)
        if (!inherits(test, "try-error")) test$p.value else NA
      } else {
        NA1:(n - 1)], y[2:n], use = "pairwise.complete.obs", method = "pearson")
      }se {
 NA
      # Spearman correlation for trend and p-value      }
      spearman <- if (n > 1) {
        test <- try(stats::cor.test(time_index, y, method = "spearman", exact = FALSE), silent = TRUE)
        if (!inherits(test, "try-error")) c(test$estimate, test$p.value) else c(NA, NA)
      } else {
        c(NA, NA)
      }

      # Anderson-Darling test for normality with p-value
      ad_test <- if (n > 1 && requireNamespace("nortest", quietly = TRUE)) {
        if (n < 7) {
          warning("Sample size too small for Anderson-Darling test (< 7 observations). Returning NA.")stats::Box.test(y, lag = 1, type = "Ljung-Box"), silent = TRUE)
          c(NA, NA)
        } else {
          test <- try(nortest::ad.test(y), silent = TRUE)
          if (!inherits(test, "try-error")) c(test$statistic, test$p.value) else c(NA, NA)
        }
      } else {pearman correlation for trend and p-value
        c(NA, NA) <- if (n > 1) {
      }ry(stats::cor.test(time_index, y, method = "spearman", exact = FALSE), silent = TRUE)
 if (!inherits(test, "try-error")) c(test$estimate, test$p.value) else c(NA, NA)
      # Chain base metrics per period      } else {
      abs_increase_chain <- if (n > 1) c(NA, diff(y)) else rep(NA, n)
      t_y_chain <- if (n > 1) c(NA, (y[2:n] / y[1:(n - 1)])) else rep(NA, n) # Coefficient form
      t_y_star_chain <- if (n > 1) c(NA, ((y[2:n] / y[1:(n - 1)]) - 1) * 100) else rep(NA, n) # Percentage form

      # Fixed base period metrics per period (using first non-NA value as base)      ad_test <- if (n > 1 && requireNamespace("nortest", quietly = TRUE)) {
      base_value <- y[1]
      abs_increase_fixed <- if (n > 0) y - base_value else rep(NA, n)e size too small for Anderson-Darling test (< 7 observations). Returning NA.")
      t_y_fixed <- if (n > 0) (y / base_value) else rep(NA, n) # Coefficient form
      t_y_star_fixed <- if (n > 0) ((y / base_value) - 1) * 100 else rep(NA, n) # Percentage form

      # Geometric mean growth rate (chain-base, coefficient form)          if (!inherits(test, "try-error")) c(test$statistic, test$p.value) else c(NA, NA)
      t_y_chain_coeff <- t_y_chain[!is.na(t_y_chain)]
      geom_mean_growth <- if (length(t_y_chain_coeff) > 0) {
        exp(mean(log(t_y_chain_coeff), na.rm = TRUE))
      } else {
        NA
      }ain base metrics per period
      geom_mean_growth_pct <- (geom_mean_growth - 1) * 100 # Percentage formbs_increase_chain <- if (n > 1) c(NA, diff(y)) else rep(NA, n)
 # Coefficient form
      # Standard deviation of chain-base growth rates (in coefficient form)      t_y_star_chain <- if (n > 1) c(NA, ((y[2:n] / y[1:(n - 1)]) - 1) * 100) else rep(NA, n) # Percentage form
      t_y_chain_sd <- if (length(t_y_chain_coeff) > 1) {
        sd(t_y_chain_coeff, na.rm = TRUE)t non-NA value as base)
      } else {
        NAease_fixed <- if (n > 0) y - base_value else rep(NA, n)
      }fixed <- if (n > 0) (y / base_value) else rep(NA, n) # Coefficient form
      t_y_chain_sd_pct <- t_y_chain_sd * 100 # Percentage form_y_star_fixed <- if (n > 0) ((y / base_value) - 1) * 100 else rep(NA, n) # Percentage form

      # Mean rate of increase (T̅′) based on geometric mean growth rate      # Geometric mean growth rate (chain-base, coefficient form)
      mean_incrate <- geom_mean_growth - 1 # Coefficient form
      mean_incrate_pct <- geom_mean_growth_pct # Percentage form

      # Compute means for chain base metrics      } else {
      abs_increase_chain_mean <- if (!all(is.na(abs_increase_chain))) mean(abs_increase_chain, na.rm = TRUE) else NA
      t_y_chain_mean <- if (!all(is.na(t_y_chain))) mean(t_y_chain, na.rm = TRUE) else NA
      t_y_star_chain_mean <- if (!all(is.na(t_y_star_chain))) mean(t_y_star_chain, na.rm = TRUE) else NA
      abs_increase_fixed_mean <- if (!all(is.na(abs_increase_fixed))) mean(abs_increase_fixed, na.rm = TRUE) else NA
      t_y_fixed_mean <- if (!all(is.na(t_y_fixed))) mean(t_y_fixed, na.rm = TRUE) else NA      # Standard deviation of chain-base growth rates (in coefficient form)
      t_y_star_fixed_mean <- if (!all(is.na(t_y_star_fixed))) mean(t_y_star_fixed, na.rm = TRUE) else NAoeff) > 1) {

      # Compute standard deviations for per-period metrics
      abs_increase_chain_sd <- if (!all(is.na(abs_increase_chain))) sd(abs_increase_chain, na.rm = TRUE) else NA
      t_y_chain_sd <- if (!all(is.na(t_y_chain))) sd(t_y_chain, na.rm = TRUE) else NA      }
      t_y_star_chain_sd <- if (!all(is.na(t_y_star_chain))) sd(t_y_star_chain, na.rm = TRUE) else NAform
      abs_increase_fixed_sd <- if (!all(is.na(abs_increase_fixed))) sd(abs_increase_fixed, na.rm = TRUE) else NA
      t_y_fixed_sd <- if (!all(is.na(t_y_fixed))) sd(t_y_fixed, na.rm = TRUE) else NA
      t_y_star_fixed_sd <- if (!all(is.na(t_y_star_fixed))) sd(t_y_star_fixed, na.rm = TRUE) else NA

      # Create per-period tibbles with period as date
      abs_increase_chain_df <- dplyr::tibble(period = period, value = round(abs_increase_chain, 4))
      t_y_chain_df <- dplyr::tibble(period = period, value = round(t_y_chain, 4))      abs_increase_chain_mean <- if (!all(is.na(abs_increase_chain))) mean(abs_increase_chain, na.rm = TRUE) else NA
      t_y_star_chain_df <- dplyr::tibble(period = period, value = round(t_y_star_chain, 4))ean(t_y_chain, na.rm = TRUE) else NA
      abs_increase_fixed_df <- dplyr::tibble(period = period, value = round(abs_increase_fixed, 4))se NA
      t_y_fixed_df <- dplyr::tibble(period = period, value = round(t_y_fixed, 4))
      t_y_star_fixed_df <- dplyr::tibble(period = period, value = round(t_y_star_fixed, 4))
m = TRUE) else NA
      # Create a wide-format tibble with simplified lowercase names else NA
      descriptive <- dplyr::tibble(TRUE) else NA
        autocorr_r1 = r1,
        durbinwatson_qbp = q_bp,
        ljungbox_pval = ljung_box_p,!all(is.na(abs_increase_chain))) sd(abs_increase_chain, na.rm = TRUE) else NA
        spearman_corr = spearman[1],(!all(is.na(t_y_chain))) sd(t_y_chain, na.rm = TRUE) else NA
        spearman_pval = spearman[2],all(is.na(t_y_star_chain))) sd(t_y_star_chain, na.rm = TRUE) else NA
        anderson_stat = ad_test[1],all(is.na(abs_increase_fixed))) sd(abs_increase_fixed, na.rm = TRUE) else NA
        anderson_pval = ad_test[2],(t_y_fixed))) sd(t_y_fixed, na.rm = TRUE) else NA
        mean_absinc_chain = abs_increase_chain_mean,is.na(t_y_star_fixed))) sd(t_y_star_fixed, na.rm = TRUE) else NA
        mean_devrate_chain = t_y_chain_mean,
        mean_incrate_chain = t_y_star_chain_mean,ith period as date
        mean_absinc_fixed = abs_increase_fixed_mean,= period, value = round(abs_increase_chain, 4))
        mean_devrate_fixed = t_y_fixed_mean, period, value = round(t_y_chain, 4))
        mean_incrate_fixed = t_y_star_fixed_mean, period, value = round(t_y_star_chain, 4))
        sd_absinc_chain = abs_increase_chain_sd,= period, value = round(abs_increase_fixed, 4))
        sd_devrate_chain = t_y_chain_sd, period, value = round(t_y_fixed, 4))
        sd_incrate_chain = t_y_star_chain_sd, period, value = round(t_y_star_fixed, 4))
        sd_absinc_fixed = abs_increase_fixed_sd,
        sd_devrate_fixed = t_y_fixed_sd, simplified lowercase names
        sd_incrate_fixed = t_y_star_fixed_sd,
        geom_mean_growth = geom_mean_growth,
        geom_mean_growth_pct = geom_mean_growth_pct,
        mean_incrate = mean_incrate,
        mean_incrate_pct = mean_incrate_pct,
        per_absinc_chain = list(abs_increase_chain_df),
        per_devrate_chain = list(t_y_chain_df),
        per_incrate_chain = list(t_y_star_chain_df),
        per_absinc_fixed = list(abs_increase_fixed_df),
        per_devrate_fixed = list(t_y_fixed_df),
        per_incrate_fixed = list(t_y_star_fixed_df)
      ) %>%
        dplyr::mutate_if(is.numeric, ~ round(., 4))

      return(descriptive)absinc_chain = abs_increase_chain_sd,
    }, .keep = TRUE) %>%
    dplyr::bind_rows(.id = "group")        sd_incrate_chain = t_y_star_chain_sd,
 abs_increase_fixed_sd,
  return(result) = t_y_fixed_sd,
}_fixed_sd,
        per_absinc_chain = list(abs_increase_chain_df),
        per_devrate_chain = list(t_y_chain_df),
        per_incrate_chain = list(t_y_star_chain_df),
        per_absinc_fixed = list(abs_increase_fixed_df),
        per_devrate_fixed = list(t_y_fixed_df),
        per_incrate_fixed = list(t_y_star_fixed_df)
      ) %>%
        dplyr::mutate_if(is.numeric, ~ round(., 4))

      return(descriptive)
    }, .keep = TRUE) %>%
    dplyr::bind_rows(.id = "group")

  return(result)
}
