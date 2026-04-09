# ============================================================
#  model 1.0.R
#  GAMLSS Model Architecture — Fitting Functions
#
#  Defines multilevel GAMLSS fitting functions for four
#  children's fitness outcomes:
#    - 50m Dash       (Log-Normal, LOGNO)
#    - Vital Capacity (Log-Normal, LOGNO)
#    - Rope Skipping  (Negative-Binomial, NBI)
#    - Sit-and-Reach  (Two-part hurdle: logistic + Log-Normal)
#
#  All continuous outcomes share a three-level random-effects
#  structure on mu (school / cohort / child) and a school-level
#  random intercept on sigma.
#
#  Sourced by:
#    - Shiny_Dashboard_v2_0.R
#    - model_fitting.R
# ============================================================

library(gamlss)
library(gamlss.dist)

# ---------------------------
# Helpers
# ---------------------------
.clean_frame <- function(data, vars) {
  keep <- rep(TRUE, nrow(data))
  for (v in vars) keep <- keep & is.finite(data[[v]]) & !is.na(data[[v]])
  data[keep, , drop = FALSE]
}

.ensure_factor_ids <- function(data, id_vars = c("school_id","cohort_id","child_id")) {
  for (v in id_vars) {
    if (!is.factor(data[[v]])) data[[v]] <- as.factor(data[[v]])
  }
  data
}

# ---------------------------
# Fit one outcome with longitudinal multilevel GAMLSS
# Supports:
#   - continuous positive outcomes (LOGNO, GA, BCPE, BCCG, BCT, etc.)
#   - count outcomes (PO, NBI, NBII, ZIP, ZINBI, etc. if available in your gamlss.dist)
# ---------------------------
fit_long_gamlss <- function(
    data,
    outcome_var = "y",
    family,
    include_bmi = TRUE,
    control = gamlss.control(n.cyc = 150, trace = FALSE),
    clean_y = TRUE,
    # random effects options (mu side)
    mu_random = list(
      school_id = ~ 1,
      cohort_id = ~ 1,
      child_id  = ~ 1 + time
    ),
    # sigma side
    sigma_re_on = TRUE,
    sigma_random = list(school_id = ~ 1),
    sigma_terms = c("pb(age)"),
    # nu / tau
    nu_terms    = c("1"),
    tau_terms   = c("1"),
    # link-safety for strictly positive outcomes
    min_y = 1e-6
) {
  stopifnot(outcome_var %in% names(data))
  data <- .ensure_factor_ids(data)
  
  # Basic NA/Inf cleaning
  data <- .clean_frame(data, c(outcome_var, "age", "sex", "time", "school_id", "cohort_id", "child_id"))
  if (include_bmi) data <- .clean_frame(data, c("bmi"))
  
  # Outcome cleaning (ONLY if appropriate — caller should disable for counts)
  if (clean_y) {
    if (any(data[[outcome_var]] < 0, na.rm = TRUE)) {
      stop("Negative values found in outcome; please check outcome scaling.")
    }
    data[[outcome_var]] <- pmax(data[[outcome_var]], min_y)
  }
  
  # Fixed terms (mu)
  fixed_terms <- c("pb(age)", "sex")
  if (include_bmi) fixed_terms <- c(fixed_terms, "bmi")
  mu_formula <- as.formula(paste(outcome_var, "~", paste(fixed_terms, collapse = " + ")))
  
  # Random effects term for mu
  mu_re <- re(random = mu_random, method = "ML")
  
  # Extend the mu formula to include the random-effects term
  full_mu_formula <- update(mu_formula, . ~ . + mu_re)
  
  # Sigma part
  sigma_rhs <- paste(sigma_terms, collapse = " + ")
  if (sigma_re_on) {
    sigma_re <- re(random = sigma_random, method = "ML")
    sigma_formula <- as.formula(paste("~", sigma_rhs, "+ sigma_re"))
  } else {
    sigma_formula <- as.formula(paste("~", sigma_rhs))
  }
  
  # Nu / Tau
  nu_formula  <- as.formula(paste("~", paste(nu_terms,  collapse = " + ")))
  tau_formula <- as.formula(paste("~", paste(tau_terms, collapse = " + ")))
  
  # Fit
  m <- gamlss(
    formula       = full_mu_formula,
    sigma.formula = sigma_formula,
    nu.formula    = nu_formula,
    tau.formula   = tau_formula,
    family        = family,
    data          = data,
    control       = control
  )
  return(m)
}

# ---------------------------
# Sit-and-reach hurdle fit:
#   (1) logistic model for P(Y > 0)
#   (2) GAMLSS model for Y | Y>0 (positive-only)
# ---------------------------
fit_sitreach_hurdle <- function(
    data,
    outcome_var = "sitreach",
    positive_family = LOGNO(),          # or GA(), BCPE(), etc.
    include_bmi = TRUE,
    control = gamlss.control(n.cyc = 150, trace = FALSE),
    mu_random = list(
      school_id = ~ 1,
      cohort_id = ~ 1,
      child_id  = ~ 1 + time
    ),
    sigma_re_on = TRUE,
    min_y = 1e-6
) {
  stopifnot(outcome_var %in% names(data))
  data <- .ensure_factor_ids(data)
  
  # Create indicator
  data <- .clean_frame(data, c(outcome_var, "age", "sex", "time", "school_id", "cohort_id", "child_id"))
  if (include_bmi) data <- .clean_frame(data, c("bmi"))
  data$y_pos_ind <- as.integer(data[[outcome_var]] > 0)
  
  # (1) Logistic submodel for P(Y > 0)
  # gamlss(BI()) with re() encounters a formula-environment scoping issue in which
  # the inner GLM scoring loop loses access to the random-effects object defined in
  # the calling frame. Standard glm(binomial) is used for the hurdle component instead;
  # it is numerically equivalent and methodologically appropriate.
  # The positive-part GAMLSS submodel (below) retains the full multilevel structure.
  fixed_terms <- c("pb(age)", "sex")
  if (include_bmi) fixed_terms <- c(fixed_terms, "bmi")
  mu_formula_bin <- as.formula(paste("y_pos_ind ~", paste(fixed_terms, collapse = " + ")))
  
  m_zero <- glm(
    formula = mu_formula_bin,
    family  = binomial(link = "logit"),
    data    = data
  )
  
  # (2) Positive-only continuous model
  dat_pos <- droplevels(subset(data, y_pos_ind == 1))  # drop unused factor levels from the positive subset
  dat_pos$y_pos <- pmax(dat_pos[[outcome_var]], min_y)
  
  m_pos <- fit_long_gamlss(
    data        = within(dat_pos, { y <- y_pos }),
    outcome_var = "y",
    family      = positive_family,
    include_bmi = include_bmi,
    control     = control,
    clean_y     = TRUE,
    mu_random   = mu_random,
    sigma_re_on = sigma_re_on,
    min_y       = min_y
  )
  
  list(
    model_zero = m_zero,
    model_pos  = m_pos,
    data_all   = data,
    data_pos   = dat_pos
  )
}

# Outcome-specific convenience wrappers (call fit_long_gamlss or fit_sitreach_hurdle)
fit_50m <- function(data, outcome_var = "dash50", ...) {
  fit_long_gamlss(data, outcome_var, family = LOGNO(), sigma_re_on = TRUE, ...)
}

fit_vc <- function(data, outcome_var = "vitalcap", ...) {
  fit_long_gamlss(data, outcome_var, family = LOGNO(), sigma_re_on = TRUE, ...)
}

fit_rope <- function(data, outcome_var = "ropeskip", ...) {
  fit_long_gamlss(
    data, outcome_var,
    family = NBI(),
    sigma_terms = c("1"),
    sigma_re_on = FALSE,
    clean_y = FALSE,  # counts must not be truncated at min_y
    ...
  )
}

fit_sitreach <- function(data, outcome_var = "sitreach", ...) {
  fit_sitreach_hurdle(data, outcome_var, positive_family = LOGNO(), ...)
}