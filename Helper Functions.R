# ============================================================
#  Helper Functions
#
#  Provides functions for extracting model diagnostics after
#  GAMLSS fitting, including:
#    - Convergence and stability checks
#    - Quantile residual summaries
#    - Prediction interval coverage
#    - Centile calibration
#    - RMSE against simulation truth (Monte Carlo pipeline only)
#    - Data quality metrics (floor rates, heaping, missingness)
#    - Full single-dataset fitting (fit_all_tests_one_dataset)
#    - Monte Carlo replicate runners (run_one_replicate,
#      run_mc_for_scenario)
#
#  Sourced by:
#    - Shiny Dashboard.R  (diagnostic helpers)
#
# ============================================================

suppressPackageStartupMessages({
  library(gamlss)
  library(dplyr)
})

# ============================================================
# 1) DIAGNOSTICS HELPERS (compact tables; no model storage needed)
# ============================================================

.safe <- function(expr) tryCatch(expr, error = function(e) NA)

stability_flag <- function(m) {
  if (is.null(m)) return(FALSE)
  # glm objects (sitreach hurdle zero-part) use $converged rather than $iter/$control
  if (inherits(m, "glm") && !inherits(m, "gamlss")) {
    ok_conv <- isTRUE(m$converged)
    ok_aic  <- is.finite(.safe(AIC(m)))
    ok_coef <- all(is.finite(coef(m)))
    return(isTRUE(ok_conv && ok_aic && ok_coef))
  }
  ok_iter <- is.null(m$iter) || is.null(m$control$n.cyc) || (m$iter < m$control$n.cyc) 
  ok_aic  <- is.finite(.safe(AIC(m)))
  ok_coef <- TRUE
  cf <- tryCatch(coefAll(m), error = function(e) NULL)
  if (!is.null(cf)) ok_coef <- all(is.finite(unlist(cf)))
  isTRUE(ok_iter && ok_aic && ok_coef)
}

get_qresid <- function(m) {
  if (is.null(m)) return(NULL)
  tryCatch(as.numeric(resid(m, type="quantile")), error = function(e) NULL)
}

resid_summary <- function(r) {
  if (is.null(r)) return(NULL)
  r <- r[is.finite(r)]
  if (length(r) < 10) return(NULL)
  data.frame(
    mean = mean(r), sd = sd(r),
    q05 = unname(quantile(r, 0.05)),
    q50 = unname(quantile(r, 0.50)),
    q95 = unname(quantile(r, 0.95))
  )
}

make_resid_row <- function(outcome, m, scenario, seed) {
  r <- get_qresid(m)
  s <- resid_summary(r)
  if (is.null(s) || nrow(s) == 0) return(NULL)
  out <- cbind(data.frame(outcome = outcome, stringsAsFactors = FALSE), s)
  out$scenario <- scenario
  out$seed <- seed
  out
}

tag_tbl <- function(tbl, scenario, seed, outcome) {
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  tbl$scenario <- scenario
  tbl$seed <- seed
  tbl$outcome <- outcome
  tbl
}

# ---- PI coverage via qFAM if available ----
get_pi_coverage <- function(m, y, levels = c(0.80, 0.90, 0.95)) {
  if (is.null(m)) return(NULL)
  fam <- m$family[1]
  qfun <- get0(paste0("q", fam), mode="function")
  if (is.null(qfun)) return(NULL)
  
  mu  <- fitted(m, "mu")
  sig <- fitted(m, "sigma")
  nu  <- tryCatch(fitted(m, "nu"),  error=function(e) NULL)
  tau <- tryCatch(fitted(m, "tau"), error=function(e) NULL)
  
  out <- lapply(levels, function(lev) {
    a <- (1-lev)/2
    lo_p <- a; hi_p <- 1-a
    
    lo <- if (!is.null(tau) && !is.null(nu)) {
      tryCatch(qfun(lo_p, mu=mu, sigma=sig, nu=nu, tau=tau), error=function(e) NULL)
    } else if (!is.null(nu)) {
      tryCatch(qfun(lo_p, mu=mu, sigma=sig, nu=nu), error=function(e) NULL)
    } else {
      tryCatch(qfun(lo_p, mu=mu, sigma=sig), error=function(e) NULL)
    }
    
    hi <- if (!is.null(tau) && !is.null(nu)) {
      tryCatch(qfun(hi_p, mu=mu, sigma=sig, nu=nu, tau=tau), error=function(e) NULL)
    } else if (!is.null(nu)) {
      tryCatch(qfun(hi_p, mu=mu, sigma=sig, nu=nu), error=function(e) NULL)
    } else {
      tryCatch(qfun(hi_p, mu=mu, sigma=sig), error=function(e) NULL)
    }
    
    if (is.null(lo) || is.null(hi)) return(NULL)
    data.frame(
      level = lev,
      coverage = mean(y >= lo & y <= hi, na.rm=TRUE),
      avg_width = mean(hi - lo, na.rm=TRUE)
    )
  })
  do.call(rbind, out)
}

# ---- PI coverage for LOGNO (lognormal): avoids qFAM lookup issues ----
get_pi_coverage_lognormal <- function(mu_hat, sigma_hat, y, levels = c(0.80, 0.90, 0.95)) {
  if (is.null(mu_hat) || is.null(sigma_hat) || is.null(y)) return(NULL)
  if (length(mu_hat) != length(sigma_hat) || length(mu_hat) != length(y)) return(NULL)
  
  out <- lapply(levels, function(lev) {
    a <- (1 - lev) / 2
    z_lo <- qnorm(a)
    z_hi <- qnorm(1 - a)
    
    lo <- exp(mu_hat + z_lo * sigma_hat)
    hi <- exp(mu_hat + z_hi * sigma_hat)
    
    data.frame(
      level     = lev,
      coverage  = mean(y >= lo & y <= hi, na.rm = TRUE),
      avg_width = mean(hi - lo, na.rm = TRUE)
    )
  })
  do.call(rbind, out)
}

# Centile calibration for the sitreach positive-part lognormal model.
# Returns (centile, target, empirical, diff) — same structure as get_centile_calibration()
# so the two can be row-bound into mc_diag_centile without column mismatches.
get_centile_calibration_lognormal <- function(mu_hat, sigma_hat, y,
                                              centiles = c(5, 10, 25, 50, 75, 90, 95)) {
  if (is.null(mu_hat) || is.null(sigma_hat) || is.null(y)) return(NULL)
  if (length(mu_hat) != length(sigma_hat) || length(mu_hat) != length(y)) return(NULL)
  out <- lapply(centiles, function(cen) {
    p    <- cen / 100
    qhat <- exp(mu_hat + qnorm(p) * sigma_hat)
    emp  <- mean(y <= qhat, na.rm = TRUE)
    data.frame(centile = cen, target = p, empirical = emp, diff = emp - p)
  })
  do.call(rbind, out)
}

get_centile_calibration <- function(m, y, centiles = c(5,10,25,50,75,90,95)) {
  if (is.null(m)) return(NULL)
  fam <- m$family[1]
  qfun <- get0(paste0("q", fam), mode="function")
  if (is.null(qfun)) return(NULL)
  
  mu  <- fitted(m, "mu")
  sig <- fitted(m, "sigma")
  nu  <- tryCatch(fitted(m, "nu"),  error=function(e) NULL)
  tau <- tryCatch(fitted(m, "tau"), error=function(e) NULL)
  
  out <- lapply(centiles, function(cen) {
    p <- cen/100
    qhat <- if (!is.null(tau) && !is.null(nu)) {
      tryCatch(qfun(p, mu=mu, sigma=sig, nu=nu, tau=tau), error=function(e) NULL)
    } else if (!is.null(nu)) {
      tryCatch(qfun(p, mu=mu, sigma=sig, nu=nu), error=function(e) NULL)
    } else {
      tryCatch(qfun(p, mu=mu, sigma=sig), error=function(e) NULL)
    }
    if (is.null(qhat)) return(NULL)
    emp <- mean(y <= qhat, na.rm=TRUE)
    data.frame(centile=cen, target=p, empirical=emp, diff=emp-p)
  })
  do.call(rbind, out)
}

# ---- Scenario sensitivity metrics ----
heaping_score <- function(x, step) {
  x <- x[is.finite(x)]
  if (length(x)==0) return(NA_real_)
  mean(abs(x/step - round(x/step)) < 1e-8)
}

quality_metrics <- function(df) {
  data.frame(
    n_rows = nrow(df),
    n_children = length(unique(df$child_id)),
    sitreach_floor = mean(df$sitreach == 0, na.rm=TRUE),
    rope_zero = mean(df$ropeskip == 0, na.rm=TRUE),
    dash50_heaping_0.1 = heaping_score(df$dash50, 0.1),
    vc_heaping_50 = heaping_score(df$vitalcap, 50)
  )
}

missingness_profile <- function(df) {
  tab <- with(df, table(child_id, time))
  prop <- colMeans(tab > 0)
  data.frame(
    time = as.integer(names(prop)),
    prop_observed_children = as.numeric(prop)
  )
}

# ============================================================
# 2) RMSE (requires truth columns saved in simulated data)
# ============================================================

rmse <- function(est, tru) {
  ok <- is.finite(est) & is.finite(tru)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((est[ok] - tru[ok])^2))
}

get_rmse_table <- function(df, fit_res, scenario, seed) {
  out <- list()
  
  # dash50 (LOGNO): mu=meanlog, sigma=sdlog
  m <- fit_res$dash50$model
  if (!is.null(m) && "dash50_mu_log_true" %in% names(df)) {
    out[[length(out)+1]] <- data.frame(
      scenario=scenario, seed=seed, outcome="dash50", param="mu",
      rmse=rmse(fitted(m, "mu"), df$dash50_mu_log_true)
    )
  }
  if (!is.null(m) && "dash50_sigma_log_true" %in% names(df)) {
    out[[length(out)+1]] <- data.frame(
      scenario=scenario, seed=seed, outcome="dash50", param="sigma",
      rmse=rmse(fitted(m, "sigma"), df$dash50_sigma_log_true)
    )
  }
  
  # vitalcap (LOGNO)
  m <- fit_res$vitalcap$model
  if (!is.null(m) && "vitalcap_mu_log_true" %in% names(df)) {
    out[[length(out)+1]] <- data.frame(
      scenario=scenario, seed=seed, outcome="vitalcap", param="mu",
      rmse=rmse(fitted(m, "mu"), df$vitalcap_mu_log_true)
    )
  }
  if (!is.null(m) && "vitalcap_sigma_log_true" %in% names(df)) {
    out[[length(out)+1]] <- data.frame(
      scenario=scenario, seed=seed, outcome="vitalcap", param="sigma",
      rmse=rmse(fitted(m, "sigma"), df$vitalcap_sigma_log_true)
    )
  }
  
  # ropeskip (NBI): mu + sigma
  m <- fit_res$ropeskip$model
  if (!is.null(m) && "ropeskip_mu_true" %in% names(df)) {
    out[[length(out)+1]] <- data.frame(
      scenario=scenario, seed=seed, outcome="ropeskip", param="mu",
      rmse=rmse(fitted(m, "mu"), df$ropeskip_mu_true)
    )
  }
  if (!is.null(m) && "ropeskip_sigma_true" %in% names(df)) {
    out[[length(out)+1]] <- data.frame(
      scenario=scenario, seed=seed, outcome="ropeskip", param="sigma",
      rmse=rmse(fitted(m, "sigma"), df$ropeskip_sigma_true)
    )
  }
  
  # sitreach hurdle:
  # zero part: compare fitted P(Y>0) to p_pos_true
  # m_zero is now a glm object; use predict(type="response") not fitted(m0, "mu")
  m0 <- fit_res$sitreach$model$model_zero
  if (!is.null(m0) && "sitreach_p_pos_true" %in% names(df)) {
    p_hat <- if (inherits(m0, "glm") && !inherits(m0, "gamlss")) {
      tryCatch(predict(m0, type = "response"), error = function(e) NULL)
    } else {
      tryCatch(fitted(m0, "mu"), error = function(e) NULL)
    }
    if (!is.null(p_hat)) {
      out[[length(out)+1]] <- data.frame(
        scenario=scenario, seed=seed, outcome="sitreach_zero", param="p_pos",
        rmse=rmse(p_hat, df$sitreach_p_pos_true)
      )
    }
  }
  
  # positive part: mu/sigma on log scale (LOGNO)
  mp <- fit_res$sitreach$model$model_pos
  dp <- fit_res$sitreach$model$data_pos
  
  if (!is.null(mp)) {
    
    # MUST use dp (positive-only subset), not df.
    # fitted(mp, "mu") has length == nrow(dp), while df has all rows including zeros.
    # Using df causes a silent length mismatch that makes rmse() return NA.
    mu_true  <- if (!is.null(dp) && "sitreach_mu_pos_log_true" %in% names(dp)) dp$sitreach_mu_pos_log_true else NULL
    sig_true <- if (!is.null(dp) && "sitreach_sigma_pos_true"  %in% names(dp)) dp$sitreach_sigma_pos_true  else NULL
    
    if (!is.null(mu_true)) {
      out[[length(out)+1]] <- data.frame(
        scenario=scenario, seed=seed, outcome="sitreach_pos", param="mu",
        rmse=rmse(fitted(mp, "mu"), mu_true)
      )
    }
    
    if (!is.null(sig_true)) {
      out[[length(out)+1]] <- data.frame(
        scenario=scenario, seed=seed, outcome="sitreach_pos", param="sigma",
        rmse=rmse(fitted(mp, "sigma"), sig_true)
      )
    }
  }
  
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

# ============================================================
# 3) EXTRACT ALL DIAGNOSTICS FOR ONE FIT
# ============================================================

extract_diagnostics <- function(df, fit_res, scenario, seed) {
  
  # stability (5 rows)
  stab <- rbind(
    data.frame(outcome="dash50",        stable=stability_flag(fit_res$dash50$model),   AIC=.safe(AIC(fit_res$dash50$model))),
    data.frame(outcome="vitalcap",      stable=stability_flag(fit_res$vitalcap$model), AIC=.safe(AIC(fit_res$vitalcap$model))),
    data.frame(outcome="ropeskip",      stable=stability_flag(fit_res$ropeskip$model), AIC=.safe(AIC(fit_res$ropeskip$model))),
    data.frame(outcome="sitreach_zero", stable=stability_flag(fit_res$sitreach$model$model_zero), AIC=.safe(AIC(fit_res$sitreach$model$model_zero))),
    data.frame(outcome="sitreach_pos",  stable=stability_flag(fit_res$sitreach$model$model_pos),  AIC=.safe(AIC(fit_res$sitreach$model$model_pos)))
  )
  stab$scenario <- scenario
  stab$seed <- seed
  
  # residual summaries (4 rows; sitreach uses pos-part)
  rs <- do.call(rbind, Filter(Negate(is.null), list(
    make_resid_row("dash50",       fit_res$dash50$model, scenario, seed),
    make_resid_row("vitalcap",     fit_res$vitalcap$model, scenario, seed),
    make_resid_row("ropeskip",     fit_res$ropeskip$model, scenario, seed),
    make_resid_row("sitreach_pos", fit_res$sitreach$model$model_pos, scenario, seed)
  )))
  
  # PI coverage
  pi <- do.call(rbind, Filter(Negate(is.null), list(
    tag_tbl(get_pi_coverage(fit_res$dash50$model, df$dash50), scenario, seed, "dash50"),
    tag_tbl(get_pi_coverage(fit_res$vitalcap$model, df$vitalcap), scenario, seed, "vitalcap"),
    tag_tbl(get_pi_coverage(fit_res$ropeskip$model, df$ropeskip), scenario, seed, "ropeskip"),
    tag_tbl(
      get_pi_coverage_lognormal(
        mu_hat    = fitted(fit_res$sitreach$model$model_pos, "mu"),
        sigma_hat = fitted(fit_res$sitreach$model$model_pos, "sigma"),
        y         = fit_res$sitreach$model$data_pos$y_pos
      ),
      scenario, seed, "sitreach_pos"
    )
  )))
  
  # centile calibration
  cent <- do.call(rbind, Filter(Negate(is.null), list(
    tag_tbl(get_centile_calibration(fit_res$dash50$model, df$dash50), scenario, seed, "dash50"),
    tag_tbl(get_centile_calibration(fit_res$vitalcap$model, df$vitalcap), scenario, seed, "vitalcap"),
    tag_tbl(get_centile_calibration(fit_res$ropeskip$model, df$ropeskip), scenario, seed, "ropeskip"),
    tag_tbl(
      get_centile_calibration_lognormal(
        mu_hat    = fitted(fit_res$sitreach$model$model_pos, "mu"),
        sigma_hat = fitted(fit_res$sitreach$model$model_pos, "sigma"),
        y         = fit_res$sitreach$model$data_pos$y_pos
      ),
      scenario, seed, "sitreach_pos"
    )
  )))
  
  # data quality + missingness
  dq <- quality_metrics(df); dq$scenario <- scenario; dq$seed <- seed
  miss <- missingness_profile(df); miss$scenario <- scenario; miss$seed <- seed
  
  # RMSE
  rmse_tbl <- get_rmse_table(df, fit_res, scenario, seed)
  
  list(
    stability = stab,
    resid     = rs,
    pi        = pi,
    centile   = cent,
    quality   = dq,
    missing   = miss,
    rmse      = rmse_tbl
  )
}

# ============================================================
# 4) FIT ALL OUTCOMES FOR ONE DATASET
# ============================================================

fit_all_tests_one_dataset <- function(df,
                                      include_bmi = TRUE,
                                      control = gamlss.control(n.cyc = 150, trace = FALSE)) {
  out <- list()
  
  # helper to safely fit and capture warnings/errors + runtime
  safe_fit <- function(expr) {
    warn <- character(0)
    err <- NULL
    t0 <- proc.time()[3]
    
    fit <- withCallingHandlers(
      tryCatch(expr, error = function(e) { err <<- conditionMessage(e); NULL }),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    t1 <- proc.time()[3]
    list(model = fit, error = err, warnings = warn, runtime_sec = (t1 - t0))
  }
  
  out$dash50   <- safe_fit(fit_50m(df, outcome_var = "dash50",
                                   include_bmi = include_bmi, control = control))
  
  out$vitalcap <- safe_fit(fit_vc(df, outcome_var = "vitalcap",
                                  include_bmi = include_bmi, control = control))
  
  out$ropeskip <- safe_fit(fit_rope(df, outcome_var = "ropeskip",
                                    include_bmi = include_bmi, control = control))
  
  out$sitreach <- safe_fit(fit_sitreach(df, outcome_var = "sitreach",
                                        include_bmi = include_bmi, control = control))
  
  # pragmatic "convergence" check
  conv_flag <- function(fit_obj) {
    if (is.null(fit_obj)) return(NA)
    # glm objects (sitreach zero part) use $converged
    if (inherits(fit_obj, "glm") && !inherits(fit_obj, "gamlss")) {
      return(isTRUE(fit_obj$converged))
    }
    if (!is.null(fit_obj$iter) && !is.null(fit_obj$control$n.cyc)) {
      return(fit_obj$iter < fit_obj$control$n.cyc)
    }
    NA
  }
  
  out$summary <- data.frame(
    outcome = c("dash50","vitalcap","ropeskip","sitreach_zero","sitreach_pos"),
    converged = c(
      conv_flag(out$dash50$model),
      conv_flag(out$vitalcap$model),
      conv_flag(out$ropeskip$model),
      conv_flag(out$sitreach$model$model_zero),
      conv_flag(out$sitreach$model$model_pos)
    ),
    runtime_sec = c(
      out$dash50$runtime_sec,
      out$vitalcap$runtime_sec,
      out$ropeskip$runtime_sec,
      out$sitreach$runtime_sec,
      out$sitreach$runtime_sec
    ),
    n = c(
      nrow(df),
      nrow(df),
      nrow(df),
      nrow(df),
      sum(df$sitreach > 0, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  out
}

# ============================================================
# 5) ONE REPLICATE: simulate -> scenario -> fit -> diagnostics
# ============================================================

run_one_replicate <- function(seed,
                              scenario = "S1",
                              design = list(J=10, H=60, n=50, T=3),
                              ICC = list(ICC_school=0.01, ICC_cohort=0.05, ICC_child=0.50),
                              total_var_eta = 1.0,
                              rope_zero_infl_baseline = FALSE,
                              control = gamlss.control(n.cyc = 150, trace = FALSE)) {
  
  df <- simulate_fitness_baseline(
    J = design$J, H = design$H, n = design$n, T = design$T,
    ICC_school = ICC$ICC_school, ICC_cohort = ICC$ICC_cohort, ICC_child = ICC$ICC_child,
    total_var_eta = total_var_eta,
    seed = seed,
    rope_zero_infl = rope_zero_infl_baseline
  )
  
  df <- apply_scenarios_fitness(df, scenario = scenario, seed = seed)
  
  fit_res <- fit_all_tests_one_dataset(df, control = control)
  
  fit_res$meta <- list(seed = seed, scenario = scenario,
                       design = design, ICC = ICC, total_var_eta = total_var_eta)
  
  fit_res$diag <- extract_diagnostics(df, fit_res, scenario = scenario, seed = seed)
  
  # optionally free the big dataset explicitly (not required but helps when iterating)
  rm(df); gc(verbose = FALSE)
  
  fit_res
}

# ============================================================
# 6) MC RUNNERS (scenario + all scenarios) — SUMMARY ONLY
# ============================================================

run_mc_for_scenario <- function(scenario, R=20, seed0=1000,
                                control=gamlss.control(n.cyc=150, trace=FALSE)) {
  
  summary_rows <- vector("list", R)
  stab_rows <- vector("list", R)
  resid_rows <- vector("list", R)
  pi_rows <- vector("list", R)
  cent_rows <- vector("list", R)
  dq_rows <- vector("list", R)
  miss_rows <- vector("list", R)
  rmse_rows <- vector("list", R)
  
  for (r in seq_len(R)) {
    seed <- seed0 + r
    res <- run_one_replicate(seed = seed, scenario = scenario, control = control)
    
    s <- res$summary
    s$replicate <- r; s$seed <- seed; s$scenario <- scenario
    summary_rows[[r]] <- s
    
    stab_rows[[r]]  <- res$diag$stability
    resid_rows[[r]] <- res$diag$resid
    pi_rows[[r]]    <- res$diag$pi
    cent_rows[[r]]  <- res$diag$centile
    dq_rows[[r]]    <- res$diag$quality
    miss_rows[[r]]  <- res$diag$missing
    rmse_rows[[r]]  <- res$diag$rmse
  }
  
  list(
    scenario = scenario,
    summary = do.call(rbind, summary_rows),
    diag_stability = do.call(rbind, stab_rows),
    diag_resid = do.call(rbind, Filter(Negate(is.null), resid_rows)),
    diag_pi = do.call(rbind, Filter(Negate(is.null), pi_rows)),
    diag_centile = do.call(rbind, Filter(Negate(is.null), cent_rows)),
    diag_quality = do.call(rbind, dq_rows),
    diag_missing = do.call(rbind, miss_rows),
    diag_rmse = do.call(rbind, Filter(Negate(is.null), rmse_rows))
  )
}