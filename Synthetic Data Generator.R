# ============================================================
#  Synthetic Data Generator — Baseline Design and Stress-Test Scenarios
#
#  Generates longitudinal children's fitness data under a three-level
#  hierarchical structure (school / cohort / child) with a shared latent
#  fitness predictor. Seven stress-test scenarios perturb the baseline
#  to probe model robustness under realistic data-quality violations.
#
#  Default design:
#    J = 10 schools,  H = 60 cohorts per school,
#    n = 50 children per cohort,  T = 3 annual measurement waves
#
# ============================================================

library(gamlss.dist)

# ============================================================
# 1) Design + covariates
# ============================================================

make_design <- function(J=10, H=60, n=50, T=3, seed=1) {
  set.seed(seed)
  
  school_id <- rep(seq_len(J), each = H * n)
  cohort_in_school <- rep(seq_len(H), times = J, each = n)
  cohort_id <- interaction(school_id, cohort_in_school, drop = TRUE)
  
  child_in_cohort <- rep(seq_len(n), times = J * H)
  child_id <- interaction(cohort_id, child_in_cohort, drop = TRUE)
  
  base <- data.frame(
    school_id = factor(school_id),
    cohort_id = factor(cohort_id),
    child_id  = factor(child_id)
  )
  
  long <- base[rep(seq_len(nrow(base)), each = T), ]
  long$time <- rep(0:(T-1), times = nrow(base))
  long
}

simulate_covariates <- function(dat, seed=1) {
  set.seed(seed + 10)
  sex_child <- rbinom(nlevels(dat$child_id), 1, 0.5)
  dat$sex <- sex_child[as.integer(dat$child_id)]
  
  set.seed(seed + 11)
  age0_child <- runif(nlevels(dat$child_id), min=6, max=10)
  dat$age <- age0_child[as.integer(dat$child_id)] + dat$time
  
  dat
}

# BMI simulator (quadratic + clustering)
simulate_bmi_quadratic <- function(dat,
                                   mu_9 = 17.0,
                                   beta_age = 0.45,
                                   beta_age2 = 0.08,
                                   beta_sex = 0.35,
                                   sd_school = 0.20,
                                   sd_cohort = 0.12,
                                   sd_child  = 0.85,
                                   sd_eps_9 = 1.10,
                                   sd_eps_age = 0.05,
                                   clamp = TRUE,
                                   clamp_min = 10,
                                   clamp_max = 35,
                                   seed = 1) {
  
  stopifnot(all(c("age","sex","school_id","cohort_id","child_id") %in% names(dat)))
  age_c <- dat$age - 9
  
  set.seed(seed + 1)
  u_school <- rnorm(nlevels(dat$school_id), 0, sd_school)
  set.seed(seed + 2)
  v_cohort <- rnorm(nlevels(dat$cohort_id), 0, sd_cohort)
  set.seed(seed + 3)
  b_child <- rnorm(nlevels(dat$child_id), 0, sd_child)
  
  sd_eps <- sd_eps_9 * (1 + sd_eps_age * abs(age_c))
  set.seed(seed + 4)
  eps <- rnorm(nrow(dat), 0, sd_eps)
  
  mu_fixed <- mu_9 + beta_age * age_c + beta_age2 * (age_c^2) + beta_sex * dat$sex
  
  bmi <- mu_fixed +
    u_school[as.integer(dat$school_id)] +
    v_cohort[as.integer(dat$cohort_id)] +
    b_child[as.integer(dat$child_id)] +
    eps
  
  if (clamp) bmi <- pmin(pmax(bmi, clamp_min), clamp_max)
  as.numeric(bmi)
}

# ============================================================
# 2) Shared latent structure (hierarchical longitudinal predictor)
# NOTE: ICC targets refer to baseline (time=0) latent scale partitioning.
# ============================================================

simulate_latent_predictors <- function(
    dat,
    ICC_school=0.01, ICC_cohort=0.05, ICC_child=0.50,
    total_var_eta = 1.0,          # variance of latent eta at baseline (time=0)
    var_s = 0.10, cor_bs = -0.20, # child slope variability
    seed=1
) {
  var_u <- ICC_school * total_var_eta
  var_v <- (ICC_cohort * total_var_eta) - var_u
  var_b <- (ICC_child  * total_var_eta) - var_u - var_v
  var_e <- total_var_eta - var_u - var_v - var_b
  if (any(c(var_u, var_v, var_b, var_e) < 0)) {
    stop("ICC targets incompatible with total_var_eta; negative component.")
  }
  
  cov_bs <- cor_bs * sqrt(var_b) * sqrt(var_s)
  
  set.seed(seed + 20)
  u_school <- rnorm(nlevels(dat$school_id), 0, sqrt(var_u))
  v_cohort <- rnorm(nlevels(dat$cohort_id), 0, sqrt(var_v))
  
  set.seed(seed + 21)
  Sigma_bs <- matrix(c(var_b, cov_bs,
                       cov_bs, var_s), 2, 2)
  L <- chol(Sigma_bs + 1e-10 * diag(2))
  Z <- matrix(rnorm(2 * nlevels(dat$child_id)), ncol=2)
  bs <- Z %*% L
  b_child <- bs[,1]
  s_child <- bs[,2]
  
  dat$u_school <- u_school[as.integer(dat$school_id)]
  dat$v_cohort <- v_cohort[as.integer(dat$cohort_id)]
  dat$b_child  <- b_child[as.integer(dat$child_id)]
  dat$s_child  <- s_child[as.integer(dat$child_id)]
  
  # Smooth-ish age effect (proxy for pb(age))
  f_age <- function(a) 0.15*(a-8) - 0.02*(a-8)^2
  
  # fixed effects (shared latent; outcomes will use scaled versions)
  beta0    <- 2.6
  beta_sex <- 0.20
  beta_bmi <- -0.03
  
  set.seed(seed + 22)
  eps_eta <- rnorm(nrow(dat), 0, sqrt(var_e))
  
  dat$eta_mu_true <- beta0 + f_age(dat$age) + beta_sex*dat$sex + beta_bmi*dat$bmi +
    dat$u_school + dat$v_cohort + dat$b_child + dat$s_child*dat$time + eps_eta
  
  # sigma predictor (mild heteroskedasticity + school RE)
  set.seed(seed + 30)
  u_school_sigma <- rnorm(nlevels(dat$school_id), 0, sqrt(0.02))
  dat$u_school_sigma <- u_school_sigma[as.integer(dat$school_id)]
  dat$eta_sigma_true <- -0.3 + 0.03*(dat$age - 8) + dat$u_school_sigma
  
  dat
}

# ============================================================
# 3) Outcome generators
# ============================================================

# 3A) Sit-and-reach: hurdle (P(Y>0) + positive LogNormal)
simulate_sitreach <- function(dat,
                              gamma0 = 0.5, gamma_age = 0.8, gamma_sex = 0.2, gamma_bmi = -0.05,
                              a0 = 2.0, a_age = 0.10, a_sex = 0.05, a_bmi = -0.01,
                              sigma_pos = 0.35,
                              re_scale_zero = 0.5,
                              seed=1) {
  set.seed(seed + 100)
  age_c <- (dat$age - 8)
  
  eta0 <- gamma0 + gamma_age*age_c + gamma_sex*dat$sex + gamma_bmi*dat$bmi +
    re_scale_zero*(dat$u_school + dat$v_cohort + dat$b_child)
  
  p_pos <- plogis(eta0)
  z_pos <- rbinom(nrow(dat), 1, p_pos)
  
  mu_log <- a0 + a_age*age_c + a_sex*dat$sex + a_bmi*dat$bmi + 0.25*dat$eta_mu_true
  y_pos <- rlnorm(nrow(dat), meanlog = mu_log, sdlog = sigma_pos)
  
  y <- ifelse(z_pos == 1, y_pos, 0)
  
  list(
    y = y,
    p_pos_true = p_pos,
    mu_pos_log_true = mu_log,
    sigma_pos_true = rep(sigma_pos, nrow(dat))
  )
}

# 3B) Rope skipping: NegBin (overdispersed) + optional structural zeros
# Using rNBI(mu, sigma) from gamlss.dist
simulate_ropeskip <- function(dat,
                              b0 = 4.5, b_age = 0.18, b_sex = 0.08, b_bmi = -0.02,
                              sigma_nb = 0.8,
                              zero_infl = FALSE,
                              omega0 = -2.5, omega_age = -0.6, omega_sex = 0.1, omega_bmi = 0.05,
                              # If provided, we calibrate omega0 so that mean(omega) ~= target_omega_mean
                              # where omega = P(structural zero). This makes the scenario interpretable.
                              target_omega_mean = NULL,
                              cohort_zero_scale = 1.0,
                              seed=1) {
  
  set.seed(seed + 200)
  age_c <- (dat$age - 8)
  
  eta <- b0 + b_age*age_c + b_sex*dat$sex + b_bmi*dat$bmi + 0.15*dat$eta_mu_true
  mu <- exp(eta)
  
  y_nb <- rNBI(nrow(dat), mu = mu, sigma = rep(sigma_nb, nrow(dat)))
  
  if (!zero_infl) {
    return(list(
      y = y_nb,
      mu_count_true = mu,
      sigma_nb_true = rep(sigma_nb, nrow(dat)),
      omega_true = rep(0, nrow(dat))
    ))
  }
  
  set.seed(seed + 201)
  cohort_shock <- rnorm(nlevels(dat$cohort_id), 0, 0.7)
  lin_zi_no_int <- omega_age*age_c + omega_sex*dat$sex + omega_bmi*dat$bmi +
    cohort_zero_scale * cohort_shock[as.integer(dat$cohort_id)]
  
  # Optional calibration: choose omega0 so the average structural-zero rate hits target.
  if (!is.null(target_omega_mean)) {
    f_root <- function(w0) mean(plogis(w0 + lin_zi_no_int)) - target_omega_mean
    omega0 <- uniroot(f_root, interval = c(-12, 12))$root
  }
  
  eta_zi <- omega0 + lin_zi_no_int
  omega <- plogis(eta_zi)
  
  z0 <- rbinom(nrow(dat), 1, omega)
  y <- ifelse(z0 == 1, 0, y_nb)
  
  list(
    y = y,
    mu_count_true = mu,
    sigma_nb_true = rep(sigma_nb, nrow(dat)),
    omega_true = omega,
    omega0_used = omega0
  )
}

# 3C) 50m dash time: LogNormal
simulate_50m <- function(dat,
                         c0 = 2.6, c_age = -0.05, c_sex = -0.03, c_bmi = 0.01,
                         sdlog_base = 0.12, sdlog_age = 0.01,
                         seed=1) {
  set.seed(seed + 300)
  age_c <- (dat$age - 8)
  
  mu_log <- c0 + c_age*age_c + c_sex*dat$sex + c_bmi*dat$bmi + 0.08*dat$eta_mu_true
  sdlog <- pmax(0.05, sdlog_base + sdlog_age*abs(age_c))
  y <- rlnorm(nrow(dat), meanlog = mu_log, sdlog = sdlog)
  
  list(y = y, mu_log_true = mu_log, sdlog_true = sdlog)
}

# 3D) Vital capacity: LogNormal
simulate_vitalcap <- function(dat,
                              d0 = 7.0, d_age = 0.10, d_sex = 0.06, d_bmi = 0.03,
                              sdlog_base = 0.18, sdlog_age = 0.01,
                              seed=1) {
  set.seed(seed + 400)
  age_c <- (dat$age - 8)
  
  mu_log <- d0 + d_age*age_c + d_sex*dat$sex + d_bmi*dat$bmi + 0.10*dat$eta_mu_true
  sdlog <- pmax(0.08, sdlog_base + sdlog_age*abs(age_c))
  y <- rlnorm(nrow(dat), meanlog = mu_log, sdlog = sdlog)
  
  list(y = y, mu_log_true = mu_log, sdlog_true = sdlog)
}

# ============================================================
# 4) Baseline generator: ALL tests
# ============================================================

simulate_fitness_baseline <- function(
    J=10, H=60, n=50, T=3,
    ICC_school=0.01, ICC_cohort=0.05, ICC_child=0.50,
    total_var_eta = 1.0,
    seed=1,
    rope_zero_infl = FALSE
) {
  dat <- make_design(J,H,n,T,seed=seed)
  dat <- simulate_covariates(dat, seed=seed)
  dat$bmi <- simulate_bmi_quadratic(dat, seed=seed+12)
  
  dat <- simulate_latent_predictors(
    dat,
    ICC_school=ICC_school, ICC_cohort=ICC_cohort, ICC_child=ICC_child,
    total_var_eta = total_var_eta,
    seed=seed
  )
  
  sr  <- simulate_sitreach(dat, seed=seed)
  rp  <- simulate_ropeskip(dat, zero_infl = rope_zero_infl, seed=seed)
  d50 <- simulate_50m(dat, seed=seed)
  vc  <- simulate_vitalcap(dat, seed=seed)
  
  dat$sitreach <- sr$y
  dat$ropeskip <- rp$y
  dat$dash50   <- d50$y
  dat$vitalcap <- vc$y
  
  # Store some truth (optional)
  # Store truth for RMSE (mu, sigma, and hurdle pieces)
  dat$sitreach_p_pos_true      <- sr$p_pos_true
  dat$sitreach_mu_pos_log_true <- sr$mu_pos_log_true
  dat$sitreach_sigma_pos_true  <- sr$sigma_pos_true
  
  dat$ropeskip_mu_true         <- rp$mu_count_true
  dat$ropeskip_sigma_true      <- rp$sigma_nb_true
  dat$ropeskip_omega_true      <- rp$omega_true
  
  dat$dash50_mu_log_true       <- d50$mu_log_true
  dat$dash50_sigma_log_true    <- d50$sdlog_true
  
  dat$vitalcap_mu_log_true     <- vc$mu_log_true
  dat$vitalcap_sigma_log_true  <- vc$sdlog_true
  
  dat
}

# ============================================================
# 5) Scenario modifiers
# ============================================================

apply_heaping <- function(x, prob = 0.30, step = 1, type = c("round","floor","ceil"), seed=1) {
  type <- match.arg(type)
  set.seed(seed)
  flag <- runif(length(x)) < prob
  y <- x
  if (type == "round") y[flag] <- round(y[flag] / step) * step
  if (type == "floor") y[flag] <- floor(y[flag] / step) * step
  if (type == "ceil")  y[flag] <- ceiling(y[flag] / step) * step
  y
}

apply_scenarios_fitness <- function(dat,
                                    scenario = c("S1","S2","S3","S4","S5","S6","S7"),
                                    seed=1) {
  scenario <- match.arg(scenario)
  
  if (scenario == "S1") return(dat)
  
  # S2: Strong floor effects in sit-and-reach
  if (scenario == "S2") {
    sr <- simulate_sitreach(
      dat,
      gamma0 = 0.1, gamma_age = 0.5, gamma_sex = 0.15, gamma_bmi = -0.10,
      re_scale_zero = 0.7,
      seed = seed + 10
    )
    dat$sitreach                  <- sr$y               # update outcome column with scenario values
    dat$sitreach_p_pos_true      <- sr$p_pos_true
    dat$sitreach_mu_pos_log_true <- sr$mu_pos_log_true
    dat$sitreach_sigma_pos_true  <- sr$sigma_pos_true
    return(dat)
  }
  
  # S3: Rope skipping zero inflation + higher overdispersion
  if (scenario == "S3") {
    rp <- simulate_ropeskip(
      dat,
      sigma_nb = 1.2,
      zero_infl = TRUE,
      # Target ~5–10% structural zeros (overall zeros will be close to this,
      # since NB zeros are typically tiny when mean counts are large)
      target_omega_mean = 0.07,
      omega_age = -0.4, omega_bmi = 0.08,
      cohort_zero_scale = 1.0,
      seed = seed + 20
    )
    dat$ropeskip <- rp$y
    dat$ropeskip_mu_true    <- rp$mu_count_true
    dat$ropeskip_sigma_true <- rp$sigma_nb_true
    dat$ropeskip_omega_true <- rp$omega_true
    return(dat)
  }
  
  # S4: Measurement heaping (partial rounding to realistic reporting units)
  # dash50 heaped to nearest 0.5s (~50% of observations); vitalcap to nearest 50ml (~50%)
  # This is lighter than S7 (which applies 100% heaping AND school calibration bias).
  if (scenario == "S4") {
    dat$dash50   <- apply_heaping(dat$dash50,   prob = 0.50, step = 0.5,  type = "round", seed = seed + 41)
    dat$vitalcap <- apply_heaping(dat$vitalcap, prob = 0.50, step = 50,   type = "round", seed = seed + 42)
    return(dat)
  }
  
  # S5: Informative dropout based on baseline (time==0) composite low performance
  if (scenario == "S5") {
    base <- dat[dat$time == 0, c("child_id","sitreach","ropeskip","dash50","vitalcap")]
    z <- function(x) as.numeric(scale(x))
    base$low_index <- (-z(base$sitreach)) + (-z(base$ropeskip)) + (z(base$dash50)) + (-z(base$vitalcap))
    
    dat <- merge(dat, base[, c("child_id","low_index")], by="child_id", all.x=TRUE)
    dat <- dat[order(dat$child_id, dat$time), ]
    
    set.seed(seed + 50)
    p_drop <- plogis(-1.0 + 0.9*dat$low_index)
    
    dat$observed <- TRUE
    dat$observed[dat$time >= 1] <- runif(sum(dat$time >= 1)) > p_drop[dat$time >= 1]
    
    # monotone dropout per child
    keep <- logical(nrow(dat))
    current_child <- dat$child_id[1]
    still_in <- TRUE
    
    for (k in seq_len(nrow(dat))) {
      if (dat$child_id[k] != current_child) {
        current_child <- dat$child_id[k]
        still_in <- TRUE
      }
      if (!still_in) {
        keep[k] <- FALSE
      } else {
        keep[k] <- dat$observed[k]
        if (!keep[k]) still_in <- FALSE
      }
    }
    
    dat <- dat[keep, ]
    dat$low_index <- NULL
    dat$observed <- NULL
    return(dat)
  }
  
  # S6: Cohort-level intervention shock
  if (scenario == "S6") {
    set.seed(seed + 60)
    
    cohorts <- levels(dat$cohort_id)
    treated <- sample(cohorts, size = max(1, floor(0.10*length(cohorts))), replace = FALSE)
    dat$treated <- as.integer(dat$cohort_id %in% treated)
    post <- as.integer(dat$time >= 1)
    
    mult_sr  <- 1 + 0.10*(dat$treated*post)
    mult_rp  <- 1 + 0.12*(dat$treated*post)
    mult_d50 <- 1 - 0.03*(dat$treated*post)
    mult_vc  <- 1 + 0.08*(dat$treated*post)
    
    dat$sitreach <- dat$sitreach * mult_sr
    dat$ropeskip <- round(dat$ropeskip * mult_rp)
    dat$dash50   <- dat$dash50   * mult_d50
    dat$vitalcap <- dat$vitalcap * mult_vc
    
    return(dat)
  }
  
  # S7: Dash/VC device/administrator calibration bias + rounding/heaping
  if (scenario == "S7") {
    set.seed(seed + 80)
    
    # log-scale SD for multiplicative calibration differences across schools
    cal_dash_sd <- 0.03  # mild (~3% SD multiplicative)
    cal_vc_sd   <- 0.06  # moderate (~6% SD multiplicative)
    
    c_dash <- rnorm(nlevels(dat$school_id), 0, cal_dash_sd)
    c_vc   <- rnorm(nlevels(dat$school_id), 0, cal_vc_sd)
    
    dat$dash50   <- dat$dash50   * exp(c_dash[as.integer(dat$school_id)])
    dat$vitalcap <- dat$vitalcap * exp(c_vc[as.integer(dat$school_id)])
    
    # Always-rounded recording (common in practice); adjust prob if desired
    dat$dash50   <- apply_heaping(dat$dash50,   prob=1.00, step=0.1, type="round", seed=seed+81)
    dat$vitalcap <- apply_heaping(dat$vitalcap, prob=1.00, step=50,  type="round", seed=seed+82)
    
    return(dat)
  }
  
  dat
}