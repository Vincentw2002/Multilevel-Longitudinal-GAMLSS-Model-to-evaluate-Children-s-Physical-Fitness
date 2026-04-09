# ============================================================
#  Shiny Dashboard.R  -  Children's Fitness Assessment Dashboard 
#
#  Place in the same folder as:
#    GAMLSS Model Architecture.R
#    Helper Functions.R
#
#  Required CSV columns:
#    school_id, cohort_id, child_id, time, age, sex
#    + at least one of: sitreach, ropeskip, dash50, vitalcap
#  Optional: bmi
# ============================================================

# ── 1. Libraries ─────────────────────────────────────────────
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(gamlss)
  library(gamlss.dist)
  library(dplyr)
  library(ggplot2)
  library(DT)
  library(tidyr)
  library(scales)
})

# ── 2. Source model scripts ───────────────────────────────────
source("GAMLSS Model Architecture.R")
local({
  interactive <- function() FALSE
  source("Helper Functions.R", local = FALSE)
})

# ── 3. Constants ──────────────────────────────────────────────
OUTCOME_LABELS <- c(
  sitreach = "Sit-and-Reach (cm)",
  ropeskip = "Rope Skipping (count)",
  dash50   = "50m Dash (seconds)",
  vitalcap = "Vital Capacity (mL)"
)

OUTCOME_DIRECTION <- c(
  sitreach = "higher = better flexibility",
  ropeskip = "higher = better coordination",
  dash50   = "lower  = faster runner",
  vitalcap = "higher = stronger lungs"
)

# For dash50, higher z-score = worse (slower), so flag upper tail
OUTCOME_LOWER_BETTER <- c(dash50 = TRUE)

PALETTE_SCHOOLS <- c(
  "#2980b9","#e74c3c","#27ae60","#f39c12",
  "#8e44ad","#16a085","#d35400","#2c3e50",
  "#c0392b","#1abc9c","#e91e63","#607d8b"
)

PALETTE_COHORTS <- c(
  "#3498db","#e74c3c","#2ecc71","#f39c12","#9b59b6",
  "#1abc9c","#e67e22","#34495e","#e91e63","#00bcd4",
  "#8bc34a","#ff5722","#607d8b","#795548","#ff9800"
)

# ── 4. Helper functions ───────────────────────────────────────

validate_data <- function(df) {
  problems <- character(0)
  required <- c("school_id","cohort_id","child_id","time","age","sex")
  missing  <- required[!required %in% names(df)]
  if (length(missing) > 0)
    problems <- c(problems, paste("Missing required columns:", paste(missing, collapse=", ")))
  if (length(names(OUTCOME_LABELS)[names(OUTCOME_LABELS) %in% names(df)]) == 0)
    problems <- c(problems, "No fitness outcome columns found.")
  problems
}

parse_time_column <- function(x) {
  raw <- suppressWarnings(as.integer(x))
  if (all(is.na(raw[!is.na(x)])))
    raw <- suppressWarnings(
      as.integer(gsub("^[^0-9]*([0-9]+).*$", "\\1", as.character(x))))
  mn <- min(raw, na.rm = TRUE)
  if (is.finite(mn) && mn > 0) raw <- raw - mn
  raw
}

clean_data <- function(df) {
  for (v in c("school_id","cohort_id","child_id"))
    if (v %in% names(df)) df[[v]] <- as.factor(df[[v]])
  if ("sex" %in% names(df)) {
    sx <- tolower(trimws(as.character(df$sex)))
    df$sex <- as.numeric(ifelse(sx %in% c("m","male","boy","1","men","man"), 1L, 0L))
  }
  for (v in c("age","bmi","sitreach","dash50","vitalcap"))
    if (v %in% names(df)) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
  if ("ropeskip" %in% names(df))
    df$ropeskip <- suppressWarnings(as.integer(round(as.numeric(df$ropeskip))))
  if ("time" %in% names(df))
    df$time <- parse_time_column(df$time)
  df
}

compute_variance_decomposition <- function(data, outcomes) {
  rows <- lapply(outcomes, function(oc) {
    if (!oc %in% names(data)) return(NULL)
    y_all <- data[[oc]]; ok <- is.finite(y_all)
    if (sum(ok) < 30) return(NULL)
    y <- y_all[ok]; sid <- data$school_id[ok]; cid <- data$cohort_id[ok]
    N <- length(y); grand <- mean(y); tv <- var(y)
    if (tv < 1e-12) return(NULL)
    sn <- as.numeric(table(sid)); smns <- tapply(y, sid, mean)
    vbs <- max(0, sum(sn * (smns - grand)^2) / (N-1))
    cn  <- as.numeric(table(cid)); cmns <- tapply(y, cid, mean)
    vbc <- max(0, sum(cn * (cmns - grand)^2) / (N-1) - vbs)
    vres <- max(0, tv - vbs - vbc)
    data.frame(
      outcome          = OUTCOME_LABELS[oc],
      pct_school       = 100*vbs/tv,
      pct_cohort       = 100*vbc/tv,
      pct_residual     = 100*vres/tv,
      icc_school       = vbs/tv,
      icc_school_cohort = (vbs+vbc)/tv
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}


safe_diagnostics <- function(df, fit_res) {
  
  .safe <- function(expr) tryCatch(expr, error = function(e) NA)
  
  # Helper: compute quantile residuals for a fitted GAMLSS model.
  get_qr <- function(m, y_obs = NULL) {
    if (is.null(m)) return(NULL)
    
    r <- tryCatch(as.numeric(resid(m, type = "quantile")), error = function(e) NULL)
    if (!is.null(r) && sum(is.finite(r)) > 10) return(r[is.finite(r)])
    
    if (is.null(y_obs) || length(y_obs) == 0) return(NULL)
    if (!inherits(m, "gamlss")) return(NULL)
    
    mu  <- tryCatch(fitted(m, "mu"),    error = function(e) NULL)
    sig <- tryCatch(fitted(m, "sigma"), error = function(e) NULL)
    if (is.null(mu) || is.null(sig)) return(NULL)
    
    n_fit <- length(mu)
    if (length(y_obs) != n_fit) return(NULL)
    
    fam  <- m$family[1]
    
    r <- if (fam == "LOGNO") {
      y_pos <- pmax(y_obs, 1e-9)
      (log(y_pos) - mu) / pmax(sig, 1e-9)
    } else {
      
      pfun <- get0(paste0("p", fam), mode = "function")
      if (is.null(pfun)) return(NULL)
      nu  <- tryCatch(fitted(m, "nu"),  error = function(e) NULL)
      tau <- tryCatch(fitted(m, "tau"), error = function(e) NULL)
      p_vals <- tryCatch(
        if (!is.null(tau) && !is.null(nu))
          pfun(y_obs, mu=mu, sigma=sig, nu=nu, tau=tau)
        else if (!is.null(nu))
          pfun(y_obs, mu=mu, sigma=sig, nu=nu)
        else
          pfun(y_obs, mu=mu, sigma=sig),
        error = function(e) NULL)
      if (is.null(p_vals)) return(NULL)
      # Clamp to avoid Inf from qnorm(0) or qnorm(1)
      p_vals <- pmin(pmax(p_vals, 1e-6), 1 - 1e-6)
      qnorm(p_vals)
    }
    
    r <- r[is.finite(r)]
    if (length(r) < 10) return(NULL)
    r
  }
  
  resid_row <- function(label, m, y_obs = NULL) {
    r <- get_qr(m, y_obs)
    if (is.null(r) || length(r) < 10) return(NULL)
    data.frame(outcome = label,
               mean = mean(r), sd = sd(r),
               q05  = unname(quantile(r, .05)),
               q50  = unname(quantile(r, .50)),
               q95  = unname(quantile(r, .95)),
               stringsAsFactors = FALSE)
  }
  
  pi_rows <- function(label, m, y_vec, levels = c(.80, .90, .95)) {
    if (is.null(m) || !inherits(m, "gamlss")) return(NULL)
    fam  <- m$family[1]
    qfun <- get0(paste0("q", fam), mode = "function")
    if (is.null(qfun)) return(NULL)
    mu  <- fitted(m, "mu")
    sig <- fitted(m, "sigma")
    nu  <- tryCatch(fitted(m, "nu"),  error = function(e) NULL)
    tau <- tryCatch(fitted(m, "tau"), error = function(e) NULL)
    do.call(rbind, lapply(levels, function(lev) {
      a  <- (1 - lev) / 2
      lo <- tryCatch(
        if (!is.null(tau) && !is.null(nu)) qfun(a,     mu=mu, sigma=sig, nu=nu, tau=tau)
        else if (!is.null(nu))             qfun(a,     mu=mu, sigma=sig, nu=nu)
        else                               qfun(a,     mu=mu, sigma=sig),
        error = function(e) NULL)
      hi <- tryCatch(
        if (!is.null(tau) && !is.null(nu)) qfun(1-a,   mu=mu, sigma=sig, nu=nu, tau=tau)
        else if (!is.null(nu))             qfun(1-a,   mu=mu, sigma=sig, nu=nu)
        else                               qfun(1-a,   mu=mu, sigma=sig),
        error = function(e) NULL)
      if (is.null(lo) || is.null(hi)) return(NULL)
      data.frame(outcome   = label,
                 level     = lev,
                 coverage  = mean(y_vec >= lo & y_vec <= hi, na.rm = TRUE),
                 avg_width = mean(hi - lo, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }
  
  pi_rows_logno <- function(label, m, y_vec, levels = c(.80, .90, .95)) {
    if (is.null(m)) return(NULL)
    mu  <- tryCatch(fitted(m, "mu"),    error = function(e) NULL)
    sig <- tryCatch(fitted(m, "sigma"), error = function(e) NULL)
    if (is.null(mu) || is.null(sig)) return(NULL)
    do.call(rbind, lapply(levels, function(lev) {
      a  <- (1 - lev) / 2
      lo <- exp(mu + qnorm(a)   * sig)
      hi <- exp(mu + qnorm(1-a) * sig)
      data.frame(outcome   = label,
                 level     = lev,
                 coverage  = mean(y_vec >= lo & y_vec <= hi, na.rm = TRUE),
                 avg_width = mean(hi - lo, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }
  
  cent_rows_logno <- function(label, m, y_vec,
                              centiles = c(5,10,25,50,75,90,95)) {
    if (is.null(m)) return(NULL)
    mu  <- tryCatch(fitted(m, "mu"),    error = function(e) NULL)
    sig <- tryCatch(fitted(m, "sigma"), error = function(e) NULL)
    if (is.null(mu) || is.null(sig)) return(NULL)
    do.call(rbind, lapply(centiles, function(cen) {
      p    <- cen / 100
      qhat <- exp(mu + qnorm(p) * sig)
      data.frame(outcome   = label,
                 centile   = cen,
                 target    = p,
                 empirical = mean(y_vec <= qhat, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }
  
  cent_rows_fam <- function(label, m, y_vec,
                            centiles = c(5,10,25,50,75,90,95)) {
    if (is.null(m) || !inherits(m, "gamlss")) return(NULL)
    fam  <- m$family[1]
    qfun <- get0(paste0("q", fam), mode = "function")
    if (is.null(qfun)) return(NULL)
    mu  <- fitted(m, "mu")
    sig <- fitted(m, "sigma")
    nu  <- tryCatch(fitted(m, "nu"),  error = function(e) NULL)
    tau <- tryCatch(fitted(m, "tau"), error = function(e) NULL)
    do.call(rbind, lapply(centiles, function(cen) {
      p    <- cen / 100
      qhat <- tryCatch(
        if (!is.null(tau) && !is.null(nu)) qfun(p, mu=mu, sigma=sig, nu=nu, tau=tau)
        else if (!is.null(nu))             qfun(p, mu=mu, sigma=sig, nu=nu)
        else                               qfun(p, mu=mu, sigma=sig),
        error = function(e) NULL)
      if (is.null(qhat)) return(NULL)
      data.frame(outcome   = label,
                 centile   = cen,
                 target    = p,
                 empirical = mean(y_vec <= qhat, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }
  
  # ── Assemble per-outcome ─────────────────────────────────
  m_d50 <- fit_res$dash50$model
  m_vc  <- fit_res$vitalcap$model
  m_rp  <- fit_res$ropeskip$model
  m_sr0 <- if (!is.null(fit_res$sitreach$model)) fit_res$sitreach$model$model_zero else NULL
  m_srp <- if (!is.null(fit_res$sitreach$model)) fit_res$sitreach$model$model_pos  else NULL
  dp    <- if (!is.null(fit_res$sitreach$model)) fit_res$sitreach$model$data_pos   else NULL
  
  y_d50 <- if ("dash50"   %in% names(df)) df$dash50[is.finite(df$dash50)]     else numeric(0)
  y_vc  <- if ("vitalcap" %in% names(df)) df$vitalcap[is.finite(df$vitalcap)] else numeric(0)
  y_rp  <- if ("ropeskip" %in% names(df)) df$ropeskip[is.finite(df$ropeskip)] else numeric(0)
  y_srp <- if (!is.null(dp) && "y_pos" %in% names(dp)) dp$y_pos else numeric(0)
  
  # Residual summaries
  resid <- do.call(rbind, Filter(Negate(is.null), list(
    resid_row("dash50",       m_d50, y_d50),
    resid_row("vitalcap",     m_vc,  y_vc),
    resid_row("ropeskip",     m_rp,  y_rp),
    resid_row("sitreach_pos", m_srp, y_srp)
  )))
  
  # PI coverage 
  pi <- do.call(rbind, Filter(Negate(is.null), list(
    if (length(y_d50) > 0) pi_rows_logno("dash50",       m_d50, y_d50) else NULL,
    if (length(y_vc)  > 0) pi_rows_logno("vitalcap",     m_vc,  y_vc)  else NULL,
    if (length(y_rp)  > 0) pi_rows("ropeskip",           m_rp,  y_rp)  else NULL,
    if (length(y_srp) > 0) pi_rows_logno("sitreach_pos", m_srp, y_srp) else NULL
  )))
  
  # Centile calibration
  centile <- do.call(rbind, Filter(Negate(is.null), list(
    if (length(y_d50) > 0) cent_rows_logno("dash50",       m_d50, y_d50) else NULL,
    if (length(y_vc)  > 0) cent_rows_logno("vitalcap",     m_vc,  y_vc)  else NULL,
    if (length(y_rp)  > 0) cent_rows_fam("ropeskip",       m_rp,  y_rp)  else NULL,
    if (length(y_srp) > 0) cent_rows_logno("sitreach_pos", m_srp, y_srp) else NULL
  )))
  
  list(resid = resid, pi = pi, centile = centile,
       stability = NULL, quality = NULL, missing = NULL, rmse = NULL)
}

# ── Flag at-risk children (>1.5 SD from cohort mean) ─────────
flag_atrisk <- function(data, outcomes, school_filter = NULL) {
  if (!is.null(school_filter))
    data <- data[data$school_id == school_filter, ]
  rows <- lapply(outcomes, function(oc) {
    if (!oc %in% names(data)) return(NULL)
    lower_better <- isTRUE(OUTCOME_LOWER_BETTER[oc])
    data %>%
      filter(is.finite(.data[[oc]])) %>%
      group_by(school_id, cohort_id) %>%
      mutate(
        cohort_mean = mean(.data[[oc]], na.rm=TRUE),
        cohort_sd   = sd(.data[[oc]], na.rm=TRUE),
        z_score     = (.data[[oc]] - cohort_mean) / pmax(cohort_sd, 1e-6),
        at_risk     = if (lower_better) z_score > 1.5 else z_score < -1.5
      ) %>%
      ungroup() %>%
      filter(at_risk) %>%
      transmute(
        School       = as.character(school_id),
        Class        = as.character(cohort_id),
        Child        = as.character(child_id),
        Test         = OUTCOME_LABELS[oc],
        Score        = round(.data[[oc]], 1),
        `Class Mean` = round(cohort_mean, 1),
        `Class SD`   = round(cohort_sd, 1),
        `Z-score`    = round(z_score, 2),
        `Time Point` = time
      )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

# ── Data-driven PE recommendations ───────────────────────────
make_recommendations <- function(data, var_decomp, outcomes, school_filter = NULL) {
  df_s <- if (!is.null(school_filter)) data[data$school_id==school_filter,] else data
  recs <- list()
  
  if (!is.null(var_decomp) && nrow(var_decomp)>0) {
    high_het <- var_decomp %>% filter(pct_cohort+pct_school>15) %>% pull(outcome)
    if (length(high_het)>0)
      recs[["het"]] <- list(icon="⚠️", colour="#e74c3c",
                            title="High Between-Class Variation Detected",
                            body=paste0(
                              "Significant differences exist between classes for: ",
                              paste(high_het,collapse=", "),". ",
                              "Audit PE lesson frequency and quality across classes. ",
                              "Consider peer-mentoring schemes or sharing PE resources between ",
                              "high- and low-performing classes."))
  }
  
  if ("sitreach" %in% names(df_s)) {
    pz <- mean(df_s$sitreach==0, na.rm=TRUE)
    if (pz>0.10)
      recs[["flex"]] <- list(icon="🧘", colour="#f39c12",
                             title=sprintf("Flexibility: %.0f%% of children cannot reach past their toes", 100*pz),
                             body=paste(
                               "Incorporate 5–10 minutes of structured stretching into every PE warm-up.",
                               "Focus on hamstring and lower-back flexibility (seated toe touches, cat-cow).",
                               "Daily classroom movement breaks (stretch-and-stand every 30 minutes) can supplement PE time."))
  }
  
  if ("ropeskip" %in% names(df_s)) {
    med_rp <- median(df_s$ropeskip, na.rm=TRUE)
    if (med_rp<30)
      recs[["coord"]] <- list(icon="🪢", colour="#8e44ad",
                              title=sprintf("Coordination: Median rope-skipping count is low (%.0f jumps)", med_rp),
                              body=paste(
                                "Introduce progressive rope-skipping: start with jumping in place (no rope),",
                                "progress to single jumps, then continuous skipping.",
                                "Short daily practice (5 minutes) yields faster improvement than weekly sessions alone."))
  }
  
  if ("dash50" %in% names(df_s)) {
    p90_dash <- quantile(df_s$dash50, 0.90, na.rm=TRUE)
    if (p90_dash>12)
      recs[["speed"]] <- list(icon="🏃", colour="#2980b9",
                              title="Running Speed: Slowest 10% take over 12 seconds for 50m",
                              body=paste(
                                "Include sprint-interval drills in PE (e.g. 3×30m sprints with 60s recovery).",
                                "Tag and relay games build speed in an engaging way.",
                                "Check for musculoskeletal issues in the slowest-performing children."))
  }
  
  if ("vitalcap" %in% names(df_s)) {
    p10_vc <- quantile(df_s$vitalcap, 0.10, na.rm=TRUE)
    if (p10_vc<1200)
      recs[["lung"]] <- list(icon="🫁", colour="#16a085",
                             title="Lung Capacity: Lowest 10% show reduced vital capacity",
                             body=paste(
                               "Aerobic activities (swimming, sustained running, cycling) develop lung capacity most effectively.",
                               "Refer children with persistently very low vital capacity to a school nurse or GP",
                               "to rule out conditions such as asthma."))
  }
  
  recs[["monitor"]] <- list(icon="📅", colour="#27ae60",
                            title="Track Progress Longitudinally",
                            body=paste(
                              "Repeat assessments each academic year using identical protocols.",
                              "A consistent decline in z-score across time points is a stronger signal than a single low score.",
                              "Share anonymised class-level summaries with parents to increase home support for physical activity."))
  
  recs[["curriculum"]] <- list(icon="📚", colour="#2c3e50",
                               title="Curriculum Design — Ensure Coverage of All Four Fitness Domains",
                               body=paste(
                                 "(1) Flexibility — yoga/stretching units;",
                                 "(2) Coordination — gymnastics, rope skills;",
                                 "(3) Speed/Power — athletics units;",
                                 "(4) Cardiorespiratory Fitness — sustained aerobic activities.",
                                 "Aim for ≥150 minutes of moderate-to-vigorous physical activity per week",
                                 "(WHO guidelines for children aged 5–17)."))
  
  recs
}


# ── 5. UI ─────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title    = tags$span("🏃 Fitness Dashboard"),
    titleWidth = 275
  ),
  
  dashboardSidebar(
    width = 275,
    tags$div(style="padding:14px;",
             
             tags$h4("📁 Upload Data", style="color:#b8c7ce;"),
             fileInput("csv_files", label=NULL, multiple=TRUE, accept=".csv",
                       buttonLabel="Browse CSV…", placeholder="No file selected"),
             tags$small("Multiple CSVs will be row-bound.", style="color:#8aa4af;"),
             
             tags$hr(style="border-color:#374850;"),
             
             tags$h4("⚙️ Tests to Analyse", style="color:#b8c7ce;"),
             checkboxGroupInput("outcomes_sel", label=NULL,
                                choices  = c("Sit-and-Reach"="sitreach","Rope Skipping"="ropeskip",
                                             "50m Dash"="dash50","Vital Capacity"="vitalcap"),
                                selected = c("sitreach","ropeskip","dash50","vitalcap")),
             
             tags$hr(style="border-color:#374850;"),
             
             tags$h4("🔬 Model Options", style="color:#b8c7ce;"),
             checkboxInput("include_bmi", "Include BMI as covariate", value=TRUE),
             sliderInput("n_cyc","Max GAMLSS cycles", min=30, max=200, value=100, step=10),
             tags$small("⚠️ Fitting may take several minutes for large datasets.",
                        style="color:#e67e22;font-size:11px;"),
             br(), br(),
             actionBttn("fit_btn","▶  Fit Models", style="fill", color="primary",
                        size="sm", block=TRUE),
             
             tags$hr(style="border-color:#374850;"),
             
             tags$h4("📥 Export", style="color:#b8c7ce;"),
             downloadButton("dl_template","CSV Template",
                            style="width:100%;margin-bottom:6px;background:#2c3e50;color:#ccc;border:1px solid #374850;"),
             downloadButton("dl_report","Download Report",
                            style="width:100%;background:#2c3e50;color:#ccc;border:1px solid #374850;")
    )
  ),
  
  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color:#f4f6f9; }
      .box { border-top:3px solid #3c8dbc; }
      .nav-tabs-custom > .nav-tabs > li.active { border-top-color:#3c8dbc; }
      .small-box h3 { font-size:32px; }
      .shiny-notification { font-size:14px; }
      .dataTables_wrapper .dataTables_filter input {
        border:1px solid #ccc; border-radius:4px; padding:3px 6px; }
      .rec-card { border-radius:8px; padding:14px 18px; margin-bottom:12px;
                  border-left:5px solid #ccc; background:white;
                  box-shadow:0 1px 4px rgba(0,0,0,.06); }
      .stat-chip { background:#f4f6f9; padding:6px 12px; border-radius:4px;
                   text-align:center; display:inline-block; margin:4px; }
    "))),
    
    fluidRow(
      valueBoxOutput("vbox_children", width=3),
      valueBoxOutput("vbox_schools",  width=3),
      valueBoxOutput("vbox_cohorts",  width=3),
      valueBoxOutput("vbox_records",  width=3)
    ),
    
    uiOutput("upload_errors_ui"),
    
    tabBox(width=12, id="main_tabs",
           
           # ══════════════════════════════════════════
           # TAB 1 — Data Overview
           # ══════════════════════════════════════════
           tabPanel(title=tagList(icon("table")," Data Overview"),
                    
                    fluidRow(
                      box(title="Data Preview (first 200 rows)", width=12,
                          solidHeader=TRUE, status="primary", DTOutput("tbl_preview"))
                    ),
                    fluidRow(
                      box(title="Data Quality & Summary Statistics", width=6,
                          solidHeader=TRUE, status="info", DTOutput("tbl_quality")),
                      box(title="Children Observed per Time Point", width=6,
                          solidHeader=TRUE, status="info",
                          plotOutput("plt_missing", height="260px"))
                    ),
                    fluidRow(
                      box(title="Score Distributions", width=12,
                          solidHeader=TRUE, status="primary",
                          plotOutput("plt_distributions", height="380px"))
                    ),
                    fluidRow(
                      box(title="Age & Sex Profile", width=4,
                          solidHeader=TRUE, status="success",
                          plotOutput("plt_age_sex", height="280px")),
                      box(title="Score vs Age (by Sex) — Filter & Zoom", width=8,
                          solidHeader=TRUE, status="success",
                          tags$p("Scatter = individual measurements. Curves = loess smooth ± 95% CI by sex.",
                                 style="color:#666;font-size:12px;margin-bottom:6px;"),
                          fluidRow(
                            column(3, selectInput("age_plot_outcome","Fitness test:",
                                                  choices=setNames(names(OUTCOME_LABELS), unname(OUTCOME_LABELS)),
                                                  width="100%")),
                            column(3, uiOutput("age_plot_school_ui")),
                            column(3, uiOutput("age_plot_cohort_ui")),
                            column(3, checkboxInput("age_plot_raw","Show raw points", value=TRUE))
                          ),
                          plotOutput("plt_age_score", height="320px"))
                    )
           ),
           
           # ══════════════════════════════════════════
           # TAB 2 — Model Fit Metrics
           # ══════════════════════════════════════════
           tabPanel(title=tagList(icon("chart-line")," Model Fit Metrics"),
                    
                    uiOutput("fit_status_ui"),
                    fluidRow(
                      box(title="Convergence & AIC", width=5,
                          solidHeader=TRUE, status="primary", DTOutput("tbl_convergence")),
                      box(title="Quantile Residual Summary",
                          footer=paste0("Mean ≈ 0 and SD ≈ 1 are expected for Log-Normal models (50m Dash, Vital Capacity, Sit-and-Reach) "," because the GAMLSS simultaneously fits both location and scale — confirming the distribution "," family is appropriate. For count outcomes (Rope Skipping, NBI), values closer to 0 / 1 indicate "," better fit. Use the Centile Calibration plot for a more discriminating check."),
                          width=7, solidHeader=TRUE, status="info", DTOutput("tbl_resid"))
                    ),
                    fluidRow(
                      box(title="Prediction Interval Coverage",
                          footer="Actual coverage should match the nominal level (90% PI → ~90% inside).",
                          width=6, solidHeader=TRUE, status="success", DTOutput("tbl_pi")),
                      box(title="Centile Calibration Plot",
                          footer="Points on the diagonal = well-calibrated centiles.",
                          width=6, solidHeader=TRUE, status="warning",
                          plotOutput("plt_centile", height="300px"))
                    ),
                    fluidRow(
                      box(title="Quantile Residual Distributions (vs Standard Normal)",
                          width=12, solidHeader=TRUE, status="primary",
                          plotOutput("plt_resid_hist", height="380px"))
                    )
           ),
           
           # ══════════════════════════════════════════
           # TAB 3 — Heterogeneity
           # ══════════════════════════════════════════
           tabPanel(title=tagList(icon("search")," Heterogeneity"),
                    
                    fluidRow(
                      box(title="Distribution by Cohort / Class", width=12,
                          solidHeader=TRUE, status="primary",
                          tags$p(
                            "Showing the ", tags$strong("highest and lowest scoring classes per school"),
                            " — highlights within-school inequality. Dashed = school average.",
                            style="color:#666;font-size:12px;margin-bottom:8px;"),
                          fluidRow(
                            column(4, selectInput("ht_outcome","Fitness test:",
                                                  choices=setNames(names(OUTCOME_LABELS), unname(OUTCOME_LABELS)))),
                            column(4, sliderInput("ht_n_extreme","Extreme classes per school:",
                                                  min=1, max=5, value=2, step=1)),
                            column(4, checkboxInput("ht_show_points","Overlay data points", value=FALSE))
                          ),
                          plotOutput("plt_cohort_box", height="460px"))
                    ),
                    
                    fluidRow(
                      box(title="School Mean Comparison (± 95% CI)", width=6,
                          solidHeader=TRUE, status="info",
                          tags$p("Dashed = grand average.", style="color:#666;font-size:12px;"),
                          plotOutput("plt_school_means", height="360px")),
                      
                      box(title="Variance Decomposition — Between-Level Focus",
                          footer="Y-axis zoomed to 80–100%. The remaining ~80–95% (individual child variance) lies below the visible range.",
                          width=6, solidHeader=TRUE, status="success",
                          plotOutput("plt_icc_bar", height="360px"))
                    ),
                    
                    fluidRow(
                      box(title="ICC Summary Table", width=12,
                          solidHeader=TRUE, status="primary",
                          tags$p("ICC = Intraclass Correlation Coefficient (empirical, unadjusted for covariates).",
                                 style="color:#666;font-size:12px;"),
                          DTOutput("tbl_icc"))
                    )
           ),
           
           # ══════════════════════════════════════════
           # TAB 4 — Teacher Summary (per-school)
           # ══════════════════════════════════════════
           tabPanel(title=tagList(icon("file-alt")," Summary for Teachers"),
                    
                    tags$div(style="max-width:960px;margin:0 auto;padding:10px;",
                             
                             # School selector + model quality badge
                             fluidRow(
                               column(5, uiOutput("ts_school_selector_ui")),
                               column(7, tags$div(style="padding-top:26px;",
                                                  uiOutput("ts_model_quality_ui")))
                             ),
                             
                             tags$hr(),
                             
                             # Dynamic header + stat cards
                             uiOutput("ts_header_ui"),
                             
                             # Score vs Age plot 
                             box(title="📈 Score vs Age (by Sex) — Selected School / Class",
                                 footer="Scatter = individual measurements. Curves = loess smooth ± 95% CI. Connect = per-child trajectory.",
                                 width=12, solidHeader=TRUE, status="primary",
                                 fluidRow(
                                   column(4, selectInput("ts_traj_outcome","Fitness test:",
                                                         choices=setNames(names(OUTCOME_LABELS), unname(OUTCOME_LABELS)),
                                                         width="100%")),
                                   column(4, uiOutput("ts_traj_cohort_ui")),
                                   column(4, checkboxInput("ts_traj_lines","Connect each child's measurements", value=FALSE))
                                 ),
                                 plotOutput("ts_traj_plot", height="400px")),
                             
                             # Student lookup
                             box(title="🔍 Individual Student Lookup",
                                 footer="Enter a Child ID to highlight that student's measurements across all time points.",
                                 width=12, solidHeader=TRUE, status="warning",
                                 fluidRow(
                                   column(4,
                                          selectizeInput("ts_child_id","Child ID:",
                                                         choices=NULL, options=list(placeholder="Type or select a child ID..."),
                                                         width="100%")),
                                   column(4, selectInput("ts_child_outcome","Fitness test:",
                                                         choices=setNames(names(OUTCOME_LABELS), unname(OUTCOME_LABELS)),
                                                         width="100%")),
                                   column(4, tags$div(style="padding-top:25px;",
                                                      actionButton("ts_child_clear","Clear", class="btn-sm btn-default")))
                                 ),
                                 plotOutput("ts_child_plot", height="340px")),
                             
                             # Within-school cohort comparison
                             box(title="📦 Class-Level Score Distributions",
                                 width=12, solidHeader=TRUE, status="info",
                                 tags$p("Boxplots for all classes, sorted by median. Dashed = school-wide average.",
                                        style="color:#666;font-size:12px;"),
                                 selectInput("ts_cohort_outcome","Fitness test:",
                                             choices=setNames(names(OUTCOME_LABELS), unname(OUTCOME_LABELS)),
                                             width="280px"),
                                 plotOutput("ts_cohort_plot", height="360px")),
                             
                             # At-risk student table
                             box(title="🚨 Children Needing Additional Support",
                                 footer="Children whose score is >1.5 SD below their class mean (or above for 50m Dash).",
                                 width=12, solidHeader=TRUE, status="danger",
                                 DTOutput("ts_atrisk_tbl")),
                             
                             # Data-driven recommendations
                             uiOutput("ts_recs_ui")
                    )
           )
    )
  )
)


# ── 6. Server ─────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── 6a. Load & merge CSVs ────────────────────────────────
  raw_data <- reactive({
    req(input$csv_files)
    dfs <- lapply(input$csv_files$datapath, function(fp)
      tryCatch(read.csv(fp, stringsAsFactors=FALSE), error=function(e) NULL))
    dfs <- Filter(Negate(is.null), dfs)
    if (length(dfs)==0) return(NULL)
    tryCatch(do.call(rbind, dfs),
             error=function(e){ showNotification("Could not merge CSVs.", type="error"); NULL })
  })
  
  # ── 6b. Clean & validate ─────────────────────────────────
  dat <- reactive({
    df <- raw_data(); req(df)
    if ("time" %in% names(df)) {
      ri <- suppressWarnings(as.integer(df$time))
      if (all(is.na(ri[!is.na(df$time)]))) {
        sv <- paste(head(unique(as.character(df$time)), 3), collapse=", ")
        showNotification(
          paste0("Time column auto-converted from labels (", sv, ") to 0-based integers."),
          type="message", duration=8)
      }
    }
    clean_data(df)
  })
  
  # ── 6c. Upload error banner ──────────────────────────────
  output$upload_errors_ui <- renderUI({
    df <- dat(); if (is.null(df)) return(NULL)
    probs <- validate_data(df); if (length(probs)==0) return(NULL)
    tags$div(
      style="background:#f8d7da;border:1px solid #f5c6cb;border-radius:6px;padding:12px 18px;margin:0 15px 15px;color:#721c24;",
      tags$strong("⚠️ Data issues:"), tags$ul(lapply(probs, tags$li)))
  })
  
  # ── 6d. Value boxes ──────────────────────────────────────
  output$vbox_children <- renderValueBox({
    df<-dat()
    v<-if(!is.null(df)&&"child_id"%in%names(df)) length(unique(df$child_id)) else "—"
    valueBox(v,"Children",icon=icon("child"),color="blue") })
  output$vbox_schools <- renderValueBox({
    df<-dat()
    v<-if(!is.null(df)&&"school_id"%in%names(df)) length(unique(df$school_id)) else "—"
    valueBox(v,"Schools",icon=icon("school"),color="green") })
  output$vbox_cohorts <- renderValueBox({
    df<-dat()
    v<-if(!is.null(df)&&"cohort_id"%in%names(df)) length(unique(df$cohort_id)) else "—"
    valueBox(v,"Cohorts / Classes",icon=icon("users"),color="yellow") })
  output$vbox_records <- renderValueBox({
    df<-dat()
    v<-if(!is.null(df)) nrow(df) else "—"
    valueBox(v,"Records",icon=icon("database"),color="red") })
  
  # ── 6e. Tab 1: Data Overview ─────────────────────────────
  output$tbl_preview <- renderDT({
    df<-dat(); req(df)
    show<-head(df,200)
    num<-sapply(show,is.numeric)
    show[,num]<-lapply(show[,num,drop=FALSE],round,2)
    datatable(show,rownames=FALSE,options=list(scrollX=TRUE,pageLength=10,dom="tip"))
  })
  
  output$tbl_quality <- renderDT({
    df<-dat(); req(df)
    ocs<-names(OUTCOME_LABELS)[names(OUTCOME_LABELS)%in%names(df)]; req(length(ocs)>0)
    rows<-lapply(ocs,function(oc){
      y<-df[[oc]]; ok<-is.finite(y)
      data.frame(Test=OUTCOME_LABELS[oc],`N valid`=sum(ok),Missing=sum(!ok),
                 `Mean ± SD`=sprintf("%.1f ± %.1f",mean(y,na.rm=T),sd(y,na.rm=T)),
                 Median=round(median(y,na.rm=T),1),
                 `Min – Max`=sprintf("%.1f – %.1f",min(y,na.rm=T),max(y,na.rm=T)),
                 `% Zero`=sprintf("%.1f%%",100*mean(y==0,na.rm=T)),check.names=FALSE)})
    datatable(do.call(rbind,rows),rownames=FALSE,options=list(dom="t",paging=FALSE))
  })
  
  output$plt_missing <- renderPlot({
    df<-dat(); req(df,"time"%in%names(df),"child_id"%in%names(df))
    miss<-df%>%group_by(time)%>%summarise(n=n_distinct(child_id),.groups="drop")
    ggplot(miss,aes(x=factor(time),y=n))+
      geom_col(fill="#3c8dbc",alpha=.85,width=.6)+
      geom_text(aes(label=n),vjust=-0.4,fontface="bold",size=4.5)+
      labs(x="Time Point",y="# Children",title="Children per Time Point")+
      theme_minimal(base_size=13)+
      theme(plot.title=element_text(face="bold"),panel.grid.major.x=element_blank())
  })
  
  output$plt_distributions <- renderPlot({
    df<-dat(); req(df)
    ocs<-names(OUTCOME_LABELS)[names(OUTCOME_LABELS)%in%names(df)]; req(length(ocs)>0)
    df_long<-df%>%select(all_of(ocs))%>%
      pivot_longer(everything(),names_to="test",values_to="score")%>%
      mutate(test=factor(test,levels=ocs,labels=OUTCOME_LABELS[ocs]))%>%
      filter(is.finite(score))
    ggplot(df_long,aes(x=score))+
      geom_histogram(aes(y=after_stat(density)),bins=40,fill="#3c8dbc",colour="white",alpha=.75)+
      geom_density(colour="#e74c3c",linewidth=1)+
      facet_wrap(~test,scales="free",ncol=2)+
      labs(x="Score",y="Density",title="Score Distributions")+
      theme_minimal(base_size=13)+
      theme(plot.title=element_text(face="bold"),strip.text=element_text(face="bold",size=11))
  })
  
  output$plt_age_sex <- renderPlot({
    df<-dat(); req(df,"age"%in%names(df),"sex"%in%names(df))
    pd<-df%>%mutate(sex_label=ifelse(sex==1,"Male","Female"))%>%filter(is.finite(age))
    ggplot(pd,aes(x=age,fill=sex_label))+
      geom_histogram(binwidth=.5,alpha=.75,colour="white",position="identity")+
      scale_fill_manual(values=c(Male="#3c8dbc",Female="#e74c3c"))+
      labs(x="Age (years)",y="Count",fill="Sex",title="Age & Sex Profile")+
      theme_minimal(base_size=13)+
      theme(plot.title=element_text(face="bold"),legend.position="bottom")
  })
  
  # ── Dynamic filter UIs for the age/score plot ────────────
  output$age_plot_school_ui <- renderUI({
    df <- dat(); req(df, "school_id" %in% names(df))
    schools <- sort(unique(as.character(df$school_id)))
    selectInput("age_plot_school", "School:",
                choices = c("All Schools" = "__all__", schools), width = "100%")
  })
  
  output$age_plot_cohort_ui <- renderUI({
    df  <- dat(); req(df)
    sel <- input$age_plot_school
    df2 <- if (!is.null(sel) && sel != "__all__")
      df[df$school_id == sel, ] else df
    cohorts <- sort(unique(as.character(df2$cohort_id)))
    selectInput("age_plot_cohort", "Class:",
                choices = c("All Classes" = "__all__", cohorts), width = "100%")
  })
  
  output$plt_age_score <- renderPlot({
    df  <- dat(); req(df)
    oc  <- input$age_plot_outcome
    req(oc %in% names(df), "age" %in% names(df), "sex" %in% names(df))
    
    # Apply school / cohort filters
    sel_sc <- input$age_plot_school
    sel_co <- input$age_plot_cohort
    pd <- df
    if (!is.null(sel_sc) && sel_sc != "__all__")
      pd <- pd[pd$school_id == sel_sc, ]
    if (!is.null(sel_co) && sel_co != "__all__")
      pd <- pd[pd$cohort_id == sel_co, ]
    
    pd <- pd %>%
      mutate(sex_label = ifelse(sex == 1, "Male", "Female")) %>%
      filter(is.finite(.data[[oc]]), is.finite(age))
    
    req(nrow(pd) > 5)
    
    # Subtitle describing current filter
    sc_lbl <- if (!is.null(sel_sc) && sel_sc != "__all__") sel_sc else "All Schools"
    co_lbl <- if (!is.null(sel_co) && sel_co != "__all__") sel_co else "All Classes"
    subtitle <- sprintf("Showing: %s  |  %s  |  n = %d measurements", sc_lbl, co_lbl, nrow(pd))
    
    p <- ggplot(pd, aes(x = age, y = .data[[oc]], colour = sex_label))
    
    if (isTRUE(input$age_plot_raw))
      p <- p + geom_point(alpha = .12, size = 1.0)
    
    p <- p +
      geom_line(aes(group = child_id), colour = "grey70", alpha = .18,
                linewidth = .35, show.legend = FALSE) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.4,
                  aes(fill = sex_label), alpha = .15) +
      scale_colour_manual(values = c(Male = "#2980b9", Female = "#e74c3c")) +
      scale_fill_manual(values   = c(Male = "#2980b9", Female = "#e74c3c")) +
      labs(x = "Age (years)", y = OUTCOME_LABELS[oc],
           colour = "Sex", fill = "Sex",
           title   = paste(OUTCOME_LABELS[oc], "vs Age — by Sex"),
           subtitle = subtitle) +
      theme_minimal(base_size = 13) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(colour = "#666", size = 10),
            legend.position = "bottom",
            panel.grid.minor = element_blank())
    p
  })
  
  # ── 6f. Model fitting ────────────────────────────────────
  fit_res    <- reactiveVal(NULL)
  is_fitting <- reactiveVal(FALSE)
  
  observeEvent(input$fit_btn, {
    df<-dat()
    if(is.null(df)){showNotification("Upload a CSV first.",type="warning");return()}
    probs<-validate_data(df)
    if(length(probs)>0){showNotification(probs[1],type="error");return()}
    selected_ocs<-input$outcomes_sel
    if(length(selected_ocs)==0){showNotification("Select at least one test.",type="warning");return()}
    is_fitting(TRUE); on.exit(is_fitting(FALSE))
    ctrl    <- gamlss.control(n.cyc=input$n_cyc,trace=FALSE)
    inc_bmi <- input$include_bmi && "bmi"%in%names(df)
    n_sel   <- length(selected_ocs)
    
    safe_wrap<-function(expr){
      warn<-character(0);err<-NULL;t0<-proc.time()[3]
      fit<-withCallingHandlers(
        tryCatch(expr,error=function(e){err<<-conditionMessage(e);NULL}),
        warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart("muffleWarning")})
      list(model=fit,error=err,warnings=warn,runtime_sec=proc.time()[3]-t0)}
    
    out<-list()
    withProgress(message="Fitting models — please wait…",value=0,{
      if("dash50"%in%selected_ocs&&"dash50"%in%names(df)){
        incProgress(1/n_sel,message="50m Dash…")
        out$dash50<-safe_wrap(fit_50m(df,outcome_var="dash50",include_bmi=inc_bmi,control=ctrl))}
      if("vitalcap"%in%selected_ocs&&"vitalcap"%in%names(df)){
        incProgress(1/n_sel,message="Vital Capacity…")
        out$vitalcap<-safe_wrap(fit_vc(df,outcome_var="vitalcap",include_bmi=inc_bmi,control=ctrl))}
      if("ropeskip"%in%selected_ocs&&"ropeskip"%in%names(df)){
        incProgress(1/n_sel,message="Rope Skipping…")
        out$ropeskip<-safe_wrap(fit_rope(df,outcome_var="ropeskip",include_bmi=inc_bmi,control=ctrl))}
      if("sitreach"%in%selected_ocs&&"sitreach"%in%names(df)){
        incProgress(1/n_sel,message="Sit-and-Reach…")
        out$sitreach<-safe_wrap(fit_sitreach(df,outcome_var="sitreach",include_bmi=inc_bmi,control=ctrl))}
      incProgress(0,message="Diagnostics…")
      conv_flag<-function(m){
        if(is.null(m)) return(NA)
        if(inherits(m,"glm")&&!inherits(m,"gamlss")) return(isTRUE(m$converged))
        if(!is.null(m$iter)&&!is.null(m$control$n.cyc)) return(m$iter<m$control$n.cyc); NA}
      mlist<-list(out$dash50$model,out$vitalcap$model,out$ropeskip$model,
                  if(!is.null(out$sitreach$model)) out$sitreach$model$model_zero else NULL,
                  if(!is.null(out$sitreach$model)) out$sitreach$model$model_pos  else NULL)
      out$summary<-data.frame(
        outcome=c("dash50","vitalcap","ropeskip","sitreach_zero","sitreach_pos"),
        converged=sapply(mlist,conv_flag),
        AIC=sapply(mlist,function(m) tryCatch(round(AIC(m),1),error=function(e) NA_real_)),
        runtime_sec=c(
          if(!is.null(out$dash50))   out$dash50$runtime_sec   else NA,
          if(!is.null(out$vitalcap)) out$vitalcap$runtime_sec else NA,
          if(!is.null(out$ropeskip)) out$ropeskip$runtime_sec else NA,
          if(!is.null(out$sitreach)) out$sitreach$runtime_sec else NA,
          if(!is.null(out$sitreach)) out$sitreach$runtime_sec else NA),
        stringsAsFactors=FALSE)
      out$diag<-safe_diagnostics(df,out)
    })
    fit_res(out)
    showNotification("✅ Model fitting complete!", type="message", duration=6)
  })
  
  # ── 6g. Tab 2: Model Fit outputs ─────────────────────────
  output$fit_status_ui <- renderUI({
    if(is.null(fit_res()))
      return(tags$div(
        style="background:#fff3cd;border:1px solid #ffc107;border-radius:6px;padding:14px 18px;margin-bottom:10px;",
        tags$strong("ℹ️ No models fitted yet. "),
        "Upload data then click ",tags$strong("'Fit Models'"),"." ))
    s<-fit_res()$summary
    nc<-sum(s$converged==TRUE,na.rm=TRUE);nt<-sum(!is.na(s$converged))
    col<-if(nc==nt)"#d4edda" else "#fff3cd"; bc<-if(nc==nt)"#28a745" else "#ffc107"
    tags$div(style=sprintf("background:%s;border-left:5px solid %s;border-radius:0 6px 6px 0;padding:10px 16px;margin-bottom:10px;",col,bc),
             sprintf("Models fitted: %d / %d converged.", nc, nt))
  })
  
  output$tbl_convergence <- renderDT({
    req(fit_res()); s<-fit_res()$summary
    s$Model<-c("50m Dash","Vital Capacity","Rope Skipping","Sit-Reach (zero)","Sit-Reach (pos.)")
    s$Converged<-ifelse(is.na(s$converged),"❓",ifelse(s$converged,"✅ Yes","❌ No"))
    s$`Runtime (s)`<-round(s$runtime_sec,1)
    datatable(s[,c("Model","Converged","AIC","Runtime (s)")],rownames=FALSE,escape=FALSE,
              options=list(dom="t",paging=FALSE))
  })
  
  output$tbl_resid <- renderDT({
    req(fit_res()); rs<-fit_res()$diag$resid
    if(is.null(rs)||nrow(rs)==0) return(datatable(data.frame(Message="No residuals available.")))
    # Extended label map: diagnostic outcomes use sub-names not in OUTCOME_LABELS
    diag_labels <- c(OUTCOME_LABELS,
                     sitreach_pos  = "Sit-and-Reach (pos. part)",
                     sitreach_zero = "Sit-and-Reach (hurdle)")
    rs$Test <- diag_labels[rs$outcome]
    disp<-rs[,c("Test","mean","sd","q05","q50","q95")]
    disp[,-1]<-lapply(disp[,-1],round,3)
    names(disp)<-c("Test","Mean","SD","5th %ile","Median","95th %ile")
    datatable(disp,rownames=FALSE,options=list(dom="t",paging=FALSE))%>%
      formatStyle("Mean",backgroundColor=styleInterval(c(-0.3,0.3),c("#f8d7da","white","#f8d7da")))
  })
  
  output$tbl_pi <- renderDT({
    req(fit_res()); pi<-fit_res()$diag$pi
    if(is.null(pi)||nrow(pi)==0) return(datatable(data.frame(Message="No PI data available.")))
    pi$Test<-OUTCOME_LABELS[pi$outcome]
    pi$`PI Level`<-paste0(round(100*pi$level),"%")
    pi$Actual<-sprintf("%.1f%%",100*pi$coverage)
    pi$`Avg Width`<-round(pi$avg_width,2)
    datatable(pi[,c("Test","PI Level","PI Level","Actual","Avg Width")],rownames=FALSE,
              options=list(dom="t",paging=FALSE))
  })
  
  output$plt_centile <- renderPlot({
    req(fit_res()); cent<-fit_res()$diag$centile
    if(is.null(cent)||nrow(cent)==0){plot.new();text(.5,.5,"No centile data.");return()}
    cent$Test<-OUTCOME_LABELS[cent$outcome]
    ggplot(cent,aes(x=target,y=empirical,colour=Test))+
      geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey50",linewidth=1)+
      geom_line(alpha=.6)+geom_point(size=3)+
      scale_x_continuous(labels=percent)+scale_y_continuous(labels=percent)+
      labs(x="Target centile",y="Empirical centile",title="Centile Calibration",colour=NULL)+
      theme_minimal(base_size=13)+
      theme(plot.title=element_text(face="bold"),legend.position="bottom")
  })
  
  output$plt_resid_hist <- renderPlot({
    req(fit_res()); fr<-fit_res(); df<-dat(); req(df)
    
    # Retrieve models and matching observed y vectors
    mods_y <- list(
      list(nm="50m Dash",        m=fr$dash50$model,
           y=if("dash50"   %in%names(df)) df$dash50[is.finite(df$dash50)]     else NULL),
      list(nm="Vital Capacity",  m=fr$vitalcap$model,
           y=if("vitalcap" %in%names(df)) df$vitalcap[is.finite(df$vitalcap)] else NULL),
      list(nm="Rope Skipping",   m=fr$ropeskip$model,
           y=if("ropeskip" %in%names(df)) df$ropeskip[is.finite(df$ropeskip)] else NULL),
      list(nm="Sit-Reach (pos.)",
           m=if(!is.null(fr$sitreach$model)) fr$sitreach$model$model_pos else NULL,
           y={dp<-if(!is.null(fr$sitreach$model)) fr$sitreach$model$data_pos else NULL
           if(!is.null(dp)&&"y_pos"%in%names(dp)) dp$y_pos else NULL})
    )
    
    rdfs <- lapply(mods_y, function(x) {
      m <- x$m; if(is.null(m)) return(NULL)
      # Try built-in first, then manual fallback
      r <- tryCatch(as.numeric(resid(m, type="quantile")), error=function(e) NULL)
      if(is.null(r) || sum(is.finite(r)) < 10) {
        # Manual: LOGNO shortcut or generic pFAM
        mu  <- tryCatch(fitted(m,"mu"),    error=function(e) NULL)
        sig <- tryCatch(fitted(m,"sigma"), error=function(e) NULL)
        y   <- x$y
        if(!is.null(mu) && !is.null(sig) && !is.null(y) && length(y)==length(mu)) {
          fam  <- tryCatch(m$family[1], error=function(e) "")
          r <- if(fam=="LOGNO") {
            (log(pmax(y,1e-9)) - mu) / pmax(sig, 1e-9)
          } else {
            pfun <- get0(paste0("p",fam), mode="function")
            if(!is.null(pfun)) {
              nu  <- tryCatch(fitted(m,"nu"),  error=function(e) NULL)
              tau <- tryCatch(fitted(m,"tau"), error=function(e) NULL)
              pv  <- tryCatch(
                if(!is.null(tau)&&!is.null(nu)) pfun(y,mu=mu,sigma=sig,nu=nu,tau=tau)
                else if(!is.null(nu))           pfun(y,mu=mu,sigma=sig,nu=nu)
                else                            pfun(y,mu=mu,sigma=sig),
                error=function(e) NULL)
              if(!is.null(pv)) qnorm(pmin(pmax(pv,1e-6),1-1e-6)) else NULL
            }
          }
        }
      }
      if(is.null(r)) return(NULL)
      r <- r[is.finite(r)]; if(length(r)<5) return(NULL)
      data.frame(Test=x$nm, residual=r)
    })
    rdf<-do.call(rbind,Filter(Negate(is.null),rdfs))
    if(is.null(rdf)||nrow(rdf)==0){plot.new();text(.5,.5,"Residuals not available.");return()}
    ggplot(rdf,aes(x=residual))+
      geom_histogram(aes(y=after_stat(density)),bins=35,fill="#3c8dbc",colour="white",alpha=.75)+
      stat_function(fun=dnorm,colour="#e74c3c",linewidth=1.1)+
      facet_wrap(~Test,ncol=2,scales="free_y")+
      labs(x="Quantile residual",y="Density",title="Quantile Residual Distributions",
           subtitle="Red = standard Normal — residuals should follow this curve if model is correct")+
      theme_minimal(base_size=13)+
      theme(plot.title=element_text(face="bold"),strip.text=element_text(face="bold"))
  })
  
  # ── 6h. Variance decomposition (reactive) ────────────────
  var_decomp <- reactive({
    df<-dat(); req(df)
    ocs<-names(OUTCOME_LABELS)[names(OUTCOME_LABELS)%in%names(df)]
    compute_variance_decomposition(df, ocs)
  })
  
  # ── 6i. Tab 3: Heterogeneity ─────────────────────────────
  output$plt_cohort_box <- renderPlot({
    df<-dat(); req(df)
    oc<-input$ht_outcome
    req(oc%in%names(df),"cohort_id"%in%names(df),"school_id"%in%names(df))
    n_extreme<-input$ht_n_extreme
    co_stats<-df%>%filter(is.finite(.data[[oc]]))%>%group_by(school_id,cohort_id)%>%
      summarise(co_med=median(.data[[oc]],na.rm=TRUE),co_n=n(),.groups="drop")
    sc_stats<-co_stats%>%group_by(school_id)%>%summarise(sc_mean=mean(co_med),.groups="drop")
    extreme_ids<-co_stats%>%group_by(school_id)%>%arrange(co_med,.by_group=TRUE)%>%
      filter(row_number()<=n_extreme|row_number()>(n()-n_extreme))%>%
      mutate(rank_label=case_when(
        row_number()<=n_extreme              ~paste0("Bottom ",n_extreme),
        row_number()>(n()-n_extreme)         ~paste0("Top ",n_extreme),
        TRUE                                 ~"Middle"),
        cohort_label=paste0(as.character(school_id),"\n",as.character(cohort_id)))%>%
      ungroup()
    pd<-df%>%filter(is.finite(.data[[oc]]),cohort_id%in%extreme_ids$cohort_id)%>%
      left_join(extreme_ids[,c("cohort_id","school_id","co_med","rank_label","cohort_label")],
                by=c("cohort_id","school_id"))%>%
      left_join(sc_stats,by="school_id")%>%arrange(school_id,co_med)
    req(nrow(pd)>0)
    pd$cohort_label<-factor(pd$cohort_label,levels=unique(pd$cohort_label))
    n_schools<-length(unique(pd$school_id))
    pal<-if(n_schools<=length(PALETTE_SCHOOLS))
      setNames(PALETTE_SCHOOLS[seq_len(n_schools)],levels(factor(pd$school_id)))
    else setNames(scales::hue_pal()(n_schools),levels(factor(pd$school_id)))
    rank_alpha<-c("Bottom 1"=.55,"Bottom 2"=.65,"Bottom 3"=.70,"Bottom 4"=.75,"Bottom 5"=.80,
                  "Top 1"=.95,"Top 2"=.90,"Top 3"=.85,"Top 4"=.80,"Top 5"=.75)
    p<-ggplot(pd,aes(x=cohort_label,y=.data[[oc]]))+
      geom_hline(data=sc_stats%>%filter(school_id%in%unique(pd$school_id)),
                 aes(yintercept=sc_mean),linetype="dashed",colour="#7f8c8d",linewidth=.6)+
      geom_boxplot(aes(fill=school_id,alpha=rank_label),
                   outlier.size=1,outlier.alpha=.4,linewidth=.4,colour="#2c3e50")+
      facet_wrap(~school_id,scales="free_x",nrow=1)+
      scale_fill_manual(values=pal,guide="none")+
      scale_alpha_manual(values=rank_alpha,name="Rank within school",
                         breaks=c(paste0("Bottom ",n_extreme),paste0("Top ",n_extreme)))+
      labs(x=NULL,y=OUTCOME_LABELS[oc],
           title=paste("Most Divergent Classes per School —",OUTCOME_LABELS[oc]),
           subtitle=paste0(n_extreme," lowest (lighter) and ",n_extreme,
                           " highest (darker) classes | Dashed = school average"))+
      theme_minimal(base_size=12)+
      theme(strip.text=element_text(face="bold",size=11,colour="white"),
            strip.background=element_rect(fill="#2c3e50",colour=NA),
            axis.text.x=element_text(size=7.5,lineheight=.85,colour="#444"),
            axis.ticks.x=element_blank(),panel.grid.major.x=element_blank(),
            panel.spacing=unit(1.2,"lines"),legend.position="bottom",
            plot.title=element_text(face="bold",size=13),
            plot.subtitle=element_text(colour="#555",size=10))
    if(input$ht_show_points)
      p<-p+geom_jitter(aes(colour=school_id),width=.18,alpha=.2,size=.7,show.legend=FALSE)+
      scale_colour_manual(values=pal)
    p
  })
  
  output$plt_school_means <- renderPlot({
    df<-dat(); req(df)
    oc<-input$ht_outcome; req(oc%in%names(df),"school_id"%in%names(df))
    gm<-mean(df[[oc]],na.rm=TRUE)
    sm<-df%>%filter(is.finite(.data[[oc]]))%>%group_by(school_id)%>%
      summarise(mn=mean(.data[[oc]]),se=sd(.data[[oc]])/sqrt(n()),.groups="drop")
    ggplot(sm,aes(x=reorder(school_id,mn),y=mn))+
      geom_hline(yintercept=gm,linetype="dashed",colour="grey50",linewidth=1)+
      geom_errorbar(aes(ymin=mn-1.96*se,ymax=mn+1.96*se),width=.2,colour="#7f8c8d")+
      geom_point(size=4.5,colour="#2980b9")+coord_flip()+
      labs(x="School",y=OUTCOME_LABELS[oc],title="School Mean Scores (± 95% CI)",
           subtitle="Dashed = overall average")+
      theme_minimal(base_size=13)+theme(plot.title=element_text(face="bold"))
  })
  
  output$plt_icc_bar <- renderPlot({
    vd <- var_decomp(); req(vd)
    
    pd <- vd %>%
      select(outcome, pct_school, pct_cohort, pct_residual) %>%
      pivot_longer(c(pct_school, pct_cohort),
                   names_to  = "level",
                   values_to = "pct") %>%
      mutate(level = factor(level,
                            levels = c("pct_school", "pct_cohort"),
                            labels = c("School", "Cohort / Class")))
    
    # Ceiling for y-axis: max of school+cohort combined, rounded up to next whole %
    y_max <- max(vd$pct_school + vd$pct_cohort, na.rm = TRUE)
    y_max <- max(ceiling(y_max) + 1, 5)   # at least 5% so thin bars are visible
    
    annot <- vd %>%
      mutate(label = sprintf("Individual: %.1f%%", pct_residual),
             x_pos = outcome,
             y_pos = pct_school + pct_cohort + y_max * 0.04)
    
    ggplot(pd, aes(x = outcome, y = pct, fill = level)) +
      geom_col(position = "stack", alpha = 0.88, width = 0.6) +

      geom_text(aes(label = sprintf("%.1f%%", pct)),
                position = position_stack(vjust = 0.5),
                colour = "white", fontface = "bold", size = 3.8) +

      geom_text(data = annot,
                aes(x = x_pos, y = y_pos, label = label),
                inherit.aes = FALSE,
                colour = "#888", size = 3.2, fontface = "italic") +
      scale_fill_manual(values = c(
        "School"         = "#e74c3c",
        "Cohort / Class" = "#f39c12")) +
      scale_y_continuous(
        labels = function(x) paste0(x, "%"),
        limits = c(0, y_max),
        expand = c(0, 0)) +
      labs(x = NULL, y = "% of Total Variance", fill = NULL,
           title = "Between-Level Variance Decomposition",
           subtitle = paste0(
             "Bars show school-level (red) and class-level (orange) variance only. ",
             "Grey text = remaining individual-child variance.")) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "#666", size = 10),
        legend.position = "bottom",
        axis.text.x   = element_text(angle = 15, hjust = 1, size = 10),
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_blank())
  })
  
  output$tbl_icc <- renderDT({
    vd<-var_decomp(); req(vd)
    disp<-vd%>%
      mutate(across(c(pct_school,pct_cohort,pct_residual),round,1),
             across(c(icc_school,icc_school_cohort),round,3))%>%
      select(outcome,icc_school,icc_school_cohort,pct_school,pct_cohort,pct_residual)
    names(disp)<-c("Fitness Test","ICC (School)","ICC (School+Class)",
                   "% Var: School","% Var: Class","% Var: Individual")
    datatable(disp,rownames=FALSE,options=list(dom="t",paging=FALSE))
  })
  
  # ── 6j. Tab 4: Teacher Summary ───────────────────────────
  
  # Populate school selector when data loads
  observe({
    df<-dat(); req(df,"school_id"%in%names(df))
    schools<-sort(unique(as.character(df$school_id)))
    output$ts_school_selector_ui <- renderUI({
      selectInput("ts_school","🏫 Select School:",
                  choices=c("All Schools"="__all__",schools),
                  width="100%")
    })
  })
  
  # Filtered dataset for the selected school
  dat_school <- reactive({
    df<-dat(); req(df)
    sel<-input$ts_school
    if(is.null(sel)||sel=="__all__") return(df)
    df[df$school_id==sel,]
  })
  
  output$ts_model_quality_ui <- renderUI({
    if(is.null(fit_res())) return(tags$div(
      style="background:#f4f6f9;border:1px dashed #bdc3c7;border-radius:6px;padding:8px 14px;color:#7f8c8d;font-size:13px;",
      "ℹ️ Fit models for statistical diagnostics."))
    s<-fit_res()$summary
    nc<-sum(s$converged==TRUE,na.rm=TRUE);nt<-sum(!is.na(s$converged))
    col<-if(nc==nt)"#d4edda" else "#fff3cd";bc<-if(nc==nt)"#28a745" else "#ffc107"
    icon_<-if(nc==nt)"✅" else "⚠️"
    tags$div(style=sprintf("background:%s;border-left:4px solid %s;border-radius:0 6px 6px 0;padding:8px 14px;font-size:13px;",col,bc),
             sprintf("%s Models: %d / %d converged.", icon_, nc, nt))
  })
  
  # ── Header + stat cards ───────────────────────────────────
  output$ts_header_ui <- renderUI({
    df  <- dat_school(); req(df)
    sel <- input$ts_school
    label <- if(is.null(sel)||sel=="__all__") "All Schools" else sel
    ocs   <- names(OUTCOME_LABELS)[names(OUTCOME_LABELS)%in%names(df)]
    n_ch  <- length(unique(df$child_id))
    n_co  <- length(unique(df$cohort_id))
    n_sc  <- length(unique(df$school_id))
    
    chips <- lapply(ocs, function(oc){
      y  <- df[[oc]][is.finite(df[[oc]])]
      mn <- mean(y); sdev <- sd(y); pz <- 100*mean(y==0)
      p10 <- quantile(y,.10); p90 <- quantile(y,.90)
      tags$div(
        style="background:white;border:1px solid #dee2e6;border-radius:8px;padding:14px;margin-bottom:12px;",
        tags$h5(OUTCOME_LABELS[oc],
                style="color:#2980b9;margin-top:0;border-bottom:2px solid #ecf0f1;padding-bottom:6px;"),
        tags$div(style="display:flex;gap:10px;flex-wrap:wrap;",
                 lapply(list(c("Mean",sprintf("%.1f",mn)),c("SD",sprintf("%.1f",sdev)),
                             c("10th %ile",sprintf("%.1f",p10)),c("90th %ile",sprintf("%.1f",p90))),
                        function(kv) tags$div(class="stat-chip",
                                              tags$div(kv[2],style="font-size:17px;font-weight:bold;color:#2c3e50;"),
                                              tags$div(kv[1],style="font-size:11px;color:#7f8c8d;")))),
        tags$small(paste("ℹ️",OUTCOME_DIRECTION[oc]),style="color:#95a5a6;"),
        if(oc=="sitreach"&&pz>5) tags$div(
          style="background:#fff3cd;border-left:4px solid #ffc107;padding:6px 10px;margin-top:8px;border-radius:0 4px 4px 0;font-size:12px;",
          sprintf("%.0f%% of children scored zero (cannot reach past their toes).",pz)))
    })
    
    tagList(
      tags$div(
        style="background:linear-gradient(135deg,#1a6496 0%,#1abc9c 100%);color:white;padding:18px 22px;border-radius:10px;margin-bottom:16px;",
        tags$h3(paste("📋 Fitness Assessment Summary —", label),
                style="margin:0 0 6px 0;font-size:19px;"),
        tags$p(sprintf("%s — %d children | %d classes | %s",
                       format(Sys.Date(),"%d %B %Y"),n_ch,n_co,
                       if(is.null(sel)||sel=="__all__") paste(n_sc,"schools") else label),
               style="margin:0;font-size:13px;opacity:.9;")),
      tagList(chips)
    )
  })
  
  # ── Cohort filter UI for teacher trajectory ──────────────
  output$ts_traj_cohort_ui <- renderUI({
    df <- dat_school(); req(df, "cohort_id" %in% names(df))
    cohorts <- sort(unique(as.character(df$cohort_id)))
    selectInput("ts_traj_cohort", "Class (optional):",
                choices = c("All Classes" = "__all__", cohorts), width = "100%")
  })
  
  # ── Teacher score-vs-age plot (replaces time-point trajectory) ──
  output$ts_traj_plot <- renderPlot({
    df  <- dat_school(); req(df)
    oc  <- input$ts_traj_outcome
    req(oc %in% names(df), "age" %in% names(df), "sex" %in% names(df))
    
    sel_co <- input$ts_traj_cohort
    if (!is.null(sel_co) && sel_co != "__all__")
      df <- df[df$cohort_id == sel_co, ]
    
    pd <- df %>%
      mutate(sex_label = ifelse(sex == 1, "Male", "Female")) %>%
      filter(is.finite(.data[[oc]]), is.finite(age))
    req(nrow(pd) > 5)
    
    co_lbl   <- if (!is.null(sel_co) && sel_co != "__all__") sel_co else "All Classes"
    sel_sc   <- input$ts_school
    sc_lbl   <- if (!is.null(sel_sc) && sel_sc != "__all__") sel_sc else "All Schools"
    subtitle <- sprintf("%s  |  %s  |  n = %d measurements", sc_lbl, co_lbl, nrow(pd))
    
    p <- ggplot(pd, aes(x = age, y = .data[[oc]], colour = sex_label))
    
    if (isTRUE(input$ts_traj_lines))
      p <- p + geom_line(aes(group = child_id), colour = "grey65",
                         alpha = .2, linewidth = .3, show.legend = FALSE)
    
    p <- p +
      geom_point(alpha = .12, size = 1.0) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1.4,
                  aes(fill = sex_label), alpha = .15) +
      scale_colour_manual(values = c(Male = "#2980b9", Female = "#e74c3c")) +
      scale_fill_manual(values   = c(Male = "#2980b9", Female = "#e74c3c")) +
      labs(x = "Age (years)", y = OUTCOME_LABELS[oc],
           colour = "Sex", fill = "Sex",
           title   = paste(OUTCOME_LABELS[oc], "vs Age — by Sex"),
           subtitle = subtitle) +
      theme_minimal(base_size = 13) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(colour = "#666", size = 10),
            legend.position = "bottom",
            panel.grid.minor = element_blank())
    p
  })
  
  # ── Student lookup: populate selectize from school filter ─
  observe({
    df  <- dat_school(); req(df, "child_id" %in% names(df))
    ids <- sort(unique(as.character(df$child_id)))
    updateSelectizeInput(session, "ts_child_id",
                         choices  = ids,
                         server   = TRUE)
  })
  
  observeEvent(input$ts_child_clear, {
    updateSelectizeInput(session, "ts_child_id", selected = "")
  })
  
  output$ts_child_plot <- renderPlot({
    df  <- dat_school(); req(df)
    oc  <- input$ts_child_outcome
    req(oc %in% names(df), "age" %in% names(df), "sex" %in% names(df))
    
    child_sel <- input$ts_child_id
    has_child <- !is.null(child_sel) && nchar(trimws(child_sel)) > 0
    
    pd_all <- df %>%
      mutate(sex_label = ifelse(sex == 1, "Male", "Female")) %>%
      filter(is.finite(.data[[oc]]), is.finite(age))
    req(nrow(pd_all) > 5)
    
    sel_sc <- input$ts_school
    sc_lbl <- if (!is.null(sel_sc) && sel_sc != "__all__") sel_sc else "All Schools"
    
    p <- ggplot(pd_all, aes(x = age, y = .data[[oc]])) +
      # Background cloud of all children (grey)
      geom_point(colour = "grey75", alpha = .12, size = 0.9) +
      # Loess by sex (context)
      geom_smooth(aes(colour = sex_label, fill = sex_label),
                  method = "loess", se = TRUE, linewidth = 1.1, alpha = .12) +
      scale_colour_manual(values = c(Male = "#2980b9", Female = "#e74c3c"),
                          aesthetics = c("colour","fill"))
    
    if (has_child) {
      pd_child <- df %>%
        filter(as.character(child_id) == child_sel,
               is.finite(.data[[oc]]), is.finite(age)) %>%
        mutate(sex_label = ifelse(sex == 1, "Male", "Female")) %>%
        arrange(age)
      
      if (nrow(pd_child) > 0) {
        child_sex_col <- if (pd_child$sex_label[1] == "Male") "#1a5276" else "#922b21"
        p <- p +
          # Child trajectory line
          geom_line(data = pd_child,
                    aes(x = age, y = .data[[oc]], group = child_id),
                    colour = child_sex_col, linewidth = 1.3, linetype = "solid") +
          # Child measurement points (large, dark)
          geom_point(data = pd_child,
                     aes(x = age, y = .data[[oc]]),
                     colour = child_sex_col, size = 4.5, shape = 21,
                     fill = child_sex_col, stroke = 1.5) +
          # Age labels for each measurement
          geom_text(data = pd_child,
                    aes(x = age, y = .data[[oc]],
                        label = sprintf("Age %.1f
%.1f", age, .data[[oc]])),
                    hjust = -0.15, size = 3.2, colour = child_sex_col, fontface = "bold")
        subtitle_child <- sprintf("Highlighted: Child %s  |  %d measurement(s)  |  %s",
                                  child_sel, nrow(pd_child), sc_lbl)
      } else {
        subtitle_child <- sprintf("Child '%s' not found in the selected school/class.", child_sel)
      }
    } else {
      subtitle_child <- sprintf("No child selected — showing all children in %s", sc_lbl)
    }
    
    p +
      labs(x = "Age (years)", y = OUTCOME_LABELS[oc],
           colour = "Sex", fill = "Sex",
           title   = paste("Individual Student Lookup —", OUTCOME_LABELS[oc]),
           subtitle = subtitle_child) +
      theme_minimal(base_size = 13) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(colour = "#555", size = 10),
            legend.position = "bottom",
            panel.grid.minor = element_blank())
  })
  
  # ── Within-school cohort boxplot ──────────────────────────
  output$ts_cohort_plot <- renderPlot({
    df  <- dat_school(); req(df)
    oc  <- input$ts_cohort_outcome; req(oc%in%names(df))
    pd  <- df%>%filter(is.finite(.data[[oc]])); req(nrow(pd)>0)
    n_cohorts <- length(unique(pd$cohort_id))
    pal_co <- if(n_cohorts<=length(PALETTE_COHORTS))
      PALETTE_COHORTS[seq_len(n_cohorts)] else scales::hue_pal()(n_cohorts)
    grand_mean <- mean(pd[[oc]],na.rm=TRUE)
    ggplot(pd,aes(x=reorder(as.character(cohort_id),.data[[oc]],median),
                  y=.data[[oc]],fill=cohort_id))+
      geom_hline(yintercept=grand_mean,linetype="dashed",colour="#7f8c8d",linewidth=.8)+
      geom_boxplot(outlier.size=1,outlier.alpha=.4,alpha=.8,linewidth=.4)+
      annotate("text",x=Inf,y=grand_mean,
               label=sprintf("Mean: %.1f",grand_mean),
               hjust=1.1,vjust=-0.4,size=3.5,colour="#7f8c8d")+
      scale_fill_manual(values=pal_co,guide="none")+
      labs(x="Class (sorted by median)",y=OUTCOME_LABELS[oc],
           title=paste("Class Score Distributions:",OUTCOME_LABELS[oc]),
           subtitle="Dashed = school-wide average")+
      theme_minimal(base_size=13)+
      theme(plot.title=element_text(face="bold"),
            axis.text.x=element_text(angle=25,hjust=1,size=9),
            panel.grid.major.x=element_blank())
  })
  
  # ── At-risk students table ────────────────────────────────
  output$ts_atrisk_tbl <- renderDT({
    df  <- dat(); req(df)
    sel <- input$ts_school
    ocs <- names(OUTCOME_LABELS)[names(OUTCOME_LABELS)%in%names(df)]
    school_arg <- if(is.null(sel)||sel=="__all__") NULL else sel
    ar <- flag_atrisk(df, ocs, school_filter=school_arg)
    if(is.null(ar)||nrow(ar)==0)
      return(datatable(
        data.frame(Message="No at-risk children identified — great news! 🎉"),
        rownames=FALSE,options=list(dom="t",paging=FALSE)))
    datatable(ar,rownames=FALSE,filter="top",
              options=list(pageLength=15,scrollX=TRUE,dom="tip",
                           order=list(list(7,"asc"))))%>%
      formatStyle("Z-score",
                  backgroundColor=styleInterval(c(-2,-1.5),c("#f8d7da","#fff3cd","white")))
  })
  
  # ── Data-driven recommendations ───────────────────────────
  output$ts_recs_ui <- renderUI({
    df  <- dat(); req(df)
    sel <- input$ts_school
    school_arg <- if(is.null(sel)||sel=="__all__") NULL else sel
    ocs <- names(OUTCOME_LABELS)[names(OUTCOME_LABELS)%in%names(df)]
    recs <- make_recommendations(df, var_decomp(), ocs, school_filter=school_arg)
    
    cards <- lapply(recs, function(r)
      tags$div(class="rec-card",
               style=sprintf("border-left-color:%s;",r$colour),
               tags$div(style="display:flex;align-items:center;gap:8px;margin-bottom:6px;",
                        tags$span(r$icon,style="font-size:20px;"),
                        tags$strong(r$title,style=sprintf("color:%s;font-size:14px;",r$colour))),
               tags$p(r$body,style="margin:0;font-size:13px;color:#444;line-height:1.7;")))
    
    tagList(
      tags$hr(),
      tags$h4("📌 Recommendations for PE Teachers & Curriculum Planners",
              style="margin-bottom:14px;"),
      tagList(cards))
  })
  
  # ── 6k. Downloads ────────────────────────────────────────
  output$dl_template <- downloadHandler(
    filename="fitness_data_template.csv",
    content=function(file) write.csv(data.frame(
      school_id=c("School_A","School_A","School_B"),
      cohort_id=c("A_Class1","A_Class1","B_Class1"),
      child_id=c("C001","C001","C101"),
      time=c(0L,1L,0L),age=c(8.0,9.0,7.5),sex=c(1L,1L,0L),
      bmi=c(16.5,17.0,14.9),sitreach=c(12.5,14.0,8.0),
      ropeskip=c(45L,52L,38L),dash50=c(10.2,9.8,12.1),
      vitalcap=c(1500,1650,1300),stringsAsFactors=FALSE),file,row.names=FALSE))
  
  output$dl_report <- downloadHandler(
    filename=function() paste0("fitness_report_",Sys.Date(),".html"),
    content=function(file) {
      df<-dat()
      if(is.null(df)){writeLines("No data uploaded.",file);return()}
      writeLines(c("<!DOCTYPE html><html><head><meta charset='utf-8'>",
                   "<title>Fitness Assessment Report</title>",
                   "<style>body{font-family:Arial,sans-serif;background:#f4f6f9;max-width:900px;margin:auto;padding:30px}</style>",
                   "</head><body>",
                   "<h1>Children's Fitness Assessment Report</h1>",
                   paste0("<p>Generated: ",format(Sys.Date(),"%d %B %Y"),"</p>"),
                   "<p>Open the dashboard for full interactive plots and class-level analysis.</p>",
                   "</body></html>"),file)
    })
}

# ── 7. Launch ─────────────────────────────────────────────────
shinyApp(ui, server)