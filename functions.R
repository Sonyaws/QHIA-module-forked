
# ----- Assign relative risks -----

assignPA_RRs <- function(df, diseases, mmet_col = "mmet_ref") {
  # mmet_col: name of the column containing mMET-hours/week values.
  #           Use "mmet_ref" for the reference scenario (default) or
  #           "mmet_scen" for the counterfactual scenario.
  
  for (i in seq_along(diseases)) {
    
    disease <- diseases[i]
    new_colname <- paste0("rr_PHYSICAL_ACTIVITY_", disease)
    
    df[[new_colname]] <- case_when(
      
      # 0 for male breast/endometrial cancer
      # sex is character "male"/"female" in the Bogotá dataset
      df$sex == "male" & disease %in% c("breast-cancer", "endometrial-cancer") ~ 0,
      
      # 0 for dementia if age < 60
      df$age < 60 & disease == "all-cause-dementia" ~ 0,
      
      # 1 for age < 18, except depression
      df$age < 18 & disease != "depression" ~ 1,
      
      # otherwise - calculate the dose response
      TRUE ~ drpa::dose_response(
        cause = disease,
        outcome_type = ifelse(disease == "diabetes", "non-fatal", "fatal-and-non-fatal"),
        dose = df[[mmet_col]],   # use the column name supplied by the caller
        quantile = 0.5,          # deterministic (median)
        censor_method = "75thPercentile",
        confidence_intervals = FALSE
      )$rr
    )
  }
  
  return(df)
}


# ----- Calculate PIFs -----

# Computes the Population Impact Fraction (PIF) for all-cause mortality,
# comparing a reference PA distribution to a counterfactual PA distribution.
#
# PIF = (E[RR(ref)] - E[RR(scen)]) / E[RR(ref)]
#
# Arguments:
#   pop_ref  — output of assignPA_RRs(..., mmet_col = "mmet_ref")
#   pop_scen — output of assignPA_RRs(..., mmet_col = "mmet_scen")
#
# Both data frames must have a column named rr_PHYSICAL_ACTIVITY_all-cause-mortality,
# a participant_id column, and age_cat and sex columns.
#
# Returns a list:
#   $pif_table   — PIF by sex and age_cat for all-cause mortality
#   $pop_combined — individual-level ref and scen RRs joined (for inspection)

calculate_pif <- function(pop_ref, pop_scen) {
  
  # Pivot PA RR columns to long form, keeping only all-cause-mortality
  prep_pop <- function(pop) {
    pop |>
      dplyr::select(participant_id, age, sex, age_cat,
                    starts_with("rr_PHYSICAL_ACTIVITY_")) |>
      pivot_longer(
        cols      = starts_with("rr_PHYSICAL_ACTIVITY_"),
        names_to  = "cause",
        values_to = "rr_pa"
      ) |>
      mutate(
        cause = sub("rr_PHYSICAL_ACTIVITY_", "", cause),
        cause = chartr("-", "_", cause))   # e.g. all-cause-mortality -> all_cause_mortality
       }
  
  # pop_ref <- bogota_ref
  # pop_scen <- bogota_scen
  
  ref_long  <- prep_pop(pop_ref)
  scen_long <- prep_pop(pop_scen)
  
  # Join reference and scenario RRs by individual
  pop_combined <- ref_long |>
    dplyr::rename(rr_ref = rr_pa) |>
    dplyr::left_join(
      scen_long |> dplyr::select(participant_id, cause, rr_scen = rr_pa),
      by = c("participant_id", "cause")
    )
  
  # PIF by sex and age_cat
  # Formula: (sum[RR_ref] - sum[RR_scen]) / sum[RR_ref]

  pif_table <- pop_combined |>
    group_by(sex, age_cat, cause) |>
    summarise(
      n            = n(),
      sum_rr_ref  = sum(rr_ref,  na.rm = TRUE),
      sum_rr_scen = sum(rr_scen, na.rm = TRUE),
      pif_pa       = (sum_rr_ref - sum_rr_scen) / sum_rr_ref,
      .groups = "drop"
    ) |>
    # Guard against NaN / Inf (e.g. if mean_rr_ref is 0 or NA)
    mutate(pif_pa = replace(pif_pa, !is.finite(pif_pa), 0))
  
  return(list(pif_table = pif_table, pop_combined = pop_combined))
}

# ----- Run life table -----

### Function to create a life table for a age and sex cohort

RunLifeTable <- function(in_data, in_sex, in_mid_age, mx_trend = NA, scenario = 0) {
  # in_data    = mslt_general
  # in_sex     = "male" or "female"
  # in_mid_age = cohort entry age (e.g. 17, 22, ..., 97)
  # mx_trend   = annual % change in mortality across simulation years
  #              NA  → no trend, use static mx
  #              2   → 2% reduction per year  (improving survival)
  #             -1   → 1% increase per year   (worsening mortality)
  # scenario   = one-off proportional reduction in mx (e.g. PIF); 0 = no change
  
  # ── Base life table data ─────────────────────────────────────────────────────
  lf_df <- in_data %>%
    dplyr::filter(age >= in_mid_age & sex == in_sex) %>%
    dplyr::select(sex, age, pyld_rate, mx)
  
  num_row <- nrow(lf_df)
  
  # ── Apply mortality trend (compounding across simulation years) ───────────────
  # Year 1 (cohort entry) is unchanged; each subsequent year compounds.
  # Positive mx_trend = reduction; negative = increase.
  if (!is.na(mx_trend)) {
    trend_multiplier <- (1 - mx_trend / 100) ^ (seq_len(num_row) - 1)
    lf_df <- lf_df %>%
      mutate(mx = mx * trend_multiplier)
  }
  
  # ── Apply scenario reduction (e.g. from PIF) on top of any trend ─────────────
  lf_df <- lf_df %>%
    mutate(mx = mx * (1 - scenario))
  
  # ── Life table calculations ──────────────────────────────────────────────────
  
  # probability of dying
  qx <- ifelse(lf_df$age < 100, 1 - exp(-lf_df$mx), 1)
  
  # number of survivors — seed from population at cohort entry age
  lx <- rep(0, num_row)
  lx[1] <- as.numeric(
    in_data$population[in_data$age == in_mid_age & in_data$sex == in_sex][1]
  )
  
  # deaths and survivors
  dx <- rep(0, num_row)
  dx[1] <- lx[1] * qx[1]
  for (i in 2:num_row) {
    lx[i] <- lx[i - 1] - dx[i - 1]
    dx[i] <- lx[i] * qx[i]
  }
  
  # person-years lived
  Lx <- rep(0, num_row)
  for (i in 1:(num_row - 1))
    Lx[i] <- (lx[i] + lx[i + 1]) / 2
  # terminal age: guard against mx == 0
  Lx[num_row] <- lx[num_row] * qx[num_row]
  
  # life expectancy — guard against lx == 0
  ex <- rep(0, num_row)
  for (i in 1:num_row)
    ex[i] <- ifelse(lx[i] > 0, sum(Lx[i:num_row]) / lx[i], 0)
  
  # health-adjusted person-years
  Lwx <- Lx * (1 - lf_df$pyld_rate)
  
  # health-adjusted life expectancy
  ewx <- rep(0, num_row)
  for (i in 1:num_row)
    ewx[i] <- ifelse(lx[i] > 0, sum(Lwx[i:num_row]) / lx[i], 0)
  
  lf_df$qx  <- qx
  lf_df$lx  <- lx
  lf_df$dx  <- dx
  lf_df$Lx  <- Lx
  lf_df$ex  <- ex
  lf_df$Lwx <- Lwx
  lf_df$ewx <- ewx
  lf_df
}

# ----- Run disease process -----
## Function to generate disease process for an age and sex cohort
## Based on formulas in Barendregt JJ, Oortmarssen GJ, van, Vos T, Murray CJL. A generic model for the assessment of
## disease epidemiology: the computational basis of DisMod II. Population Health Metrics.2003;1(1):4.

RunDisease <- function(in_data,
                       in_sex,
                       in_mid_age,
                       in_disease,
                       inc_trend    = NA,   # % annual change in incidence (compounding)
                       cf_trend     = NA,   # % annual change in case fatality (compounding)
                       scenario_inc = 0,    # one-off proportional reduction in incidence
                       scenario_cf  = 0     # one-off proportional reduction in case fatality
) {
  
  # ── 1. Filter data to the requested sex and cohort entry age ───────────────
  in_data_f <- in_data %>%
    filter(sex == in_sex, age >= in_mid_age)
  
  num_rows <- nrow(in_data_f)
  
  # ── 2. Extract base disease rates from the data ────────────────────────────
  incidence_base     <- in_data_f[[paste0("incidence_",     in_disease)]]
  case_fatality_base <- in_data_f[[paste0("case_fatality_", in_disease)]]
  dw_disease         <- in_data_f[[paste0("dw_adj_",        in_disease)]]
  pyld_rate          <- in_data_f$pyld_rate
  
  # ── 3. Apply compounding trends (mirrors mx_trend in RunLifeTable) ─────────
  # t = 1 at cohort entry → multiplier (1 − trend/100)^0 = 1 (no change)
  # t = 2 one year later  → multiplier (1 − trend/100)^1, etc.
  t_seq <- seq_len(num_rows) - 1L  # 0, 1, 2, ..., n−1
  
  if (!is.na(inc_trend) && inc_trend != 0) {
    incidence_disease <- incidence_base * (1 - inc_trend / 100) ^ t_seq
  } else {
    incidence_disease <- incidence_base
  }
  
  if (!is.na(cf_trend) && cf_trend != 0) {
    case_fatality_disease <- case_fatality_base * (1 - cf_trend / 100) ^ t_seq
  } else {
    case_fatality_disease <- case_fatality_base
  }
  
  # ── 4. Apply one-off scenario reductions ──────────────────────────────────
  incidence_disease     <- incidence_disease     * (1 - scenario_inc)
  case_fatality_disease <- case_fatality_disease * (1 - scenario_cf)
  
  # ── 5. Initialise state vectors ───────────────────────────────────────────
  # Starting cohort of 1,000 — all healthy at entry age
  Sx <- numeric(num_rows)
  Cx <- numeric(num_rows)
  Dx <- numeric(num_rows)
  
  Sx[1] <- 1000
  Cx[1] <- 0
  Dx[1] <- 0
  
  # ── 6. Propagate the three-state disease model ────────────────────────────
  # Based on Barendregt et al. (1998) and DisMod II (Barendregt et al. 2003)
  # Within each 1-year interval, hazards are assumed constant.
  #
  #  Healthy → Diseased  governed by incidence hazard γ
  #  Diseased → Dead     governed by case fatality hazard φ
  #  Sx, Cx, Dx are proportions of 1,000
  #
  # The double-transition probability r_a (healthy → diseased → dead within
  # the same year) follows Barendregt 1998 eq. 15:
  #   r_a = [φ(1−e^{−γ}) − γ(1−e^{−φ})] / (φ − γ)  when γ ≠ φ
  #   r_a = 1 − e^{−γ} − γ·e^{−γ}                   when γ = φ
  
  for (i in seq_len(num_rows - 1L)) {
    gamma_i <- incidence_disease[i]
    phi_i   <- case_fatality_disease[i]
    
    # Probability of incidence within the year
    prob_inc <- 1 - exp(-gamma_i)
    
    # Double-transition probability (healthy → diseased → dead within 1 yr)
    if (abs(phi_i - gamma_i) > 1e-10) {
      r_a <- (phi_i * (1 - exp(-gamma_i)) - gamma_i * (1 - exp(-phi_i))) /
        (phi_i - gamma_i)
    } else {
      r_a <- 1 - exp(-gamma_i) - gamma_i * exp(-gamma_i)
    }
    r_a <- max(r_a, 0)
    
    # Probability of dying from disease within the year (for current prevalent)
    prob_cf <- 1 - exp(-phi_i)
    
    # Update states (following Barendregt 1998 eqs 16−18)
    Sx[i + 1] <- Sx[i] * exp(-gamma_i)
    Cx[i + 1] <- Sx[i] * (1 - exp(-gamma_i) - r_a) + Cx[i] * exp(-phi_i)
    Dx[i + 1] <- Sx[i] * r_a + Cx[i] * (1 - exp(-phi_i)) + Dx[i]
  }
  
  # Ensure non-negative states (rounding protection)
  Sx <- pmax(Sx, 0)
  Cx <- pmax(Cx, 0)
  
  Tx <- Sx + Cx + Dx  # should always equal 1,000
  
  # ── 7. Derive disease epidemiology outputs ────────────────────────────────
  # Person-years at risk (half-cycle correction)
  PYx <- 0.5 * (Sx + Cx + c(Sx[-1], Sx[num_rows]) + c(Cx[-1], Cx[num_rows]))
  
  # Prevalence: proportion diseased among alive
  alive <- Sx + Cx
  px    <- ifelse(alive > 0, Cx / alive, 0)
  
  # Disease-specific mortality rate (deaths per person-year alive)
  # Numerator: expected deaths from disease in interval
  deaths_disease <- Cx * (1 - exp(-case_fatality_disease))
  mx_disease     <- ifelse(PYx > 0, deaths_disease / PYx, 0)
  
  # Incidence numbers (cases per person-year among non-prevalent)
  inc_num <- Sx * (1 - exp(-incidence_disease))
  
  # ── 8. Return as data frame ───────────────────────────────────────────────
  out <- data.frame(
    age                   = in_data_f$age,
    sex                   = in_sex,
    disease               = in_disease,
    Sx                    = Sx,
    Cx                    = Cx,
    Dx                    = Dx,
    Tx                    = Tx,
    PYx                   = PYx,
    px                    = px,
    mx                    = mx_disease,
    inc_num               = inc_num,
    incidence_disease     = incidence_disease,
    case_fatality_disease = case_fatality_disease,
    dw                    = dw_disease
  )
  
  return(out)
}


# ----- Run Model -----
CalculationModel <- function(
    in_data,
    in_sex,
    in_mid_age,
    in_disease,              # character vector of disease short-names
    mx_trend       = NA,     # % annual change in all-cause mx (compounding)
    inc_trend      = NA,     # % annual change in disease incidence (compounding)
    cf_trend       = NA,     # % annual change in case fatality (compounding)
    scenario       = 0,      # direct PIF on all-cause mx (like RunLifeTable)
    scenario_inc   = 0,      # scalar OR named vector of incidence PIFs per disease
    scenario_cf    = 0,      # scalar OR named vector of case-fatality PIFs per disease
    use_disease_mx = TRUE    # if TRUE  (default): disease-specific mortality changes
                             #   from Step 4 feed into general life table mx.
                             # if FALSE: disease mortality feedback is suppressed —
                             #   mx in the general LT is modified only by 'scenario'
                             #   (direct all-cause PIF). Disease sections still run
                             #   with scenario_inc so case counts are produced, but
                             #   the life-year effect enters through the all-cause
                             #   PIF rather than the disease-mortality pathway.
                             #   pyld is still updated from disease prevalence changes.
) {
  
  # ── Helper: resolve per-disease scenario values ────────────────────────────
  # If scenario_inc / scenario_cf is a single scalar, replicate it for every
  # disease. If it is a named vector, validate that all disease names are present.
  resolve_disease_scenario <- function(sc_val, diseases) {
    if (length(sc_val) == 1) {
      # scalar → same value for every disease
      out <- setNames(rep(sc_val, length(diseases)), diseases)
    } else {
      # named vector — check all diseases are covered
      missing <- setdiff(diseases, names(sc_val))
      if (length(missing) > 0) {
        stop("scenario_inc/scenario_cf: missing disease(s) in named vector: ",
             paste(missing, collapse = ", "),
             "\nProvide either a single scalar or a named vector covering all diseases: ",
             paste(diseases, collapse = ", "))
      }
      out <- sc_val[diseases]
    }
    out
  }
  
  sc_inc_vec <- resolve_disease_scenario(scenario_inc, in_disease)
  sc_cf_vec  <- resolve_disease_scenario(scenario_cf,  in_disease)
  
  # ── Step 1: Baseline general life table ───────────────────────────────────
  general_lt_bl <- RunLifeTable(
    in_data    = in_data,
    in_sex     = in_sex,
    in_mid_age = in_mid_age,
    mx_trend   = mx_trend,
    scenario   = 0          # no direct scenario — baseline
  )
  message("Step 1 complete: baseline general life table")
  
  # ── Step 2: Baseline disease life tables ──────────────────────────────────
  disease_lt_list_bl <- setNames(
    lapply(in_disease, function(dis) {
      RunDisease(
        in_data      = in_data,
        in_sex       = in_sex,
        in_mid_age   = in_mid_age,
        in_disease   = dis,
        inc_trend    = inc_trend,
        cf_trend     = cf_trend,
        scenario_inc = 0,
        scenario_cf  = 0
      )
    }),
    in_disease
  )
  message("Step 2 complete: baseline disease life tables")
  
  # ── Step 3: Scenario disease life tables ──────────────────────────────────
  disease_lt_list_sc <- setNames(
    lapply(in_disease, function(dis) {
      RunDisease(
        in_data      = in_data,
        in_sex       = in_sex,
        in_mid_age   = in_mid_age,
        in_disease   = dis,
        inc_trend    = inc_trend,
        cf_trend     = cf_trend,
        scenario_inc = sc_inc_vec[[dis]],   # disease-specific PIF (or shared scalar)
        scenario_cf  = sc_cf_vec[[dis]]
      ) %>%
        mutate(
          diff_mort_disease  = mx - disease_lt_list_bl[[dis]]$mx,
          diff_pylds_disease = (px - disease_lt_list_bl[[dis]]$px) * dw
        )
    }),
    in_disease
  )
  message("Step 3 complete: scenario disease life tables")
  
  # ── Step 4: Aggregate disease-level changes in mx and pyld ────────────────
  disease_lt_sc_all <- bind_rows(disease_lt_list_sc)
  
  mx_pylds_sc_summary <- disease_lt_sc_all %>%
    group_by(age) %>%
    summarise(
      mortality_sum = sum(diff_mort_disease,  na.rm = TRUE),
      pylds_sum     = sum(diff_pylds_disease, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ── Step 5: Recalculate general life table with modified mx and pyld ───────
  #
  # use_disease_mx = TRUE  (default):
  #   Disease-specific mortality changes from Step 4 are added to baseline mx,
  #   then the direct all-cause 'scenario' PIF is applied on top.
  #   mx_modified = (mx_base + Σ Δmx_disease) × (1 − scenario)
  #
  # use_disease_mx = FALSE:
  #   Disease mortality feedback is suppressed — the general LT mx is modified
  #   only by the direct 'scenario' PIF. Disease incidence changes (scenario_inc)
  #   still run and produce case counts, but their mortality effect does not
  #   propagate to the general LT. pyld is still updated from disease prevalence.
  #   mx_modified = mx_base × (1 − scenario)
  #
  # Comparing use_disease_mx = TRUE vs FALSE (with identical scenario_inc and
  # scenario = pif_acm) isolates whether routing PA benefit through the disease
  # pathway changes life-year estimates compared to applying the all-cause PIF.

  if (use_disease_mx) {
    in_data_sc <- in_data %>%
      filter(sex == in_sex) %>%
      left_join(mx_pylds_sc_summary, by = "age") %>%
      mutate(
        mortality_sum = replace_na(mortality_sum, 0),
        pylds_sum     = replace_na(pylds_sum,     0),
        mx        = (mx + mortality_sum) * (1 - scenario),
        pyld_rate = pmax(pyld_rate + pylds_sum, 0)
      ) %>%
      select(-mortality_sum, -pylds_sum)

  } else {
    # Disease mortality feedback suppressed: mx modified by direct PIF only.
    # pyld still updated from disease prevalence changes.
    in_data_sc <- in_data %>%
      filter(sex == in_sex) %>%
      left_join(mx_pylds_sc_summary, by = "age") %>%
      mutate(
        pylds_sum = replace_na(pylds_sum, 0),
        mx        = mx * (1 - scenario),
        pyld_rate = pmax(pyld_rate + pylds_sum, 0)
      ) %>%
      select(-mortality_sum, -pylds_sum)
  }

  general_lt_sc <- RunLifeTable(
    in_data    = in_data_sc,
    in_sex     = in_sex,
    in_mid_age = in_mid_age,
    mx_trend   = mx_trend,
    scenario   = 0    # scenario already baked into in_data_sc
  )
  message("Step 5 complete: scenario general life table")
  
  # ── Step 6: Build output data frame ───────────────────────────────────────
  
  # 6a. Disease-level outputs (wide on disease)
  disease_bl_df <- bind_rows(disease_lt_list_bl) %>%
    select(sex, age, disease, incidence_disease, mx, px)
  
  disease_sc_df <- bind_rows(disease_lt_list_sc) %>%
    select(sex, age, disease, incidence_disease, mx, px)
  
  disease_wide <- inner_join(
    disease_sc_df %>% rename_with(~ paste0(., "_sc"), -c(sex, age, disease)),
    disease_bl_df %>% rename_with(~ paste0(., "_bl"), -c(sex, age, disease)),
    by = c("sex", "age", "disease")
  ) %>%
    inner_join(general_lt_sc %>% select(sex, age, Lx, Lwx), by = c("sex", "age")) %>%
    rename(Lx_sc = Lx, Lwx_sc = Lwx) %>%
    inner_join(general_lt_bl %>% select(sex, age, Lx, Lwx), by = c("sex", "age")) %>%
    rename(Lx_bl = Lx, Lwx_bl = Lwx) %>%
    mutate(
      # Incidence numbers = incidence rate × (1 − prevalence) × person-years
      inc_num_bl   = incidence_disease_bl * (1 - px_bl) * Lx_bl,
      inc_num_sc   = incidence_disease_sc * (1 - px_sc) * Lx_sc,
      inc_num_diff = inc_num_sc - inc_num_bl,
      # Disease-specific deaths = mortality rate × person-years
      mx_num_bl    = mx_bl * Lx_bl,
      mx_num_sc    = mx_sc * Lx_sc,
      mx_num_diff  = mx_num_sc - mx_num_bl
    ) %>%
    select(sex, age, disease,
           inc_num_diff, mx_num_diff,
           incidence_disease_sc, incidence_disease_bl,
           mx_sc, mx_bl, px_sc, px_bl) %>%
    # Pivot disease to columns so each row = one age × sex
    pivot_wider(
      names_from  = disease,
      values_from = c(inc_num_diff, mx_num_diff,
                      incidence_disease_sc, incidence_disease_bl,
                      mx_sc, mx_bl, px_sc, px_bl)
    )
  
  # 6b. General life table differences
  general_lf <- inner_join(
    general_lt_sc %>% select(sex, age, Lx, ex, Lwx, ewx) %>%
      rename_with(~ paste0(., "_sc"), -c(sex, age)),
    general_lt_bl %>% select(sex, age, Lx, ex, Lwx, ewx) %>%
      rename_with(~ paste0(., "_bl"), -c(sex, age)),
    by = c("sex", "age")
  ) %>%
    mutate(
      Lx_diff  = Lx_sc  - Lx_bl,
      Lwx_diff = Lwx_sc - Lwx_bl,
      ex_diff  = ex_sc  - ex_bl,
      ewx_diff = ewx_sc - ewx_bl
    )
  
  # 6c. Combined output
  output_df <- inner_join(disease_wide, general_lf, by = c("sex", "age"))
  
  output_df
}
