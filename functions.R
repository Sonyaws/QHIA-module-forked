
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
      pif       = (sum_rr_ref - sum_rr_scen) / sum_rr_ref,
      .groups = "drop"
    ) |>
    # Guard against NaN / Inf (e.g. if mean_rr_ref is 0 or NA)
    mutate(pif = replace(pif, !is.finite(pif), 0))
  
  return(list(pif_table = pif_table, pop_combined = pop_combined))
}

# ----- Run life table -----

### Function to create a life table for a age and sex cohort

RunLifeTable <- function(in_data, in_sex, in_mid_age, mx_trend = NA,
                         mx_trend_years = NA, scenario = 0) {
  # in_data         = mslt_general
  # in_sex          = "male" or "female"
  # in_mid_age      = cohort entry age (e.g. 17, 22, ..., 97)
  # mx_trend        = annual % change in mortality across simulation years
  #                   NA  → no trend, use static mx
  #                   2   → 2% reduction per year  (improving survival)
  #                  -1   → 1% increase per year   (worsening mortality)
  # mx_trend_years  = number of years to apply mx_trend (NA = no cap, trend
  #                   applies for entire simulation)
  # scenario        = one-off proportional reduction in mx. Can be:
  #                   - scalar: same PIF for all ages (e.g. 0.05)
  #                   - data.frame: age-specific PIFs with columns:
  #                     sex, age_cat, and a PIF column (e.g. pif_pa).
  #                     The cohort uses PIF for their current age/sex as they age.
  
  # ── Base life table data ─────────────────────────────────────────────────────
  lf_df <- in_data %>%
    dplyr::filter(age >= in_mid_age & sex == in_sex) %>%
    dplyr::select(sex, age, pyld_rate, mx)
  
  num_row <- nrow(lf_df)
  
  # ── Apply mortality trend (compounding across simulation years) ───────────────
  # Year 1 (cohort entry) is unchanged; each subsequent year compounds.
  # Positive mx_trend = reduction; negative = increase.
  # mx_trend_years caps how many years the trend applies (NA = no cap).
  if (!is.na(mx_trend)) {
    t_seq <- seq_len(num_row) - 1
    if (!is.na(mx_trend_years)) {
      t_seq <- pmin(t_seq, mx_trend_years) #  caps the time index so the trend stops increasing after that many years.
    }
    trend_multiplier <- (1 - mx_trend / 100) ^ t_seq
    lf_df <- lf_df %>%
      mutate(mx = mx * trend_multiplier)
  }
  
  # ── Apply scenario reduction (e.g. from PIF) on top of any trend ─────────────
  # scenario can be:
  #   - scalar: single PIF applied to all ages (e.g. 0.05)
  #   - data.frame: age-specific PIFs with columns sex, age_cat, and a PIF column
  if (is.data.frame(scenario)) {
    # lf_df <- lt_ref
    pif_df <- scenario
    pif_col <- setdiff(names(pif_df), c("sex", "age_cat", "cause", "n", "sum_rr_ref", "sum_rr_scen"))
    if (length(pif_col) == 0) {
      stop("scenario data.frame must contain a PIF column (e.g. pif_pa)")
    }
    pif_col <- pif_col[1]
    lf_df <- lf_df %>%
      mutate(age_cat = paste0(floor(age / 5) * 5, "-", floor(age / 5) * 5 + 4)) %>%
      left_join(
        pif_df %>% filter(sex == in_sex) %>% select(age_cat, pif = all_of(pif_col)),
        by = "age_cat"
      ) %>%
      mutate(pif = as.numeric(ifelse(is.na(pif), 0, pif))) %>%
      mutate(mx = mx * (1 - pif)) %>%
      select(-age_cat)
  } else {
    lf_df <- lf_df %>%
      mutate(mx = mx * (1 - scenario))
  }
  
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
    Lx[i] = (lx[i] + lx[i + 1]) / 2
  # terminal age: use GitHub formula
  Lx[num_row] = lx[num_row] / lf_df$mx[num_row]
  
  # life expectancy — guard against lx == 0
  ex <- rep(0, num_row)
  for (i in 1:num_row)
    ex[i] = sum(Lx[i:num_row]) / lx[i]
  
  # health-adjusted person-years
  Lwx <- Lx * (1 - lf_df$pyld_rate)
  
  # health-adjusted life expectancy
  ewx <- rep(0, num_row)
  for (i in 1:num_row)
    ewx[i] = sum(Lwx[i:num_row]) / lx[i]
  
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
## Based on formulas in Barendregt JJ, Oortmarssen GJ, van, Vos T, Murray CJL. 
## A generic model for the assessment of disease epidemiology: the computational basis of DisMod II. 
## Population Health Metrics.2003;1(1):4.

RunDisease <- function(in_data,
                       in_sex,
                       in_mid_age,
                       in_disease,
                       inc_trend       = NA,  # % annual change in incidence (compounding)
                       cf_trend        = NA,  # % annual change in case fatality (compounding)
                       rem_trend       = NA,  # % annual change in remission (compounding)
                       inc_trend_years = NA,  # years to apply inc_trend (NA = no cap)
                       cf_trend_years  = NA,  # years to apply cf_trend  (NA = no cap)
                       rem_trend_years = NA,  # years to apply rem_trend (NA = no cap)
                       scenario_inc    = 0,   # one-off proportional reduction in incidence
                       scenario_cf     = 0,   # one-off proportional reduction in case fatality
                       scenario_rem    = 0    # one-off proportional reduction in remission
) {
  
  # ── 1. Filter data to the requested sex and cohort entry age ───────────────
  in_data_f <- in_data %>%
    filter(sex == in_sex, age >= in_mid_age)
  
  num_rows <- nrow(in_data_f)
  
  # ── 2. Extract base disease rates from the data ────────────────────────────
  dw_disease         <- in_data_f[[paste0("dw_adj_",        in_disease)]]
  incidence_disease  <- in_data_f[[paste0("incidence_",     in_disease)]]
  case_fatality_base <- in_data_f[[paste0("case_fatality_", in_disease)]]
  remission_base     <- in_data_f[[paste0("remission_",      in_disease)]]
  remission_base[is.na(remission_base)] <- 0
  
  # ── 3. Apply compounding trends ─────────────────────────────────────────────
  t_seq <- seq_len(num_rows) - 1L
  
  if (!is.na(inc_trend) && inc_trend != 0) {
    t_inc <- t_seq
    if (!is.na(inc_trend_years)) t_inc <- pmin(t_inc, inc_trend_years)
    incidence_disease <- incidence_disease * (1 - inc_trend / 100) ^ t_inc
  }
  
  if (!is.na(cf_trend) && cf_trend != 0) {
    t_cf <- t_seq
    if (!is.na(cf_trend_years)) t_cf <- pmin(t_cf, cf_trend_years)
    case_fatality_base <- case_fatality_base * (1 - cf_trend / 100) ^ t_cf
  }
  
  if (!is.na(rem_trend) && rem_trend != 0) {
    t_rem <- t_seq
    if (!is.na(rem_trend_years)) t_rem <- pmin(t_rem, rem_trend_years)
    remission_base <- remission_base * (1 - rem_trend / 100) ^ t_rem
  }
  
  # ── 4. Apply scenario reductions ────────────────────────────────────────────
  incidence_disease    <- incidence_disease  * (1 - scenario_inc)
  case_fatality_disease <- case_fatality_base * (1 - scenario_cf)
  remission_disease    <- remission_base     * (1 - scenario_rem)
  
  # ── 5. Barendregt 1998 formulas with remission ────────────────────────────────
  # lx = incidence + case_fatality + remission (total hazard leaving diseased state)
  # qx = sqrt((incidence + remission - case_fatality)^2 + 4*incidence*remission)
  lx <- incidence_disease + case_fatality_disease + remission_disease
  qx <- sqrt((incidence_disease + remission_disease - case_fatality_disease)^2 + 
             4 * incidence_disease * remission_disease)
  wx <- exp(-1 * (lx + qx) / 2)
  vx <- exp(-1 * (lx - qx) / 2)
  
  # ── 6. Initialise state vectors ───────────────────────────────────────────
  number_of_ages <- num_rows
  Sx <- Cx <- Dx <- Tx <- Ax <- PYx <- px <- mx <- rep(0, number_of_ages)
  cfds <- case_fatality_disease
  
  Sx[1] <- 1000
  Ax[1] <- 1000
  
  # ── 7. Propagate three-state model using Barendregt 1998 ───────────────────
  for (i in 2:(number_of_ages - 1)) {
    if (qx[i - 1] > 0) {
      vxmwx <- vx[i - 1] - wx[i - 1]
      SxpCx <- Sx[i - 1] + Cx[i - 1]
      dqx <- 2 * qx[i - 1]
      qxmlx <- qx[i - 1] - lx[i - 1]
      qxplx <- qx[i - 1] + lx[i - 1]
      
      Sx[i] <- Sx[i - 1] * (2 * vxmwx * cfds[i - 1] + 
                            (vx[i - 1] * qxmlx + wx[i - 1] * qxplx)) / dqx
      Cx[i] <- -1 * (vxmwx * (2 * (cfds[i - 1] * SxpCx - lx[i - 1] * Sx[i - 1]) - 
                              Cx[i - 1] * lx[i - 1]) - 
                     Cx[i - 1] * qx[i - 1] * (vx[i - 1] + wx[i - 1])) / dqx
      Dx[i] <- (vxmwx * (2 * cfds[i - 1] * Cx[i - 1] - lx[i - 1] * SxpCx) - 
                qx[i - 1] * SxpCx * (vx[i - 1] + wx[i - 1]) + 
                dqx * (SxpCx + Dx[i - 1])) / dqx
    } else {
      Sx[i] <- Sx[i - 1]
      Cx[i] <- Cx[i - 1]
      Dx[i] <- Dx[i - 1]
    }
  }
  
  # GitHub leaves last element at initial 0 (no special handling needed)
  # Sx, Cx, Dx were initialized to 0, so last element stays 0
  
  Tx   <- Sx + Cx + Dx
  Ax <- Sx + Cx
  
  # ── 8. Derive PYx, px, mx ──────────────────────────────────────────────────
  first_indices <- 1:(number_of_ages - 1)
  last_indices <- 2:number_of_ages
  PYx <- (Ax[first_indices] + Ax[last_indices]) / 2
  mx[first_indices] <- (Dx[last_indices] - Dx[first_indices]) / PYx[first_indices]
  mx[mx < 0] <- 0
  px[first_indices] <- (Cx[last_indices] + Cx[first_indices]) / 2 / PYx[first_indices]
  
  # Pad px, mx, and PYx to match full length (add NA for last age)
  px[number_of_ages] <- NA
  mx[number_of_ages] <- NA
  PYx[number_of_ages] <- NA
  
  # ── 9. Return as data frame ────────────────────────────────────────────────
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
    mx                    = mx,
    incidence_disease     = incidence_disease,
    case_fatality_disease = case_fatality_disease,
    remission_disease     = remission_disease,
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
    mx_trend_years  = NA,    # years to apply mx_trend  (NA = no cap)
    inc_trend_years = NA,    # years to apply inc_trend (NA = no cap)
    cf_trend_years  = NA,    # years to apply cf_trend  (NA = no cap)
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
    in_data         = in_data,
    in_sex          = in_sex,
    in_mid_age      = in_mid_age,
    mx_trend        = mx_trend,
    mx_trend_years  = mx_trend_years,
    scenario        = 0          # no direct scenario — baseline
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
        rem_trend    = NA,
        scenario_inc = 0,
        scenario_cf  = 0,
        scenario_rem = 0
      )
    }),
    in_disease
  )
  message("Step 2 complete: baseline disease life tables")
  
  # ── Step 3: Scenario disease life tables ──────────────────────────────────
  disease_lt_list_sc <- setNames(
    lapply(in_disease, function(dis) {
      RunDisease(
        in_data         = in_data,
        in_sex          = in_sex,
        in_mid_age      = in_mid_age,
        in_disease      = dis,
        inc_trend       = inc_trend,
        cf_trend        = cf_trend,
        rem_trend       = NA,
        inc_trend_years = inc_trend_years,
        cf_trend_years  = cf_trend_years,
        scenario_inc    = sc_inc_vec[[dis]],   # disease-specific PIF (or shared scalar)
        scenario_cf     = sc_cf_vec[[dis]],
        scenario_rem    = 0
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
