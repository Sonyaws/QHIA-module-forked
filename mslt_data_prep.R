# =============================================================================
# MSLT Data Preparation — Colombia, GBD 2023
# =============================================================================
# Produces a dataset with one row per age (0–100) × sex containing:
#   - all-cause mortality rate (mx) and population
#   - per-disease: deaths_rate, ylds_rate, dw_adj, incidence, prevalence,
#                  remission, case_fatality
#
# Pipeline:
#   1. Load & reshape GBD data
#   2. Compute disability weights (pyld_rate, dw_adj)
#   3. Prepare disbayes inputs (numerators / denominators from CIs)
#   4. Run disbayes (Bayesian estimation of case fatality, incidence, remission)
#   5. Interpolate 5-year → 1-year age groups
#   6. Join everything and write output
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(disbayes)
library(tempdisagg)
library(doSNOW)
library(parallel)
library(foreach)

# =============================================================================
# 0. Configuration
# =============================================================================

GBD_FILE    <- "data/lecture_5_6/gbd_mslt_2023_colombia.csv"   # <-- set path as needed
OUTPUT_FILE <- "data/lecture_5_6/mslt_colombia.csv"

# Diseases present in this GBD extract and their short names.
DISEASE_MAP <- tibble(
  gbd_name = c(
    "Ischemic heart disease",
    "Chronic obstructive pulmonary disease",
    "Diabetes mellitus type 2",
    "All causes"
  ),
  sname = c("ishd", "copd", "dmt2", "allc")
)


# =============================================================================
# 1. Load & reshape GBD data
# =============================================================================

gbd_raw <- read_csv(GBD_FILE)

# Standardise column names: strip "_name" and "_id" suffixes
names(gbd_raw) <- gsub("_name$", "", names(gbd_raw))
gbd <- gbd_raw %>%
  select(-contains("_id"))

# Keep only the diseases we want and attach short names
gbd <- gbd %>%
  inner_join(DISEASE_MAP, by = c("cause" = "gbd_name")) %>%
  mutate(
    cause   = tolower(cause),
    sex     = tolower(sex),
    measure = case_when(
      measure == "YLDs (Years Lived with Disability)" ~ "ylds",
      TRUE ~ measure
    ),
    metric  = tolower(metric)
  )

# Pivot Rate / Number into columns, convert rate to per-1 (GBD gives per 100,000)
gbd_wide <- gbd %>%
  select(-upper, -lower) %>%
  pivot_wider(names_from = "metric", values_from = "val") %>%
  mutate(
    rate = rate / 100000,
    pop  = number / rate     # implicit population
  )

# Derive population from all-cause deaths (most stable denominator)
gbd_pop <- gbd_wide %>%
  filter(sname == "allc", measure == "Deaths") %>%
  select(age, sex, pop)

# Full rate table joined to consistent population 
gbd_rate <- gbd_wide %>%
  select(-pop) %>%
  left_join(gbd_pop, by = c("age", "sex")) %>%
  # Drop younger ages (MSLT models adults)
  filter(!age %in% c("<5 years", "5-9 years", "10-14 years")) %>%
  mutate(
    age = case_when(
      age == "95+ years" ~ "95-99",
      TRUE ~ str_remove(age, " years")
    )
  ) %>%
  rowwise() %>%
  separate(age, into = c("from_age", "to_age"), sep = "-", remove = FALSE) %>%
  mutate(
    from_age = as.integer(from_age),
    to_age   = as.integer(to_age),
    age_cat  = from_age + 2L        # midpoint label used for joining
  ) %>%
  ungroup()

# =============================================================================
# 2. Compute disability weights
# =============================================================================

# Pivot to wide: one row per age × sex × disease, columns = measure × (rate/number)
gbd_wider <- gbd_rate %>%
  mutate(age_sex = paste(age_cat, sex, sep = "_")) %>%
  pivot_wider(
    id_cols     = c(age, age_cat, age_sex, sex, pop, from_age, to_age),
    names_from  = c(measure, sname),
    values_from = c(rate, number),
    names_glue  = "{measure}_{.value}_{sname}"
  ) %>%
  mutate(across(everything(), ~ replace(., is.infinite(.) | is.na(.), 0)))

# Residual all-cause YLD rate not explained by modelled diseases
all_disease_yld_count <- gbd_wider %>%
  select(matches("^ylds_number_(?!allc)", perl = TRUE)) %>%
  rowSums(na.rm = TRUE)

gbd_wider <- gbd_wider %>%
  mutate(
    pyld_rate = (ylds_number_allc - all_disease_yld_count) / pop,
    pyld_rate = pmax(pyld_rate, 0)   # clip at zero (rounding artefacts)
  )

# =============================================================================
# 3. Prepare disbayes inputs
# =============================================================================
# Disbayes needs counts (numerator/denominator) estimated from GBD's CIs.
# This mirrors the ci2num approach in the original disbayes.R.

gbd_raw_db <- read_csv(GBD_FILE)
names(gbd_raw_db) <- gsub("_name$", "", names(gbd_raw_db))

gbd_db_prep <- gbd_raw_db %>%
  select(-contains("_id")) %>%
  inner_join(DISEASE_MAP, by = c("cause" = "gbd_name")) %>%
  filter(sname != "allc") %>%
  mutate(
    cause   = tolower(cause),
    sex     = tolower(sex),
    metric  = tolower(metric)
  ) %>%
  filter(metric == "rate") %>%
  mutate(
    val   = val   / 100000,
    lower = lower / 100000,
    upper = upper / 100000,
    # Guard against boundary values that break the beta CI fitting
    lower = if_else(lower <= 0, pmin(val / 2, 1e-5), lower),
    upper = if_else(upper >= 1, pmax((1 + val) / 2, 0.99999), upper),
    age   = case_when(
      age == "95+ years" ~ "95-99",
      age == "<5 years"  ~ "0-4",
      TRUE ~ str_remove(age, " years")
    ),
    measure = case_when(
      measure == "YLDs (Years Lived with Disability)" ~ "YLDs",
      TRUE ~ measure
    )
  ) %>%
  filter(measure %in% c("Deaths", "Incidence", "Prevalence")) %>%
  rowwise() %>%
  separate(age, into = c("from_age", "to_age"), sep = "-", remove = FALSE) %>%
  mutate(
    from_age = as.integer(from_age),
    to_age   = as.integer(to_age),
    agediff  = to_age - from_age + 1L
  ) %>%
  ungroup()

# Population size per age-sex-cause (for capping denominators)
gbd_num_db <- gbd_raw_db %>%
  select(-contains("_id")) %>%
  inner_join(DISEASE_MAP, by = c("cause" = "gbd_name")) %>%
  filter(sname != "allc") %>%
  mutate(
    cause   = tolower(cause),
    sex     = tolower(sex),
    metric  = tolower(metric),
    age     = case_when(
      age == "95+ years" ~ "95-99",
      TRUE ~ str_remove(age, " years")
    )
  ) %>%
  filter(metric == "number",
         measure %in% c("Deaths", "Incidence", "Prevalence")) %>%
  select(measure, sex, age, cause, sname, Number = val)

gbd_db_prep <- gbd_db_prep %>%
  left_join(
    gbd_num_db, by = c("measure", "sex", "age", "cause", "sname")
  ) %>%
  mutate(pop_actual = Number / val)

# Compute effective sample sizes using disbayes::ci2num (parallel)
message("Computing effective sample sizes from GBD confidence intervals...")
library(doParallel)
registerDoParallel(max(1L, detectCores() - 1L))

n_rows   <- nrow(gbd_db_prep)
numdenom <- foreach(i = seq_len(n_rows), .combine = rbind,
                    .packages = "disbayes") %dopar% {
                      row <- gbd_db_prep[i, ]
                      if (!is.na(row$val) && row$val > row$lower && row$val < row$upper) {
                        counts <- tryCatch(
                          disbayes:::ci2num(row$val, row$lower, row$upper, denom0 = row$pop_actual),
                          error = function(e) list(num = NA_real_, denom = NA_real_)
                        )
                        c(rowid = i, num = counts$num, denom = counts$denom)
                      } else {
                        c(rowid = i, num = NA_real_, denom = NA_real_)
                      }
                    }
stopImplicitCluster()

numdenom    <- as.data.frame(numdenom)
gbd_db_prep <- gbd_db_prep %>%
  mutate(rowid = row_number()) %>%
  left_join(numdenom, by = "rowid") %>%
  mutate(
    pop_actual = replace_na(pop_actual, 5000),
    denom      = if_else(is.na(denom) | denom > pop_actual, pop_actual, denom),
    num        = if_else(is.na(num), round(val * denom), num)
  ) %>%
  select(-rowid)

# Disaggregate 5-year groups to 1-year using tempdisagg
message("Disaggregating 5-year age groups to 1-year...")

gbd_db_prep <- gbd_db_prep %>%
  mutate(
    num1yr   = round(num   / agediff),
    denom1yr = round(denom / agediff)
  )

gbd_grp <- gbd_db_prep %>%
  group_by(measure, sex, sname) %>%
  arrange(measure, sex, sname, from_age)

disagg_fn <- function(dat, key) {
  res <- with(dat, {
    data.frame(
      from_age = rep(from_age, agediff) + sequence(agediff) - 1,
      num      = predict(td(num   ~ 1, to = 5, method = "fast")),
      denom    = predict(td(denom ~ 1, to = 5, method = "fast"))
    )
  })
  res
}

gbd_disagg_smooth <- group_modify(gbd_grp, disagg_fn) %>% ungroup()

# Fall back to simple /5 for any groups that produced negatives
neg_num   <- gbd_disagg_smooth %>% filter(num   < 0) %>%
  select(measure, sex, sname, agegroup = from_age) %>% distinct() %>%
  mutate(fix_num = TRUE)
neg_denom <- gbd_disagg_smooth %>% filter(denom < 0) %>%
  select(measure, sex, sname, agegroup = from_age) %>% distinct() %>%
  mutate(fix_denom = TRUE)

gbd_disagg <- gbd_disagg_smooth %>%
  rename(ageyr = from_age) %>%
  left_join(
    gbd_db_prep %>% select(measure, sex, sname, from_age, num1yr, denom1yr),
    by = c("measure", "sex", "sname", "ageyr" = "from_age")
  ) %>%
  mutate(
    num   = if_else(num   < 0, num1yr,   round(num)),
    denom = if_else(denom < 0, denom1yr, round(denom))
  ) %>%
  select(measure, sex, sname, age = ageyr, num, denom)

# -----------------------------------------------------------------------------
# FIX 1: Enforce num/denom consistency before handing data to disbayes.
#
# Stroke (and some other diseases) has zero or near-zero denominators at young
# ages, but nonzero numerators can appear after disaggregation — an impossible
# probability that causes Stan to crash. Three rules enforced:
#   (a) floor both at 0  (tempdisagg smoothing can leave tiny negatives)
#   (b) num cannot exceed denom  (rate > 1 is impossible)
#   (c) if denom is 0, num must also be 0
# -----------------------------------------------------------------------------
gbd_disagg <- gbd_disagg %>%
  mutate(
    denom = pmax(denom, 0L),
    num   = pmax(num,   0L),
    num   = pmin(num, denom),
    num   = if_else(denom == 0L, 0L, as.integer(num))
  )

# Format for disbayes: wide on measure
gbddb <- gbd_disagg %>%
  filter(measure %in% c("Deaths", "Incidence", "Prevalence")) %>%
  mutate(measure = recode(measure,
                          Deaths     = "mort",
                          Incidence  = "inc",
                          Prevalence = "prev")) %>%
  pivot_wider(
    names_from  = measure,
    values_from = c(num, denom),
    values_fill = 0
  ) %>%
  rename(
    inc_num    = num_inc,   inc_denom  = denom_inc,
    prev_num   = num_prev,  prev_denom = denom_prev,
    mort_num   = num_mort,  mort_denom = denom_mort
  ) %>%
  arrange(sname, sex, age) %>%
  mutate_all(~ replace_na(., 0))

# -----------------------------------------------------------------------------
# Disease-specific remission assumptions (recovery duration in years)
# Higher = less remission. 999 = effectively no remission (but causes disbayes issues)
# Use 10000 for very low remission that still works with disbayes
# -----------------------------------------------------------------------------
RECOVERY_YEARS <- c(
  ishd = 10000,   # IHD: very low remission (chronic heart disease) - use 10000 not 999
  copd = 20,      # COPD: some recovery possible with treatment
  dmt2 = 100      # DMT2: higher = lower remission = higher prevalence
)

message("Disease-specific recovery assumptions:")
for (d in names(RECOVERY_YEARS)) {
  if (RECOVERY_YEARS[d] >= 999) {
    message("  ", d, ": no remission")
  } else {
    message("  ", d, ": ", RECOVERY_YEARS[d], " years (rate = ", round(1/RECOVERY_YEARS[d], 4), ")")
  }
}

saveRDS(gbddb, "data/lecture_5_6/disbayes_input_colombia.rds")
message("Saved disbayes input to data/lecture_5_6/disbayes_input_colombia.rds")

# Quick diagnostic — flag any remaining impossible rows
impossible <- gbddb %>%
  filter(inc_num > inc_denom | prev_num > prev_denom | mort_num > mort_denom)
if (nrow(impossible) > 0) {
  message("WARNING: ", nrow(impossible),
          " rows still have num > denom after cleaning — inspect before running disbayes:")
  print(impossible)
} else {
  message("All num <= denom checks passed.")
}

# =============================================================================
# 4. Run disbayes
# =============================================================================

## incidence, prevalence and disease mortality from results tool. 
## remission based on documentation of GBD indicating after how many years 
## recovery is assumed. 

diseases_db  <- unique(gbddb$sname) %>% sort()
genders_db   <- unique(gbddb$sex)   %>% sort()

combinations <- crossing(disease = diseases_db, gender = genders_db) %>%
  mutate(combo = paste(disease, gender, sep = "_")) %>%
  pull(combo)

message(paste("Running disbayes for", length(combinations), "disease × sex combinations..."))

gbddb <- gbddb %>%
  rowwise() %>%
  mutate(
    recovery_yrs = RECOVERY_YEARS[sname],
    rem_num = if (is.na(recovery_yrs) || recovery_yrs >= 50000) {
      0L
    } else {
      floor(prev_num / recovery_yrs)
    },
    rem_denom = prev_num
  ) %>%
  ungroup() %>%
  mutate(
    rem_num = if_else(rem_denom > 0 & rem_num > rem_denom, rem_denom, rem_num),
    rem_num = if_else(is.na(rem_num), 0L, rem_num)
  ) %>%
  select(-recovery_yrs)

message("Synthetic remission data added:")
for (d in names(RECOVERY_YEARS)) {
  rec_yrs <- RECOVERY_YEARS[d]
  if (rec_yrs >= 999) {
    message("  ", d, ": no remission")
  } else {
    message("  ", d, ": ", rec_yrs, " years (rate = ", round(1/rec_yrs, 4), ")")
  }
}

# -----------------------------------------------------------------------------
# Add GBD-style value priors for DMT2 only
# - Remission: max 0.01 for ages 15+
# - Case fatality (excess mortality): max 0.15
# - Incidence: max 0.0008 for ages 1-15, max 0.1 for ages 15+
# -----------------------------------------------------------------------------
gbddb <- gbddb %>%
  mutate(
    # Constrain incidence rates for dmt2 only
    inc_rate = ifelse(inc_denom > 0, inc_num / inc_denom, 0),
    inc_rate = ifelse(sname == "dmt2",
      case_when(
        age < 1 ~ 0,
        age < 15 ~ pmin(inc_rate, 0.0008),
        TRUE ~ pmin(inc_rate, 0.2)  # increased from 0.1 to allow higher incidence
      ),
      inc_rate
    ),
    inc_num = ifelse(inc_denom > 0, floor(inc_rate * inc_denom), 0L),
    
    # Constrain case fatality rates for dmt2 only (max 0.15)
    cf_rate = ifelse(inc_denom > 0, mort_num / inc_denom, 0),
    cf_rate = ifelse(sname == "dmt2", pmin(cf_rate, 0.15), cf_rate),
    mort_num = ifelse(inc_denom > 0, floor(cf_rate * inc_denom), 0L),
    
    # Constrain remission rates for dmt2 only (max 0.01 for ages 15+)
    rem_rate = ifelse(rem_denom > 0, rem_num / rem_denom, 0),
    rem_rate = ifelse(sname == "dmt2" & age < 15, 0, rem_rate),
    rem_rate = ifelse(sname == "dmt2" & age >= 15, pmin(rem_rate, 0.01), rem_rate),
    rem_num = ifelse(rem_denom > 0, floor(rem_rate * rem_denom), 0L)
  ) %>%
  select(-inc_rate, -cf_rate, -rem_rate)


message("GBD-style value priors applied to DMT2 only:")
message("  Incidence: max 0.0008 (ages 1-15), max 0.1 (ages 15+)")
message("  Case fatality: max 0.15")
message("  Remission: max 0.01 (ages 15+)")

disbayes_results <- vector("list", length(combinations))

for (i in seq_along(combinations)) {
  
 
  sel_disease <- strsplit(combinations[i], "_")[[1]][1]
  sel_gender  <- strsplit(combinations[i], "_")[[1]][2]
  
  dat <- gbddb %>% filter(sname == sel_disease, sex == sel_gender)
  
  # Run disbayes - for IHD skip remission (causes numerical issues)
  if (sel_disease == "ishd") {
    dbres <- disbayes(
      data       = dat,        age        = "age",
      inc_num    = "inc_num",  inc_denom  = "inc_denom",
      prev_num   = "prev_num", prev_denom = "prev_denom",
      mort_num   = "mort_num", mort_denom = "mort_denom",
      method     = "opt", iter = 10000, eqage = 0
    )
  } else {
    dbres <- disbayes(
      data       = dat,        age        = "age",
      inc_num    = "inc_num",  inc_denom  = "inc_denom",
      prev_num   = "prev_num", prev_denom = "prev_denom",
      mort_num   = "mort_num", mort_denom = "mort_denom",
      rem_num    = "rem_num",  rem_denom  = "rem_denom",
      method     = "opt", iter = 10000, eqage = 0,
      rem_prior  = c(1.1, 1)
    )
  }
  
  
  summ <- disbayes::tidy(dbres) 
  
  
  message("disbayes tidy columns: ", paste(names(summ), collapse = ", "))
  message("disbayes tidy vars: ", paste(unique(summ$var), collapse = ", "))
  
  # Find the median column - could be "50%" or "mode" depending on method
  median_col <- if ("50%" %in% names(summ)) "50%" else "mode"
  
  cf_col  <- paste0("case_fatality_", sel_disease)
  inc_col <- paste0("incidence_",     sel_disease)
  rem_col <- paste0("remission_",     sel_disease)
  
  cf_df <- summ %>%
    filter(var == "cf", between(age, 0, 100)) %>%
    transmute(
      sex_age_cat = paste0(sel_gender, "_", age),
      !!cf_col    := .data[[median_col]]
    )

  inc_df <- summ %>%
    filter(var == "inc", between(age, 0, 100)) %>%
    transmute(
      sex_age_cat = paste0(sel_gender, "_", age),
      !!inc_col   := .data[[median_col]]
    )

  rem_df <- summ %>%
    filter(var == "rem", between(age, 0, 100)) %>%
    transmute(
      sex_age_cat = paste0(sel_gender, "_", age),
      !!rem_col   := .data[[median_col]]
    )

  if (nrow(rem_df) == 0) {
    message("WARNING: No remission from disbayes for ", sel_disease, "/", sel_gender, " - using 0")
    rem_df <- data.frame(
      sex_age_cat = cf_df$sex_age_cat,
      !!rem_col   := 0
    )
  }

  out_df <- cf_df %>%
    inner_join(inc_df, by = "sex_age_cat") %>%
    inner_join(rem_df, by = "sex_age_cat")

  message(sel_disease, "/", sel_gender, " - remission range: ",
          round(range(out_df[[rem_col]], na.rm = TRUE), 6))
  
  disbayes_results[[i]] <- list(df = out_df, dbres = dbres)
}

names(disbayes_results) <- combinations


# Assemble results ----

all_long <- disbayes_results %>%
  purrr::keep(~ is.list(.) && "df" %in% names(.)) %>%
  purrr::map("df") %>%
  purrr::map(function(df) {
    # Extract disease name from the case_fatality column
    dis <- sub("case_fatality_", "",
               grep("^case_fatality_", names(df), value = TRUE))
    # Standardise column names to generic names before stacking
    df %>%
      dplyr::rename(
        case_fatality = !!paste0("case_fatality_", dis),
        incidence     = !!paste0("incidence_",     dis),
        remission     = !!paste0("remission_",     dis)
      ) %>%
      dplyr::mutate(disease = dis)
  }) %>%
  dplyr::bind_rows()
# all_long now has: sex_age_cat | case_fatality | incidence | remission | disease
# Each row is one age × sex × disease — no column conflicts possible

disbayes_output <- all_long %>%
  tidyr::pivot_wider(
    id_cols     = sex_age_cat,
    names_from  = disease,
    values_from = c(case_fatality, incidence, remission),
    names_glue  = "{.value}_{disease}"   # → case_fatality_copd, incidence_copd, ...
  )


write_csv(disbayes_output, "data/lecture_5_6/disbayes_output_colombia.csv")
message("Saved disbayes output to data/lecture_5_6/disbayes_output_colombia.csv")


# =============================================================================
# 5. Interpolate 5-year → 1-year and assemble final dataset
# =============================================================================

# Scaffold: every age 0–100 × sex
mslt_scaffold <- data.frame(
  age = rep(0:100, 2),
  sex = c(rep("male", 101), rep("female", 101))
) %>%
  mutate(age_sex = paste(age, sex, sep = "_"))

# 5-year midpoints map to age_cat labels (e.g., 0-4 → age_cat = 2)
age_cat_labels <- tibble(
  from_age = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
               55, 60, 65, 70, 75, 80, 85, 90, 95),
  age_cat  = from_age + 2L
)

# Spline interpolation on log scale (avoids negative interpolated rates)
interpolate_log <- function(values_at_cats) {
  age_full <- 0:100
  keep     <- !is.na(values_at_cats) & values_at_cats > 0
  if (sum(keep) < 2) return(rep(0, 101))
  f <- stats::splinefun(
    x      = age_full[keep],
    y      = log(values_at_cats[keep]),
    method = "monoH.FC",
    ties   = mean
  )
  exp(f(age_full))
}

# Join gbd_wider to scaffold via age_cat
gbd_for_interp <- gbd_wider %>%
  select(age_cat, sex, pop, pyld_rate,
         matches("^(deaths|ylds)_(rate|number)_"),
         matches("^prevalence_")) %>%
  mutate(sex = tolower(sex))

# Build per-age 0-100 data frame with interpolated columns
interp_diseases <- setdiff(DISEASE_MAP$sname, "allc")

result_list <- map(c("male", "female"), function(sx) {

  sub <- gbd_for_interp %>% filter(sex == sx) %>% arrange(age_cat)

  # Helper: get a column by age_cat index for this sex
  get_vec <- function(col) {
    v <- rep(NA_real_, 101)
    idx <- sub$age_cat + 1L  # age_cat 2 → index 3, etc.
    idx <- pmin(pmax(idx, 1L), 101L)
    v[idx] <- as.numeric(sub[[col]])
    v
  }

  out <- mslt_scaffold %>% filter(sex == sx) %>% select(age, sex, age_sex)

  # Population: use value at midpoint; NA elsewhere
  pop_vec  <- get_vec("pop")
  out$population <- pop_vec

  # All-cause mortality rate → mx
  mx_col <- "Deaths_rate_allc"
  if (mx_col %in% names(sub)) {
    out$mx <- interpolate_log(get_vec(mx_col))
  } else {
    out$mx <- NA_real_
  }

  # pyld_rate
  out$pyld_rate <- interpolate_log(get_vec("pyld_rate"))

  # Per-disease columns
  for (dis in interp_diseases) {

    dr_col  <- paste0("deaths_rate_", dis)
    yr_col  <- paste0("ylds_rate_",   dis)
    yn_col  <- paste0("ylds_number_", dis)
    pn_col  <- paste0("prevalence_number_", dis)  # note: name from pivot

    # deaths_rate
    out[[paste0("deaths_rate_", dis)]] <-
      if (dr_col %in% names(sub)) interpolate_log(get_vec(dr_col)) else 0

    # ylds_rate
    out[[paste0("ylds_rate_", dis)]] <-
      if (yr_col %in% names(sub)) interpolate_log(get_vec(yr_col)) else 0

    # dw_adj (disability weight adjusted)
    if (yn_col %in% names(sub) && pn_col %in% names(sub)) {
      pyld_vec <- get_vec("pyld_rate")
      raw_dw   <- get_vec(yn_col) / pmax(get_vec(pn_col), 1e-10) /
                  pmax(1 - pyld_vec, 1e-10)
      raw_dw[is.nan(raw_dw) | is.infinite(raw_dw)] <- 0
      out[[paste0("dw_adj_", dis)]] <- interpolate_log(raw_dw)
    } else {
      out[[paste0("dw_adj_", dis)]] <- 0
    }

    # prevalence_rate (needed for some downstream MSLT models)
    prev_r_col <- paste0("prevalence_rate_", dis)
    out[[paste0("prevalence_", dis)]] <-
      if (prev_r_col %in% names(sub)) interpolate_log(get_vec(prev_r_col)) else 0
  }

  out
})

mslt_base <- bind_rows(result_list)

# Join disbayes outputs (case fatality, incidence, remission)
disbayes_tidy <- disbayes_output %>%
  separate(sex_age_cat, into = c("sex", "age"), sep = "_") %>%
  mutate(
    age     = as.integer(age),
    age_sex = paste(age, sex, sep = "_")
  )

mslt_final <- mslt_base %>%
  left_join(disbayes_tidy %>% select(-sex, -age), by = "age_sex") %>%
  mutate(across(everything(),
                ~ replace(., is.infinite(.) | is.na(.), 0))) %>%
  arrange(sex, age)

# =============================================================================
# 6. Write output
# =============================================================================

write_csv(mslt_final, OUTPUT_FILE)
message("Done. Output written to: ", OUTPUT_FILE)
message("Rows: ", nrow(mslt_final), "  Columns: ", ncol(mslt_final))
message("Columns: ", paste(names(mslt_final), collapse = ", "))
