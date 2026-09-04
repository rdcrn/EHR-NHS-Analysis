#!/usr/bin/env Rscript
# =============================================================================
# 04_growth_ehr_table3.R
#
# Table 3 - EHR growth records per site, split at 20 years. For each measure:
#           the value, observations per participant, and duration of follow-up.
# =============================================================================

#Load shared configuration and helpers. Works whether this script is run from
#the repository root or from inside the R/ folder.
local({
  p <- c("R/00_utils.R", "00_utils.R")
  hit <- p[file.exists(p)]
  if (!length(hit)) {
    stop("Cannot find 00_utils.R. Open EHR-NHS-Analysis.Rproj, or set the ",
         "working directory to the repository root or its R/ folder.")
  }
  source(hit[1], local = FALSE)
})

suppressPackageStartupMessages(library(cdcanthro))

#Get vitals data from REDcap with data pulled from FHIR
vitals <- redcap_export(redcap_token_ehr,
                        c(redcap_fields(c('record_id', 'mrn')),
                          list("forms[0]" = "vital_signs")))
vitals <- fill_mrn(vitals)
vitals$local_id <- vitals$mrn

#Site membership comes from the export the measurements come from, so it does
#not depend on the data access group being set on the demographics form too
ehr_sites <- vitals %>%
  filter(!is.na(local_id)) %>%
  distinct(local_id, redcap_data_access_group)

#Get demographics data from REDcap with data pulled from FHIR
demo_raw <- redcap_export(redcap_token_ehr, list("forms[0]" = "demographics_5a25")) %>%
  rename(local_id = mrn)
demo <- demo_raw %>% select(local_id, dob, sex)

#The three Table 3 rows for one measure
table3_rows <- function(df, col, label, scale = 1) {
  per_participant <- df %>%
    group_by(local_id) %>%
    summarise(n_obs = dplyr::n(),
              n_years = n_distinct(substr(vitals_time, 1, 4)),
              .groups = 'drop')
  data.frame(
    measure = c(label, "No. of observations per participant",
                "Duration of follow-up, years"),
    median_iqr = c(desc_continuous(df[[col]], scale = scale),
                   desc_continuous(per_participant$n_obs),
                   desc_continuous(per_participant$n_years))
  )
}

print_block <- function(df, label, include_z) {
  cat("\n--- Table 3: EHR growth,", label, "---\n")
  cat("Participant-days with height and weight measured: ", nrow(df),
      " (N=", n_distinct(df$local_id), ")\n", sep = "")
  if (!nrow(df)) return(invisible(NULL))

  #trim on the raw measure, then on its z-score where z-scores exist
  trim <- function(col, zcol) {
    out <- remove_outliers(df, col)
    if (include_z) out <- remove_outliers(out, zcol)
    out
  }

  #Row order follows Table 3: BMI, BMI z, Weight, Weight z, Height, Height z
  rows <- rbind(
    table3_rows(trim("mean_BMI", "mean_bmiz"), "mean_BMI", "BMI, kg/m2"),
    if (include_z) table3_rows(remove_outliers(df, "mean_bmiz"), "mean_bmiz", "BMI, z-score"),
    table3_rows(trim("mean_weight", "mean_weightz"), "mean_weight", "Weight, kg"),
    if (include_z) table3_rows(remove_outliers(df, "mean_weightz"), "mean_weightz", "Weight for age, z-score"),
    #height is recorded in cm and reported in meters
    table3_rows(trim("mean_height", "mean_heightz"), "mean_height", "Height, m", scale = 1/100),
    if (include_z) table3_rows(remove_outliers(df, "mean_heightz"), "mean_heightz", "Height for age, z-score")
  )
  print(rows, row.names = FALSE)
}

analyse_site <- function(dag) {
  cat("\n\n==================== ", dag, " ====================\n")
  site_ids <- site_participants(ehr_sites, dag)
  site_vitals <- semi_join(vitals, site_ids, by = 'local_id')
  if (!nrow(site_vitals)) {
    cat("No EHR vital signs for this site.\n")
    return(invisible(NULL))
  }

  #ehr_vitals_measures() keeps only days on which both height and weight were
  #measured, so every measure below rests on the same set of days
  measures <- ehr_vitals_measures(site_vitals)
  if (!nrow(measures)) {
    cat("No participant-days with both height and weight for this site.\n")
    return(invisible(NULL))
  }

  ages <- ehr_growth_ages(measures, demo)

  #under 20 uses the z-score frame; 20-and-over is summarised from the
  #measurements before cdcanthro, which does not compute z-scores at that age
  daily_z <- ehr_growth_daily(ehr_growth_zscores(ages))
  daily_raw <- ages %>%
    group_by(local_id, vitals_time) %>%
    summarize(mean_height = mean(height, na.rm = TRUE),
              mean_weight = mean(weight, na.rm = TRUE),
              mean_BMI    = mean(BMI, na.rm = TRUE),
              mean_age    = mean(age, na.rm = TRUE),
              .groups = 'drop')

  print_block(filter(daily_z, mean_age < 240), "age 2 to <20 years", include_z = TRUE)
  print_block(filter(daily_raw, mean_age >= 240), "age >=20 years", include_z = FALSE)

  invisible(daily_z)
}

invisible(lapply(ehr_data_access_groups(ehr_sites), analyse_site))
