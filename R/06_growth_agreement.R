#!/usr/bin/env Rscript
# =============================================================================
# 06_growth_agreement.R
#
# Table 2b - growth agreement between EHR and NHS per site, as bias with 95%
#            limits of agreement (Bland-Altman).
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

suppressPackageStartupMessages({
  library(cdcanthro)
  library(SimplyAgree)
})

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

#Get vitals data from REDcap with data pulled from NHS
nhs_raw <- redcap_export(redcap_token_nhs, redcap_fields(nhs_growth_export_fields),
                         export_dags = 'false')
nhs <- nhs_growth(nhs_raw)
nhs_z <- cdcanthro(nhs, age = age, wt = weight, ht = height, bmi = BMI, all = TRUE)

#Average NHS measurements taken on the same date for the same participant
nhs_daily <- nhs_z %>%
  group_by(local_id, vitals_time) %>%
  summarise(across(c(height, weight, BMI, haz, waz, bmiz), ~ mean(.x, na.rm = TRUE)),
            .groups = 'drop')

analyse_site <- function(dag) {
  cat("\n\n==================== ", dag, " ====================\n")
  site_ids <- site_participants(ehr_sites, dag)
  site_vitals <- semi_join(vitals, site_ids, by = 'local_id')
  if (!nrow(site_vitals)) {
    cat("No EHR vital signs for this site.\n")
    return(invisible(NULL))
  }

  #EHR measurements, averaged within participant and date
  ehr_daily <- ehr_growth_daily(ehr_growth_zscores(ehr_growth_ages(
    ehr_vitals_measures(site_vitals), demo)))

  #Match EHR and NHS measurements on participant and date, taking the closest
  #EHR observation when the dates do not coincide exactly
  paired <- match_nearest_date(nhs_daily, ehr_daily,
                               nhs_date = 'vitals_time', ehr_date = 'vitals_time')
  if (!nrow(paired)) {
    cat("No EHR/NHS growth pairs within", match_window_days, "day(s) for this site.\n")
    return(invisible(NULL))
  }
  cat("Date matching: ", describe_match(paired), "\n", sep = "")

  #agree_test(nhs, ehr): the NHS value is the first argument throughout, so the
  #sign of the bias is consistent across scripts
  fit <- function(nhs_col, ehr_col, scale = 1) {
    d <- paired %>% filter(!is.na(.data[[nhs_col]]), !is.na(.data[[ehr_col]]))
    if (!nrow(d)) return(NULL)
    list(fit = agree_test(d[[nhs_col]] * scale, d[[ehr_col]] * scale,
                          agree.level = 0.95, conf.level = 0.95),
         n = nrow(d))
  }

  #Row order follows Table 2b: BMI, BMI z, Weight, Weight z, Height, Height z
  specs <- list(
    list("BMI, kg/m2",              "BMI",  "mean_BMI",     1),
    list("BMI, z-score",            "bmiz", "mean_bmiz",    1),
    list("Weight, kg",              "weight", "mean_weight", 1),
    list("Weight for age, z-score", "waz",  "mean_weightz", 1),
    #height is recorded in cm and reported in meters
    list("Height, m",               "height", "mean_height", 1/100),
    list("Height for age, z-score", "haz",  "mean_heightz", 1)
  )

  rows <- lapply(specs, function(s) {
    r <- fit(s[[2]], s[[3]], s[[4]])
    if (is.null(r)) return(NULL)
    data.frame(measure = s[[1]], bias_loa = loa_fmt(r$fit), n_pairs = r$n)
  })

  cat("\n--- Table 2b: growth agreement (EHR vs NHS), bias [95% limits of agreement] ---\n")
  print(bind_rows(rows), row.names = FALSE)

  invisible(paired)
}

invisible(lapply(ehr_data_access_groups(ehr_sites), analyse_site))
