#!/usr/bin/env Rscript
# =============================================================================
# 02_growth_nhs_table1.R
#
# Table 1 - NHS growth measures (height, weight, BMI and z-scores), per site,
#           restricted to ages 2 to <20 years.
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

#Participant list and site membership come from the EHR project
ehr <- redcap_export(redcap_token_ehr, list("forms[0]" = "demographics_5a25")) %>%
  rename(local_id = mrn)

#Get vitals data from REDcap with data pulled from NHS
nhs_raw <- redcap_export(redcap_token_nhs, redcap_fields(nhs_growth_export_fields),
                         export_dags = 'false')
nhs <- nhs_growth(nhs_raw)

#Calculate z-scores and percentiles for height, weight, and BMI
nhs_z <- cdcanthro(nhs, age = age, wt = weight, ht = height, bmi = BMI, all = TRUE)

print_site <- function(dag) {
  site_ids <- site_participants(ehr, dag)
  df <- nhs_for_site(nhs_z, site_ids)

  #Table 1 growth rows cover ages 2 to <20 years
  df <- filter(df, age >= 24 & age < 240)

  cat("\n\n==================== ", dag, " ====================\n")
  cat("\n--- Table 1: NHS growth measures (age 2 to <20 years) ---\n")
  cat("Participant-visits with height and weight recorded: ", nrow(df),
      " (N=", n_distinct(df$local_id), ")\n", sep = "")

  #Row order follows Table 1: BMI, BMI z, Weight, Weight z, Height, Height z
  print(data.frame(
    measure = c("BMI, kg/m2", "BMI age-adjusted z-score",
                "Weight, kg", "Weight for age, z-score",
                "Height, m", "Height for age, z-score"),
    median_iqr = c(desc_continuous(df$BMI),
                   desc_continuous(df$bmiz),
                   desc_continuous(df$weight),
                   desc_continuous(df$waz),
                   #height is recorded in cm and reported in meters
                   desc_continuous(df$height, scale = 1/100),
                   desc_continuous(df$haz)),
    n = c(sum(!is.na(df$BMI)), sum(!is.na(df$bmiz)),
          sum(!is.na(df$weight)), sum(!is.na(df$waz)),
          sum(!is.na(df$height)), sum(!is.na(df$haz)))
  ), row.names = FALSE)

  invisible(df)
}

invisible(lapply(ehr_data_access_groups(ehr), print_site))
