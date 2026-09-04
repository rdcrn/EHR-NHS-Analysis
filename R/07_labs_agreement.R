#!/usr/bin/env Rscript
# =============================================================================
# 07_labs_agreement.R
#
# Agreement of absolute eosinophil counts between EHR and NHS per site, as bias
# with 95% limits of agreement (Bland-Altman). Repeat values on the same date
# are averaged within each source before pairing.
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

suppressPackageStartupMessages(library(SimplyAgree))

#Get lab data from REDcap with data pulled from FHIR
labs <- redcap_export(redcap_token_ehr,
                      redcap_fields(c('record_id', 'mrn', 'labs_label',
                                      'labs_loinc_code', 'labs_time', 'labs_value')))
labs <- fill_mrn(labs)
labs$local_id <- labs$mrn

ehr_sites <- labs %>%
  filter(!is.na(local_id)) %>%
  distinct(local_id, redcap_data_access_group)

#Get lab data from REDCap with data from NHS
nhs_raw <- redcap_export(redcap_token_nhs, redcap_fields(nhs_aec_export_fields))
nhs_daily <- nhs_aec_daily(nhs_aec(nhs_raw))

analyse_site <- function(dag) {
  cat("\n\n==================== ", dag, " ====================\n")
  site_ids <- site_participants(ehr_sites, dag)
  site_labs <- semi_join(labs, site_ids, by = 'local_id')
  if (!nrow(site_labs)) {
    cat("No EHR labs for this site.\n")
    return(invisible(NULL))
  }

  ehr_daily <- ehr_aec_daily(ehr_aec(site_labs))

  #Match on participant and draw date, taking the closest EHR observation when
  #the dates do not coincide exactly
  paired <- match_nearest_date(nhs_for_site(nhs_daily, site_ids), ehr_daily)
  if (!nrow(paired)) {
    cat("No EHR/NHS AEC pairs within", match_window_days, "day(s) for this site.\n")
    return(invisible(NULL))
  }
  cat("Date matching: ", describe_match(paired), "\n", sep = "")

  fit <- agree_test(paired$aec_nhs, paired$aec_ehr,
                    agree.level = 0.95, conf.level = 0.95)

  cat("\n--- AEC agreement (EHR vs NHS), bias [95% limits of agreement] ---\n")
  print(data.frame(
    measure = "Absolute eosinophil counts, x10^9/L",
    bias_loa = loa_fmt(fit),
    n_pairs = nrow(paired),
    n_participants = n_distinct(paired$local_id)
  ), row.names = FALSE)

  invisible(paired)
}

invisible(lapply(ehr_data_access_groups(ehr_sites), analyse_site))
