#!/usr/bin/env Rscript
# =============================================================================
# 05_labs_ehr_table3.R
#
# Table 3 - EHR absolute eosinophil counts per site: the value, observations
#           per participant, and duration of follow-up.
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

#Get lab data from REDcap with data pulled from FHIR
labs <- redcap_export(redcap_token_ehr,
                      redcap_fields(c('record_id', 'mrn', 'labs_label',
                                      'labs_loinc_code', 'labs_time', 'labs_value')))
labs <- fill_mrn(labs)
labs$local_id <- labs$mrn

#Site membership comes from the EHR project
ehr_sites <- labs %>%
  filter(!is.na(local_id)) %>%
  distinct(local_id, redcap_data_access_group)

analyse_site <- function(dag) {
  cat("\n\n==================== ", dag, " ====================\n")
  site_ids <- site_participants(ehr_sites, dag)
  site_labs <- semi_join(labs, site_ids, by = 'local_id')
  if (!nrow(site_labs)) {
    cat("No EHR labs for this site.\n")
    return(invisible(NULL))
  }

  daily <- ehr_aec_daily(ehr_aec(site_labs))

  #Per-participant summaries
  median_abs <- daily %>% group_by(local_id) %>%
    summarise(median_abs = median(aec_ehr, na.rm = TRUE), .groups = 'drop')
  count_obs <- daily %>% group_by(local_id) %>% tally()
  follow_up <- daily %>% group_by(local_id) %>%
    summarise(follow_up = n_distinct(year), .groups = 'drop')

  cat("\n--- Table 3: EHR labs ---\n")
  cat("Total EHR records (AEC): ", nrow(daily),
      " (N=", n_distinct(daily$local_id), ")\n", sep = "")
  print(data.frame(
    measure = c("Absolute eosinophil counts, x10^9/L",
                "No. of observations per participant",
                "Duration of follow-up, years"),
    median_iqr = c(desc_continuous(median_abs$median_abs),
                   desc_continuous(count_obs$n),
                   desc_continuous(follow_up$follow_up))
  ), row.names = FALSE)

  invisible(daily)
}

invisible(lapply(ehr_data_access_groups(ehr_sites), analyse_site))
