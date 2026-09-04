#!/usr/bin/env Rscript
# =============================================================================
# 03_labs_nhs_table1.R
#
# Table 1 - NHS absolute eosinophil counts (AEC), per site.
# Also leaves one AEC dataframe per site in the environment (nhs_aec_by_site).
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

#Participant list and site membership come from the EHR project
ehr <- redcap_export(redcap_token_ehr, list("forms[0]" = "demographics_5a25")) %>%
  rename(local_id = mrn)

#Get lab data from REDCap with data from NHS
nhs_raw <- redcap_export(redcap_token_nhs, redcap_fields(nhs_aec_export_fields))
aec_all <- nhs_aec(nhs_raw)

build_site <- function(dag) {
  site_ids <- site_participants(ehr, dag)
  df <- nhs_for_site(aec_all, site_ids)

  cat("\n\n==================== ", dag, " ====================\n")
  cat("\n--- Table 1: NHS lab values ---\n")
  cat("Total NHS records (AEC): ", nrow(df),
      " (N=", n_distinct(df$local_id), ")\n", sep = "")
  #The AEC row is the median across all records, not across per-participant medians
  cat("Absolute eosinophil counts, x10^9/L: ", desc_continuous(df$value), "\n", sep = "")

  df
}

dags <- ehr_data_access_groups(ehr)
nhs_aec_by_site <- lapply(dags, build_site)
names(nhs_aec_by_site) <- dags
