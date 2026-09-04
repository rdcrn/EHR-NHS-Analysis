#!/usr/bin/env Rscript
# =============================================================================
# 08_problem_list_agreement.R
#
# Table 2a - simple agreement of problem lists between EHR and NHS, per site:
#            EoE alone, any of the three conditions, and all three.
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

suppressPackageStartupMessages(library(binom))

#Get problem list data from REDCap with data pulled from FHIR
prob_ehr <- redcap_export(redcap_token_ehr,
                          c(redcap_fields(c('record_id', 'mrn')),
                            list("forms[0]" = "problem_list")))
prob_ehr <- fill_mrn(prob_ehr)
prob_ehr$local_id <- prob_ehr$mrn

ehr_sites <- prob_ehr %>%
  filter(!is.na(local_id)) %>%
  distinct(local_id, redcap_data_access_group)

#Get problem list data from REDCap with data from NHS
nhs_raw <- redcap_export(redcap_token_nhs, redcap_fields(nhs_problem_export_fields))
data_nhs <- nhs_problem_list(nhs_raw)

analyse_site <- function(dag) {
  cat("\n\n==================== ", dag, " ====================\n")
  site_ids <- site_participants(ehr_sites, dag)
  site_prob <- semi_join(prob_ehr, site_ids, by = 'local_id')

  data_fhir <- ehr_problem_list(site_prob)

  #Every participant at this site, with whichever problem lists exist
  merged <- site_ids %>%
    left_join(data_fhir, by = 'local_id') %>%
    left_join(data_nhs, by = 'local_id') %>%
    filter(!if_all(c(FHIR_problem, NHS_problem), is.na))

  if (!nrow(merged)) {
    cat("No problem list data for this site.\n")
    return(invisible(NULL))
  }

  #All three conditions: the proportion of codes the two sources share
  overlap_score <- mapply(function(a, b) {
    s1 <- strsplit(ifelse(is.na(a), "", a), ",")[[1]]
    s2 <- strsplit(ifelse(is.na(b), "", b), ",")[[1]]
    denom <- max(length(s1), length(s2))
    if (denom == 0) return(NA_real_)
    length(intersect(s1, s2)) / denom
  }, merged$FHIR_problem, merged$NHS_problem)

  #Any of the three conditions: recorded in one source at all, or not
  any_ehr <- as.integer(!is.na(merged$FHIR_problem) & merged$FHIR_problem != "")
  any_nhs <- as.integer(!is.na(merged$NHS_problem) & merged$NHS_problem != "")
  any_score <- as.integer(any_ehr == any_nhs)

  #EoE alone, dropping participants with no EoE in either source
  eoe_ehr <- as.integer(grepl('K20.0', merged$FHIR_problem, fixed = TRUE))
  eoe_nhs <- as.integer(grepl('K20.0', merged$NHS_problem, fixed = TRUE))
  keep <- !(eoe_ehr == 0 & eoe_nhs == 0)
  eoe_score <- as.integer(eoe_ehr[keep] == eoe_nhs[keep])

  #Row order follows Table 2a
  cat("\n--- Table 2a: simple agreement of problem lists (EHR vs NHS) ---\n")
  print(data.frame(
    variable = c("Only EoE", "Any of EoE, EoGN, EoC", "All three (EoE, EoGN, EoC)"),
    agreement_ci = c(simple_agreement(eoe_score),
                     simple_agreement(any_score),
                     simple_agreement(na.omit(overlap_score))),
    n = c(length(eoe_score), length(any_score), sum(!is.na(overlap_score)))
  ), row.names = FALSE)

  invisible(merged)
}

invisible(lapply(ehr_data_access_groups(ehr_sites), analyse_site))
