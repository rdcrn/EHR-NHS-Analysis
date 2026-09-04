#!/usr/bin/env Rscript
# =============================================================================
# 09_bmiz_aec_correlation.R
#
# Supplementary temporal coherence check: Spearman correlation between
# participant-level mean BMI z-score and participant-level mean AEC, per site,
# plus a formal test of whether the two site estimates differ.
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

#Get lab data from REDcap with data pulled from FHIR
labs <- redcap_export(redcap_token_ehr,
                      redcap_fields(c('record_id', 'mrn', 'labs_label',
                                      'labs_loinc_code', 'labs_time', 'labs_value')))
labs <- fill_mrn(labs)
labs$local_id <- labs$mrn

#Participant-level means for one site
site_pairs <- function(dag) {
  site_ids <- site_participants(ehr_sites, dag)

  site_vitals <- semi_join(vitals, site_ids, by = 'local_id')
  site_labs <- semi_join(labs, site_ids, by = 'local_id')
  if (!nrow(site_vitals) || !nrow(site_labs)) return(NULL)

  #BMI z-score exists for ages 2 to <20 years only
  bmiz <- ehr_growth_zscores(ehr_growth_ages(ehr_vitals_measures(site_vitals), demo)) %>%
    filter(!is.na(bmiz)) %>%
    group_by(local_id) %>%
    summarise(mean_bmiz = mean(bmiz, na.rm = TRUE),
              n_bmiz = dplyr::n(), .groups = 'drop')

  aec <- ehr_aec_daily(ehr_aec(site_labs)) %>%
    group_by(local_id) %>%
    summarise(mean_aec = mean(aec_ehr, na.rm = TRUE),
              n_aec = dplyr::n(), .groups = 'drop')

  inner_join(bmiz, aec, by = 'local_id') %>% mutate(site = dag)
}

dags <- ehr_data_access_groups(ehr_sites)
paired <- bind_rows(lapply(dags, site_pairs))

#Spearman's rho with a two-sided p-value. exact = FALSE because ties are
#certain in this data and the exact method cannot be used with them.
spearman_row <- function(d) {
  fit <- suppressWarnings(cor.test(d$mean_bmiz, d$mean_aec, method = 'spearman',
                                  alternative = 'two.sided', exact = FALSE))
  data.frame(site = d$site[1],
             n = nrow(d),
             rho = formatC(round(unname(fit$estimate), 2), format = 'f', digits = 2),
             p_value = formatC(round(fit$p.value, 3), format = 'f', digits = 3))
}

cat("\n--- Spearman correlation: participant-level mean BMIz vs mean AEC ---\n")
by_site <- split(paired, paired$site)
print(bind_rows(lapply(by_site, spearman_row)), row.names = FALSE)

#Fisher z comparison of two Spearman estimates. 1.06 is the standard variance
#inflation factor for rank correlations.
if (length(by_site) == 2) {
  rho_n <- function(d) {
    fit <- suppressWarnings(cor.test(d$mean_bmiz, d$mean_aec,
                                     method = 'spearman', exact = FALSE))
    c(rho = unname(fit$estimate), n = nrow(d))
  }
  a <- rho_n(by_site[[1]]); b <- rho_n(by_site[[2]])
  se <- sqrt(1.06/(a['n'] - 3) + 1.06/(b['n'] - 3))
  stat <- (atanh(a['rho']) - atanh(b['rho'])) / se
  p <- 2 * pnorm(-abs(stat))
  cat("\nBetween-site comparison (Fisher z on Spearman rho):\n")
  cat("  ", names(by_site)[1], " vs ", names(by_site)[2],
      ": z = ", formatC(unname(stat), format = 'f', digits = 2),
      ", p = ", formatC(unname(p), format = 'f', digits = 3), "\n", sep = "")
} else {
  cat("\nBetween-site comparison is only reported when exactly two sites are present.\n")
}
