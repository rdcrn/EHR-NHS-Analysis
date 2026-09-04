#!/usr/bin/env Rscript
# =============================================================================
# 01_demographics_agreement.R
#
# Table 1  - NHS cohort demographics, per site
# Table 2a - Cohen's kappa for sex, race, ethnicity, per site
# Table S1 - simple agreement for the same variables, per site
# Also runs a Bland-Altman analysis of age.
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
  library(collapse)
  library(psych)
  library(SimplyAgree)
  library(binom)
})

# --- Category labels and display order (follows Table 1) ---------------------

sex_labels <- c("1" = "Male", "2" = "Female")
sex_order  <- c("Male", "Female", "Missing")

ethnicity_labels <- c("1" = "Hispanic", "2" = "Not Hispanic or Latino")
ethnicity_order  <- c("Hispanic", "Not Hispanic or Latino", "Missing")

race_labels <- c("1" = "American Indian or Alaskan Native", "2" = "Asian",
                 "3" = "Native Hawaiian or Pacific Islander",
                 "4" = "Black or African American", "5" = "White",
                 "6" = "More than one race", "9999" = "Other race",
                 "3333" = "Missing/unknown")
race_order <- c("American Indian or Alaskan Native", "Asian",
                "Black or African American", "White", "More than one race",
                "Other race", "Missing/unknown")

# --- EHR (FHIR-derived) demographics ----------------------------------------

ehr_raw <- redcap_export(redcap_token_ehr, list("forms[0]" = "demographics_5a25"))

ehr <- ehr_raw %>%
  mutate(
    sex_suppl = case_when(
      is.na(sex) ~ NA_character_,
      sex == 'M' ~ '1',
      sex == 'F' ~ '2',
      TRUE ~ NA_character_
    ),
    ethnicity_suppl = case_when(
      is.na(ethnicity) ~ NA_character_,
      ethnicity == '2135-2' ~ '1',
      ethnicity == '2186-5' ~ '2',
      TRUE ~ NA_character_
    ),
    race_suppl = case_when(
      is.na(race) ~ NA_character_,
      race == '2106-3' ~ '5',
      race == '1002-5' ~ '1',
      race == '2028-9' ~ '2',
      race == '2076-8' ~ '3',
      race == '2054-5' ~ '4',
      race == '2131-1' ~ '9999',
      race == 'UNK' ~ '3333',
      TRUE ~ NA_character_
    ),
    dob_suppl = dob,
    age_suppl_years = age_in_months(Sys.Date(), dob) / 12
  ) %>%
  rename(local_id = mrn) %>%
  select(local_id, redcap_data_access_group, sex_suppl, ethnicity_suppl,
         race_suppl, dob_suppl, age_suppl_years)

# --- NHS demographics --------------------------------------------------------

nhs_raw <- redcap_export(
  redcap_token_nhs,
  redcap_fields(c('local_id', 'dateofbirth', 'gender', 'race', 'ethnicity'))
)
nhs_raw <- nhs_raw %>% distinct(local_id, .keep_all = TRUE)

race_cols <- c("race___1", "race___2", "race___3", "race___4", "race___5",
               "race___3333", "race___4444", "race___9999")
race_map <- c(race___1 = "1", race___2 = "2", race___3 = "4", race___4 = "3",
              race___5 = "5", race___3333 = "3333", race___4444 = "4444",
              race___9999 = "9999")

nhs <- nhs_raw %>%
  rowwise() %>%
  mutate(
    #more than one box ticked reports as "more than one race"
    race_7801 = {
      selected <- race_cols[which(c_across(all_of(race_cols)) == '1')]
      if (length(selected) == 0) NA_character_
      else if (length(selected) > 1) '6'
      else unname(race_map[selected])
    }
  ) %>%
  ungroup() %>%
  mutate(
    sex_7801 = gender,
    ethnicity_7801 = ethnicity,
    dob_7801 = dateofbirth,
    age_7801_years = age_in_months(Sys.Date(), dateofbirth) / 12
  ) %>%
  select(local_id, sex_7801, ethnicity_7801, race_7801, dob_7801, age_7801_years)

# --- Per-site analysis -------------------------------------------------------

analyse_site <- function(dag) {
  site_ids <- site_participants(ehr, dag)
  merged <- ehr %>%
    semi_join(site_ids, by = 'local_id') %>%
    left_join(nhs, by = 'local_id')

  cat("\n\n==================== ", dag, " ====================\n")

  # Table 1: NHS cohort demographics. Row order follows Table 1.
  cat("\n--- Table 1: NHS cohort descriptive statistics ---\n")
  cat("Sex, n (%)\n")
  print(desc_categorical(merged$sex_7801, labels = sex_labels,
                         category_order = sex_order), row.names = FALSE)
  cat("Age, years, median (IQR): ", desc_continuous(merged$age_7801_years), "\n")
  cat("Race, n (%)\n")
  print(desc_categorical(merged$race_7801, labels = race_labels,
                         merge_missing_into = "Missing/unknown",
                         category_order = race_order), row.names = FALSE)
  cat("Ethnicity, n (%)\n")
  print(desc_categorical(merged$ethnicity_7801, labels = ethnicity_labels,
                         category_order = ethnicity_order), row.names = FALSE)

  # Paired subsets, one per variable
  paired_sex  <- na_omit(merged, cols = c("sex_suppl", "sex_7801"))
  paired_eth  <- na_omit(merged, cols = c("ethnicity_suppl", "ethnicity_7801"))
  paired_race <- filter(merged, !is.na(race_7801), race_7801 != 'NA')

  cat("\n--- Concordance: sex ---\n")
  print(table(paired_sex$sex_7801, paired_sex$sex_suppl))
  cat("\n--- Concordance: ethnicity ---\n")
  print(table(paired_eth$ethnicity_7801, paired_eth$ethnicity_suppl))
  cat("\n--- Concordance: race ---\n")
  print(table(paired_race$race_7801, paired_race$race_suppl))

  # Table 2a: Cohen's kappa. Row order follows Table 2a.
  agree_sex  <- cohen.kappa(x = cbind(paired_sex$sex_7801, paired_sex$sex_suppl))
  agree_race <- cohen.kappa(x = cbind(paired_race$race_7801, paired_race$race_suppl))
  agree_eth  <- cohen.kappa(x = cbind(paired_eth$ethnicity_7801, paired_eth$ethnicity_suppl))

  cat("\n--- Table 2a: Cohen's kappa (95% CI), EHR vs NHS ---\n")
  print(data.frame(
    variable = c("Sex", "Race", "Ethnicity"),
    kappa_ci = c(kappa_fmt(agree_sex), kappa_fmt(agree_race), kappa_fmt(agree_eth))
  ), row.names = FALSE)

  # Table S1: simple agreement
  score <- function(df, a, b) ifelse(df[[a]] == df[[b]], 1, 0)
  cat("\n--- Table S1: simple agreement (95% CI), EHR vs NHS ---\n")
  print(data.frame(
    variable = c("Sex", "Race", "Ethnicity"),
    agreement_ci = c(simple_agreement(score(paired_sex, "sex_7801", "sex_suppl")),
                     simple_agreement(score(paired_race, "race_7801", "race_suppl")),
                     simple_agreement(score(paired_eth, "ethnicity_7801", "ethnicity_suppl"))),
    n = c(nrow(paired_sex), nrow(paired_race), nrow(paired_eth))
  ), row.names = FALSE)

  # Bland-Altman analysis of age (years)
  paired_age <- na_omit(merged, cols = c("age_suppl_years", "age_7801_years"))
  fit_age <- agree_test(paired_age$age_7801_years, paired_age$age_suppl_years,
                        agree.level = 0.95, conf.level = 0.95)
  cat("\n--- Age (years) agreement, bias [95% limits of agreement] ---\n")
  cat(loa_fmt(fit_age), "\n")
  cat("Age (years) summary:\n")
  print(round(summary(paired_age$age_7801_years), n_digits))

  invisible(list(dag = dag, merged = merged, age_fit = fit_age))
}

results <- lapply(ehr_data_access_groups(ehr), analyse_site)
names(results) <- ehr_data_access_groups(ehr)
