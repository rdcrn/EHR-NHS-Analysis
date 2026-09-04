#!/usr/bin/env Rscript
# =============================================================================
# 00_utils.R - shared configuration, REDCap access, extraction and formatting
#
# Sourced by every analysis script. Nothing here runs an analysis; it only
# defines configuration and functions.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(rlist)
  library(stringi)
})

# --- Configuration -----------------------------------------------------------

redcap_url <- "https://rc.rarediseasesnetwork.org/api/"
redcap_token_ehr <- Sys.getenv("REDCAP_API_TOKEN_EHR")
redcap_token_nhs <- Sys.getenv("REDCAP_API_TOKEN_NHS")

#Where figures and exports are written.
#
#Resolved once here so that the repository root and the R/ folder both write to
#<repo>/output rather than to <repo>/R/output. That matters now that the two
#Figure 2 scripts share their participant selection through a file in this
#folder: if they resolved it differently they would not see each other's cache.
output_dir <- if (file.exists("00_utils.R") && !file.exists("R/00_utils.R")) {
  file.path("..", "output")   # working directory is the R/ folder
} else {
  "output"                    # working directory is the repository root
}

#Local IDs excluded during QC
excluded_local_ids <- c('11879663')

#The NHS lab fields record when a value was entered (visit_date) separately from
#when the blood was drawn (date), so values with old draw dates keep entering the
#export as new visits occur. Set this to a date to reproduce an earlier export
#snapshot; NA uses all currently available data.
visit_date_cutoff <- NA

#Methods: "When EHR and NHS entries were not recorded on the same date, the
#closest EHR observation was used." This is the window, in days, within which an
#EHR observation may be matched to an NHS entry. 0 means exact-date matching
#only; increase it if a site's EHR and NHS dates are systematically offset.
match_window_days <- 1

#Decimal places for all reported statistics
n_digits <- 1

#LOINC codes. Only height and weight are read; BMI is derived from them.
loinc_height <- '8302-2'
loinc_weight <- '29463-7'
loinc_aec    <- c('6690-2', '26450-7', '713-8')  # WBC, eosinophil %, and related

#ICD-10 codes for the three EGID diagnoses
icd_egid <- c('K20.0', 'K52.81', 'K52.82')

#Participants drawn at random for the trajectory figures
plot_n_participants <- 5
plot_seed <- 20240101

# --- REDCap access -----------------------------------------------------------

#Pull a flat, raw-coded export. Fails fast with a readable message if REDCap
#returns an error object instead of data.
redcap_export <- function(token, extra_fields = list(), export_dags = 'true',
                          url = redcap_url) {
  if (identical(token, "")) {
    stop("REDCap token is empty. Set REDCAP_API_TOKEN_EHR and REDCAP_API_TOKEN_NHS.")
  }
  formData <- c(
    list(token = token,
         content = 'record',
         action = 'export',
         format = 'csv',
         type = 'flat',
         csvDelimiter = '',
         rawOrLabel = 'raw',
         rawOrLabelHeaders = 'raw',
         exportCheckboxLabel = 'false',
         exportSurveyFields = 'false',
         exportDataAccessGroups = export_dags,
         returnFormat = 'json'),
    extra_fields
  )
  response <- httr::POST(url, body = formData, encode = "form")
  result <- httr::content(response)
  if (is.list(result) && !is.data.frame(result)) {
    stop("REDCap API error: ", paste(unlist(result), collapse = " "))
  }
  as.data.frame(result)
}

#Build the fields[n] parameter list REDCap expects from a plain vector
redcap_fields <- function(fields) {
  setNames(as.list(fields), sprintf("fields[%d]", seq_along(fields) - 1))
}

#mrn is only populated on the first row of each record, so carry it forward
fill_mrn <- function(result) {
  for (i in 2:nrow(result))
  {
    row <- result[i,]
    id<-row$record_id
    row_b<-result[i-1,]
    id_b=row_b$record_id
    if(id==id_b)
    {
      row$mrn<-row_b$mrn
      result[i,]$mrn<-row$mrn
    }
  }
  result
}

# --- Data access groups ------------------------------------------------------

#The list of sites is read from the EHR project rather than hardcoded, so adding
#a site to REDCap is enough for every script to pick it up.
ehr_data_access_groups <- function(ehr_df) {
  dags <- unique(ehr_df$redcap_data_access_group)
  dags <- dags[!is.na(dags) & dags != ""]
  sort(as.character(dags))
}

#Participants belonging to one site, from an EHR export carrying local_id and
#redcap_data_access_group
site_participants <- function(ehr_df, dag) {
  ehr_df %>%
    filter(redcap_data_access_group == dag,
           !local_id %in% excluded_local_ids,
           !is.na(local_id)) %>%
    distinct(local_id)
}

#The NHS project uses different data access group names from the EHR project, so
#NHS records are attributed to a site by participant membership instead.
nhs_for_site <- function(nhs_df, site_ids) {
  semi_join(nhs_df, site_ids, by = 'local_id')
}

# --- Age ---------------------------------------------------------------------

#Age in months between two dates.
#
#difftime is given explicit units. Subtracting POSIXlt values uses units="auto",
#which picks a single unit for the whole vector based on its smallest element, so
#one same-day record makes every age come back in hours instead of days.
age_in_months <- function(date, dob) {
  as.numeric(difftime(as.Date(substr(as.character(date), 1, 10), format = '%Y-%m-%d'),
                      as.Date(substr(as.character(dob),  1, 10), format = '%Y-%m-%d'),
                      units = "days")) / 30.4
}

# --- Extraction: EHR vital signs --------------------------------------------

#One row per participant and measurement date, with height, weight and BMI.
#
#Only height and weight are read from the EHR; BMI is derived from them, so the
#stored BMI observation is not exported.
#
#Height and weight are paired on the calendar date. Sites differ in whether the
#two share a timestamp, so pairing on the timestamp would separate measures
#taken minutes apart. Repeat measures on a date are averaged. Nothing is carried
#across dates or participants.
#
#With require_both = TRUE (the default) a growth record is a day on which both
#height and weight were measured, so every reported measure - height, weight,
#BMI and the z-scores - rests on the same set of days. Days holding only one of
#the two are dropped rather than being completed from another visit. Set it to
#FALSE to keep those days, with NA for the measure that is absent.
#
#The date is taken with substr because sites differ in format: some FHIR servers
#return "YYYY-MM-DD HH:MM:SS" and others ISO 8601 "YYYY-MM-DDTHH:MM:SS".
#
#Columns are selected by name, so the result does not depend on the column order
#of the export.
ehr_vitals_measures <- function(vitals_df, require_both = TRUE) {
  wide <- vitals_df %>%
    filter(vitals_loinc_code %in% c(loinc_height, loinc_weight)) %>%
    select(local_id, vitals_loinc_code, vitals_time, vitals_value) %>%
    mutate(vitals_time  = substr(as.character(vitals_time), 1, 10),
           vitals_value = suppressWarnings(as.numeric(vitals_value))) %>%
    filter(!is.na(local_id), !is.na(vitals_value)) %>%
    group_by(local_id, vitals_time, vitals_loinc_code) %>%
    summarise(vitals_value = mean(vitals_value), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = vitals_loinc_code, values_from = vitals_value)

  #a site may record only one of the two measures
  if (!loinc_height %in% names(wide)) wide[[loinc_height]] <- NA_real_
  if (!loinc_weight %in% names(wide)) wide[[loinc_weight]] <- NA_real_

  out <- wide %>%
    transmute(local_id,
              vitals_time,
              height = .data[[loinc_height]],
              weight = .data[[loinc_weight]],
              BMI = weight / ((height / 100)^2)) %>%
    arrange(local_id, vitals_time) %>%
    as.data.frame()

  if (require_both) out <- out[!is.na(out$height) & !is.na(out$weight), , drop = FALSE]
  out
}

#Attach date of birth and age in months
ehr_growth_ages <- function(measures, demo) {
  df10 <- left_join(measures, demo, by = 'local_id') %>%
    mutate(age = age_in_months(vitals_time, dob)) %>%
    as.data.frame()  #cdcanthro uses data.table internally
  df10$age <- as.numeric(as.character(df10$age))
  df10
}

#z-scores from cdcanthro, which covers ages 2 to <20 years only. Rows outside
#that range come back without z-scores, so any 20-and-over summary must be built
#from the ehr_growth_ages() frame instead.
ehr_growth_zscores <- function(df10) {
  cdcanthro::cdcanthro(df10, age = age, wt = weight, ht = height, bmi = BMI, all = FALSE)
}

#Average repeat measurements taken on the same date
ehr_growth_daily <- function(zscores) {
  zscores %>%
    group_by(local_id, vitals_time) %>%
    summarize(mean_height = mean(height, na.rm = TRUE),
              mean_weight = mean(weight, na.rm = TRUE),
              mean_BMI    = mean(BMI, na.rm = TRUE),
              mean_age    = mean(age, na.rm = TRUE),
              mean_heightz = mean(haz, na.rm = TRUE),
              mean_weightz = mean(waz, na.rm = TRUE),
              mean_bmiz    = mean(bmiz, na.rm = TRUE),
              .groups = 'drop')
}

# --- Extraction: EHR labs ----------------------------------------------------

#Absolute eosinophil count, derived from the white cell count and the
#eosinophil percentage recorded at the same timestamp. A count can only be
#derived where both components are present, so timestamps with a single
#component are dropped.
#
#Columns are referenced by name. Renaming them by position is unsafe here: the
#export order depends on how local_id was added to the frame.
ehr_aec <- function(labs_df) {
  labs_df %>%
    filter(labs_loinc_code %in% loinc_aec) %>%
    select(local_id, record_id, labs_loinc_code, labs_time, labs_value) %>%
    mutate(labs_value = suppressWarnings(as.numeric(labs_value))) %>%
    group_by(local_id, record_id, labs_time) %>%
    summarise(n_components = dplyr::n(),
              absolute_eosonophil_count = if (dplyr::n() == 2)
                prod(labs_value) / 100 else NA_real_,
              .groups = 'drop') %>%
    filter(!is.na(absolute_eosonophil_count)) %>%
    select(local_id, record_id, labs_time, absolute_eosonophil_count) %>%
    distinct()
}

#One AEC value per participant per day
ehr_aec_daily <- function(aec_df) {
  aec_df %>%
    mutate(date = as.Date(substr(labs_time, 1, 10), format = '%Y-%m-%d'),
           aec  = suppressWarnings(as.numeric(absolute_eosonophil_count))) %>%
    filter(!is.na(local_id), !is.na(date), !is.na(aec)) %>%
    group_by(local_id, date) %>%
    summarise(aec_ehr = mean(aec), n_obs = dplyr::n(), .groups = 'drop') %>%
    mutate(year = as.integer(format(date, '%Y'))) %>%
    arrange(local_id, date)
}

# --- Extraction: NHS labs ----------------------------------------------------

#AEC is captured in two places on the medical history form: the count drawn at
#the current visit, and the count carried forward from the previous visit. Each
#has its own month/day/year date fields, so both are pulled and combined.
nhs_aec_fields <- list(
  current = c(value = 'mhlababscnt2',
              m = 'mhlababscntmdt2', d = 'mhlababscntddt2', y = 'mhlababscntydt2'),
  historical = c(value = 'mhlabhistabscnt2',
                 m = 'mhlabhistabscntmdt2', d = 'mhlabhistabscntddt2', y = 'mhlabhistabscntydt2')
)

nhs_aec_export_fields <- c('local_id', 'visit_date',
                           unlist(lapply(nhs_aec_fields, unname)))

nhs_aec <- function(nhs_df, apply_cutoff = TRUE) {
  extract_one <- function(f) {
    nhs_df %>%
      transmute(
        local_id,
        #an explicit format returns NA for blank/partial dates instead of erroring
        visit_date = as.Date(as.character(visit_date), format = '%Y-%m-%d'),
        value = suppressWarnings(as.numeric(.data[[f['value']]])),
        date = as.Date(sprintf('%s-%s-%s',
                               .data[[f['y']]], .data[[f['m']]], .data[[f['d']]]),
                       format = '%Y-%m-%d')
      ) %>%
      filter(!is.na(value))
  }
  out <- bind_rows(lapply(nhs_aec_fields, extract_one))
  if (apply_cutoff && !is.na(visit_date_cutoff)) {
    out <- filter(out, is.na(visit_date) | visit_date <= visit_date_cutoff)
  }
  #the same count appears at the visit it was drawn and again in the next
  #visit's carry-forward field, so duplicates are dropped
  out %>%
    filter(!is.na(local_id), !is.na(date)) %>%
    distinct(local_id, date, value) %>%
    arrange(local_id, date)
}

#One NHS AEC value per participant per day
nhs_aec_daily <- function(aec_df) {
  aec_df %>%
    group_by(local_id, date) %>%
    summarise(aec_nhs = mean(value), .groups = 'drop')
}

# --- Extraction: NHS growth --------------------------------------------------

nhs_growth_export_fields <- c('local_id', 'visit_date_mh', 'mhheight', 'mhweight',
                              'dateofbirth', 'gender')

#Rename the NHS growth export to the shared column names.
#
#As on the EHR side, require_both = TRUE keeps only visits at which both height
#and weight were recorded, so both sources are treated the same way.
nhs_growth <- function(nhs_df, require_both = TRUE) {
  x <- c('local_id','redcap_event_name','redcap_repeat_instrument',
         'redcap_repeat_instance','dob','sex','vitals_time','height','weight')
  stopifnot(ncol(nhs_df) == length(x))
  colnames(nhs_df) <- x
  nhs_df$height <- as.numeric(as.character(nhs_df$height))
  nhs_df$weight <- as.numeric(as.character(nhs_df$weight))
  nhs_df$BMI    <- nhs_df$weight / ((nhs_df$height / 100)^2)
  nhs_df$age    <- age_in_months(nhs_df$vitals_time, nhs_df$dob)
  #normalise to a plain calendar date so it can be matched to the EHR side
  nhs_df$vitals_time <- substr(as.character(nhs_df$vitals_time), 1, 10)
  nhs_df <- as.data.frame(nhs_df)
  if (require_both) {
    nhs_df <- nhs_df[!is.na(nhs_df$height) & !is.na(nhs_df$weight), , drop = FALSE]
  }
  nhs_df
}

# --- Extraction: problem lists ----------------------------------------------

#One comma-separated list of EGID ICD-10 codes per participant, from the EHR
ehr_problem_list <- function(prob_df) {
  df3p <- prob_df %>%
    filter(problem_icd10_code %in% icd_egid) %>%
    select(local_id, problem_icd10_code) %>%
    group_by(local_id) %>%
    summarise(problem_icd10_code = paste(problem_icd10_code, collapse = ","),
              .groups = 'drop')
  supp1 = data.frame()
  for (i in 1:nrow(df3p))
  {
    row <- df3p[i,]
    l1 <- as.list(strsplit(row$problem_icd10_code, ",")[[1]])
    l2 <- unique(unlist(l1))
    row$problem_icd10_code <- stri_paste(l2, collapse = ',')
    supp1 <- rbind(supp1, row)
  }
  supp1 %>% rename(FHIR_problem = problem_icd10_code)
}

nhs_problem_export_fields <- c('local_id', 'mheoediagnosed', 'mhegdiagnosed',
                               'mhecdiagnosed', 'mhegediagnosed', 'selectcondition',
                               'visit_date_mh')

#Map the NHS diagnosis questions onto the same ICD-10 codes. A condition counts
#as present if it was answered yes, or if the answer is missing/unknown but the
#participant was enrolled under that condition.
nhs_problem_list <- function(nhs_df) {
  result <- nhs_df[with(nhs_df, {
    !(is.na(selectcondition) &
        is.na(mheoediagnosed) &
        is.na(mhegdiagnosed) &
        is.na(mhecdiagnosed) &
        is.na(mhegediagnosed))
  }), ]
  df2p<-result
  df2p[, 'problem_icd10_code'] = "NA"
  df2p$mheoediagnosed[is.na(df2p$mheoediagnosed)] <- 9999
  df2p$mhegdiagnosed[is.na(df2p$mhegdiagnosed)] <- 9999
  df2p$mhecdiagnosed[is.na(df2p$mhecdiagnosed)] <- 9999
  df2p$mhegediagnosed[is.na(df2p$mhegediagnosed)] <- 9999
  df2p$selectcondition[is.na(df2p$selectcondition)] <- 9999
  df4p = data.frame()
  for (i in 1:nrow(df2p))
  {
    lst = list()
    row <- df2p[i,]
    if (row$mheoediagnosed == '1')
    {
      lst <- list.append(lst, 'K20.0')
    }
    if ((row$mheoediagnosed == '9999' | row$mheoediagnosed == '3333' | row$mheoediagnosed == '0') & row$selectcondition == '1')
    {
      lst <- list.append(lst, 'K20.0')
    }
    if (row$mhegdiagnosed == '1')
    {
      lst <- list.append(lst, 'K52.81')
    }
    if ((row$mhegdiagnosed == '9999' | row$mhegdiagnosed == '3333' | row$mhegdiagnosed == '0') & (row$selectcondition == '2' | row$selectcondition == '5'))
    {
      lst <- list.append(lst, 'K52.81')
    }
    if (row$mhecdiagnosed == '1')
    {
      lst <- list.append(lst, 'K52.82')
    }
    if ((row$mhecdiagnosed == '9999' | row$mhecdiagnosed == '3333' | row$mhecdiagnosed == '0') & row$selectcondition == '3')
    {
      lst <- list.append(lst, 'K52.82')
    }
    if (row$mhegediagnosed == '1')
    {
      lst <- list.append(lst, 'K52.81')
    }
    if ((row$mhegediagnosed == '9999' | row$mhegediagnosed == '3333' | row$mhegediagnosed == '0') & (row$selectcondition == '2' | row$selectcondition == '5'))
    {
      lst <- list.append(lst, 'K52.81')
    }
    c1 = unique(unlist(lst))
    c2 <- stri_paste(c1, collapse = ',')
    row$problem_icd10_code <- c2
    df4p <- rbind(df4p, row)
  }
  supp_nhs <- df4p %>% group_by(local_id) %>%
    summarise(problem_icd10_code = paste(problem_icd10_code, collapse = ","),
              .groups = 'drop')
  supp_nhs_2 = data.frame()
  for (i in 1:nrow(supp_nhs))
  {
    row <- supp_nhs[i,]
    l1 <- as.list(strsplit(row$problem_icd10_code, ",")[[1]])
    l2 <- unique(unlist(l1))
    row$problem_icd10_code <- stri_paste(l2, collapse = ',')
    supp_nhs_2 <- rbind(supp_nhs_2, row)
  }
  supp_nhs_2 %>%
    select(local_id, problem_icd10_code) %>%
    rename(NHS_problem = problem_icd10_code)
}

# --- Matching NHS entries to EHR observations -------------------------------

#Pair each NHS entry with the closest EHR observation for the same participant,
#within match_window_days. An exact same-date match always wins, since its
#offset is zero. Ties are broken toward the earlier EHR date so the result does
#not depend on row order.
#
#Both frames need an id column and a date column; date columns may be Date or
#"YYYY-MM-DD" character. Returns one row per matched NHS entry, with the offset
#in days so it can be reported.
match_nearest_date <- function(nhs_df, ehr_df,
                               nhs_date = 'date', ehr_date = 'date',
                               id = 'local_id',
                               max_days = match_window_days) {
  as_date <- function(x) {
    if (inherits(x, "Date")) return(x)
    as.Date(substr(as.character(x), 1, 10), format = '%Y-%m-%d')
  }

  nhs_df <- as.data.frame(nhs_df)
  ehr_df <- as.data.frame(ehr_df)
  nhs_df$.nhs_date <- as_date(nhs_df[[nhs_date]])
  ehr_df$.ehr_date <- as_date(ehr_df[[ehr_date]])
  nhs_df <- nhs_df[!is.na(nhs_df$.nhs_date), , drop = FALSE]
  ehr_df <- ehr_df[!is.na(ehr_df$.ehr_date), , drop = FALSE]
  if (!nrow(nhs_df) || !nrow(ehr_df)) return(nhs_df[0, , drop = FALSE])

  #every NHS/EHR combination for the same participant, then keep the closest
  m <- merge(nhs_df, ehr_df, by = id, suffixes = c('.nhs', '.ehr'))
  if (!nrow(m)) return(m)

  m$offset_days <- abs(as.numeric(m$.nhs_date - m$.ehr_date))
  m <- m[m$offset_days <= max_days, , drop = FALSE]
  if (!nrow(m)) return(m)

  #one row per NHS entry: smallest offset wins, ties broken to the earlier date
  key <- paste(m[[id]], m$.nhs_date)
  m <- m[order(key, m$offset_days, m$.ehr_date), , drop = FALSE]
  m[!duplicated(paste(m[[id]], m$.nhs_date)), , drop = FALSE]
}

#One-line description of how well the dates lined up, for the console output
describe_match <- function(paired) {
  if (!nrow(paired)) return("no pairs")
  sprintf("%d pairs (%d exact date, %d within %d day(s), max offset %g)",
          nrow(paired), sum(paired$offset_days == 0),
          sum(paired$offset_days > 0), match_window_days, max(paired$offset_days))
}

#Report every stage of the pairing for one site, to locate where rows are lost.
#Pass the two daily frames, e.g.
#   diagnose_pairing(nhs_daily, ehr_daily, "colorado", ehr_sites)
diagnose_pairing <- function(nhs_daily, ehr_daily, dag, sites_df,
                             nhs_date = 'date', ehr_date = 'date') {
  ids <- site_participants(sites_df, dag)
  nd <- nhs_for_site(nhs_daily, ids)
  ed <- semi_join(as.data.frame(ehr_daily), ids, by = 'local_id')

  cat("\n--- pairing diagnostic:", dag, "---\n")
  cat("site participants          :", nrow(ids), "\n")
  cat("NHS rows for site          :", nrow(nd), "\n")
  cat("EHR rows for site          :", nrow(ed), "\n")
  cat("participants in both       :", length(intersect(nd$local_id, ed$local_id)), "\n")
  cat("NHS date column class      :", class(nd[[nhs_date]]), "\n")
  cat("EHR date column class      :", class(ed[[ehr_date]]), "\n")
  if (nrow(nd)) cat("NHS example dates          :", paste(utils::head(as.character(nd[[nhs_date]]), 3), collapse = ", "), "\n")
  if (nrow(ed)) cat("EHR example dates          :", paste(utils::head(as.character(ed[[ehr_date]]), 3), collapse = ", "), "\n")

  m <- match_nearest_date(nd, ed, nhs_date = nhs_date, ehr_date = ehr_date)
  cat("pairs within", match_window_days, "day(s)      :", nrow(m), "\n")

  #the smallest gap actually available, ignoring the window
  wide <- match_nearest_date(nd, ed, nhs_date = nhs_date, ehr_date = ehr_date,
                             max_days = Inf)
  if (nrow(wide)) {
    cat("nearest-offset distribution (days, no window):\n")
    print(table(cut(wide$offset_days, c(-1, 0, 1, 3, 7, 14, 30, 90, 365, Inf),
                    labels = c("0","1","2-3","4-7","8-14","15-30","31-90","91-365",">365"))))
  } else {
    cat("no participant has data in both sources at this site\n")
  }
  invisible(m)
}

# --- Formatting --------------------------------------------------------------

#Statistics are computed at full precision and rounded once, here, for display.

#median (IQR)
desc_continuous <- function(x, scale = 1) {
  x <- as.numeric(x[!is.na(x)]) * scale
  if (!length(x)) return(NA_character_)
  fmt <- function(v) formatC(round(v, n_digits), format = "f", digits = n_digits)
  sprintf("%s (%s - %s)", fmt(median(x)), fmt(quantile(x, 0.25)), fmt(quantile(x, 0.75)))
}

#n (%), with optional labels, a category to fold missing values into, and a
#fixed display order (listed categories always appear, as 0 if unobserved)
desc_categorical <- function(x, labels = NULL, merge_missing_into = NULL,
                             category_order = NULL) {
  total <- length(x)
  x_lab <- if (!is.null(labels)) ifelse(x %in% names(labels), labels[x], x) else x
  x_lab[is.na(x_lab) | x_lab == 'NA'] <- "Missing"
  if (!is.null(merge_missing_into)) {
    x_lab[x_lab == "Missing"] <- merge_missing_into
  }
  tab <- table(x_lab)
  if (!is.null(category_order)) {
    missing_cats <- setdiff(category_order, names(tab))
    if (length(missing_cats) > 0) {
      tab <- c(tab, setNames(rep(0, length(missing_cats)), missing_cats))
    }
    tab <- tab[c(category_order, setdiff(names(tab), category_order))]
  }
  pct <- 100 * as.vector(tab) / total
  data.frame(category = names(tab),
             n_pct = sprintf("%d (%s%%)", as.vector(tab),
                             formatC(round(pct, n_digits), format = "f", digits = n_digits)))
}

#estimate [lower - upper]
ci_fmt <- function(est, lower, upper) {
  fmt <- function(v) formatC(round(v, n_digits), format = "f", digits = n_digits)
  sprintf("%s [%s - %s]", fmt(est), fmt(lower), fmt(upper))
}

#Cohen's kappa with 95% CI
kappa_fmt <- function(fit) ci_fmt(fit$kappa, fit$confid[1, 1], fit$confid[1, 3])

#Bland-Altman bias with 95% limits of agreement
loa_fmt <- function(fit) ci_fmt(fit$loa$estimate[1], fit$loa$estimate[2], fit$loa$estimate[3])

#Simple agreement with a 95% Wilson interval
simple_agreement <- function(score) {
  ci <- binom::binom.confint(round(sum(score)), length(score), conf.level = 0.95)
  ci <- ci[ci$method == "wilson", ]
  ci_fmt(mean(score), ci$lower, ci$upper)
}

#Trim values beyond 1.5*IQR
remove_outliers <- function(df, col) {
  Q <- quantile(df[[col]], probs = c(.25, .75), na.rm = TRUE)
  iqr <- IQR(df[[col]], na.rm = TRUE)
  subset(df, df[[col]] > (Q[1] - 1.5 * iqr) & df[[col]] < (Q[2] + 1.5 * iqr))
}

# --- Figures -----------------------------------------------------------------

#Colours and marker shapes only distinguish individuals; they carry no meaning
plot_palette <- c("#D2B48C",  # tan
                  "#0000FF",  # blue
                  "#228B22",  # forestgreen
                  "#FFFF00",  # yellow
                  "#696969")  # dimgray
plot_shapes  <- c(15, 11, 13, 5, 1)

#Draw participants at random for a figure. Sampling is restricted to
#participants who have a value in both sources, so every selected line has at
#least one red NHS marker to compare against. The seed makes it reproducible.
sample_plot_ids <- function(candidate_ids, n = plot_n_participants, seed = plot_seed) {
  candidate_ids <- sort(unique(as.character(candidate_ids)))
  if (!length(candidate_ids)) return(character(0))
  set.seed(seed)
  sample(candidate_ids, min(n, length(candidate_ids)))
}

# --- Problem-list legend labels ---------------------------------------------

#Legend entries name each participant's diagnoses and where they were recorded,
#so the figures identify no one. Used by both Figure 2 panels, which is why this
#lives here rather than in either script.
#
#(EHR and NHS) = recorded in both sources; (EHR only) / (NHS only) = one source.
problem_code_map <- c("K20.0" = "EoE", "K52.81" = "EoGN", "K52.82" = "EoC")

problem_as_cats <- function(codes) {
  if (is.na(codes) || !nzchar(codes)) return(character(0))
  unname(problem_code_map[trimws(strsplit(codes, ",")[[1]])])
}

problem_label <- function(fhir, nhs) {
  f <- problem_as_cats(fhir); n <- problem_as_cats(nhs)
  both <- f[f %in% n]; fhir_only <- f[!f %in% n]; nhs_only <- n[!n %in% f]
  parts <- character(0)
  if (length(both))      parts <- c(parts, paste0(paste(both, collapse = ", "), " (EHR and NHS)"))
  if (length(fhir_only)) parts <- c(parts, paste0(paste(fhir_only, collapse = ", "), " (EHR only)"))
  if (length(nhs_only))  parts <- c(parts, paste0(paste(nhs_only,  collapse = ", "), " (NHS only)"))
  paste(parts, collapse = "; ")
}

#Named vector of legend labels, keyed by local_id in factor-level order, for a
#frame already carrying FHIR_problem and NHS_problem
problem_legend_labels <- function(d) {
  key <- unique(as.data.frame(d)[, c("local_id", "FHIR_problem", "NHS_problem")])
  key <- key[order(key$local_id), ]
  setNames(mapply(problem_label, key$FHIR_problem, key$NHS_problem),
           as.character(key$local_id))
}

# --- Figure 2: one participant sample shared by both panels ------------------

#Figures 2a and 2b must show the same people, so the sample cannot be drawn
#inside either script from its own data: height pairs and AEC pairs are
#different sets of participants, and sampling each separately gives two
#different fives even with one seed.
#
#Instead the sample is drawn once, at each site, from participants who have at
#least one same-date EHR/NHS pair in BOTH measures. That intersection is the
#only pool for which every plotted line carries a red NHS marker in both panels.
#
#The selection is cached to output/figure2_participants.csv. Whichever script
#runs first does the work; the second reads the file. Delete it, or call with
#refresh = TRUE, to draw a new sample.

figure2_ids_file <- function() file.path(output_dir, "figure2_participants.csv")

#The Figure 2a series for one site: same-date EHR/NHS height pairs in meters,
#after the same age limit and outlier trimming the figure applies. Shared with
#the figure script so the candidate pool and the plotted data cannot drift.
figure2a_site_data <- function(site_vitals, demo, nhs_height) {
  daily <- ehr_growth_daily(ehr_growth_zscores(ehr_growth_ages(
    ehr_vitals_measures(site_vitals), demo)))
  under_20 <- filter(daily, mean_age < 240)
  #Trim height outliers on the measure, then on its z-score, and convert to meters
  trimmed <- remove_outliers(remove_outliers(under_20, "mean_height"), "mean_heightz")
  trimmed <- mutate(trimmed, mean_height = mean_height / 100)
  trimmed %>% left_join(nhs_height, by = c('local_id', 'vitals_time'))
}

#The Figure 2b series for one site
figure2b_site_data <- function(site_labs, nhs_daily, site_ids, demo) {
  ehr_aec_daily(ehr_aec(site_labs)) %>%
    left_join(nhs_for_site(nhs_daily, site_ids), by = c('local_id', 'date')) %>%
    left_join(demo, by = 'local_id') %>%
    mutate(age = age_in_months(date, dob))
}

#Participants with at least one NHS value to compare against, i.e. at least one
#red marker on their line
paired_participants <- function(df, nhs_col) {
  df %>%
    filter(!is.na(.data[[nhs_col]])) %>%
    distinct(local_id) %>%
    pull(local_id) %>%
    as.character()
}

#NHS heights in meters, averaged per participant-date, as both panels need them
nhs_height_daily <- function(nhs_df) {
  nhs_df %>%
    mutate(height_nhs = height / 100) %>%
    filter(!is.na(local_id), !is.na(vitals_time), !is.na(height_nhs)) %>%
    group_by(local_id, vitals_time) %>%
    summarise(height_nhs = mean(height_nhs), .groups = 'drop')
}

#Named list, one character vector of local_ids per site
figure2_participant_ids <- function(refresh = FALSE, path = figure2_ids_file()) {
  if (!refresh && file.exists(path)) {
    cached <- utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character")
    cat("Figure 2 participants read from ", path, "\n", sep = "")
    return(split(cached$local_id, cached$dag))
  }

  cat("Selecting Figure 2 participants; this pulls both vitals and labs so that\n",
      "both panels can show the same people.\n", sep = "")

  vitals <- fill_mrn(redcap_export(
    redcap_token_ehr,
    c(redcap_fields(c('record_id', 'mrn')), list("forms[0]" = "vital_signs"))))
  vitals$local_id <- vitals$mrn

  labs <- fill_mrn(redcap_export(
    redcap_token_ehr,
    redcap_fields(c('record_id', 'mrn', 'labs_label',
                    'labs_loinc_code', 'labs_time', 'labs_value'))))
  labs$local_id <- labs$mrn

  demo <- redcap_export(redcap_token_ehr, list("forms[0]" = "demographics_5a25")) %>%
    rename(local_id = mrn) %>%
    select(local_id, dob, sex)

  nhs_h <- nhs_height_daily(nhs_growth(redcap_export(
    redcap_token_nhs, redcap_fields(nhs_growth_export_fields), export_dags = 'false')))

  nhs_a <- nhs_aec_daily(nhs_aec(redcap_export(
    redcap_token_nhs, redcap_fields(nhs_aec_export_fields))))

  #Site membership from both exports; script 10 reads it from vitals and script
  #11 from labs, and a participant must appear in both to be eligible anyway
  ehr_sites <- bind_rows(
    vitals %>% filter(!is.na(local_id)) %>% distinct(local_id, redcap_data_access_group),
    labs   %>% filter(!is.na(local_id)) %>% distinct(local_id, redcap_data_access_group)
  ) %>% distinct()

  dags <- ehr_data_access_groups(ehr_sites)
  chosen <- list()
  for (dag in dags) {
    site_ids <- site_participants(ehr_sites, dag)
    sv <- semi_join(vitals, site_ids, by = 'local_id')
    sl <- semi_join(labs,   site_ids, by = 'local_id')
    if (!nrow(sv) || !nrow(sl)) {
      cat("  ", dag, ": no vitals or no labs at this site\n", sep = "")
      chosen[[dag]] <- character(0)
      next
    }
    with_height <- paired_participants(
      figure2a_site_data(sv, demo, nhs_h), "height_nhs")
    with_aec <- paired_participants(
      figure2b_site_data(sl, nhs_a, site_ids, select(demo, local_id, dob)), "aec_nhs")
    both <- intersect(with_height, with_aec)
    chosen[[dag]] <- sample_plot_ids(both)
    cat(sprintf("  %s: %d with height pairs, %d with AEC pairs, %d with both -> %s\n",
                dag, length(with_height), length(with_aec), length(both),
                if (length(chosen[[dag]])) paste(chosen[[dag]], collapse = ", ") else "none"))
  }

  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(
    data.frame(dag = rep(names(chosen), lengths(chosen)),
               local_id = as.character(unlist(chosen, use.names = FALSE)),
               stringsAsFactors = FALSE),
    path, row.names = FALSE)
  cat("Figure 2 participants written to ", path, "\n", sep = "")
  chosen
}
