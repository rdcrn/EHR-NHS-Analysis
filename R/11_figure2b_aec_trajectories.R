#!/usr/bin/env Rscript
# =============================================================================
# 11_figure2b_aec_trajectories.R
#
# Figure 2b - absolute eosinophil count trajectories for a random sample of
#             participants at each site. Black markers are EHR observations,
#             red markers the NHS observations recorded on the same date.
#             Legend entries give each participant's problem list and the source
#             of each diagnosis, with no participant identifiers shown.
#
# The participants are the same ones Figure 2a shows: the sample is drawn once
# by figure2_participant_ids() in 00_utils.R, from those with a same-date
# EHR/NHS pair in BOTH height and AEC, and cached in output/ so both scripts
# read the same selection whichever runs first.
#
# One figure per site, returned in a named list and also assigned to p1, p2, ...
# in site order, so a single site's figure can be worked with directly.
# =============================================================================

#Load shared configuration and helpers. Works from the repository root or from
#inside the R/ folder, the two working directories the README supports.
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
  library(ggplot2)
  library(ggstar)
})

#The participants both Figure 2 panels show. Pass refresh = TRUE, or delete
#output/figure2_participants.csv, to draw a new sample.
figure2_ids <- figure2_participant_ids()

#Get lab data from REDcap with data pulled from FHIR
labs <- redcap_export(redcap_token_ehr,
                      redcap_fields(c('record_id', 'mrn', 'labs_label',
                                      'labs_loinc_code', 'labs_time', 'labs_value')))
labs <- fill_mrn(labs)
labs$local_id <- labs$mrn

ehr_sites <- labs %>%
  filter(!is.na(local_id)) %>%
  distinct(local_id, redcap_data_access_group)

#Get demographics data from REDcap with data pulled from FHIR
demo <- redcap_export(redcap_token_ehr, list("forms[0]" = "demographics_5a25")) %>%
  rename(local_id = mrn) %>%
  select(local_id, dob)

#Get problem list data from REDCap with data pulled from FHIR
prob_ehr <- redcap_export(redcap_token_ehr,
                          c(redcap_fields(c('record_id', 'mrn')),
                            list("forms[0]" = "problem_list")))
prob_ehr <- fill_mrn(prob_ehr)
prob_ehr$local_id <- prob_ehr$mrn
data_fhir <- ehr_problem_list(prob_ehr)

#Get problem list data from REDCap with data from NHS
data_nhs <- nhs_problem_list(
  redcap_export(redcap_token_nhs, redcap_fields(nhs_problem_export_fields))
)

#Get lab data from REDCap with data from NHS
nhs_daily <- nhs_aec_daily(nhs_aec(
  redcap_export(redcap_token_nhs, redcap_fields(nhs_aec_export_fields))
))

# --- One figure per site -----------------------------------------------------

build_figure <- function(dag) {
  cat("\n\n==================== ", dag, " ====================\n")
  site_ids <- site_participants(ehr_sites, dag)
  site_labs <- semi_join(labs, site_ids, by = 'local_id')
  if (!nrow(site_labs)) {
    cat("No EHR labs for this site.\n")
    return(NULL)
  }

  #Same pipeline the participant selection used, so the plotted series and the
  #candidate pool cannot drift apart
  all_aec <- figure2b_site_data(site_labs, nhs_daily, site_ids, demo)

  chosen <- figure2_ids[[dag]]
  if (is.null(chosen) || !length(chosen)) {
    cat("No shared Figure 2 participants at this site.\n")
    return(NULL)
  }
  #a selected participant should always have rows here, but never plot an empty
  #factor level: it would add a legend entry with no line
  absent <- setdiff(chosen, unique(as.character(all_aec$local_id)))
  if (length(absent)) {
    warning("no EHR AEC at ", dag, " for: ", paste(absent, collapse = ", "))
    chosen <- setdiff(chosen, absent)
  }
  if (!length(chosen)) {
    cat("None of the shared participants have EHR AEC at this site.\n")
    return(NULL)
  }
  cat("Participants plotted: ", paste(chosen, collapse = ", "), "\n", sep = "")

  d <- all_aec %>%
    filter(local_id %in% chosen) %>%
    left_join(data_fhir, by = 'local_id') %>%
    left_join(data_nhs, by = 'local_id')
  d$local_id <- factor(as.character(d$local_id), levels = sort(chosen))

  legend_labels <- problem_legend_labels(d)
  print(legend_labels)

  p <- ggplot(d) +
    geom_line(aes(x = age, y = aec_ehr, color = local_id, group = local_id),
              linewidth = 1, key_glyph = "path") +
    geom_star(aes(x = age, y = aec_ehr, starshape = local_id),
              size = 2, fill = "black", color = "black") +
    geom_star(aes(x = age, y = aec_nhs, starshape = local_id),
              size = 4, fill = "red", color = "red") +
    scale_color_manual(name = NULL, values = plot_palette, labels = legend_labels) +
    scale_starshape_manual(values = plot_shapes, guide = "none") +
    guides(color = guide_legend(ncol = 1, order = 1)) +
    labs(x = "Age in months",
         y = expression(bold("Absolute eosinophil count (10"^"9"*"/L)"))) +
    #headroom so the inset legend does not sit on the data
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    theme_minimal() +
    theme(
      legend.position             = "inside",
      legend.position.inside      = c(0.99, 0.99),
      legend.justification.inside = c(1, 1),
      # for ggplot2 < 3.5.0 replace the three lines above with:
      # legend.position      = c(0.99, 0.99),
      # legend.justification = c(1, 1),
      legend.direction  = "vertical",
      legend.background = element_rect(fill = scales::alpha("white", 0.85), color = NA),
      legend.margin     = margin(3, 5, 3, 3),
      legend.text       = element_text(size = 8.5, color = "black"),
      legend.key.width  = unit(1.4, "lines"),
      legend.key.height = unit(0.95, "lines"),
      legend.spacing.y  = unit(0.05, "lines"),
      axis.title = element_text(size = 16, face = "bold", color = "black"),
      axis.text  = element_text(size = 14, face = "bold", color = "black")
    )
  print(p)
  p
}

dags <- ehr_data_access_groups(ehr_sites)
figures_2b <- lapply(dags, build_figure)
names(figures_2b) <- dags

#Also expose each site's figure as p1, p2, ... in site order, so a single plot
#can be printed or edited on its own: p1 is the first site in dags, p2 the second
for (i in seq_along(figures_2b)) assign(paste0("p", i), figures_2b[[i]])
cat("\nfigures available as: ", paste0("p", seq_along(dags), " = ", dags, collapse = ", "), "\n", sep = "")

#Save at publication resolution. 600 dpi suits line art and small markers;
#cairo_pdf embeds fonts properly for the superscript axis label.
fig_dir <- output_dir
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
for (i in seq_along(figures_2b)) {
  if (is.null(figures_2b[[i]])) next
  stem <- file.path(fig_dir, paste0("Figure2b_", dags[i]))
  ggsave(paste0(stem, ".png"), plot = figures_2b[[i]],
         width = 7.5, height = 5.5, units = "in", dpi = 600)
  ggsave(paste0(stem, ".pdf"), plot = figures_2b[[i]],
         width = 7.5, height = 5.5, units = "in", device = cairo_pdf)
  cat("saved ", stem, ".png and .pdf\n", sep = "")
}
