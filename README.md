# EHR-NHS-Analysis

Agreement and descriptive analyses comparing FHIR-extracted electronic health
record (EHR) data against manually entered Natural History Study (NHS) data.

Analysis code for the manuscript *A middleware framework for secure, multi-site
EHR data ingestion into REDCap via HL7 FHIR: augmenting a rare disease natural
history study*.

Every script reports results **for each site automatically**. The list of sites
is read from the EHR project's data access groups at runtime, so no site name is
hardcoded and adding a site in REDCap needs no code change.

## Layout

```
R/
  00_utils.R                          shared config, REDCap access, extraction, formatting
  01_demographics_agreement.R         Table 1 (demographics), Table 2a (kappa), Table S1
  02_growth_nhs_table1.R              Table 1 (NHS growth)
  03_labs_nhs_table1.R                Table 1 (NHS labs, AEC)
  04_growth_ehr_table3.R              Table 3 (EHR growth)
  05_labs_ehr_table3.R                Table 3 (EHR labs, AEC)
  06_growth_agreement.R               Table 2b (growth, Bland-Altman)
  07_labs_agreement.R                 AEC agreement (Bland-Altman)
  08_problem_list_agreement.R         Table 2a (problem lists, simple agreement)
  09_bmiz_aec_correlation.R           supplementary BMIz vs AEC coherence check
  10_figure2a_height_trajectories.R   Figure 2a (height; no legend)
  11_figure2b_aec_trajectories.R      Figure 2b (AEC; problem list in the legend)
output/                               gitignored; figures and exports stay local
EHR-NHS-Analysis.Rproj                RStudio project; opening it sets the working directory
```

Every script is self-contained: it sources `R/00_utils.R` itself and depends on
none of the others. The numeric prefixes only fix the order in a file listing.

## Requirements

R >= 4.4 and:

```r
install.packages(c("dplyr", "tidyr", "collapse", "psych", "SimplyAgree", "binom",
                   "rlist", "stringi", "httr", "ggplot2", "ggstar", "scales"))
install.packages('https://raw.github.com/CDC-DNPAO/CDCAnthro/master/cdcanthro_0.1.3.tar.gz',
                 type = 'source', repos = NULL)
```

## Configuration

Credentials come from environment variables. No token is stored in this
repository, and none should be added to it.

| Variable | Description |
|---|---|
| `REDCAP_API_TOKEN_EHR` | Token for the EHR/FHIR project (export rights) |
| `REDCAP_API_TOKEN_NHS` | Token for the NHS project (export rights) |

The simplest permanent setup is a `~/.Renviron` file, which is gitignored:

```
REDCAP_API_TOKEN_EHR=yourEhrToken
REDCAP_API_TOKEN_NHS=yourNhsToken
```

Run `usethis::edit_r_environ()` to open it, save, then **restart R** — `.Renviron`
is only read at startup. Confirm with:

```r
nchar(Sys.getenv("REDCAP_API_TOKEN_EHR"))   # should be 32, not 0
```

An empty token produces the same "You do not have permissions to use the API"
message as an unauthorised one, so this check is worth doing first.

## Running

Open `EHR-NHS-Analysis.Rproj` in RStudio, which sets the working directory,
then open a script and run it. Either the repository root or the `R/` folder
works as the working directory.

```bash
Rscript R/01_demographics_agreement.R
Rscript R/06_growth_agreement.R
```

Results print to the console under a heading per site. Figures 2a and 2b are
saved to `output/` as PNG at 600 dpi and PDF, and are also left in the session as
`p1`, `p2`, ... in site order.

Other settings, all at the top of `R/00_utils.R`:

| Setting | Purpose |
|---|---|
| `excluded_local_ids` | participants dropped during QC |
| `visit_date_cutoff` | `NA` uses all current data; set a date to reproduce an earlier export snapshot |
| `match_window_days` | how far apart an EHR and NHS observation may be and still be paired (agreement scripts only) |
| `n_digits` | decimal places for reported statistics |
| `plot_n_participants`, `plot_seed` | sample size and seed for the figures |

## Notes on the data

These are the decisions that most affect the numbers. Several were established
by checking the raw exports directly.

**Sites are matched by participant, not by group name.** The EHR and NHS projects
use different data access group names for the same site, so NHS records are
attributed to a site by whether the participant appears in that site's EHR data
access group. Site membership is taken from the export the measurements come
from, not from the demographics form.

**A growth record is a day with both height and weight.** `ehr_vitals_measures()`
and `nhs_growth()` default to `require_both = TRUE`, so height, weight, BMI and
all three z-scores rest on the same set of participant-days. Days holding only
one of the two are dropped rather than completed from another visit. This matters:
about a quarter of EHR days have a weight but no height.

**Height and weight are paired on the calendar date, not the timestamp.** Sites
differ in whether the two share a timestamp, and pairing on the timestamp would
separate measures taken minutes apart. Repeat measures on a day are averaged.

**BMI is derived, not exported.** Only the height and weight LOINC codes are
read; BMI is computed as `weight / (height/100)^2`.

**AEC is derived too.** In the EHR it is the white cell count multiplied by the
eosinophil percentage recorded at the same timestamp, divided by 100. In the NHS
project it appears in two places — the count taken at the current visit and the
count carried forward from the previous visit — and both are needed: roughly half
the NHS values appear only in the carry-forward field.

**Timestamp formats differ by site.** Some FHIR servers return
`YYYY-MM-DD HH:MM:SS` and others ISO 8601 `YYYY-MM-DDTHH:MM:SS`, so dates are
taken with `substr(x, 1, 10)` rather than by splitting on a space.

**Ages use explicit day units.** `age_in_months()` calls
`difftime(..., units = "days")`. Subtracting `POSIXlt` values uses `units = "auto"`,
which picks one unit for the whole vector from its smallest element, so a single
same-day record can make every age come back in hours.

**Growth z-scores cover ages 2 to <20 years.** `cdcanthro` returns nothing
outside that range, so the 20-and-over rows of Table 3 are summarised from the
measurements before the z-score step, and records under 2 years appear in neither
age band.

**`visit_date_cutoff` reproduces a snapshot.** The NHS lab fields record the draw
date separately from the visit at which the value was entered, so values with old
draw dates keep entering the export as new visits occur. Left at `NA`, results
reflect current data and will exceed the counts in the published tables.

**Rounding happens once, at display.** Statistics are computed at full precision;
`R/00_utils.R` rounds only when formatting.

**Figure participants are sampled at random** from those with a same-date
EHR/NHS pair in *both* height and AEC at that site, so every plotted line carries
an NHS marker to compare against in both panels. `plot_seed` makes the selection
reproducible. The sample is drawn once by `figure2_participant_ids()` and cached
in `output/figure2_participants.csv`, so Figures 2a and 2b show the same
participants whichever script runs first. Delete that file, or call
`figure2_participant_ids(refresh = TRUE)`, to draw a new sample.

## Diagnostics

`diagnose_pairing()` in `R/00_utils.R` reports each stage of an EHR/NHS pairing
for one site — row counts, shared participants, date column classes, example
dates, and the distribution of nearest available offsets. Useful when a site
returns fewer pairs than expected:

```r
diagnose_pairing(nhs_daily, ehr_aec_daily(ehr_aec(labs)), "cincinnati", ehr_sites)
```

## Data handling

Only de-identified, aggregate results should leave the secure analysis
environment. `output/` is gitignored because figures and intermediate exports are
derived from participant-level data. Do not commit tokens, participant
identifiers, or REDCap exports.
