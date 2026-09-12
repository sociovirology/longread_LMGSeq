#Code modified by Claude based on Sam's pairwise_infections.R from 2023 paper

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(gtools)

# ============================================================================
# 1. Importing Data and Initial Preparation of Data Frame ----
# ============================================================================
# ONT/usearch strain assignment output format (one row per file x segment):
#   V1 file            - b6 filename (e.g., plate02_well01_90_merged.b6)
#   V2 seg             - assigned segment (e.g., "NA", "HA", "NONE")
#   V3 match_for_seg   - top-hit strain_segment (e.g., "PAN99_NA", "NONE")
#   V4 max_for_seg     - reads assigned to the top-hit strain at that segment
#   V5 total_seg       - total reads aligning to that segment (any strain)
#   V6 total_reads_well - total reads for the WELL, independent of strain
#                         assignment (constant across all segment rows for a
#                         given well). NOTE: "_90_merged" in the filename is
#                         leftover Illumina-era naming; there is no read
#                         merging step in this ONT pipeline.
#   V7 unsorted        - ambiguous/multi-locus hits; near-always 0
#
# No multi_hit column in this file version.

file_names <- list.files(
  "outputs/",   # adjust to your working directory
  pattern = "*_strain_assignment_output_all_samples.txt",
  full.names = TRUE
)

strain_assignment_raw <- NULL
for (file in file_names) {
  df <- read_tsv(
    file,
    col_names = c("file", "seg", "match_for_seg", "max_for_seg",
                  "total_seg", "total_reads_well", "unsorted"),
    na = character(),   # <-- ADD THIS: disables readr's default "NA"-string-as-missing
    #     behavior, since the flu NA (neuraminidase) segment is a
    #     literal string "NA" in this data, not a missing value.
    #     "NONE" remains the only sentinel for "nothing called."
    col_types = cols(
      file             = col_character(),
      seg              = col_character(),
      match_for_seg    = col_character(),
      max_for_seg      = col_double(),
      total_seg        = col_double(),
      total_reads_well = col_double(),
      unsorted         = col_double()
    )
  )
  strain_assignment_raw <- bind_rows(strain_assignment_raw, df)
}

nrow(strain_assignment_raw)

# Strip the (legacy) usearch file suffix, and build the plate_well key up front
# -- used both for the locus table and for joining sample metadata below.
strain_assignment_raw <- strain_assignment_raw %>%
  mutate(file = str_remove(file, "_90_merged\\.b6$"))

# Keep a WELL-level table before dropping anything -- this preserves all 96
# wells (including complete assignment failures) for QC purposes, e.g.
# tying back to read_counts_demultiplexed_*.txt to compute what fraction of
# raw reads made it through alignment. One row per well.
well_summary <- strain_assignment_raw %>%
  distinct(file, total_reads_well)

nrow(well_summary)  # should be 96 for a full plate

# Drop rows where NO segment was called at all for the well (seg == "NONE").
# These represent wells with zero usable aligned reads -- a well-level QC
# failure, not a locus-level miss -- and would break the string-splitting below.
n_none <- sum(strain_assignment_raw$seg == "NONE")
message(sprintf(
  "%d well(s) with no segments assigned at all (seg == 'NONE'); excluded from locus-level table (still present in well_summary).",
  n_none
))
strain_assignment <- strain_assignment_raw %>% filter(seg != "NONE")

# Recreate the original column names/semantics used throughout the rest of the script:
#   cross_sample_locus, total_reads, majority_assigned_reads, majority_strain_locus
two_strain_database17_98_locus <- strain_assignment %>%
  mutate(
    cross_sample_locus = paste(file, seg, sep = "_"),
    total_reads = total_seg,
    majority_assigned_reads = max_for_seg,
    majority_strain_locus = match_for_seg
  ) %>%
  select(cross_sample_locus, total_reads, majority_assigned_reads,
         majority_strain_locus, total_reads_well, unsorted, file, seg)

# Proportion of reads assigned to the majority strain at each locus (same as original)
two_strain_database17_98_locus <- mutate(
  two_strain_database17_98_locus,
  proportion_assigned = majority_assigned_reads / total_reads
)

# QC plot -- same as original script
qplot(two_strain_database17_98_locus$proportion_assigned, two_strain_database17_98_locus$total_reads)

# Split file into cross ("plateNN") and sample ("wellNN")
two_strain_database17_98_locus <- separate(
  two_strain_database17_98_locus, file, c("cross", "sample"), sep = "_", remove = FALSE
)
two_strain_database17_98_locus <- rename(two_strain_database17_98_locus, locus = seg)
two_strain_database17_98_locus <- unite(
  two_strain_database17_98_locus, cross_sample, c("cross", "sample"), sep = "_", remove = FALSE
)

# Split majority strain and majority locus (matches original)
two_strain_database17_98_locus <- separate(
  two_strain_database17_98_locus, majority_strain_locus,
  c("majority_strain", "majority_locus"), sep = "_", remove = FALSE
)

nrow(two_strain_database17_98_locus)

# ----------------------------------------------------------------------------
# Join per-well sample metadata (replaces cross_data_runA.csv join)
# ----------------------------------------------------------------------------
# The header row has a stray trailing space on "plate " -- trimming column
# names defensively rather than assuming the CSV is clean.
sample_list <- read_csv("sample_list_CA09xPAN99.csv") %>%
  rename_with(str_trim) %>%
  mutate(cross_sample = paste(plate, well, sep = "_"))

# Sanity check: every locus-table cross_sample should have exactly one metadata match.
# (For now this data set is plate02 only, all "coinfection" type CA09 x PAN99.)
unmatched <- anti_join(two_strain_database17_98_locus, sample_list, by = "cross_sample")
if (nrow(unmatched) > 0) {
  message(sprintf("%d row(s) in the locus table have no matching sample metadata -- check plate/well naming.", nrow(unmatched)))
}

two_strain_database17_98_locus <- left_join(two_strain_database17_98_locus, sample_list, by = "cross_sample")

# Reconstruct strainAB for compatibility with downstream sections that facet/group by it
two_strain_database17_98_locus <- mutate(
  two_strain_database17_98_locus,
  strainAB = paste(parent1_label, parent2_label, sep = "_")
)

# ============================================================================
# 3. Quality measures of strain assignments and reassortment frequencies
#    (Part A: per-sample/per-cross calculations -- generalizes to any number
#    of crosses, including the current n=1. Part B -- multi-cross comparative
#    stats, strain-level means, PID correlation, subtype tests -- is deferred
#    until more plates are analyzed; see note at end of this block.)
# ============================================================================

# Compatibility alias so downstream facet_grid(strainA ~ strainB) calls,
# copied verbatim from the original script, keep working with the new
# per-well metadata column names.
two_strain_database17_98_locus <- two_strain_database17_98_locus %>%
  mutate(strainA = parent1_label, strainB = parent2_label)

# Number of genome segments -- used below instead of hardcoding 254/256.
# Influenza A has 8 segments; change this if you ever apply this pipeline
# to a virus with a different segment count.
n_segments <- 8
theoretical_free_reassortment <- (2^n_segments - 2) / 2^n_segments  # 0.9921875 for n=8

# Canonical large-to-small influenza A segment order (replaces the old
# suffix-laden locus names like 'PB2f'/'NS1d' from the two-loci-per-antigenic-
# segment Illumina scheme, which no longer exist in the ONT/usearch output).
segment_order <- c('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS')

#### QC: how good are the strain assignments? ----
mean(two_strain_database17_98_locus$proportion_assigned)
sd(two_strain_database17_98_locus$proportion_assigned)

nrow(subset(two_strain_database17_98_locus, proportion_assigned > .80)) / nrow(two_strain_database17_98_locus)
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .75)) / nrow(two_strain_database17_98_locus)
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .70)) / nrow(two_strain_database17_98_locus)
nrow(subset(two_strain_database17_98_locus, proportion_assigned > .65)) / nrow(two_strain_database17_98_locus)

ggplot(two_strain_database17_98_locus, aes(x = proportion_assigned)) + geom_histogram()

ggplot(two_strain_database17_98_locus, aes(x = locus, y = sample)) +
  geom_point(aes(col = majority_strain, alpha = proportion_assigned), shape = 15, size = 6) +
  theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~ cross)

#### Build controls_df: segments typed per sample, majority-strain counts ----
controls_df <- two_strain_database17_98_locus

controls_df_loci_numbers <- controls_df %>%
  group_by(cross_sample) %>%
  summarise(
    total_segments = length(locus),
    max_majority_strain = max(table(majority_strain))
  )

controls_df <- right_join(controls_df, controls_df_loci_numbers)

ggplot(subset(controls_df, total_segments > 7), aes(x = locus, y = reorder(sample, max_majority_strain))) +
  geom_point(aes(col = majority_strain, alpha = proportion_assigned), shape = 15, size = 6) +
  theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~ cross, scales = "free_y")

ggplot(controls_df, aes(x = factor(locus, level = segment_order), y = sample)) +
  geom_point(aes(col = majority_strain), shape = 15, size = 6) +
  xlab("Segment") + ylab("Plaque Isolate") +
  facet_grid(strainA ~ strainB, scales = "free_y") +
  theme_tufte() +
  theme(text = element_text(size = 20, family = "Helvetica")) +
  theme(axis.text.x = element_text(angle = 90, size = 17)) +
  theme(axis.text.y = element_blank()) +
  theme(strip.text = element_text(face = "bold")) +
  theme(legend.position = "none")
ggsave("outputs/figure2A.pdf", width = 8, height = 11)

#### Parentals table ----
parentals <- controls_df %>%
  group_by(cross_sample, cross, sample) %>%
  summarise(parental = max_majority_strain == total_segments) %>%
  distinct()

table(parentals$cross, parentals$parental)

#### Reassortant classification and per-cross reassortment frequency ----
controls_df_number_parents <- controls_df %>%
  group_by(cross, cross_sample) %>%
  summarise(
    number_parents = length(unique(majority_strain)),
    max_majority_strain = max(table(majority_strain))
  )

controls_df_number_parents <- mutate(controls_df_number_parents,
                                     reassortant = ifelse(number_parents > 1, yes = 1, no = 0))

cross_stats <- controls_df_number_parents %>%
  group_by(cross) %>%
  summarise(
    clones = length(number_parents),
    reassortants = sum(reassortant)
  )

cross_stats <- mutate(cross_stats, prop_reassortant = reassortants / clones)

#### Sensitivity check: does restricting to complete (8-segment) genotypes ----
#### change the reassortment estimate? Valid within a single cross, so this
#### works fine at n=1 and scales automatically as you add crosses.
controls_df_number_parents_complete <- subset(controls_df, total_segments > 7) %>%
  group_by(cross, cross_sample) %>%
  summarise(number_parents = length(unique(majority_strain)))

controls_df_number_parents_complete <- mutate(controls_df_number_parents_complete,
                                              reassortant = ifelse(number_parents > 1, yes = 1, no = 0))

cross_stats_complete <- controls_df_number_parents_complete %>%
  group_by(cross) %>%
  summarise(
    clones = length(number_parents),
    reassortants = sum(reassortant)
  )

cross_stats_complete <- mutate(cross_stats_complete, prop_reassortant = reassortants / clones)

# Loop over every cross present (was 10 hardcoded prop.test() calls indexed
# [1] through [10] in the original -- this scales to however many crosses
# you have, now or in future runs).
for (i in seq_len(nrow(cross_stats))) {
  this_cross <- cross_stats$cross[i]
  complete_row <- cross_stats_complete[cross_stats_complete$cross == this_cross, ]
  if (nrow(complete_row) == 0) {
    message(sprintf("%s: no complete (8-segment) genotypes found -- skipping comparison.", this_cross))
    next
  }
  result <- prop.test(
    c(complete_row$reassortants, cross_stats$reassortants[i]),
    c(complete_row$clones, cross_stats$clones[i])
  )
  message(sprintf("%s: complete-vs-all prop.test p = %.4f", this_cross, result$p.value))
}

#### Binomial test against theoretical free reassortment (per cross) ----
for (i in seq_len(nrow(cross_stats))) {
  result <- binom.test(
    c(cross_stats$clones[i] - cross_stats$reassortants[i], cross_stats$reassortants[i]),
    p = theoretical_free_reassortment
  )
  message(sprintf(
    "%s: observed reassortment = %.4f, theoretical = %.4f, binom.test p = %.4g",
    cross_stats$cross[i], cross_stats$prop_reassortant[i], theoretical_free_reassortment, result$p.value
  ))
}

#### Join strainAB info and plot ----
cross_stats <- right_join(cross_stats, unique(subset(controls_df, select = c(cross, strainA, strainB, strainAB))))

ggplot(cross_stats, aes(x = reorder(strainAB, prop_reassortant), y = prop_reassortant)) +
  geom_col() +
  ylab("Proportion of Reassortant Plaque Isolates") +
  xlab("Strains in Experimental Coinfection") +
  geom_hline(yintercept = 0.40, linetype = 2, color = "grey", alpha = 0.75) +
  geom_hline(yintercept = theoretical_free_reassortment, linetype = 2, color = "red", alpha = 0.75) +
  theme_tufte() +
  theme(text = element_text(size = 20, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(panel.grid.major.y = element_line(color = "lightgray", size = 0.5))
ggsave("outputs/figure2B.pdf", width = 12, height = 8.5)

controls_df <- right_join(controls_df, cross_stats)

#### Segment-level strain bias (Figure 6 style) ----
locus_strain <- group_by(controls_df, cross, locus, majority_strain) %>%
  summarise(total_samples = length(unique(sample)))

locus_strain <- right_join(locus_strain, subset(controls_df, select = c("cross", "strainA", "strainB", "strainAB")))

ggplot(locus_strain, aes(x = factor(locus, level = segment_order), y = total_samples, fill = majority_strain)) +
  geom_col(position = "fill") +
  facet_grid(strainA ~ strainB) +
  xlab("Segment") + ylab("Proportion of Plaque Isolates Assigned to Each Strain") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2) +
  theme_tufte() +
  theme(text = element_text(size = 20, family = "Helvetica")) +
  theme(axis.text.x = element_text(angle = 90, size = 17)) +
  theme(strip.text = element_text(face = "bold")) +
  theme(legend.position = "none")

#### Same, restricted to reassortants only ----
controls_df_reassortant <- right_join(controls_df, subset(controls_df_number_parents, select = c("cross_sample", "reassortant")), by = "cross_sample")
controls_df_reassortants_only <- subset(controls_df_reassortant, reassortant == 1)

locus_strain_reassortants_only <- group_by(controls_df_reassortants_only, cross, locus, majority_strain) %>%
  summarise(total_samples = length(unique(sample)))

locus_strain_reassortants_only <- right_join(locus_strain_reassortants_only, subset(controls_df_reassortants_only, select = c("cross", "strainA", "strainB", "strainAB")))

ggplot(locus_strain_reassortants_only, aes(x = factor(locus, level = segment_order), y = total_samples, fill = majority_strain)) +
  geom_col(position = "fill") +
  facet_grid(strainA ~ strainB) +
  xlab("Segment") + ylab("Proportion of Plaque Isolates Assigned to Each Strain") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)

# ============================================================================
# 4. Calculate Segment Representation ----
#    Per-cross, per-locus test of whether the majority-strain proportion at
#    each segment falls within the 95% CI expected under free (50/50)
#    assortment. This is inherently per-cross, so it generalizes fine to your
#    current n=1 and needs no changes as you add more plates.
# ============================================================================

# Binomial 95% CI around a 50-50 split, for a given sample size. Used instead
# of a single fixed CI (e.g. the "0.396-0.604 for 96 trials" example in the
# original comments) because sample size varies by locus and by cross once
# you have wells with incomplete/dropped-out segments.
get_low_CI <- function(sample_size) {
  half_integer <- round(sample_size / 2, digits = 0)
  as.numeric(binom.test(c(half_integer, half_integer), p = 0.5)$conf.int[1])
}

get_high_CI <- function(sample_size) {
  half_integer <- round(sample_size / 2, digits = 0)
  as.numeric(binom.test(c(half_integer, half_integer), p = 0.5)$conf.int[2])
}

#### 4.0 All data ----
binomial_segment <- group_by(controls_df, cross, locus, majority_strain) %>%
  summarise(count = n()) %>%
  group_by(cross, locus) %>%
  summarise(proportion = max(count) / sum(count), total = sum(count))

binomial_segment <- data.frame(
  binomial_segment,
  ci_low = sapply(binomial_segment$total, get_low_CI),
  ci_high = sapply(binomial_segment$total, get_high_CI)
)

binomial_segment <- right_join(binomial_segment, unique(subset(controls_df, select = c("cross", "locus", "strainA", "strainB", "strainAB"))))

inside_ci <- binomial_segment$proportion >= binomial_segment$ci_low & binomial_segment$proportion <= binomial_segment$ci_high
binomial_segment <- data.frame(binomial_segment, inside_ci)

ggplot(binomial_segment, aes(x = locus, y = proportion, colour = locus)) +
  geom_point(aes(shape = inside_ci), size = 2) +
  geom_point(colour = "grey90", size = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  facet_grid(strainA ~ strainB) +
  geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)

ggplot(binomial_segment, aes(x = factor(locus, level = segment_order), y = proportion, colour = inside_ci)) +
  geom_point(aes(shape = inside_ci), size = 2) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Proportion of Plaque Isolates Assigned to Each Strain") + xlab("Segment") +
  geom_point(colour = "grey90", size = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  facet_grid(strainA ~ strainB) +
  geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2) +
  theme_tufte() +
  theme(text = element_text(size = 20, family = "Helvetica")) +
  theme(axis.text.x = element_text(angle = 90, size = 17)) +
  theme(strip.text = element_text(face = "bold")) +
  theme(legend.position = "none")

ggplot(binomial_segment, aes(x = inside_ci, y = proportion)) +
  geom_boxplot() + geom_point(aes(colour = strainAB)) + facet_wrap(~ locus)

# Summary by locus (pools across whatever crosses exist -- 1 row per locus
# right now, will average across crosses once you have more)
binomial_segment %>% group_by(locus) %>% summarise(
  count = n(),
  mean_proportion_max = mean(proportion),
  sd_proportion_max = sd(proportion),
  number_inside_ci = length(which(inside_ci == TRUE))
)

# Summary by cross
binomial_segment %>% group_by(cross, strainAB) %>% summarise(
  count = n(),
  mean_proportion = mean(proportion),
  sd_proportion = sd(proportion),
  number_inside_ci = length(which(inside_ci == TRUE))
)

# Relationship between reassortment rate and number of freely-assorting segments.
# NOTE: this is a cross-level correlation -- with n=1 cross it will plot a
# single point and mean nothing statistically. Kept in place (harmless, won't
# error) so it activates automatically once you have multiple plates/crosses.
segment_assortment_reassortment <- right_join(
  binomial_segment %>% group_by(cross, strainAB) %>% summarise(
    count = n(),
    mean_proportion = mean(proportion),
    sd_proportion = sd(proportion),
    number_inside_ci = length(which(inside_ci == TRUE))
  ),
  cross_stats
)

ggplot(segment_assortment_reassortment, aes(x = prop_reassortant, y = number_inside_ci)) +
  geom_point(aes(colour = strainAB), size = 6)

# ----------------------------------------------------------------------------
# Number of segments involved in reassortment, per cross
# ----------------------------------------------------------------------------
# A segment is "involved in reassortment" in a given cross if more than one
# parental strain was ever assigned as the majority strain at that locus,
# across all wells sampled for that cross. Max possible value is 8 (segments
# in an influenza A genome). This REPLACES the manually-tallied
# `segments_reassortment <- c(3,8,8,8,8,7,6,8,7,8)` vector from the original
# script (eyeballed off the Figure 6/6B plots) with a reproducible calculation
# straight from controls_df.
segments_reassortment_df <- controls_df %>%
  group_by(cross, locus) %>%
  summarise(
    n_strains_at_locus = length(unique(majority_strain)),
    segment_reassorted = n_strains_at_locus > 1
  ) %>%
  group_by(cross) %>%
  summarise(
    segments_reassorted = sum(segment_reassorted),
    segments_typed = n()   # how many of the 8 segments had ANY calls at all in this cross
  )

segments_reassortment_df <- left_join(
  segments_reassortment_df,
  unique(subset(controls_df, select = c(cross, strainA, strainB, strainAB)))
)

segments_reassortment_df

mean(segments_reassortment_df$segments_reassorted)
sd(segments_reassortment_df$segments_reassorted)   # NA with a single cross -- resolves once more plates are added

ggplot(segments_reassortment_df, aes(x = reorder(strainAB, segments_reassorted), y = segments_reassorted)) +
  geom_col() +
  ylab("Number of Segments Involved in Reassortment (out of 8)") +
  xlab("Strains in Experimental Coinfection") +
  ylim(0, 8) +
  theme_tufte() +
  theme(text = element_text(size = 20, family = "Helvetica")) +
  theme(axis.text.x = element_text(size = 12))

#### 4.1 Same, but reassortants only ----
ggplot(locus_strain_reassortants_only, aes(x = factor(locus, level = segment_order), y = total_samples, fill = majority_strain)) +
  geom_col(position = "fill") + facet_grid(strainA ~ strainB) +
  xlab("Segment") + ylab("Proportion of Plaque Isolates Assigned to Each Strain") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)

binomial_segment_reassortants_only <- group_by(controls_df_reassortants_only, cross, locus, majority_strain) %>%
  summarise(count = n()) %>%
  group_by(cross, locus) %>%
  summarise(proportion = max(count) / sum(count), total = sum(count))

# FIX: use this data frame's own $total (sample sizes among reassortants only),
# not binomial_segment$total (all-data sample sizes) as in the original.
binomial_segment_reassortants_only <- data.frame(
  binomial_segment_reassortants_only,
  ci_low = sapply(binomial_segment_reassortants_only$total, get_low_CI),
  ci_high = sapply(binomial_segment_reassortants_only$total, get_high_CI)
)

binomial_segment_reassortants_only <- right_join(binomial_segment_reassortants_only, unique(subset(controls_df_reassortants_only, select = c("cross", "locus", "strainA", "strainB", "strainAB"))))

inside_ci_reassortants_only <- binomial_segment_reassortants_only$proportion >= binomial_segment_reassortants_only$ci_low & binomial_segment_reassortants_only$proportion <= binomial_segment_reassortants_only$ci_high
binomial_segment_reassortants_only <- data.frame(binomial_segment_reassortants_only, inside_ci_reassortants_only)

ggplot(binomial_segment_reassortants_only, aes(x = locus, y = proportion, colour = locus)) +
  geom_point(aes(shape = inside_ci_reassortants_only), size = 2) +
  geom_point(colour = "grey90", size = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  facet_grid(strainA ~ strainB) +
  geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)

ggplot(binomial_segment_reassortants_only, aes(x = factor(locus, level = segment_order), y = proportion, colour = inside_ci_reassortants_only)) +
  geom_point(aes(shape = inside_ci_reassortants_only), size = 2) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Proportion of Plaque Isolates Assigned to Each Strain") + xlab("Segment") +
  geom_point(colour = "grey90", size = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
  facet_grid(strainA ~ strainB) +
  geom_hline(yintercept = 0.50, alpha = 0.75, linetype = 2)

ggplot(binomial_segment_reassortants_only, aes(x = inside_ci_reassortants_only, y = proportion)) +
  geom_boxplot() + geom_point(aes(colour = strainAB)) + facet_wrap(~ locus)

binomial_segment_reassortants_only %>% group_by(locus) %>% summarise(
  count = n(),
  mean_proportion_max = mean(proportion),
  sd_proportion_max = sd(proportion),
  number_inside_ci = length(which(inside_ci_reassortants_only == TRUE))
)

binomial_segment_reassortants_only %>% group_by(cross, strainAB) %>% summarise(
  count = n(),
  mean_proportion = mean(proportion),
  sd_proportion = sd(proportion),
  number_inside_ci = length(which(inside_ci_reassortants_only == TRUE))
)

segment_assortment_reassortment_only <- right_join(
  binomial_segment_reassortants_only %>% group_by(cross, strainAB) %>% summarise(
    count = n(),
    mean_proportion = mean(proportion),
    sd_proportion = sd(proportion),
    number_inside_ci = length(which(inside_ci_reassortants_only == TRUE))
  ),
  cross_stats
)

ggplot(segment_assortment_reassortment_only, aes(x = prop_reassortant, y = number_inside_ci)) +
  geom_point(aes(colour = strainAB))

# ============================================================================
# 5. Linkage of Segments: Heterologous/Homologous Pairwise Genotypes ----
#    Part A only (per-cross, generalizes to n=1). Part B (cross-level ANOVA
#    and reassortment-rate correlation) is deferred until multiple plates
#    exist -- see note at the end of this block.
#
#    NOTE: the original script's first approach here (a linkage_disequilibrium()
#    function computing pairwise LD from marginal proportions) was statistically
#    incorrect and was abandoned by the original author partway through in favor
#    of the direct pairwise-genotype approach below. It has been dropped rather
#    than ported.
# ============================================================================

loci <- segment_order  # c('PB2','PB1','PA','HA','NP','NA','M','NS'), from Section 3

loci_pairings <- as.data.frame(combinations(n = 8, r = 2, v = loci))
colnames(loci_pairings) <- c("locusA", "locusB")
pairwise_locus_cols <- paste(loci_pairings$locusA, loci_pairings$locusB, sep = "_")

#### Build the per-sample, per-locus-pair genotype table ----
genotypes_table <- group_by(controls_df, cross, strainA, strainB, strainAB, cross_sample, locus) %>%
  summarise(majority_strain = majority_strain)

genotypes_table <- pivot_wider(genotypes_table, names_from = locus, values_from = majority_strain)

# Full 8-segment genotype code (also used to identify complete genotypes)
genotypes_table <- unite(genotypes_table, "genotype", all_of(loci), remove = FALSE)

# Restrict to complete (8-segment) genotypes -- incomplete genotypes can't
# contribute a valid pairwise call for any pair involving a missing locus.
genotypes_table_two_locus <- na.omit(genotypes_table)
genotypes_table_two_locus <- select(genotypes_table_two_locus, -genotype)

for (i in 1:nrow(loci_pairings)) {
  a <- as.vector(loci_pairings$locusA[i])
  b <- as.vector(loci_pairings$locusB[i])
  ab <- paste(a, b, sep = "_")
  genotypes_table_two_locus <- unite(genotypes_table_two_locus, !!ab, c(a, b), remove = FALSE)
}

# Drop the original single-locus columns, now that all 28 pairwise columns exist
genotypes_table_two_locus <- select(genotypes_table_two_locus, -all_of(loci))

#### Long format: one row per sample x locus-pair, then collapse to frequencies ----
genotypes_table_two_locus_long <- pivot_longer(
  genotypes_table_two_locus,
  cols = all_of(pairwise_locus_cols),
  names_to = "multilocus", values_to = "genotype"
) %>%
  group_by(cross, strainA, strainB, strainAB, multilocus, genotype) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))

# Classify each pairwise genotype as homologous (same parent at both loci) or heterologous
genotypes_table_two_locus_long <- separate(genotypes_table_two_locus_long, genotype,
                                           c("genotype_locusA", "genotype_locusB"), remove = FALSE)

homologous <- ifelse(
  genotypes_table_two_locus_long$genotype_locusA == genotypes_table_two_locus_long$genotype_locusB,
  "Homologous", "Heterologous"
)
genotypes_table_two_locus_long <- data.frame(genotypes_table_two_locus_long, homologous)

#### Plots ----
# Per locus-pair, per cross
ggplot(subset(genotypes_table_two_locus_long, multilocus == "HA_NA"), aes(x = genotype, y = freq)) +
  geom_col() + facet_wrap(~ cross, scales = "free_x")

# All locus-pairs within one cross
ggplot(genotypes_table_two_locus_long, aes(x = genotype, y = freq, fill = homologous)) +
  geom_col() + scale_fill_manual(values = c("orange", "darkblue")) +
  facet_wrap(~ multilocus, scales = "free_x")

# Summary bar of homologous vs. heterologous counts, per cross
ggplot(genotypes_table_two_locus_long, aes(x = homologous, y = count, fill = homologous)) +
  geom_col() + theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("orange", "darkblue")) +
  ylab("Frequency of progeny plaque isolates") +
  xlab("Strain combinations of pairwise segment combinations") +
  facet_wrap(~ strainAB, scales = "free_y", ncol = 5)

#### Overall homologous vs. heterologous frequency (pooled across locus pairs) ----
overall_homologous_proportions <- genotypes_table_two_locus_long %>% group_by(homologous) %>%
  summarise(freq = freq)
aggregate(freq ~ homologous, data = overall_homologous_proportions, mean, na.rm = TRUE)

#### Per-cross heterologous frequency, across the 28 locus pairs ----
genotypes_table_two_locus_long %>% filter(homologous == "Heterologous") %>% group_by(cross, strainAB) %>%
  summarise(
    mean_frequency = mean(freq),
    max_frequency = max(freq),
    sd_frequency = sd(freq)
  )

#### Per-cross, per-locus-pair heterologous frequency ----
multilocus_average_heterologous <- genotypes_table_two_locus_long %>% filter(homologous == "Heterologous") %>%
  group_by(cross, strainAB, multilocus) %>%
  summarise(mean_frequency = mean(freq))

#### Per-cross overall heterologous frequency (folded into cross_stats for later use) ----
cross_average_heterologous <- genotypes_table_two_locus_long %>% filter(homologous == "Heterologous") %>%
  group_by(cross, strainAB) %>%
  summarise(pairwise_heterologous_mean_frequency = mean(freq))

cross_stats <- left_join(cross_stats, cross_average_heterologous, by = c("cross", "strainAB"))
cross_stats$pairwise_heterologous_mean_frequency[is.na(cross_stats$pairwise_heterologous_mean_frequency)] <- 0

# ----------------------------------------------------------------------------
# DEFERRED (Part B -- requires multiple crosses):
#
# multilocus_average_heterologous %>% facet by same_subtype, and
#   model9 <- lm(mean_frequency ~ strainAB, multilocus_average_heterologous)
#   -- "same_subtype" doesn't exist in your metadata, and strainAB as a
#   categorical predictor needs multiple crosses to have >1 level.
#
# ggplot(cross_stats, aes(x = prop_reassortant, y = pairwise_heterologous_mean_frequency)) + ...
#   reassortment_pairwise_model <- lm(pairwise_heterologous_mean_frequency ~ prop_reassortant, cross_stats)
#   -- commented out rather than left in place: with n=1 cross this is a
#   regression on a single point, which is rank-deficient and will either
#   error or return a meaningless fit rather than just an uninteresting plot.
# ----------------------------------------------------------------------------