#Code by Claude
library(readr)
library(dplyr)

# ============================================================================
# Import
# ============================================================================
# No header in the file, so we name columns explicitly.
reads <- read_csv(
  "outputs/read_counts_demultiplexed_CA09xPAN99.txt",
  col_names = c("n_reads", "plate_well", "filename"),
  col_types = cols(
    n_reads    = col_integer(),
    plate_well = col_character(),
    filename   = col_skip()   # ignoring filename as requested
  )
)

# ============================================================================
# Parse plate/well
# ============================================================================
# plate_well looks like "plate02_well01" -> split into plate and well
reads <- reads %>%
  mutate(
    plate = sub("_well.*$", "", plate_well),
    well  = sub("^plate[0-9]+_", "", plate_well)
  )

# ============================================================================
# Reads per plate (pooling wells)
# ============================================================================
reads_per_plate <- reads %>%
  group_by(plate) %>%
  summarise(total_reads = sum(n_reads), n_wells = n(), .groups = "drop") %>%
  arrange(plate)

print(reads_per_plate)

# ============================================================================
# Plate02-specific: mean and SD per well
# ============================================================================
plate02 <- reads %>% filter(plate == "plate02")

plate02_stats <- plate02 %>%
  summarise(
    n_wells    = n(),
    mean_reads = mean(n_reads),
    sd_reads   = sd(n_reads),
    min_reads  = min(n_reads),
    max_reads  = max(n_reads)
  )

print(plate02_stats)

# ============================================================================
# Missing well check (plate02)
# ============================================================================
# Assumes a 96-well plate. Change EXPECTED_WELLS if your plates differ (e.g., 384).
EXPECTED_WELLS <- 1:96

# Extract the numeric well index from strings like "well01", "well1", etc.
plate02 <- plate02 %>%
  mutate(well_num = as.integer(sub("^well", "", well)))

observed_wells  <- sort(unique(plate02$well_num))
missing_wells   <- setdiff(EXPECTED_WELLS, observed_wells)
duplicate_wells <- plate02$well_num[duplicated(plate02$well_num)]

if (length(missing_wells) > 0) {
  message(sprintf(
    "plate02: %d well(s) missing out of %d expected: %s",
    length(missing_wells), length(EXPECTED_WELLS),
    paste(sprintf("well%02d", missing_wells), collapse = ", ")
  ))
} else {
  message("plate02: all expected wells present.")
}

if (length(duplicate_wells) > 0) {
  message(sprintf(
    "plate02: WARNING — duplicate well entries found: %s",
    paste(sprintf("well%02d", duplicate_wells), collapse = ", ")
  ))
}

# Optional: flag wells present but with suspiciously low counts (e.g., <1% of plate mean)
low_threshold <- 0.01 * plate02_stats$mean_reads
low_wells <- plate02 %>% filter(n_reads < low_threshold)

if (nrow(low_wells) > 0) {
  message(sprintf(
    "plate02: %d well(s) with reads < 1%% of plate mean (possible failed wells): %s",
    nrow(low_wells), paste(low_wells$well, collapse = ", ")
  ))
}