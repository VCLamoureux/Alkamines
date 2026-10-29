library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(patchwork)

project_dir <- "/Users/vincent/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma"
output_dir <- file.path(project_dir, "julius_standalone_results")

input_files <- c(
  feature_table = file.path(project_dir,
    "RA_dataset_15082024_all_lib-1c9129f0a76c4b3a9d58806d4729ef03-featuretable_reformated.csv"),
  annotation = file.path(project_dir,
    "1c9129f0a76c4b3a9d58806d4729ef03-merged_results_with_gnps.tsv"),
  clinical_metadata = file.path(project_dir,
    "clinical_data_20patients_updated_may20_final_response_plus_diet_score.csv"),
  drug_metadata = file.path(project_dir, "metadata_new_drug_t1.csv"))

if (any(!file.exists(input_files))) {
  stop("Missing source files:\n", paste(input_files[!file.exists(input_files)], collapse = "\n"),
       "\nUpdate project_dir or input_files at the top of this script.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_width <- 6.8
figure_height <- 4.25
figure_dpi <- 600
detection_threshold <- 1
visits_ord <- c("T1", "T2", "T3")

feature_ids <- c("5ASA_pyruvate", "5ASA_ribose", "5ASA_deoxyribose",
                 "5ASA_2hydropropanal_1", "sf_AKG", "sf_glyceraldehyde")

compound_map <- tibble::tribble(
  ~Feature,              ~Compound,                           ~Panel, ~Exponent,
  "5ASA_pyruvate",       "5-ASA_pyruvate",                     "f",    5,
  "5ASA_2hydropropanal_1", "5-ASA\n2-hydroxypropanal",           "g",    6,
  "5ASA_ribose",         "5-ASA_ribose",                       "h",    5,
  "sf_glyceraldehyde",   "Sulfapyridine_\nglyceraldehyde",      "i",    5)

group_palette <- c("non-user" = "#A6A9A9", "user" = "#BD5060")
axis_breaks <- list(c(0, 0.5, 1, 1.5, 2), c(0, 0.4, 0.8, 1.2),
                    c(0, 2, 4, 6), c(0, 1, 2, 3))

name_map <- c(
  "1688" = "5ASA_ribose",
  "3037" = "5ASA_pyruvate",
  "1710" = "5ASA_deoxyribose",
  "2518" = "5ASA_2hydropropanal_1",
  "826" = "sf_AKG",
  "4410" = "sf_glyceraldehyde")

target_feature_ids <- names(name_map)
sulfasalazine_user_ids <- c("FD003", "FD028")
specified_features <- c("5552", "2324", "4471", "2394", "2872", "7140")
sample_pattern <- "PD|almond|black|bread|chia|FD|flaxseeds|ginger|miso|oats|sesame|tahini|turmeric|vinegar"

feature_table <- read_csv(input_files["feature_table"],
                          col_types = cols(.default = col_character()),
                          na = c("", "NA", "N/A"), name_repair = "minimal")
names(feature_table) <- sub(" Peak area$", "", trimws(names(feature_table)))

if (anyDuplicated(names(feature_table))) stop("Duplicate feature-table column names.")
if (!"row ID" %in% names(feature_table)) stop("The feature table needs a row ID column.")

feature_table <- feature_table %>%
  dplyr::rename(Feature = `row ID`) %>%
  dplyr::mutate(Feature = trimws(Feature))

if (anyNA(feature_table$Feature) || any(!nzchar(feature_table$Feature)) ||
    anyDuplicated(feature_table$Feature)) {
  stop("Feature IDs must be present and unique in the feature table.")
}

area_columns <- names(feature_table)[grepl("\\.mzXML$", names(feature_table),
                                          ignore.case = TRUE)]
if (!length(area_columns)) stop("No mzXML peak-area columns were found.")

valid_area_columns <- vapply(feature_table[area_columns], function(x) {
  y <- suppressWarnings(as.numeric(x))
  all(is.na(x) | (!is.na(y) & is.finite(y) & y >= 0))
}, logical(1))

if (any(!valid_area_columns)) {
  stop("Invalid peak areas in: ", paste(area_columns[!valid_area_columns], collapse = ", "))
}

feature_table <- feature_table %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(area_columns), as.numeric))

qc_samples <- list(
  sixmix = area_columns[grepl("6Mix", area_columns)],
  blank = area_columns[grepl("Blank", area_columns)],
  sample = area_columns[grepl(sample_pattern, area_columns)])

if (any(lengths(qc_samples) == 0)) {
  stop("Missing columns for the original QC groups: ",
       paste(names(qc_samples)[lengths(qc_samples) == 0], collapse = ", "),
       ". Select the original feature table containing the QC samples.")
}

qc_summary <- function(sample_columns, prefix) {
  m <- feature_table %>%
    dplyr::select(dplyr::all_of(sample_columns)) %>%
    as.matrix()
  tibble(Feature = feature_table$Feature,
         Mean = rowMeans(m), SD = apply(m, 1, sd)) %>%
    dplyr::mutate(CV = SD / Mean) %>%
    dplyr::filter(Mean > 0) %>%
    dplyr::rename_with(~ paste0(.x, "_", prefix), -Feature)
}

qc_summaries <- purrr::imap(qc_samples, qc_summary)
sixmix_feature_info <- qc_summaries$sixmix
blank_feature_info <- qc_summaries$blank
sample_feature_info <- qc_summaries$sample

feature_to_remove_qcmix <- sixmix_feature_info %>%
  dplyr::left_join(sample_feature_info, by = "Feature") %>%
  dplyr::mutate(sample_Mix = Mean_sample / Mean_sixmix) %>%
  dplyr::filter(sample_Mix < 5 | is.na(sample_Mix)) %>%
  dplyr::transmute(Feature, reason = "sample_to_6Mix_below_5_or_missing") %>%
  dplyr::bind_rows(blank_feature_info %>%
                    dplyr::filter(Feature %in% specified_features) %>%
                    dplyr::transmute(Feature, reason = "listed_standard_present_in_blank")) %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(reason = paste(unique(reason), collapse = ";"), .groups = "drop")

write_csv(feature_to_remove_qcmix, file.path(output_dir, "qc_removed_features.csv"))
purrr::iwalk(qc_summaries, function(x, name) {
  write_csv(x, file.path(output_dir, paste0("qc_", name, "_summary.csv")))
})

if (!all(target_feature_ids %in% feature_table$Feature)) {
  stop("Target feature IDs absent from the raw table: ",
       paste(setdiff(target_feature_ids, feature_table$Feature), collapse = ", "))
}
if (any(target_feature_ids %in% feature_to_remove_qcmix$Feature)) {
  stop("Target features failed the original QC filter: ",
       paste(intersect(target_feature_ids, feature_to_remove_qcmix$Feature), collapse = ", "),
       ". Review qc_removed_features.csv before analysing them.")
}

sample_columns <- area_columns[!grepl("6Mix|Blank", area_columns)]
visit_pattern <- paste0("_(", paste(visits_ord, collapse = "|"), ")\\.mzXML$")
fecal_columns <- sample_columns[grepl("FD", sample_columns) &
                                 grepl(visit_pattern, sample_columns, ignore.case = TRUE)]
if (!length(fecal_columns)) stop("No fecal participant-visit columns were found.")

df_melted_feces <- feature_table %>%
  dplyr::filter(!Feature %in% feature_to_remove_qcmix$Feature,
                Feature %in% target_feature_ids) %>%
  dplyr::select(Feature, dplyr::all_of(fecal_columns)) %>%
  tidyr::pivot_longer(-Feature, names_to = "Sample_Name", values_to = "Peak_Area") %>%
  dplyr::mutate(Sample_Name = sub("\\.(mzXML|mzML)$", "", Sample_Name, ignore.case = TRUE))

if (any(!grepl("^FD[0-9]+_T[0-9]+$", df_melted_feces$Sample_Name))) {
  stop("Fecal sample names must follow FD028_T3 after removing the file extension.")
}
if (anyDuplicated(df_melted_feces[c("Feature", "Sample_Name")])) {
  stop("Duplicate raw feature-sample measurements after cleaning sample names.")
}

annotation <- read_tsv(input_files["annotation"],
                        col_types = cols(`#Scan#` = col_character()),
                        na = c("", "NA", "N/A"))
if (!all(c("#Scan#", "Compound_Name") %in% names(annotation))) {
  stop("The annotation table needs #Scan# and Compound_Name columns.")
}

annotation <- annotation %>%
  dplyr::rename(Feature = `#Scan#`) %>%
  dplyr::mutate(Feature = trimws(Feature)) %>%
  dplyr::filter(Feature %in% target_feature_ids) %>%
  dplyr::distinct()

if (anyDuplicated(annotation$Feature)) {
  stop("Multiple different annotation rows for a selected feature. Resolve these",
       " annotations before joining; duplicated peak areas must not be summed.")
}

clinical_metadata <- read_csv(input_files["clinical_metadata"],
                               col_types = cols(Patient = col_character(), DMARD = col_character()),
                               na = c("", "NA", "N/A"))
new_metadata <- read_csv(input_files["drug_metadata"],
                          col_types = cols(filename2 = col_character(), DMARD = col_character()),
                          na = c("", "NA", "N/A"))

if (!all(c("Patient", "DMARD", "Pain_50_improvement") %in% names(clinical_metadata))) {
  stop("Clinical metadata needs Patient, DMARD and Pain_50_improvement columns.")
}
if (!all(c("filename2", "DMARD") %in% names(new_metadata))) {
  stop("Drug metadata needs filename2 and DMARD columns.")
}
if (anyNA(clinical_metadata$Patient) || any(!nzchar(trimws(clinical_metadata$Patient)))) {
  stop("Clinical metadata contains a missing Patient identifier.")
}

clinical_metadata_renamed_2 <- clinical_metadata %>%
  dplyr::mutate(
    Patient = paste0("F", trimws(Patient)),
    Patient = sub("^FD010d-0$", "FD010d0", Patient),
    Patient = sub("^FD001D0$", "FD001d0", Patient),
    Patient = sub("d0$", "_T2", Patient),
    Patient = sub("d14$", "_T3", Patient)) %>%
  dplyr::rename(filename = Patient) %>%
  dplyr::distinct()

drug_updates <- new_metadata %>%
  dplyr::transmute(
    Sample_Name = sub("\\.(mzXML|mzML)$", "", trimws(filename2), ignore.case = TRUE),
    DMARD_new = DMARD) %>%
  dplyr::distinct()

if (anyNA(clinical_metadata_renamed_2$filename) ||
    anyDuplicated(clinical_metadata_renamed_2$filename)) {
  stop("Clinical metadata must have one distinct record per renamed sample.")
}
if (anyNA(drug_updates$Sample_Name) || anyDuplicated(drug_updates$Sample_Name)) {
  stop("Drug metadata must have one distinct DMARD update per sample.")
}

table_julius <- df_melted_feces %>%
  dplyr::left_join(clinical_metadata_renamed_2, by = c("Sample_Name" = "filename")) %>%
  dplyr::left_join(annotation, by = "Feature") %>%
  dplyr::left_join(drug_updates, by = "Sample_Name") %>%
  dplyr::mutate(
    DMARD = dplyr::coalesce(DMARD_new, DMARD),
    Participant_ID = sub("_T[0-9]+$", "", Sample_Name),
    sulfasalazine = dplyr::if_else(Participant_ID %in% sulfasalazine_user_ids, "user", "non-user"),
    Feature_ID = Feature,
    Feature = unname(name_map[Feature]),
    Feature_Compound = paste(Feature, Compound_Name, sep = "_")) %>%
  dplyr::select(-DMARD_new) %>%
  dplyr::select(Feature_Compound, Feature, Feature_ID, Sample_Name, Peak_Area,
                Pain_50_improvement, DMARD, dplyr::everything())

if (nrow(table_julius) != nrow(df_melted_feces) ||
    anyDuplicated(table_julius[c("Feature", "Sample_Name")])) {
  stop("A metadata join duplicated measurements. Resolve the join before analysis.")
}

metadata_join_audit <- table_julius %>%
  dplyr::distinct(Sample_Name, Participant_ID, sulfasalazine, DMARD) %>%
  dplyr::mutate(
    clinical_match = Sample_Name %in% clinical_metadata_renamed_2$filename,
    drug_metadata_match = Sample_Name %in% drug_updates$Sample_Name)

write_csv(table_julius, file.path(output_dir, "table_julius.csv"))
write_csv(metadata_join_audit, file.path(output_dir, "metadata_join_audit.csv"))

df_raw <- table_julius

required_columns <- c("Feature", "Sample_Name", "Peak_Area", "sulfasalazine")
if (!all(required_columns %in% names(df_raw))) {
  stop("Missing columns: ", paste(setdiff(required_columns, names(df_raw)), collapse = ", "))
}

df_long <- df_raw %>%
  dplyr::mutate(Feature = trimws(as.character(Feature))) %>%
  dplyr::filter(Feature %in% feature_ids) %>%
  dplyr::select(dplyr::all_of(required_columns)) %>%
  dplyr::mutate(
    Sample_Name = sub("\\.(mzXML|mzML)$", "", trimws(as.character(Sample_Name)),
                      ignore.case = TRUE),
    Peak_Area_original = Peak_Area,
    Peak_Area = suppressWarnings(as.numeric(Peak_Area)),
    sulfasalazine = trimws(as.character(sulfasalazine)))

if (!all(feature_ids %in% df_long$Feature)) stop("Some requested features are absent.")
if (anyNA(df_long$Sample_Name) ||
    any(!grepl("^FD[0-9]+_T[0-9]+$", df_long$Sample_Name))) {
  stop("Sample names must follow FD028_T3.")
}
if (anyDuplicated(df_long[c("Feature", "Sample_Name")])) {
  stop("Duplicate feature-sample rows: check annotation and metadata joins.")
}
if (any(!is.na(df_long$Peak_Area_original) & is.na(df_long$Peak_Area)) ||
    any(!is.na(df_long$Peak_Area) &
        (!is.finite(df_long$Peak_Area) | df_long$Peak_Area < 0))) {
  stop("Peak_Area must contain nonnegative numeric values or missing values.")
}
if (anyNA(df_long$sulfasalazine) ||
    any(!df_long$sulfasalazine %in% c("non-user", "user"))) {
  stop("Sulfasalazine status must be user or non-user.")
}

df_long <- df_long %>%
  dplyr::mutate(
    Participant_ID = sub("_T[0-9]+$", "", Sample_Name),
    Timepoint = sub("^.*_(T[0-9]+)$", "\\1", Sample_Name)) %>%
  dplyr::select(-Peak_Area_original)

if (any(!df_long$Timepoint %in% visits_ord)) stop("Unexpected visit: use T1, T2 and T3.")

exposure_check <- df_long %>%
  dplyr::distinct(Participant_ID, sulfasalazine) %>%
  dplyr::count(Participant_ID)

if (any(exposure_check$n != 1)) stop("Exposure changes within a participant.")

df_long <- df_long %>%
  dplyr::mutate(
    Timepoint = factor(Timepoint, levels = visits_ord),
    sulfasalazine = factor(sulfasalazine, levels = c("non-user", "user"))) %>%
  dplyr::arrange(Feature, sulfasalazine, Participant_ID, Timepoint)

sample_design <- df_long %>%
  dplyr::distinct(Sample_Name, Participant_ID, Timepoint, sulfasalazine)

group_counts <- sample_design %>%
  dplyr::distinct(Participant_ID, sulfasalazine) %>%
  dplyr::count(sulfasalazine, .drop = FALSE)

participant_means <- df_long %>%
  dplyr::group_by(Feature, Participant_ID, sulfasalazine) %>%
  dplyr::summarise(
    n_visits = sum(!is.na(Peak_Area)),
    complete_profile = setequal(as.character(Timepoint[!is.na(Peak_Area)]), visits_ord),
    missing_visits = paste(setdiff(visits_ord, as.character(Timepoint[!is.na(Peak_Area)])),
                            collapse = ";"),
    mean_peak_area = if (complete_profile) mean(Peak_Area) else NA_real_,
    ever_detected = if (complete_profile) any(Peak_Area > detection_threshold) else NA,
    .groups = "drop") %>%
  tidyr::complete(Feature = feature_ids, nesting(Participant_ID, sulfasalazine),
                  fill = list(n_visits = 0L, complete_profile = FALSE,
                              missing_visits = paste(visits_ord, collapse = ";")))

write_csv(participant_means, file.path(output_dir, "participant_summaries.csv"))

if (any(!participant_means$complete_profile)) {
  warning("Incomplete participant-feature profiles are excluded from the tests.")
}

stats_longitudinal <- purrr::map_dfr(feature_ids, function(feat) {
  d <- participant_means %>%
    dplyr::filter(Feature == feat, complete_profile)
  y <- d$mean_peak_area
  users <- d$sulfasalazine == "user"
  n <- length(y)
  k <- sum(users)
  if (k < 2 || n - k < 2) {
    stop("At least two complete participants per group are needed for ", feat, ".")
  }
  n_permutations <- choose(n, k)
  if (!is.finite(n_permutations) || n_permutations > 1e6) {
    stop("Too many exact permutations for ", feat, ".")
  }
  observed <- mean(y[users]) - mean(y[!users])
  permuted <- combn(seq_len(n), k, FUN = function(i) mean(y[i]) - mean(y[-i]))
  tolerance <- 1e-12 * max(1, abs(observed), abs(permuted))
  tibble(
    Feature = feat,
    n_users = k,
    n_non_users = n - k,
    mean_user = mean(y[users]),
    mean_non_user = mean(y[!users]),
    difference = observed,
    p_value = mean(abs(permuted) >= abs(observed) - tolerance),
    n_permutations = n_permutations)
})

stats_detection <- purrr::map_dfr(feature_ids, function(feat) {
  d <- participant_means %>%
    dplyr::filter(Feature == feat, complete_profile)
  tab <- table(factor(d$sulfasalazine, levels = c("non-user", "user")),
                factor(d$ever_detected, levels = c(FALSE, TRUE)))
  p <- if (any(colSums(tab) == 0)) 1 else fisher.test(tab)$p.value
  tibble(
    Feature = feat,
    detection_threshold = detection_threshold,
    users_detected = unname(tab["user", "TRUE"]),
    n_users = unname(sum(tab["user", ])),
    non_users_detected = unname(tab["non-user", "TRUE"]),
    n_non_users = unname(sum(tab["non-user", ])),
    p_value = p)
})

plot_stats <- compound_map %>%
  dplyr::left_join(stats_longitudinal, by = "Feature") %>%
  dplyr::mutate(
    p_label = paste0("P = ", sub("e([+-])0+", "e\\1",
                                 formatC(p_value, format = "e", digits = 2))),
    stars = dplyr::case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE            ~ "ns"))

plot_data <- df_long %>%
  dplyr::filter(Feature %in% compound_map$Feature) %>%
  dplyr::group_by(Feature, Participant_ID, sulfasalazine) %>%
  tidyr::complete(Timepoint = factor(visits_ord, levels = visits_ord)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(compound_map, by = "Feature") %>%
  dplyr::arrange(sulfasalazine, Participant_ID, Timepoint) %>%
  dplyr::mutate(
    Participant_ID = factor(Participant_ID, levels = unique(Participant_ID)),
    Intensity = Peak_Area / 10^Exponent)

legend_labels <- setNames(paste0(group_counts$sulfasalazine, " (n = ", group_counts$n, ")"),
                          group_counts$sulfasalazine)

paper_theme <- theme_classic(base_size = 10.5, base_family = "Helvetica") +
    theme(
      axis.text = element_text(color = "black", size = 10),
      axis.title.x = element_text(size = 11, margin = margin(t = 5)),
      axis.title.y = element_text(face = "bold", size = 10.5, margin = margin(r = 3)),
      axis.line = element_line(color = "#555555", linewidth = 0.35),
      axis.ticks = element_line(color = "#555555", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.8, "pt"),
      plot.title = element_text(size = 10.5, hjust = 0.5, lineheight = 1.05,
                                margin = margin(t = 6, b = 5)),
      plot.subtitle = element_text(size = 8.6, hjust = 0.5, lineheight = 1.2,
                                   margin = margin(b = 5)),
      plot.tag = element_text(face = "bold", size = 15),
      plot.tag.position = c(0, 1),
      plot.margin = margin(t = 3, r = 7, b = 3, l = 2),
      legend.title = element_text(size = 9.5),
      legend.text = element_text(size = 9.5),
      legend.key.width = grid::unit(16, "pt"),
      plot.background = element_rect(fill = NA, color = NA))

p_compounds <- purrr::map(seq_len(nrow(compound_map)), function(i) {
  feat <- compound_map$Feature[i]
  d <- plot_data %>% dplyr::filter(Feature == feat)
  s <- plot_stats %>% dplyr::filter(Feature == feat)
  upper <- max(d$Intensity, na.rm = TRUE)
  if (upper <= 0) upper <- 1
  breaks <- axis_breaks[[i]]
  if (upper > max(breaks) * 1.5 || upper < max(breaks) * 0.5) {
    breaks <- pretty(c(0, upper), n = 4)
    breaks <- breaks[breaks >= 0 & breaks <= upper * 1.07]
  }
  ggplot(d, aes(x = Timepoint, y = Intensity, group = Participant_ID,
                color = sulfasalazine)) +
    geom_line(linewidth = 0.55, alpha = 0.8, na.rm = TRUE) +
    geom_point(aes(fill = sulfasalazine), shape = 21, color = "#555555",
                size = 1.9, stroke = 0.3, alpha = 0.95, na.rm = TRUE) +
    scale_color_manual(values = group_palette, labels = legend_labels, drop = FALSE) +
    scale_fill_manual(values = group_palette, labels = legend_labels, drop = FALSE) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = 0.2)) +
    scale_y_continuous(breaks = breaks, limits = c(-upper * 0.04, upper * 1.07),
                       expand = expansion(mult = 0),
                       labels = function(x) format(x, trim = TRUE, scientific = FALSE)) +
    labs(x = "Time point",
         y = bquote(bold("Intensity (" * 10^.(s$Exponent) * ")")),
         title = s$Compound, subtitle = paste(s$stars, s$p_label, sep = "\n"),
         tag = s$Panel, color = "Sulfasalazine", fill = "Sulfasalazine") +
    paper_theme
})

names(p_compounds) <- compound_map$Feature

p_longitudinal <- wrap_plots(p_compounds, nrow = 1, guides = "collect",
                             axis_titles = "collect_x") +
  plot_annotation(
    caption = "Unadjusted P: between-group tests of participant means across T1-T3.",
    theme = theme(plot.caption = element_text(size = 7.5, hjust = 0.5,
                                               family = "Helvetica", color = "#444444"))) &
  theme(legend.position = "bottom")

boxplot_data <- participant_means %>%
  dplyr::filter(complete_profile, Feature %in% compound_map$Feature) %>%
  dplyr::left_join(compound_map, by = "Feature") %>%
  dplyr::mutate(Intensity = mean_peak_area / 10^Exponent)

p_boxplots <- purrr::map(seq_len(nrow(compound_map)), function(i) {
  feat <- compound_map$Feature[i]
  d <- boxplot_data %>% dplyr::filter(Feature == feat)
  s <- plot_stats %>% dplyr::filter(Feature == feat)
  upper <- max(d$Intensity, na.rm = TRUE)
  if (upper <= 0) upper <- 1
  breaks <- pretty(c(0, upper), n = 4)
  breaks <- breaks[breaks >= 0 & breaks <= upper * 1.07]
  x_labels <- c("non-user" = paste0("non-user\n(n = ", s$n_non_users, ")"),
                 "user" = paste0("user\n(n = ", s$n_users, ")"))
  ggplot(d, aes(x = sulfasalazine, y = Intensity, fill = sulfasalazine)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.35,
                  color = "#777777", linewidth = 0.35) +
    geom_point(shape = 21, color = "#555555", size = 2.2, stroke = 0.3,
                position = position_jitter(width = 0.09, height = 0, seed = 28 + i)) +
    scale_fill_manual(values = group_palette, guide = "none") +
    scale_x_discrete(labels = x_labels, expand = expansion(add = 0.45)) +
    scale_y_continuous(breaks = breaks, limits = c(-upper * 0.04, upper * 1.07),
                       expand = expansion(mult = 0),
                       labels = function(x) format(x, trim = TRUE, scientific = FALSE)) +
    labs(x = "Sulfasalazine",
         y = bquote(bold("Mean intensity (" * 10^.(s$Exponent) * ")")),
         title = s$Compound, subtitle = paste(s$stars, s$p_label, sep = "\n"),
         tag = s$Panel) +
    paper_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
          legend.position = "none")
})

names(p_boxplots) <- compound_map$Feature

p_boxplot <- wrap_plots(p_boxplots, nrow = 1, axis_titles = "collect_x") +
  plot_annotation(
    caption = "Each dot is one participant's T1-T3 mean. P values are unadjusted.",
    theme = theme(plot.caption = element_text(size = 7.5, hjust = 0.5,
                                               family = "Helvetica", color = "#444444")))

write_csv(boxplot_data, file.path(output_dir, "participant_boxplot_data.csv"))

writeLines(c(
  "f-i, Participant mean raw metabolite peak areas by sulfasalazine use.",
  "Each point is one participant's arithmetic mean across T1, T2 and T3.",
  "Group labels give the number of complete participant profiles used in each panel.",
  "Boxes show medians and interquartile ranges of participant means; whiskers extend",
  "to the most extreme observations within 1.5 times the interquartile range.",
  "All individual points are shown. For the two-person user group, the quartiles",
  "are interpolated from those two values and the box is a descriptive summary only.",
  "P values are unadjusted, two-sided exact participant-label permutation tests of",
  "the difference in group means, counting absolute differences at least as large as observed.",
  "The test is the same as that displayed on the accompanying trajectory figure.",
  "Stars: * P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001; ns P >= 0.05.",
  "Axis units are labelled powers of ten; raw areas and linear axes are used.",
  "Participant averages summarize overall abundance and do not show temporal patterns."
), file.path(output_dir, "boxplot_caption.txt"))


write_csv(df_long, file.path(output_dir, "analysis_data.csv"))
write_csv(sample_design, file.path(output_dir, "sample_design.csv"))
write_csv(stats_longitudinal, file.path(output_dir, "abundance_results.csv"))
write_csv(stats_detection, file.path(output_dir, "detection_results.csv"))
write_csv(plot_stats, file.path(output_dir, "selected_plot_statistics.csv"))

writeLines(c(
  "f-i, Individual longitudinal profiles of raw metabolite peak areas.",
  "Lines connect visits from each participant. Gray: non-users; red: sulfasalazine users.",
  "T1-T3 are categorical visits; spacing does not indicate elapsed days.",
  "P values are unadjusted, two-sided exact participant-label permutation tests",
  "of the difference in group means after averaging T1-T3 within each participant.",
  "Only complete profiles enter the tests. The null assumes exchangeable participant summaries.",
  "Stars: * P < 0.05; ** P < 0.01; *** P < 0.001; **** P < 0.0001; ns P >= 0.05.",
  "These are overall between-group abundance comparisons, not tests of change over time.",
  "Axes are linear and use the labelled display units. Raw peak areas are used in the tests.",
  paste0("Participants: ", paste(legend_labels, collapse = "; "), "."),
  paste0("Exact allocations per feature: ", paste(unique(plot_stats$n_permutations), collapse = ", "), "."),
  "Interpret group associations as exploratory given the small exposed group.",
  "Secondary detection uses one ever-detected status per participant; the threshold is not a validated LOD."
), file.path(output_dir, "figure_caption.txt"))

pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

plots_to_save <- list(
  selected_compounds_longitudinal = p_longitudinal,
  selected_compounds_participant_boxplot = p_boxplot)

purrr::iwalk(plots_to_save, function(p, name) {
  ggsave(paste0(name, ".pdf"), p, path = output_dir,
          device = pdf_device, width = figure_width, height = figure_height,
          family = "Helvetica", bg = "transparent")
  ggsave(paste0(name, ".svg"), p, path = output_dir,
          device = grDevices::svg, width = figure_width, height = figure_height,
          family = "Helvetica", bg = "transparent")
  ggsave(paste0(name, ".png"), p, path = output_dir,
          device = grDevices::png, width = figure_width, height = figure_height,
          dpi = figure_dpi, family = "Helvetica", bg = "white")
})

print(stats_longitudinal)
print(p_longitudinal)
print(p_boxplot)
