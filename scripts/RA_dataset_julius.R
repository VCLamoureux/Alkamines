setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# import packages
library(tidyverse)
library(ggrepel)

# import tables
feature_table <- read_csv("/Users/vincentlamoureux/Downloads/RA_dataset_15082024_all_lib-1c9129f0a76c4b3a9d58806d4729ef03-featuretable_reformated.csv") |>
  dplyr::rename_with(~gsub(" Peak area", "", .))

drift <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/MZmine_3_01292025/fbmn_quant_mzmine.csv")
colnames(drift)[1] <- "Feature"

metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/ReDU_Diet_Nov19_Template.csv") |>
  dplyr::select(-MassiveID)

new_metadata <- read_csv("/Users/vincentlamoureux/Downloads/metadata_new_drug_t1.csv")

# Clinical metadata
clinical_metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/clinical_data_20patients_updated_may20_final_response_plus_diet_score.csv")

clinical_metadata_renamed_2 <- clinical_metadata |>
  dplyr::mutate(Patient = paste0("F", Patient)) |>
  dplyr::mutate(
    Patient = gsub("FD010d-0", "FD010d0", Patient),
    Patient = gsub("FD001D0",  "FD001d0", Patient),
    Patient = str_replace(Patient, "d0",  "_T2"),
    Patient = str_replace(Patient, "d14", "_T3")) |>
  dplyr::rename(filename = Patient)

clinical_metadata_renamed_plasma_2 <- clinical_metadata |>
  dplyr::mutate(Patient = paste0("P", Patient)) |>
  dplyr::mutate(
    Patient = gsub("PD010d-0", "PD010d0", Patient),
    Patient = gsub("PD001D0",  "PD001d0", Patient),
    Patient = str_replace(Patient, "d0",  "_T2"),
    Patient = str_replace(Patient, "d14", "_T3")) |>
  dplyr::rename(filename = Patient)

# annotation table
annotation <- read_tsv("/Users/vincentlamoureux/Downloads/1c9129f0a76c4b3a9d58806d4729ef03-merged_results_with_gnps.tsv") |>
  dplyr::rename(Feature = `#Scan#`) |>
  dplyr::mutate(Feature = as.character(Feature))

# get info feature
info_feature <- feature_table[, 1:3] |>
  dplyr::rename(Feature = `row ID`) |>
  dplyr::mutate(Feature = as.character(Feature))
colnames(info_feature) <- c("Feature", "mz", "RT")

info_feature_name <- info_feature |>
  left_join(annotation, by = "Feature")

# transpose feature table
data_transpose <- feature_table |>
  column_to_rownames("row ID") |>
  dplyr::select(contains(".mzXML")) |>
  t() |>
  as.data.frame() |>
  rownames_to_column("filename")

# rt drift of internal standard (Feature 2394)
rt_data_cleaned <- drift |>
  dplyr::filter(Feature == 2394) |>
  dplyr::select(contains("Feature RT")) |>
  tidyr::gather(key = "filename", value = "RT") |>
  dplyr::filter(!str_detect(filename, "Blank")) |>
  dplyr::mutate(RT = as.numeric(RT))

rt_stats <- rt_data_cleaned |>
  summarise(
    Mean_RT     = mean(RT, na.rm = TRUE),
    SD_RT       = sd(RT, na.rm = TRUE),
    Variance_RT = var(RT, na.rm = TRUE))

ggplot(rt_data_cleaned, aes(x = filename, y = RT)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Filename", y = "Retention Time (RT)", title = "Retention Time Drift for Feature 2394")

# peak area drift
peakarea_data_cleaned <- drift |>
  dplyr::filter(Feature == 2394) |>
  dplyr::select(contains("Peak area")) |>
  tidyr::gather(key = "filename", value = "Peak") |>
  dplyr::filter(!str_detect(filename, "Blank"))

peakarea_stats <- peakarea_data_cleaned |>
  summarise(
    Mean_peakarea     = mean(Peak, na.rm = TRUE),
    SD_peakarea       = sd(Peak, na.rm = TRUE),
    Variance_peakarea = var(Peak, na.rm = TRUE))

ggplot(peakarea_data_cleaned, aes(x = filename, y = Peak)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Filename", y = "Peak Area", title = "Peak Area Drift for Feature 2394")

# filter for 6 mix
data_sixmix <- data_transpose |>
  dplyr::filter(str_detect(filename, "6Mix")) |>
  dplyr::mutate(filename = gsub(".mzXML", "", filename))

sixmix_feature_info <- data.frame(
  Feature      = colnames(data_sixmix)[-1],
  Mean_sixmix  = data_sixmix |> column_to_rownames("filename") |> colMeans(),
  SD_sixmix    = data_sixmix |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sixmix = SD_sixmix / Mean_sixmix) |>
  dplyr::filter(Mean_sixmix > 0) |>
  dplyr::arrange(desc(Mean_sixmix)) |>
  dplyr::mutate(Feature = as.character(Feature))

sixmix_feature_info_name <- sixmix_feature_info |>
  left_join(info_feature_name, by = "Feature") |>
  dplyr::select(Feature, Mean_sixmix, SD_sixmix, CV_sixmix, mz, RT, Compound_Name)

# filter for blanks
data_blank <- data_transpose |>
  dplyr::filter(str_detect(filename, "Blank"))

blank_feature_info <- data.frame(
  Feature    = colnames(data_blank)[-1],
  Mean_blank = data_blank |> column_to_rownames("filename") |> colMeans(),
  SD_blank   = data_blank |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_blank = SD_blank / Mean_blank) |>
  dplyr::filter(Mean_blank > 0) |>
  dplyr::arrange(desc(Mean_blank))

blank_feature_info_name <- blank_feature_info |>
  left_join(info_feature_name, by = "Feature")

# filter for samples
data_sample <- data_transpose |>
  dplyr::filter(str_detect(filename, "PD|almond|black|bread|chia|FD|flaxseeds|ginger|miso|oats|sesame|tahini|turmeric|vinegar"))

sample_feature_info <- data.frame(
  Feature     = colnames(data_sample)[-1],
  Mean_sample = data_sample |> column_to_rownames("filename") |> colMeans(),
  SD_sample   = data_sample |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample / Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |>
  dplyr::arrange(desc(Mean_sample))

sample_feature_info_name <- sample_feature_info |>
  left_join(info_feature_name, by = "Feature")

# identify features for the qc mix
specified_features <- c("5552", "2324", "4471", "2394", "2872", "7140")

feature_to_remove_qcmix <- sixmix_feature_info |>
  left_join(sample_feature_info, by = "Feature") |>
  dplyr::filter(Mean_sixmix > 0) |>
  dplyr::mutate(sample_Mix = Mean_sample / Mean_sixmix) |>
  dplyr::filter(sample_Mix < 5 | is.na(sample_Mix)) |>
  dplyr::bind_rows(blank_feature_info |> dplyr::filter(Feature %in% specified_features))

# data clean
data_clean2 <- data_transpose |>
  dplyr::select(-all_of(feature_to_remove_qcmix$Feature)) |>
  dplyr::filter(!str_detect(filename, "6Mix|Blank"))

data_clean_transpose <- data_clean2 |>
  column_to_rownames("filename") |>
  t() |>
  as.data.frame() |>
  rownames_to_column("Feature")

# Merge with annotation
merged <- data_clean_transpose |>
  left_join(annotation, by = "Feature") |>
  dplyr::mutate(Feature = as.character(Feature))

merged_mzml_mean_all <- merged |>
  dplyr::select(any_of(c("Compound_Name", "Feature", "LibraryName")), matches("_T1|_T2|_T3"))

# melt
mzml_columns <- grep("T1|T2|T3", names(merged_mzml_mean_all), value = TRUE)
df_melted <- reshape2::melt(merged_mzml_mean_all, id.vars = "Feature", measure.vars = mzml_columns, variable.name = "Sample_Name", value.name = "Peak_Area")

df_melted_feces <- df_melted |>
  dplyr::filter(str_detect(Sample_Name, "FD")) |>
  dplyr::mutate(Sample_Name = gsub(".mzXML", "", Sample_Name))

df_melted_feces_clinical <- df_melted_feces |>
  left_join(clinical_metadata_renamed_2, by = c("Sample_Name" = "filename"))

name_map <- c(
  "1688" = "5ASA_ribose",
  "3037" = "5ASA_pyruvate",
  "1710" = "5ASA_deoxyribose",
  "2518" = "5ASA_2hydropropanal_1",
  "826" = "sf_AKG", 
  "4410" = "sf_glyceraldehyde")

feature_ids <- names(name_map)

sulfasalazine_users <- c("FD003_T1", "FD028_T1", "FD003_T2", "FD028_T2", "FD003_T3", "FD028_T3")

df_melted_metadata_annotation <- df_melted_feces_clinical |>
  dplyr::mutate(Feature = trimws(as.character(Feature))) |>
  dplyr::filter(Feature %in% feature_ids) |>
  dplyr::select(Feature, Sample_Name, Peak_Area, Pain_50_improvement, DMARD, dplyr::everything()) |>
  dplyr::mutate(sulfasalazine  = dplyr::if_else(Sample_Name %in% sulfasalazine_users, "user", "non-user"), Peak_Area_Log = log10(Peak_Area + 1)) |>
  dplyr::left_join(annotation, by = "Feature") |>
  dplyr::left_join(new_metadata |> dplyr::select(filename2, DMARD) |> dplyr::rename(DMARD_new = DMARD), by = c("Sample_Name" = "filename2")) |>
  dplyr::mutate(DMARD = dplyr::coalesce(DMARD_new, DMARD)) |>
  dplyr::select(-DMARD_new) |>
  dplyr::mutate(Feature = dplyr::recode(as.character(Feature), !!!name_map, .default = as.character(Feature)), Feature_Compound = paste(Feature, Compound_Name, sep = "_")) |>
  dplyr::select(Feature_Compound, dplyr::everything())

df_sum <- df_melted_metadata_annotation |>
  dplyr::filter(Feature == "5ASA_2hydropropanal_1") |>
  dplyr::group_by(Sample_Name) |>
  dplyr::summarise(
    Peak_Area         = sum(Peak_Area, na.rm = TRUE),
    sulfasalazine     = dplyr::first(sulfasalazine),
    Pain_50_improvement = dplyr::first(Pain_50_improvement),
    DMARD             = dplyr::first(DMARD),
    .groups = "drop") |>
  dplyr::mutate(
    Feature       = "metabolite",
    #Peak_Area     = if_else(Peak_Area == 0, 1, Peak_Area),
    Peak_Area_Log = log10(Peak_Area + 1),
    point_color   = dplyr::if_else(sulfasalazine == "user", "#d1495b", "gray50"))

# plot
plot_sum <- ggplot(df_sum, aes(x = sulfasalazine, y = Peak_Area)) +
  geom_boxplot(aes(fill = sulfasalazine), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = point_color), size = 8, alpha = 0.7, shape = 16, stroke = NA, position = position_jitter(width = 0.15, height = 0)) +
  scale_color_identity() +
  scale_fill_manual(values = c("non-user" = "gray90", "user" = "#d1495b")) +
  scale_y_continuous(labels = scales::label_scientific()) +
  labs(x = "Sulfasalazine", y = "metabolite") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size = 14), axis.title  = element_text(size = 18, face = "bold"), legend.position = "none")

plot_sum

# Fisher if presence / absence (raw values of one group all 0)
# Create detection column
fisher_df <- df_sum |> dplyr::mutate(detected = ifelse(Peak_Area > 1, "detected", "not_detected"))
table(fisher_df$sulfasalazine, fisher_df$detected)
fisher.test(table(fisher_df$sulfasalazine, fisher_df$detected))

# Wilcoxon
wilcox_result <- wilcox.test(Peak_Area_Log ~ sulfasalazine, data = df_sum, exact = FALSE)
wilcox_result
#Save 
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/sf_glyceraldehyde.pdf", plot_cholyl_sum, width = 3, height = 7, dpi = 900)
#write_csv(df_melted_metadata_annotation, "df_melted_metadata_annotation_RA_COHORT.csv")
