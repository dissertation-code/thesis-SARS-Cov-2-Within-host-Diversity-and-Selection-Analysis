################################################################################
### Quality Control and Final Longitudinal Analysis Cohort Building 
# Contents:
# 1) Loading the data
# 2) Formatting raw .csv files for analysis 
# 3) Quality control 
# 4) Combined cohort 1: LUNAR_UC_baseline_analysis_cohort
# 5) Selecting patients for longitudinal cohort
# 6) Combined cohort 2: LUNAR_UC_trajectory_analysis_cohort

## Includes figures:
# 1) Figure 9: Sampling timepoints for each LUNAR patient included in the final longitudinal analysis cohort
# 2) Figure 10: Sampling timepoints for each Usual Care patient included in the final longitudinal analysis cohort
# 3) Figure 11: Comparison between sequencing quality control metrics for LUNAR and Usual Care

# Load packages
library(tidyverse)
library(patchwork)
library(grid)


################################################################################
######### Load the Data ########################################################
################################################################################

UC <- read.csv("path/to/your/data/directory/UsualCare.csv")  # Loading variant long table 
UC <- as_tibble(UC)

LUNAR <- read.csv("path/to/your/data/directory/LUNAR_data.csv")  # Loading variant long table 
LUNAR <- as_tibble(LUNAR)


################################################################################
######### Format Usual Care to Match LUNAR #####################################
################################################################################
# 1. Change names of Filter name and Filter Sendref
# 2. Convert LUNAR Names to characters 
# 3. Add Log VL
# 4. Create factors of timepoints 
# 5. Add Lineage data and quality metrics 
# 6. Add a column identifying as immunocompromised and control for both cohorts 
# 7. Combine the two cohorts (First extract exact columns from LUNAR and ensure they match)


# View structure 
str(UC)

# Renaming UC Columns
UC <- UC %>% rename(Filter_Name = Filter.Name)
UC <- UC %>% rename(Filter_SendRef = Filter.SenderRef)

# Adding Log10 values of Viral Load 
UC %>%                                                      # Checking for the value which is a character 
  filter(is.na(suppressWarnings(as.numeric(VL)))) %>%
  distinct(VL)

UC <- UC %>% mutate(VL = as.numeric(gsub(",", "", VL)))     # Converting from chr -> num

UC <- UC %>%                                                # Calculating log values
  mutate(Log10_VL = log10(VL))

# Creating factors out of timepoints 
control_timelevels <- c("DAY1", "DAY2", "DAY3", "DAY4", "DAY5", "DAY6", "DAY7", "DAY10", "DAY14")
UC <- UC  %>% mutate(Filter_SendRef = factor(Filter_SendRef,
                                 levels = control_timelevels))

# Adding Lineage data to UC
UC_lineage <- read.csv("path/to/your/data/directory/panoramic_lineage.csv")
UC_lineage <- as_tibble(UC_lineage)

UC <- UC %>%
  left_join(
    UC_lineage %>%
      select(
        SampleID,
        LINEAGE,
        Nextclade,
        pct_Mapped_reads,
        Coverage_median,
        pct_Coverage10x
      ),
    by = "SampleID"
  )

UC <- UC %>% rename(Pct_Coverage_gt_10x =  pct_Coverage10x) # Rename column to match LUNAR
UC <- UC %>% rename(VOC_VUI =  Nextclade)
UC <- UC %>% mutate(Coverage_median = as.numeric(Coverage_median))
UC <- UC %>% select(-Coverage.median)
UC <- UC %>% rename(cov10X =  Cov_10x)
UC <- UC %>% rename(Pct_Mapped_reads =  pct_Mapped_reads)

# QC of the cohort 
UC_post_QC <- UC %>% filter(ALT_FREQ >= 0.05,
                             DP >= 100,
                             AO >= 5,
                            Pct_Coverage_gt_10x >= 87)


################################################################################
######### Combining Both Cohorts  ##############################################
################################################################################
# 1. Convert LUNAR Names to characters
# 2. Add a column identifying as LUNAR and Usual Care for both cohorts
# 3. Combine the two cohorts (First extract exact columns from LUNAR and ensure they match)
# 4. Output: 2 combined cohorts: For baseline and trajectory analysis 

# Converting int -> chr
LUNAR <- LUNAR %>% mutate(Filter_Name = as.character(Filter_Name))

# QC of Lunar
LUNAR_post_QC <- LUNAR %>% filter(ALT_FREQ >= 0.05,
                            DP >= 100,
                            AO >= 5,
                            Pct_Coverage_gt_10x >= 87)

# Load LUNAR lineage data 
LUNAR_lineage <- read.csv("path/to/your/data/directory/Lunar_lineage.csv")
LUNAR_lineage <- as_tibble(LUNAR_lineage)

# Add LUNAR lineage data to post QC LUNAR
LUNAR_lineage <- LUNAR_lineage %>% mutate(Patient = as.character(Patient))

LUNAR_lineage %>%                         # Checking if each sample has only 1 lineage 
  count(Patient, TimePoint) %>%
  filter(n > 1)

LUNAR_post_QC <- LUNAR_post_QC %>%
  left_join(
    LUNAR_lineage %>%
      select(Patient, TimePoint, VOC_VUI, LINEAGE),
    by = c(
      "Filter_Name" = "Patient",
      "Filter_SendRef" = "TimePoint"
    )
  )

########################## Comined Cohort1: Baseline Analysis ##################
# Contains all samples at all timepoints after applying QC metrics 
LUNAR_post_QC <- LUNAR_post_QC %>%
  mutate(Cohort = "LUNAR")

UC_post_QC <- UC_post_QC %>%
  mutate(Cohort = "Usual Care")

LUNAR_UC_baseline_analysis_cohort <- bind_rows(LUNAR_post_QC, UC_post_QC)


######################## Comined Cohort 2: Trajectory Analysis ##################
# Loading filtered LUNAR trajectory analysis data (filtered in the same way as Usual Care is done in this script)
LUNAR_trajectory_analysis_cohort <- read.csv("path/to/your/data/directory/master_emergent_snps_and_lineage.csv")
LUNAR_trajectory_analysis_cohort <- as_tibble(LUNAR_trajectory_analysis_cohort)

# Making LUNAR data compatible for joining
LUNAR_trajectory_analysis_cohort <- LUNAR_trajectory_analysis_cohort %>% mutate(Filter_Name = as.character(Filter_Name))
LUNAR_trajectory_analysis_cohort <- LUNAR_trajectory_analysis_cohort %>% mutate(Cohort = "LUNAR")
LUNAR_trajectory_analysis_cohort <- LUNAR_trajectory_analysis_cohort %>% select(-Patient)

# Creating the UC trajectory analysis cohort
UC_trajectory_analysis_cohort <- UC_post_QC %>%
  mutate(day = str_extract(Filter_SendRef, "\\d+"))

UC_trajectory_analysis_cohort <- UC_trajectory_analysis_cohort %>%
  group_by(Filter_Name) %>%
  filter(
    day %in% c("1", "14") |
      day == case_when(
        any(day == "5") ~ "5",
        any(day == "6") ~ "6",
        any(day == "4") ~ "4",
        TRUE ~ NA_character_
      )
  ) %>%
  ungroup()



# Creating the subset with the final patients 
UC_trajectory_analysis_cohort <- UC_trajectory_analysis_cohort %>%
  group_by(Filter_Name) %>%
  filter(
    all(c("DAY1", "DAY14") %in% Filter_SendRef) &
      ("DAY5" %in% Filter_SendRef | "DAY4" %in% Filter_SendRef)
  ) %>%
  filter(
    Filter_SendRef %in% c("DAY1", "DAY14") |
      Filter_SendRef == case_when(
        "DAY5" %in% Filter_SendRef ~ "DAY5",
        "DAY4" %in% Filter_SendRef ~ "DAY4"
      )
  ) %>%
  ungroup()



set.seed(123)   # Makes the random selection reproducible

# Find eligible patients
additional_patients <- UC_post_QC %>%
  group_by(Filter_Name) %>%
  filter(
    "DAY4" %in% Filter_SendRef
  ) %>%
  ungroup() %>%
  distinct(Filter_Name) %>%
  anti_join(
    UC_trajectory_analysis_cohort %>%
      distinct(Filter_Name),
    by = "Filter_Name"
  ) %>%
  slice_sample(n = 3)

# Extract all rows for those patients
additional_data <- UC_post_QC %>%
  semi_join(additional_patients, by = "Filter_Name")  %>%
  filter(Filter_SendRef %in% c("DAY1", "DAY4", "DAY5", "DAY6", "DAY14"))

# Add them to your cohort
UC_trajectory_analysis_cohort <- bind_rows(
  UC_trajectory_analysis_cohort,
  additional_data
)

UC_trajectory_analysis_cohort %>% distinct(Filter_Name)


### Removing baseline + propagated varuabts 
UC_baseline_variants <- UC_trajectory_analysis_cohort %>%
  filter(Filter_SendRef == "DAY1") %>%
  distinct(Filter_Name, CHR, POS, REF, ALT)

UC_trajectory_analysis_cohort <- UC_trajectory_analysis_cohort %>%
  filter(Filter_SendRef != "DAY1") %>%
  anti_join(
    UC_baseline_variants,
    by = c("Filter_Name", "CHR", "POS", "REF", "ALT")
  )


# Adding Variant ID to UC 
UC_trajectory_analysis_cohort <- UC_trajectory_analysis_cohort %>%
  mutate(
    Variant_ID = paste(POS, REF, ALT, sep = "_")
  )

# Only keeping the indels 
UC_trajectory_analysis_cohort <- UC_trajectory_analysis_cohort %>%
  filter(nchar(REF) == 1 & nchar(ALT) == 1)


LUNAR_UC_trajectory_analysis_cohort <- bind_rows(LUNAR_trajectory_analysis_cohort, UC_trajectory_analysis_cohort)



###########################################################################################################
# 1) Figure 9: Sampling timepoints for each LUNAR patient included in the final longitudinal analysis cohort
###########################################################################################################
VCF <- read.csv("path/to/your/data/directory/LUNAR_data.csv")
VCF_tib <- as_tibble(VCF)

# --- 1. One row per unique sample ---
sample_info <- VCF_tib %>%
  select(Filter_Name, Filter_SendRef) %>%
  distinct() %>%
  mutate(present = TRUE)

# --- 2. Define your 4 timepoints ---
timepoint_levels <- c("BASELINE","DAY7","DAY14","DAY28")

tp_colours <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
names(tp_colours) <- timepoint_levels

BASELINE <- timepoint_levels[1]  # adjust if your baseline isn't the first timepoint

patients_to_keep <- sample_info %>%
  filter(Filter_SendRef != BASELINE) %>%
  group_by(Filter_Name) %>%
  summarise(n_post_baseline = n_distinct(Filter_SendRef)) %>%
  filter(n_post_baseline >= 2) %>%
  pull(Filter_Name)

all_combos_filtered <- all_combos %>%
  filter(Filter_Name %in% patients_to_keep)

# --- 3. Complete grid ---
all_combos_filtered <- expand.grid(
  Filter_Name    = unique(patients_to_keep),
  Filter_SendRef = timepoint_levels,
  stringsAsFactors = FALSE
) %>%
  left_join(sample_info, by = c("Filter_Name", "Filter_SendRef")) %>%
  mutate(present = replace_na(present, FALSE))

# --- 4. Plot ---
ggplot(all_combos_filtered, aes(x = factor(Filter_SendRef, levels = timepoint_levels), 
                       y = 10, 
                       fill = ifelse(present, as.character(Filter_SendRef), NA))) +
  geom_col(colour = "white", linewidth = 0.5) +
  facet_wrap(~ Filter_Name, ncol = 10) +
  scale_fill_manual(
    values       = tp_colours,
    na.value     = "grey90",
    name         = "Timepoint",
    na.translate = FALSE
  ) +
  scale_y_continuous(breaks = NULL) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(colour = "white", face = "bold", size = 7,
                                    margin = margin(3, 3, 3, 3)),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "grey97"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 16),
    legend.text      = element_text(size = 14),
    legend.key.size  = unit(1.2, "cm"),
    plot.margin      = margin(10, 10, 10, 10),
  )

ggsave("patient_timepoint_grid_filtered.png", width = 14, height = 10, dpi = 300)


###########################################################################################################
# 2) Figure 10: Sampling timepoints for each Usual Care patient included in the final longitudinal analysis cohort
###########################################################################################################
collapse_tp <- function(x) {
  case_when(
    x %in% c("DAY4","DAY5","DAY6") ~ "DAY5 ±1",
    TRUE                            ~ x
  )
}

sample_info <- LUNAR_UC_trajectory_analysis_cohort %>%
  filter(Cohort == "Usual Care") %>%
  distinct(Filter_Name, Filter_SendRef) %>%
  mutate(Filter_SendRef = collapse_tp(Filter_SendRef)) %>%
  distinct() %>%                      # merges DAY4 + DAY5 in the same patient
  mutate(present = TRUE)

# add DAY1 as before
day1_info <- UC_post_QC %>%
  filter(Filter_Name %in% sample_info$Filter_Name,
         Filter_SendRef == "DAY1") %>%
  distinct(Filter_Name, Filter_SendRef) %>%
  mutate(present = TRUE)

sample_info <- bind_rows(sample_info, day1_info)

timepoint_levels <- c("DAY1","DAY5 ±1","DAY14")
tp_colours <- c("#E69F00","#56B4E9","#009E73")
names(tp_colours) <- timepoint_levels


# --- Complete grid ---
all_combos <- expand.grid(
  Filter_Name    = unique(sample_info$Filter_Name),
  Filter_SendRef = timepoint_levels,
  stringsAsFactors = FALSE
) %>%
  left_join(sample_info, by = c("Filter_Name", "Filter_SendRef")) %>%
  mutate(present = replace_na(present, FALSE))

# --- Plot ---
ggplot(all_combos, aes(x = factor(Filter_SendRef, levels = timepoint_levels), 
                       y = 10, 
                       fill = ifelse(present, as.character(Filter_SendRef), NA))) +
  geom_col(colour = "white", linewidth = 0.5) +
  facet_wrap(~ Filter_Name, ncol = 7) +
  scale_fill_manual(
    values       = tp_colours,
    na.value     = "grey90",
    name         = "Timepoint",
    na.translate = FALSE
  ) +
  scale_y_continuous(breaks = NULL) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(colour = "white", face = "bold", size = 20,
                                    margin = margin(3, 3, 3, 3)),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "grey97"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 16),
    legend.text      = element_text(size = 14),
    legend.key.size  = unit(1.2, "cm"),
    plot.margin      = margin(10, 10, 10, 10),
  )

ggsave("patient_timepoint_grid_new.png", width = 14, height = 10, dpi = 300)




###########################################################################################################
# 3) Figure 11: Comparison between sequencing quality control metrics for LUNAR and Usual Care
###########################################################################################################
sample_qc <- LUNAR_UC_trajectory_analysis_cohort %>%
  distinct(
    Cohort,
    Filter_Name,
    Filter_SendRef,
    Coverage_median,
    Pct_Mapped_reads,
    Pct_Coverage_gt_10x
  )


qc_long <- sample_qc %>%
  select(Cohort, Coverage_median,
         Pct_Mapped_reads,
         Pct_Coverage_gt_10x) %>%
  pivot_longer(-Cohort,
               names_to = 'Metric',
               values_to = 'Value')

ggplot(qc_long, aes(Cohort, Value)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~Metric, scales = 'free_y', nrow = 1) +
  theme_classic() +
  labs(x = NULL, y = NULL)
