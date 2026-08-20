################################################################################
######### Load the Data ########################################################
################################################################################
library(tidyverse)
library(patchwork)

UC <- read.csv("path/to/your/data/directory/UsualCare.csv")  # Loading variant long table 
UC <- as_tibble(UC)

LUNAR <- read.csv("path/to/your/data/directory/LUNAR_data.csv")  # Loading variant long table 
LUNAR <- as_tibble(LUNAR)

LUNAR_UC_trajectory_analysis_cohort <- read.csv("D:/MSc Dissertation R/LUNAR_UC_trajectory_analysis_cohort.csv")
LUNAR_UC_trajectory_analysis_cohort <- as_tibble(LUNAR_UC_trajectory_analysis_cohort)

LUNAR_UC_baseline <- read.csv("D:/MSc Dissertation R/LUNAR_UC_baseline.csv")
LUNAR_UC_baseline <- as_tibble(LUNAR_UC_baseline)


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
# 2. Add a column identifying as immunocompromised and control for both cohorts
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



