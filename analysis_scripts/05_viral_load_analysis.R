################################################################################
### Comparing Viral Load and Persistence
# Contents:
# 1) Viral load trajectory comparison between LUNAR and Usual Care  
# 2) Persistence comparison between patients with escape and non-escape mutations (LUNAR only)
# 3) Viral load trajectory comparison between patients with escape and non-escape mutations (LUNAR only)

## Includes figures:
# 1) Figure 13: Comparison of viral load between LUNAR and Usual Care
# 2) Figure 26: Viral load persistence by presence of sotrovimab-associated mutation
# 3) Figure 27: Viral load trajectory by sotrovimab-escape mutation status 

# Load packages 
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggh4x)
library(ggpubr)



################################################################################################
#### 1) Viral load trajectory comparison between LUNAR and Usual Care
################################################################################################

# Refer to 01_QC_and_analysis_cohort_building.R script for LUNAR_UC_baseline_analysis_cohort

combined_df <- LUNAR_UC_baseline_analysis_cohort %>%
  mutate(
    Day = case_when(
      Cohort == "LUNAR" & Filter_SendRef == "BASELINE" ~ 0,
      Cohort == "LUNAR" & Filter_SendRef == "DAY7" ~ 7,
      Cohort == "LUNAR" & Filter_SendRef == "DAY14" ~ 14,
      Cohort == "LUNAR" & Filter_SendRef == "DAY28" ~ 28,
      
      Cohort == "Usual Care" & Filter_SendRef == "DAY1" ~ 0,
      Cohort == "Usual Care" & Filter_SendRef == "DAY2" ~ 2,
      Cohort == "Usual Care" & Filter_SendRef == "DAY3" ~ 3,
      Cohort == "Usual Care" & Filter_SendRef == "DAY4" ~ 4,
      Cohort == "Usual Care" & Filter_SendRef == "DAY5" ~ 5,
      Cohort == "Usual Care" & Filter_SendRef == "DAY6" ~ 6,
      Cohort == "Usual Care" & Filter_SendRef == "DAY7" ~ 7,
      Cohort == "Usual Care" & Filter_SendRef == "DAY10" ~ 10,
      Cohort == "Usual Care" & Filter_SendRef == "DAY14" ~ 14,
      
      TRUE ~ NA_real_
    )
  )

# Calculate median and IQR for each cohort and day
summary_df <- combined_df %>%
  group_by(Cohort, Day) %>%
  summarise(
    median_VL = median(Log10_VL, na.rm = TRUE),
    q25 = quantile(Log10_VL, 0.25, na.rm = TRUE),
    q75 = quantile(Log10_VL, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


vl_data <- combined_df %>%
  filter(!is.na(Log10_VL), !is.na(Cohort)) %>%
  distinct(Filter_Name, Day, Cohort, Log10_VL)   # one row per sample

pvals <- vl_data %>%
  group_by(Day) %>%
  filter(n_distinct(Cohort) == 2) %>%
  summarise(
    p = wilcox.test(Log10_VL ~ Cohort)$p.value,
    n1 = sum(Cohort == "LUNAR"),
    n2 = sum(Cohort == "Usual Care"),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    group1 = "LUNAR",
    group2 = "Usual Care",
    label = paste0("p = ", signif(p_adj, 2)),
    y.position = max(vl_data$Log10_VL, na.rm = TRUE) + 0.5
  )

pvals <- pvals %>%
  mutate(label = paste0("Day ", Day, ": p = ", signif(p_adj, 2)))


############## Figure 13: Comparison of viral load between LUNAR and Usual Care ####################### 
ggplot() +
  # Individual patient trajectories
  geom_line(
    data = combined_df,
    aes(x = Day,
        y = Log10_VL,
        group = Filter_Name),
    colour = "grey75",
    alpha = 0.20,
    linewidth = 0.5
  ) +
  
  # IQR ribbon
  geom_ribbon(
    data = summary_df,
    aes(x = Day,
        ymin = q25,
        ymax = q75,
        fill = Cohort,
        group = Cohort),
    alpha = 0.25
  ) +
  
  # Median line
  geom_line(
    data = summary_df,
    aes(x = Day,
        y = median_VL,
        colour = Cohort,
        group = Cohort),
    linewidth = 1.8
  ) +
  
  # Median points
  geom_point(
    data = summary_df,
    aes(x = Day,
        y = median_VL,
        colour = Cohort),
    size = 3
  ) +
    
  geom_text(
    data = pvals,
    aes(x = Day, y = y.position, label = label),
    inherit.aes = FALSE,
    size = 4,
  #  angle = 45,
    hjust = 0,
    colour = "grey20"       
  ) +
  
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 10, 14, 28)
  ) +
  
  labs(
    x = "Days After Baseline",
    y = expression("Viral Load ("~log[10]~"copies/ml )"),
    colour = "Cohort",
    fill = "Cohort"
  ) +
  
  theme_classic(base_size = 15)





################################################################################################
#### 2) Persistence comparison between patients with escape and non-escape mutations 
####### (LUNAR only)
################################################################################################

# Load raw LUNAR data without the filtering  
VCF <- read.csv("path/to/your/data/directory/LUNAR_data.csv")
VCF_tib <- as_tibble(VCF)


drug_patients <- unique(sotrovimab_mutations$Filter_Name)

patient_duration <- VCF_tib %>%
  mutate(
    day = case_when(
      Filter_SendRef == "BASELINE" ~ 0,
      Filter_SendRef == "DAY7" ~ 7,
      Filter_SendRef == "DAY14" ~ 14,
      Filter_SendRef == "DAY28" ~ 28
    )
  ) %>%
  group_by(Filter_Name) %>%
  summarise(last_day = max(day, na.rm = TRUE), .groups = "drop") %>%
  mutate(drug_mutation_patient = Filter_Name %in% drug_patients)

# Base R Wilcoxon test
w <- wilcox.test(last_day ~ drug_mutation_patient,
                 data = patient_duration)

# Format p-value
p_lab <- paste0("p = ",
                format.pval(w$p.value, digits = 3, eps = 0.001))

# Position for label
y_pos <- max(patient_duration$last_day, na.rm = TRUE) * 1.08

########### Figure 26: Viral load persistence by presence of sotrovimab-associated mutation #################
ggplot(patient_duration,
       aes(x = drug_mutation_patient,
           y = last_day,
           fill = drug_mutation_patient)) +
  geom_boxplot(width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    x = "Sotrovimab escape mutation",
    y = "Last sampled day",
    title = "Persistence by Drug Mutation Status"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  annotate("text",
           x = 1.5,
           y = y_pos,
           label = p_lab,
           size = 4)




################################################################################################
#### 3) Viral load trajectory comparison between patients with escape and non-escape mutations 
####### (LUNAR only)
################################################################################################
VCF_tib <- VCF_tib %>%
  left_join(
    patient_duration %>%
      select(Filter_Name, drug_mutation_patient),
    by = "Filter_Name"
  )

VCF_tib <- VCF_tib %>%
  mutate(
    Filter_SendRef = factor(
      Filter_SendRef,
      levels = c( "BASELINE", "DAY7", "DAY14", "DAY28")
    )
  )


vl_plot_data <- VCF_tib %>%
  filter(!is.na(Filter_SendRef), !is.na(Log10_VL),
         !is.na(drug_mutation_patient)) %>%
  distinct(Filter_Name, Filter_SendRef, Log10_VL, drug_mutation_patient)


pvals_drug <- vl_plot_data %>%
  group_by(Filter_SendRef) %>%
  filter(n_distinct(drug_mutation_patient) == 2) %>%
  summarise(
    p     = wilcox.test(Log10_VL ~ drug_mutation_patient)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj      = p.adjust(p, method = "BH"),              
    label      = paste0("p = ", signif(p_adj, 2)),
    y.position = max(vl_plot_data$Log10_VL, na.rm = TRUE) + 0.4
  )

############# Figure 27: Viral load trajectory by sotrovimab-escape mutation status #####################
ggplot(VCF_tib %>% filter(!is.na(Filter_SendRef)),
       aes(Filter_SendRef,
           Log10_VL,
           group = Filter_Name,
           colour = drug_mutation_patient)) +

  geom_line(alpha = 0.3) +

  stat_summary(
    aes(group = drug_mutation_patient),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +

  stat_summary(
    aes(group = drug_mutation_patient),
    fun = mean,
    geom = "point",
    size = 3
  ) +

  geom_text(
    data        = pvals_drug,
    aes(x = Filter_SendRef, y = y.position, label = label),
    inherit.aes = FALSE,
    size        = 6,
    colour      = "grey20"
  )  +
  
  labs(
    x = "Timepoint",
    y = expression("Viral Load ("~log[10]~"copies/ml )"),
    colour = "Drug mutation",
    title = "Viral Load Trajectory by Sotrovimab Mutation Status"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 12),  
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 16),
  axis.text  = element_text(size = 14),
  plot.title = element_text(size = 18, )
  )
