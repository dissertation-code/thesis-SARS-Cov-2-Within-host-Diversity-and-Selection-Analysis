################################################################################
### Measuring Within-host Diversity in LUNAR and Usual Care
# Contents:
# 1) Mutation burden per sample per timepoint 
# 2) Per-site nucleotide diversity (π site) (LUNAR only)
# 3) Gene-specific nuceleotide diversity (π gene) (LUNAR only)

# Load packages 
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggh4x)
library(ggpubr)

################################################################################
######## 1) Mutation burden per sample per timepoint 
############# Figure 13 & 14 
################################################################################
# Load the LUNAR + UC longitudinal analysis cohort (refer 01_QC_and_analysis_cohort_building.R)
LUNAR_UC_trajectory_analysis_cohort <- read.csv("path/to/your/data/directory/LUNAR_UC_trajectory_analysis_cohort.csv")
LUNAR_UC_trajectory_analysis_cohort <- as_tibble(LUNAR_UC_trajectory_analysis_cohort)

# Count variants per sample
variant_count_per_sample <- LUNAR_UC_trajectory_analysis_cohort %>%
  group_by(Cohort, Filter_Name, Filter_SendRef) %>%
  summarise(
    n_variants = n_distinct(paste(POS, REF, ALT)),
    .groups = "drop"
  )

# Order timepoints
variant_count_per_sample <- variant_count_per_sample %>%
  mutate(
    Filter_SendRef = factor(
      Filter_SendRef,
      levels = c("DAY4", "DAY5","DAY6", "DAY7", "DAY14", "DAY28")
    )
  )

# Median variants per cohort at each timepoint
median_variants_per_sample <- variant_count_per_sample %>%
  group_by(Cohort, Filter_SendRef) %>%
  summarise(
    median_variants = median(n_variants),
    Q1 = quantile(n_variants, 0.25, na.rm = TRUE),
    Q3 = quantile(n_variants, 0.75, na.rm = TRUE),
    IQR_variants = IQR(n_variants, na.rm = TRUE),
    .groups = "drop"
  )


### Calulating the p-value
Day14_p_data <- variant_count_per_sample %>%
  filter(Filter_SendRef == "DAY14")

wilcox.test(n_variants ~ Cohort,
            data = Day14_p_data,
            exact = FALSE)

######## Mutation burden per sample per timepoint: Trajectory Plot (Figure 13) #########
ggplot() +
  # Individual patient trajectories
  geom_line(
    data = variant_count_per_sample,
    aes(
      x = Filter_SendRef,
      y = n_variants,
      group = interaction(Cohort, Filter_Name),
      colour = Cohort
    ),
    alpha = 0.25
  ) +
  # Cohort medians
  geom_line(
    data = median_variants_per_sample,
    aes(
      x = Filter_SendRef,
      y = median_variants,
      group = Cohort,
      colour = Cohort
    ),
    linewidth = 1.5
  ) +
  geom_point(
    data = median_variants_per_sample,
    aes(
      x = Filter_SendRef,
      y = median_variants,
      colour = Cohort
    ),
    size = 3
  ) +
  scale_colour_manual(
    values = c(
      "Usual Care" = "#1f78b4",
      "LUNAR" = "#e31a1c"
    )
  ) +
  theme_bw() +
  labs(
    x = "Timepoint",
    y = "No. of variants per sample",
    colour = "Cohort",
    title = "Number of variants per sample across timepoints",
    subtitle = "Thick coloured lines show the cohort median"
  )

########### Distribution of variants across SARS-CoV-2 genes (Day 14) (Figure 14) #################
cohort_gene_counts <- LUNAR_UC_trajectory_analysis_cohort %>% filter(Filter_SendRef == "DAY14") %>%
  group_by(Cohort, GENE_f) %>%
  summarise(
    n_mutations = n(),
    .groups = "drop"
  )

cohort_gene_props <- cohort_gene_counts %>%
  group_by(Cohort) %>%
  mutate(
    proportion = n_mutations / sum(n_mutations),
    percent = proportion * 100
  ) %>%
  ungroup()

ggplot(cohort_gene_props,
       aes(x = Cohort,
           y = percent,
           fill = GENE_f)) +
  
  geom_col(width = 0.7) +
  
  labs(
    title = "Distribution of variants across SARS-CoV-2 genes (Day 14)",
    x = NULL,
    y = "Mutations (%)",
    fill = "Gene"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )




################################################################################
######## 2) Within-host per-site nucleotide diversity (π site)  
############# Figure 15
################################################################################

# Loading only LUNAR longitudinal analysis cohort 
cohort <- read.csv("path/to/your/data/directory/master_emergent_snps_and_lineage.csv") 
cohort <- as_tibble(cohort)

# Calculate heterozygosity at each site
df_snps <- cohort %>%
  mutate(
    h = 2 * ALT_FREQ * (1 - ALT_FREQ)
  )

pi_sample <- df_snps %>%
  group_by(Filter_Name, Filter_SendRef) %>%
  summarise(
    pi = sum(h, na.rm = TRUE),
    n_variants = n(),
    Pct_Coverage_gt_10x = first(Pct_Coverage_gt_10x),
    .groups = "drop"
  )  %>%
  mutate(mean_h = pi / n_variants)


GENOME_LENGTH <- 29903

pi_sample <- pi_sample %>%
  mutate(
    callable_sites =
      (Pct_Coverage_gt_10x / 100) * GENOME_LENGTH
  )

pi_sample <- pi_sample %>%
  mutate(
    pi_site = pi / callable_sites
  )


#### Trajectory Plot (Figure 15)
pi_sample <- pi_sample %>%
  mutate(
    Filter_SendRef = factor(
      Filter_SendRef,
      levels = c("DAY7", "DAY14", "DAY28")
    )
  )

summary_pi <- pi_sample %>%
  group_by(Filter_SendRef) %>%
  summarise(
    median_pi = median(pi, na.rm = TRUE),
    adjusted_median_pi = median(mean_h, na.rm = TRUE),
    median_pi_site = median(pi_site, na.rm = TRUE),
    Q1 = quantile(pi_site, 0.25, na.rm = TRUE),
    Q3 = quantile(pi_site, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot() +
  geom_line(data = pi_sample,
            aes(Filter_SendRef,
                pi_site,
                group = Filter_Name),
            alpha = 0.3,
            colour = "grey60") +
  geom_line(data = summary_pi,
            aes(Filter_SendRef,
                median_pi_site,
                group = 1),
            linewidth = 1.5,
            colour = "red") +
  geom_point(data = summary_pi,
             aes(Filter_SendRef,
                 median_pi_site),
             size = 3,
             colour = "red") +
  theme_bw() +
  labs(
    x = "Timepoint",
    y = expression(pi),
    title = "Within-host per-site diversity trajectories",
    subtitle = "Red line shows the cohort median"
  )




################################################################################
######## 3) Gene-specific nucleotide diversity (π gene)  
############# Figure 16
################################################################################

pi_gene <- df_snps %>%
  mutate(
    h = 2 * ALT_FREQ * (1 - ALT_FREQ)
  ) %>%
  group_by(Filter_Name,
           Filter_SendRef,
           GENE_f) %>%
  summarise(
    pi = sum(h, na.rm = TRUE),
    variants = n(),
    .groups = "drop"
  )



gene_lengths <- tibble(
  GENE_f = c(
    "ORF1ab","S","NSP12",
    "NSP13","N","ORF3a",
    "NSP8","ORF7a","NSP7",
    "M","ORF6","ORF8",
    "E","ORF10","ORF7b"
  ),
  gene_length = c(
    21291,3822,2796,
    1803,1260,828,
    594,366,249,
    669,186,366,
    228,117,132
  )
)


# nucleotide diversity per-gene level 
gene_pi <- cohort %>%
  mutate(
    h = 2 * ALT_FREQ * (1 - ALT_FREQ)
  ) %>%
  group_by(GENE_f) %>%
  summarise(
    pi = sum(h, na.rm = TRUE),
    n_variants = n(),
    mean_h = mean(h, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(gene_lengths,
            by = "GENE_f") %>%
  mutate(
    pi_site_gene = pi / gene_length
  ) %>%
  arrange(desc(pi_site_gene))



# nucleotide diversity per-gene level + time
gene_time_pi <- cohort %>%
  filter(
    nchar(REF) == 1,
    nchar(ALT) == 1
  ) %>%
  mutate(
    h = 2 * ALT_FREQ * (1 - ALT_FREQ)
  ) %>%
  group_by(Filter_SendRef, GENE_f) %>%
  summarise(
    pi = sum(h, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(gene_lengths,
            by = "GENE_f") %>%
  mutate(
    pi_site_gene = pi / gene_length
  )


# Setting timepoints
gene_time_pi <- gene_time_pi %>%
  mutate(
    Filter_SendRef = factor(
      Filter_SendRef,
      levels = c("DAY7", "DAY14", "DAY28")
    )
  )


# Heatmap of gene-specific nucleotide diversity trajectory (Figure 16)
ggplot(
  gene_time_pi,
  aes(
    x = Filter_SendRef,
    y = GENE_f,
    fill = pi_site_gene
  )
) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    x = "Timepoint",
    y = "Gene",
    fill = expression(pi[site]),
    title = "Gene-specific nucleotide diversity through time"
  ) +
  theme_classic()
