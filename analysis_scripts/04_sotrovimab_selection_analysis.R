################################################################################
### Sotrovimab Selection Analysis: Detection, Recurrence and Consensus-crossing events 
# Contents:
# 1) Well-documented sotrovimab escape mutations  
# 2) Novel/Potential sotrovimab escape mutations
# 3) Non-escape mutations 

# This analysis is conducted on LUNAR data only, since no escape mutations were detected in Usual Care

# Load packages 
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggh4x)
library(ggpubr)


# Load data
cohort <- read.csv("path/to/your/data/directory/master_emergent_snps_and_lineage.csv")
cohort <- as_tibble(cohort)


######################################################################################
### 1) Well-documented sotrovimab escape mutations 
######################################################################################

variants_to_remove <- c("N440K", "S373P", "ST375FA") # Lineage defining mutations 

cohort <- cohort %>%
  filter(!AA_change_f %in% variants_to_remove)

G339D <- cohort %>% filter(AA_change_f == "G339D") 

S371F <- cohort %>% filter(AA_change_f == "S371F") 

K356T <- cohort %>% filter(AA_change_f == "K356T") 

P337L <- cohort %>% filter(AA_change_f == "P337L") 

P337R <- cohort %>% filter(AA_change_f == "P337R") 

P337H <- cohort %>% filter(AA_change_f == "P337H")  

S373P <- cohort %>% filter(AA_change_f == "S373P") 

N440K <- cohort %>% filter(AA_change_f == "N440K") 

dual <- cohort %>% filter(AA_change_f %in% c("GE339DQ", "GE339DK", "GE339DA", "GE339DV", "GE339DG",
                                             "GE339HQ", "GE339HK", "GE339HG", "GE339DD") )



# Combining all mutations into one dataframe
sotrovimab_mutations <-  bind_rows(G339D, K356T, P337H, P337L, P337R, S371F, S373P, N440K, dual)  %>%    
  mutate(
    Filter_SendRef = factor(
      Filter_SendRef,
      levels = c("DAY7", "DAY14", "DAY28")
    )
  )

# Define epitope mutatiosn lying in the sotrovimab binding regions 
epitope_positions <- c(334:356)  
sotrovimab_mutations <- sotrovimab_mutations %>%
  mutate(
    aa_pos = readr::parse_number(AA_change_f),  
    epitope_status = if_else(aa_pos %in% epitope_positions,
                             "Epitope", "Non-epitope")
  )

sotrovimab_mutations <- sotrovimab_mutations %>%
  mutate(AA_change_f = factor(AA_change_f, levels = sort(unique(AA_change_f))))

strip_cols <- sotrovimab_mutations %>%
  distinct(AA_change_f, epitope_status) %>%
  arrange(AA_change_f) %>%
  mutate(fill = if_else(epitope_status == "Epitope", "#F4A582", "#D9D9D9"))

# Creating a variant summary table for variant
variant_summary <- sotrovimab_mutations %>%
  group_by(Variant_ID, AA_change_f, GENE_f) %>%
  summarise(
    n = n(),
    patients = paste0(
      unique(paste0(Filter_Name, " (", Filter_SendRef, ")")),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>% 
  arrange(desc(n))

# Escape mutation trajectories 
epitope_cols <- c("Epitope" = "#F4A582", "Non-epitope" = "#D9D9D9")

sotrovimab_mutations %>%
  mutate(Filter_Name = factor(Filter_Name)) %>%
  ggplot(aes(Filter_SendRef, ALT_FREQ,
             group = Filter_Name, colour = Filter_Name)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
   geom_point(aes(fill = epitope_status),
             shape = 22, alpha = 0, colour = NA, show.legend = TRUE) +
  geom_hline(data = data.frame(yintercept = 0.50),
             aes(yintercept = yintercept,
                 linetype = "Consensus threshold (50%)"),
             colour = "red", inherit.aes = FALSE) +
  facet_wrap2(
    ~AA_change_f, scales = "free_y",
    strip = strip_themed(background_x = elem_list_rect(fill = strip_cols$fill))
  ) +  
  scale_y_continuous(limits = c(0, 1), name = "Variant frequency") +
   scale_colour_viridis_d(option = "turbo") + 
  scale_fill_manual(values = epitope_cols, name = "Position") +
  scale_linetype_manual(values = c("Consensus threshold (50%)" = "dashed"),
                        name = NULL) +
  guides(
    colour = "none",
    fill = guide_legend(
      order = 2,
      override.aes = list(alpha = 1, size = 5, shape = 22, colour = "black")
    )
  ) +
  labs(x = "Post-treatment Sampling Timepoints", colour = "Patient ID") +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
  legend.text = element_text(size = 13),   # the key labels
  legend.title = element_text(size = 14),
  axis.title = element_text(size = 14),   # "Days post-treatment" / "Variant frequency"
  axis.text = element_text(size = 11)  )


############## Recurrent Sotrovimab Escape Mutations #################################
ggplot(
  sotrovimab_mutations %>% distinct(Filter_Name, GENE_f, AA_change_f) %>%   # one occurrence per patient
    count(GENE_f, AA_change_f, name = "n_patients") %>%
    arrange(desc(n_patients)) %>% filter(n_patients >= 2),
  aes(x = reorder(AA_change_f, n_patients),
      y = n_patients,
      fill = GENE_f)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Mutation",
    y = "Number of patients",
    title = "Recurrent Sotrovimab-Escape Mutations"
  ) +
  theme_classic()


#################### Sotrovimab Escape Mutations reaching Consensus ################
sot_fixation_table <- sotrovimab_mutations %>% filter(isRep == "NO") %>%
  group_by(Filter_Name, AA_change_f, Filter_SendRef) %>%
  summarise(
    ALT_FREQ = max(ALT_FREQ, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Filter_SendRef,
    values_from = ALT_FREQ
  )

sample_matrix <- cohort %>%
  distinct(Filter_Name, Filter_SendRef) %>%
  mutate(sample_taken = TRUE) %>%
  pivot_wider(
    names_from = Filter_SendRef,
    values_from = sample_taken,
    values_fill = FALSE
  )

sot_fixation_table <- sot_fixation_table %>%
  left_join(sample_matrix, by = "Filter_Name")

sot_fixation_table <- sot_fixation_table %>%
  mutate(
    DAY7.x  = ifelse(DAY7.y  & is.na(DAY7.x),  0, DAY7.x),
    DAY14.x = ifelse(DAY14.y & is.na(DAY14.x), 0, DAY14.x),
    DAY28.x = ifelse(DAY28.y & is.na(DAY28.x), 0, DAY28.x)
  )


sot_fixation_table <- sot_fixation_table %>%
  mutate(
    fixation_time = case_when(
      
      DAY7.x >= 0.5 &
        (is.na(DAY14.x) | DAY14.x >= 0.5) &
        (is.na(DAY28.x) | DAY28.x >= 0.5) ~ "DAY7",
      
      (is.na(DAY7.x) | DAY7.x < 0.5) &
        DAY14.x >= 0.5 &
        (is.na(DAY28.x) | DAY28.x >= 0.5) ~ "DAY14",
      
      (is.na(DAY7.x) | DAY7.x < 0.5) &
        (is.na(DAY14.x) | DAY14.x < 0.5) &
        DAY28.x >= 0.5 ~ "DAY28",
      
      TRUE ~ NA_character_
    )
  )

sot_fixation_table <- sot_fixation_table %>% filter(!is.na(fixation_time) )



sot_fixation_counts <- sot_fixation_table %>%
  count(fixation_time) %>%    
  mutate(
    fixation_time = factor(
      fixation_time,
      levels = c("DAY7", "DAY14", "DAY28")
    )
  )

p1 <- ggplot(
  sot_fixation_counts,
  aes(
    x = fixation_time,
    y = n
  )
) +
   geom_col(fill = "#2C7FB8", colour = "black") +
  theme_bw() +
  labs(
    x = "Time of consensus-crossing",
    y = "Number of mutations"
  )


sot_recurrent_fixation <- sot_fixation_table %>%
  filter(!is.na(fixation_time)) %>%
  count(AA_change_f) %>%
  arrange(desc(n))

p2 <- ggplot(
  sot_recurrent_fixation,
  aes(
    x = reorder(AA_change_f, n),
    y = n
  )
) +
  geom_col(fill = "#2C7FB8", colour = "black") +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Mutation",
    y = "Patients reaching consensus"
  )


p1 + p2 +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "A")




######################################################################################
### 2) Novel/Potential sotrovimab escape mutations 
######################################################################################

AA_pos_cohort <- cohort %>%
  mutate(AA_pos = parse_number(AA_change_f)) 

novel_sotrovimab_mutations <- AA_pos_cohort %>%
   filter(
     AA_pos >= 334,
     AA_pos <= 356,
     GENE_f == "S",
     !AA_change_f %in% sotrovimab_mutations$AA_change_f
   ) %>%    
  mutate(
    Filter_SendRef = factor(
      Filter_SendRef,
      levels = c("DAY7", "DAY14", "DAY28")
    )
  )

novel_non_epitope <- AA_pos_cohort %>%
  filter(
    GENE_f == "S",
    AA_pos %in% c(371,373,375,417,440,446,452,477,478,484,493,496,498,501,505),
    !AA_change_f %in% sotrovimab_mutations$AA_change_f,
     !AA_change_f %in% novel_sotrovimab_mutations$AA_change_f) %>%
  mutate(
    Filter_SendRef = factor(
      Filter_SendRef,
      levels = c("DAY7", "DAY14", "DAY28")
    )
  )

novel_sotrovimab_mutations <-  bind_rows(novel_sotrovimab_mutations, novel_non_epitope)  %>%    
  mutate(
    Filter_SendRef = factor(
      Filter_SendRef,
      levels = c("DAY7", "DAY14", "DAY28")
    )
  )
novel_sotrovimab_mutations <- novel_sotrovimab_mutations %>%
  mutate(
    aa_pos = readr::parse_number(AA_change_f),   # "S:E340K" -> 340
    epitope_status = if_else(aa_pos %in% epitope_positions,
                             "Epitope", "Non-epitope")
  )

novel_sotrovimab_mutations <- novel_sotrovimab_mutations %>%
  mutate(AA_change_f = factor(AA_change_f, levels = sort(unique(AA_change_f))))

novel_strip_cols <- novel_sotrovimab_mutations %>%
  distinct(AA_change_f, epitope_status) %>%
  arrange(AA_change_f) %>%
  mutate(fill = if_else(epitope_status == "Epitope", "#F4A582", "#D9D9D9"))

novel_variant_summary <- novel_sotrovimab_mutations %>%
  group_by(Variant_ID, AA_change_f, GENE_f) %>%
  summarise(
    n = n(),
    patients = paste0(
      unique(paste0(Filter_Name, " (", Filter_SendRef, ")")),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>% 
  arrange(desc(n))

novel_variant_summary



### Novel Escape Mutation Trajectories 
epitope_cols <- c("Epitope" = "#F4A582", "Non-epitope" = "#D9D9D9")

novel_sotrovimab_mutations %>%
  mutate(Filter_Name = factor(Filter_Name)) %>%
  ggplot(aes(Filter_SendRef, ALT_FREQ,
             group = Filter_Name, colour = Filter_Name)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
   geom_point(aes(fill = epitope_status),
             shape = 22, alpha = 0, colour = NA, show.legend = TRUE) +
  geom_hline(data = data.frame(yintercept = 0.50),
             aes(yintercept = yintercept,
                 linetype = "Consensus threshold (50%)"),
             colour = "red", inherit.aes = FALSE) +
  facet_wrap2(
    ~AA_change_f, scales = "free_y",
    strip = strip_themed(background_x = elem_list_rect(fill = novel_strip_cols$fill))
  ) +
  scale_y_continuous(limits = c(0, 1), name = "Variant frequency") +
   scale_colour_viridis_d(option = "turbo") +
  scale_fill_manual(values = epitope_cols, name = "Position") +
  scale_linetype_manual(values = c("Consensus threshold (50%)" = "dashed"),
                        name = NULL) +
  guides(
    colour = "none",
    fill = guide_legend(
      order = 2,
      override.aes = list(alpha = 1, size = 5, shape = 22, colour = "black")
    )
  ) +
  labs(x = "Post-treatment Sampling Timepoints", colour = "Patient ID") +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
  legend.text = element_text(size = 13),   # the key labels
  legend.title = element_text(size = 14),
  axis.title = element_text(size = 14),   # "Days post-treatment" / "Variant frequency"
  axis.text = element_text(size = 11)  )


############# Novel/Potential Escape Mutations Recurrence ###########################
ggplot(
  novel_sotrovimab_mutations %>% distinct(Filter_Name, GENE_f, AA_change_f) %>%   # one occurrence per patient
    count(GENE_f, AA_change_f, name = "n_patients") %>%
    arrange(desc(n_patients)) %>% filter(n_patients >= 2),
  aes(x = reorder(AA_change_f, n_patients),
      y = n_patients,
      fill = GENE_f)
) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 20, by = 2)) +
  labs(
    x = "Mutation",
    y = "Number of patients",
    title = "Recurrent Potential Sotrovimab-Escape Mutations"
  ) +
  theme_classic()


############# Novel/Potential Escape Mutations reaching Consensus ####################
sample_matrix <- cohort %>%
  distinct(Filter_Name, Filter_SendRef) %>%
  mutate(sample_taken = TRUE) %>%
  pivot_wider(
    names_from = Filter_SendRef,
    values_from = sample_taken,
    values_fill = FALSE
  )


novel_sot_fixation_table <- novel_sotrovimab_mutations %>% filter(isRep == "NO") %>%
  group_by(Filter_Name, AA_change_f, Filter_SendRef) %>%
  summarise(
    ALT_FREQ = max(ALT_FREQ, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Filter_SendRef,
    values_from = ALT_FREQ
  )

novel_sot_fixation_table <- novel_sot_fixation_table %>%
  left_join(sample_matrix, by = "Filter_Name")

novel_sot_fixation_table <- novel_sot_fixation_table %>%
  mutate(
    DAY7.x  = ifelse(DAY7.y  & is.na(DAY7.x),  0, DAY7.x),
    DAY14.x = ifelse(DAY14.y & is.na(DAY14.x), 0, DAY14.x),
    DAY28.x = ifelse(DAY28.y & is.na(DAY28.x), 0, DAY28.x)
  )


novel_sot_fixation_table <- novel_sot_fixation_table %>%
  mutate(
    fixation_time = case_when(
      
      DAY7.x >= 0.5 &
        (is.na(DAY14.x) | DAY14.x >= 0.5) &
        (is.na(DAY28.x) | DAY28.x >= 0.5) ~ "DAY7",
      
      (is.na(DAY7.x) | DAY7.x < 0.5) &
        DAY14.x >= 0.5 &
        (is.na(DAY28.x) | DAY28.x >= 0.5) ~ "DAY14",
      
      (is.na(DAY7.x) | DAY7.x < 0.5) &
        (is.na(DAY14.x) | DAY14.x < 0.5) &
        DAY28.x >= 0.5 ~ "DAY28",
      
      TRUE ~ NA_character_
    )
  )

novel_sot_fixation_table <- novel_sot_fixation_table %>% filter(!is.na(fixation_time) )

novel_sot_fixation_counts <- novel_sot_fixation_table %>%
  count(fixation_time) %>%    
  mutate(
    fixation_time = factor(
      fixation_time,
      levels = c("DAY7", "DAY14", "DAY28")
    )
  )

p1 <- ggplot(
  novel_sot_fixation_counts,
  aes(
    x = fixation_time,
    y = n
  )
) +
  geom_col(fill = "#2C7FB8", colour = "black") +
  theme_bw() +
  labs(
    x = "Time of consensus-crossing",
    y = "Number of mutations"
  )


novel_sot_recurrent_fixation <- novel_sot_fixation_table %>%
  filter(!is.na(fixation_time)) %>%
  count(AA_change_f) %>%
  arrange(desc(n))

p2 <- ggplot(
  novel_sot_recurrent_fixation,
  aes(
    x = reorder(AA_change_f, n),
    y = n
  )
) +
  geom_col(fill = "#2C7FB8", colour = "black") +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Mutation",
    y = "Patients reaching consensus"
  )


p1 + p2 +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "A")




######################################################################################
### 3) Non-escape mutations 
######################################################################################

non_sotrovimab_mutations <- cohort %>% filter(
   !AA_change_f %in% sotrovimab_mutations$AA_change_f,
   !AA_change_f %in% novel_sotrovimab_mutations$AA_change_f
 ) 


############### Recurrent Non-escape Mutations #######################################
ggplot(
  non_sotrovimab_mutations %>% distinct(Filter_Name, GENE_f, AA_change_f) %>%   # one occurrence per patient
    count(GENE_f, AA_change_f, name = "n_patients") %>%
    arrange(desc(n_patients)) %>% filter(n_patients >= 2),
  aes(x = reorder(AA_change_f, n_patients),
      y = n_patients,
      fill = GENE_f)
) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 20, by = 1)) +
  labs(
    x = "Mutation",
    y = "Number of patients",
    title = "Recurrent Non-Escape Mutations"
  ) +
  theme_classic()



################ Non-escape Mutations reaching Consensus #############################
non_sot_fixation_table <- non_sotrovimab_mutations %>% filter(isRep == "NO") %>%
  group_by(Filter_Name, AA_change_f, Filter_SendRef, GENE_f) %>%
  summarise(
    ALT_FREQ = max(ALT_FREQ, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Filter_SendRef,
    values_from = ALT_FREQ
  )

non_sot_fixation_table <- non_sot_fixation_table %>%
  left_join(sample_matrix, by = "Filter_Name")

non_sot_fixation_table <- non_sot_fixation_table %>%
  mutate(
    DAY7.x  = ifelse(DAY7.y  & is.na(DAY7.x),  0, DAY7.x),
    DAY14.x = ifelse(DAY14.y & is.na(DAY14.x), 0, DAY14.x),
    DAY28.x = ifelse(DAY28.y & is.na(DAY28.x), 0, DAY28.x)
  )


non_sot_fixation_table <- non_sot_fixation_table %>%
  mutate(
    fixation_time = case_when(
      
      DAY7.x >= 0.5 &
        (is.na(DAY14.x) | DAY14.x >= 0.5) &
        (is.na(DAY28.x) | DAY28.x >= 0.5) ~ "DAY7",
      
      (is.na(DAY7.x) | DAY7.x < 0.5) &
        DAY14.x >= 0.5 &
        (is.na(DAY28.x) | DAY28.x >= 0.5) ~ "DAY14",
      
      (is.na(DAY7.x) | DAY7.x < 0.5) &
        (is.na(DAY14.x) | DAY14.x < 0.5) &
        DAY28.x >= 0.5 ~ "DAY28",
      
      TRUE ~ NA_character_
    )
  )

non_sot_fixation_table <- non_sot_fixation_table %>% filter(!is.na(fixation_time) )

non_sot_fixation_counts <- non_sot_fixation_table %>%
  count(fixation_time) %>%    
  mutate(
    fixation_time = factor(
      fixation_time,
      levels = c("DAY7", "DAY14", "DAY28")
    )
  )

p3 <- ggplot(
  non_sot_fixation_counts,
  aes(
    x = fixation_time,
    y = n
  )
) +
  geom_col() +
  theme_bw() +
  labs(
    x = "Time of consensus-crossing",
    y = "Number of mutations"
  )


non_sot_recurrent_fixation <- non_sot_fixation_table %>%
  filter(!is.na(fixation_time)) %>%
  count(AA_change_f, GENE_f) %>%
  arrange(desc(n))

p4 <- ggplot(
  non_sot_recurrent_fixation,
  aes(
    x = reorder(AA_change_f, n),
    y = n,
    fill = GENE_f
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Mutation",
    y = "Patients reaching consensus",
  )


p3 + p4 +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "A")
