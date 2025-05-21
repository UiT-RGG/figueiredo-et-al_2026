# Core tidy tools
library(ggplot2)     # plotting
library(dplyr)       # data wrangling
library(tidyr)       # data reshaping
library(readxl)      # reading Excel files
library(readr)       # reading CSV/TSV files
library(purrr)       # functional mapping
library(forcats)     # factor reordering (fct_relevel)

# Statistics
library(rstatix)     # anova_test, tukey_hsd, kruskal_test, etc.
library(dunn.test)   # dunn_test (post hoc for Kruskal-Wallis)
library(effectsize)  # computing Cohen's d

# Plotting
library(ggpubr)      # ggbarplot and plot annotations
library(patchwork)   # combining plots
library(scales)      # axis scaling and percent formatting
library(MoMAColors)  # custom color palettes

# Specific Growth Rate ----
setwd('path/to/data/')

excel_sheets('sgr_data.xlsx')

immunization_sgr <- read_excel('sgr_data.xlsx', sheet = 'SGR_immunization')
immunization_sgr

challenge_sgr <- read_excel('sgr_data.xlsx', sheet = 'SGR_challenge')

immunization_sgr <- immunization_sgr %>%
  mutate(
    delta = as.numeric(difftime(sampling_date, start, units = 'days')),
    sgr = (log(sampling_w) - log(starting_w)) / delta * 100
  )

challenge_sgr <- challenge_sgr %>%
  mutate(
    delta = as.numeric(difftime(sampling_date, start, units = 'days')),
    sgr = (log(sampling_w) - log(challenge_w)) / delta * 100
  )

df1 <- immunization_sgr %>% dplyr::select(treatment, timepoint, assay, sgr)
df2 <- challenge_sgr %>% dplyr::select(treatment, timepoint, assay, sgr)

df_sgr <- bind_rows(df1, df2)

str(df_sgr)

df_sgr <- df_sgr %>%
  mutate(treatment = fct_relevel(treatment, 'conu', 'ptag', 'ivld', 'eomes', 'gata')) %>%
  mutate(timepoint = fct_relevel(timepoint, '10WPI', '4WPC', '10WPC'))

df_sgr %>%
  group_by(timepoint, treatment) %>%
  group_split() %>%
  lapply(function(df) {
    shapiro_result <- shapiro.test(df$sgr)
    data.frame(
      treatment = unique(df$treatment),
      p_value = shapiro_result$p.value
    )
  }) %>%
  bind_rows()  # data normally distributed

# Grouped One-way ANOVA test
df_sgr %>%
  group_by(timepoint) %>%
  anova_test(sgr ~ treatment)

# Tukey's HSD multiple comparisons
stat.test_sgr <- df_sgr %>%
  group_by(timepoint) %>%
  tukey_hsd(sgr ~ treatment) %>%
  add_y_position() %>%
  add_x_position(x = 'timepoint', dodge = 0.8) %>%
  filter(p.adj < 0.05)

stat.test_sgr$y.position <- c(1.45, 1.55, 1.65, 1.75, 1.20, 1.30)  # setting bracket positions manually since stat_pvalue_manual doesn't always figure this out correctly

# Plotting SGR
sgr_plot <- df_sgr %>%
  mutate(treatment = fct_relevel(treatment, 'conu', 'ptag', 'ivld', 'eomes', 'gata')) %>%
  ggbarplot(
    x = 'timepoint',
    y = 'sgr',
    add = 'mean_se',
    fill = 'treatment',
    position = position_dodge(.8),
    error.plot = 'upper_pointrange',
    title = 'Specific Growth Rate',
    ylab = 'Specific Growth Rate (%/day)'
  ) +
  scale_fill_manual(
    values = moma.colors('Rattner'),
    name = 'Treatment',
    labels = c('CONU', 'Plasmid', 'IV-LD', 'EOMES', 'GATA3')
  ) +
  geom_vline(
    xintercept = c(1.5),
    color = 'grey',
    linetype = 'dotdash',
    lwd = .2
  ) +
  geom_hline(
    yintercept = c(0),
    color = 'black',
    linetype = 'dotdash',
    lwd = .1
  ) +
  theme(
    legend.position = 'none',
    plot.margin = unit(c(0.8, 0, .3, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial')
  ) +
  labs(x = NULL) +
  stat_pvalue_manual(
    stat.test_sgr,
    # p-values dataframe
    label = 'p adj. = {p.adj}',
    # Add 'p adj. = ' before the p-value
    digits = 3,
    size = 3,
    tip.length = .005,
    vjust = -0.2,
    # increase separation between bracket and p-value
    step.increase = 0.005  # separation between overlapping brackets
  ) +
  scale_x_discrete(labels = c('Immunization - 10WPI', 'Challenge - 4WPC', 'Challenge - 10WPC')) +
  scale_y_continuous(breaks = seq(0, 2, by = .25))

# Adding n to plot
count_df <- df_sgr %>%
  count(timepoint, treatment)

sgr_plot <- sgr_plot +
  geom_text(
    data = count_df,
    aes(
      x = timepoint,
      y = 0,
      label = paste0('n = ', n),
      group = treatment
    ),
    position = position_dodge(width = 0.8),
    vjust = 2,
    size = 3,
    family = 'Arial',
    color = 'black'
  )

print(sgr_plot)

ggsave(
  filename = 'path/to/output/PhD/Papers/Paper II/data/definitive_plots/sgr_plot_FINAL.jpg',
  width = 3700,
  height = 2500,
  dpi = 300,
  units = 'px'
)

# Cohen's d ----
# get all timepoints
timepoints <- unique(df_sgr$timepoint)

# calculate all pairwise effect sizes
effect_sizes <- map_df(timepoints, function(tp) {
  df <- df_sgr %>% filter(timepoint == tp)

  # all treatment pairs
  treatment_pairs <- combn(unique(df$treatment), 2, simplify = FALSE)

  # compute Cohen's d per pair
  map_df(treatment_pairs, function(pair) {
    df_sub <- df %>% filter(treatment %in% pair)
    df_sub <- droplevels(df_sub)  # clean factor levels

    d_result <- cohens_d(sgr ~ treatment, data = df_sub, pooled_sd = TRUE)

    tibble(
      timepoint = tp,
      group1 = pair[1],
      group2 = pair[2],
      effsize = d_result$Cohens_d,
      CI_low = d_result$CI_low,
      CI_high = d_result$CI_high
    )
  })
})

# Setting effect size limits
effect_sizes <- effect_sizes %>%
  mutate(
    interpretation = case_when(
      abs(effsize) < 0.2 ~ 'negligible',
      abs(effsize) < 0.5 ~ 'small',
      abs(effsize) < 0.8 ~ 'medium',
      TRUE ~ 'large'
    )
  )

effect_sizes %>% print(n = Inf)

# Compare treatments
effect_sizes <- effect_sizes %>%
  mutate(comparison = paste(group1, 'vs', group2))

# Formatting table
effect_sizes_formatted <- effect_sizes %>%
  mutate(
    comparison = paste(group1, 'vs', group2),
    period = factor(timepoint, levels = c('10WPI', '4WPC', '10WPC')),
    interpretation = factor(interpretation, levels = c('negligible', 'small', 'medium', 'large')),
    comparison = fct_reorder(comparison, effsize)
  )


# ELISA ----
elisa_data <-
  read.csv(
    'path/to/data/ELISA_data.tsv',
    header = T,
    check.names = F,
    sep = '\t')

elisa_data <- elisa_data %>%
  mutate(Treatment = fct_relevel(Treatment, 'CONU', 'Plasmid', 'IVLD', 'EOMES', 'GATA3')) %>%
  mutate(Sampling = fct_relevel(Sampling, '10WPI', '6WPC'))

elisa_summary <- elisa_data %>%
  group_by(Sampling, Treatment) %>%
  summarise(
    Average_a = mean(Absorbance),
    SEM_a = sd(Absorbance) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(total_sem = Average_a + SEM_a)

# Test normality
elisa_data %>%
  group_by(Sampling, Treatment) %>%
  group_split() %>%
  lapply(function(df) {
    shapiro_result <- shapiro.test(df$Absorbance)
    data.frame(
      Sampling = unique(df$Sampling),
      Treatment = unique(df$Treatment),
      p_value = shapiro_result$p.value
    )
  }) %>%
  bind_rows()  # data not normally distributed

# Grouped Kruskal-Wallis test
elisa_data %>%
  group_by(Sampling) %>%
  kruskal_test(Absorbance ~ Treatment)

# Dunn's multiple comparisons
stat.test_elisa <- elisa_data %>%
  group_by(Sampling) %>%
  dunn_test(Absorbance ~ Treatment)

stat.test_elisa <-
  stat.test_elisa %>% add_xy_position(x = 'Sampling', dodge = .8, scales = 'fixed')  # adding position of significance values
stat.test_elisa <- stat.test_elisa %>% filter(p.adj < .05)  # filtering only p < .05
stat.test_elisa$p.adj <- round(stat.test_elisa$p.adj, 3)  # rounding values

elisa_summary %>% filter(Sampling == '10WPI') %>% slice_max(total_sem)
elisa_summary %>% filter(Sampling == '6WPC') %>% slice_max(total_sem)

stat.test_elisa$y.position <- c(.36, .40, .32, .45, .53, .56)

# plotting ELISA
elisa_plot <- elisa_data %>%
  mutate(Treatment = fct_relevel(Treatment, 'CONU', 'Plasmid', 'IVLD', 'EOMES', 'GATA3')) %>%
  ggbarplot(
    x = 'Sampling',
    y = 'Absorbance',
    add = 'mean_se',
    fill = 'Treatment',
    position = position_dodge(.8),
    error.plot = 'upper_pointrange',
    title = 'Enzyme-Linked Immunosorbent Assay',
    ylab = 'Absorbance at 492 nm'
  ) +
  scale_fill_manual(
    values = moma.colors('Rattner'),
    name = 'Treatment',
    labels = c('CONU', 'Plasmid', 'IV-LD', 'EOMES', 'GATA3')
  ) +
  geom_vline(
    xintercept = c(1.5),
    color = 'grey',
    linetype = 'dotdash',
    lwd = .2
  ) +
  geom_hline(
    yintercept = c(0),
    color = 'black',
    linetype = 'dotdash',
    lwd = .1
  ) +
  theme(
    legend.position = 'right',
    plot.margin = unit(c(0.8, 0, .3, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial')
  ) +
  labs(x = NULL) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1)) +
  stat_pvalue_manual(
    stat.test_elisa,
    # p-values data frame (stat.test)
    label = 'p adj. = {p.adj}',
    # Add 'p adj. = ' before the p-value
    digits = 3,
    size = 3,
    tip.length = .005,
    vjust = -0.2,
    step.increase = 0.01  
  )

count_df <- elisa_data %>%
  count(Sampling, Treatment)

elisa_plot <- elisa_plot +
  geom_text(
    data = count_df,
    aes(
      x = Sampling,
      y = 0,
      label = paste0('n = ', n),
      group = Treatment
    ),
    position = position_dodge(width = 0.8),
    vjust = 2,
    size = 3,
    family = 'Arial',
    color = 'black'
  )

ggsave(
  filename = 'path/to/output/PhD/Papers/Paper II/data/definitive_plots/ELISA_plot_FINAL.jpg',
  width = 3700,
  height = 2500,
  dpi = 300,
  units = 'px'
)



# VNA ----
vna_data <-
  read.csv(
    'path/to/data/VNA_data.tsv',
    header = T,
    check.names = F,
    sep = '\t'
  )

vna_data <- vna_data %>%
  mutate(Treatment = fct_relevel(Treatment, 'Plasmid', 'IVLD', 'EOMES', 'GATA3')) %>%
  mutate(Sampling = fct_relevel(Sampling, '10WPI' , '6WPC'))

vna_summary <- vna_data %>%
  group_by(Treatment, Sampling) %>%
  summarise(
    mean_dilution = round(mean(Dilution), 3),
    n = n(),
    sem = round(sd(Dilution) / sqrt(n), 3),
    .groups = 'drop'
  ) %>%
  mutate(total_sem = mean_dilution + sem) %>%
  arrange(Sampling)


# Test normality
vna_data %>%
  group_by(Sampling, Treatment) %>%
  group_split() %>%
  lapply(function(df) {
    shapiro_result <- shapiro.test(df$Dilution)
    data.frame(
      Sampling = unique(df$Sampling),
      Treatment = unique(df$Treatment),
      p_value = number(shapiro_result$p.value, accuracy = 0.00001)  # data not normally distributed
    )
  }) %>%
  bind_rows()

# Grouped Kruskal-Wallis test
vna_data %>%
  group_by(Sampling) %>%
  kruskal_test(Dilution ~ Treatment)

# Dunn's multiple comparisons
stat.test_vna <- vna_data %>%
  group_by(Sampling) %>%
  dunn_test(Dilution ~ Treatment)

stat.test_vna <-
  stat.test_vna %>% add_xy_position(x = 'Sampling', dodge = .8, scales = 'fixed')  # adding position of significance values
stat.test_vna <- stat.test_vna %>% filter(p.adj < .05)  # filtering only p < .05
stat.test_vna$p.adj <- round(stat.test_vna$p.adj, 3)  # rounding values

vna_summary %>% filter(Sampling == '10WPI') %>% slice_max(total_sem)
vna_summary %>% filter(Sampling == '6WPC') %>% slice_max(total_sem)

stat.test_vna$y.position <- c(85)

rattner_colors <- moma.colors('Rattner')
colors <- moma.colors('Rattner', 4)
rattner_colors_4 <- rattner_colors[c(2, 3, 4, 5)]

# plotting VNA
vna_plot <- vna_data %>%
  mutate(Treatment = fct_relevel(Treatment, 'Plasmid', 'IVLD', 'EOMES', 'GATA3')) %>%
  ggbarplot(
    x = 'Sampling',
    y = 'Dilution',
    add = 'mean_se',
    fill = 'Treatment',
    position = position_dodge(.8),
    error.plot = 'upper_pointrange',
    title = 'Virus Neutralization Assay',
    ylab = 'Dilution 1:Y'
  ) +
  scale_fill_manual(
    values = rattner_colors_4,
    name = 'Treatment',
    labels = c('Plasmid', 'IV-LD', 'EOMES', 'GATA3')
  ) +
  geom_vline(
    xintercept = c(1.5),
    color = 'grey',
    linetype = 'dotdash',
    lwd = .2
  ) +
  geom_hline(
    yintercept = c(0),
    color = 'black',
    linetype = 'dotdash',
    lwd = .1
  ) +
  theme(
    legend.position = 'right',
    plot.margin = unit(c(0.8, 0, .3, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial')
  ) +
  labs(x = NULL) +
  scale_y_continuous(breaks = seq(0, 250, by = 50)) +
  stat_pvalue_manual(
    stat.test_vna,
    label = 'p adj. = {p.adj}',
    digits = 3,
    size = 3,
    tip.length = .005,
    vjust = -0.2,
    step.increase = 0.01
  )

count_df <- vna_data %>%
  count(Sampling, Treatment)

vna_plot +
  geom_text(
    data = count_df,
    aes(
      x = Sampling,
      y = 0,
      label = paste0('n = ', n),
      group = Treatment
    ),
    position = position_dodge(width = 0.8),
    vjust = 2,
    size = 3,
    family = 'Arial',
    color = 'black'
  )

ggsave(
  filename = 'path/to/output/PhD/Papers/Paper II/data/definitive_plots/VNA_plot_FINAL.jpg',
  width = 3700,
  height = 2500,
  dpi = 300,
  units = 'px'
)
# Histoscore ----
histoscore_data <- readxl::read_xlsx('path/to/data/histoscore_data_updated.xlsx')

histoscore_data <- histoscore_data %>%
  mutate(Treatment = factor(Treatment, levels = c('CONU', 'Plasmid', 'IV-LD', 'EOMES', 'GATA3'))) %>%
  mutate(Sampling = factor(Sampling, levels = c('4WPC', '6WPC', '10WPC'))) %>%
  mutate(Tissue = factor(Tissue, levels = c('Pancreas', 'Heart')))


histoscore_summary <- histoscore_data %>%
  group_by(Tissue, Sampling, Treatment) %>%
  summarise(
    mean_histoscore = round(mean(Histoscore), 3),
    n = n(),
    sem = round(sd(Histoscore) / sqrt(n), 3),
    .groups = 'drop'
  ) %>%
  mutate(total_sem = mean_histoscore + sem) %>%
  arrange(Sampling)

# Test normality
histoscore_data %>%
  group_by(Tissue, Sampling, Treatment) %>%
  group_split() %>%
  lapply(function(df) {
    shapiro_result <- shapiro.test(df$Histoscore)
    data.frame(
      Tissue = unique(df$Tissue),
      Sampling = unique(df$Sampling),
      Treatment = unique(df$Treatment),
      p_value = number(shapiro_result$p.value, accuracy = 0.00000001)
    )
  }) %>%
  bind_rows()


# Grouped Kruskal-Wallis test
histoscore_data %>%
  filter(Sampling == '4WPC') %>%
  group_by(Sampling) %>%
  kruskal_test(Histoscore ~ Treatment)

# Dunn's multiple comparisons
stat.test_histoscore4WPC <- histoscore_data %>%
  filter(Sampling == '4WPC') %>%
  group_by(Tissue) %>%
  dunn_test(Histoscore ~ Treatment)


stat.test_histoscore4WPC <-
  stat.test_histoscore4WPC %>% add_xy_position(x = 'Tissue', dodge = .8, scales = 'fixed')  # adding position of significance values
stat.test_histoscore4WPC <- stat.test_histoscore4WPC %>% filter(p.adj < .05)  # filtering only p < .05
stat.test_histoscore4WPC$p.adj <- round(stat.test_histoscore4WPC$p.adj, 5)  # rounding values

histoscore_summary %>% filter(Tissue == 'Pancreas') %>% slice_max(total_sem)
histoscore_summary %>% filter(Tissue == 'Heart') %>% slice_max(total_sem)

stat.test_histoscore4WPC$y.position <- c(2.7, 2.85)

histoscore_plot <- histoscore_data %>%
  mutate(Treatment = fct_relevel(Treatment, 'CONU', 'Plasmid', 'IV-LD', 'EOMES', 'GATA3')) %>%
  ggbarplot(
    x = 'Tissue',
    y = 'Histoscore',
    add = 'mean_se',
    fill = 'Treatment',
    position = position_dodge(.8),
    error.plot = 'upper_pointrange',
    title = 'Histoscore at 4 weeks post-challenge',
    ylab = 'Histoscore'
  ) +
  scale_fill_manual(
    values = moma.colors('Rattner'),
    name = 'Treatment',
    labels = c('CONU', 'Plasmid', 'IV-LD', 'EOMES', 'GATA3')
  ) +
  geom_vline(
    xintercept = c(1.5),
    color = 'grey',
    linetype = 'dotdash',
    lwd = .2
  ) +
  geom_hline(
    yintercept = c(0),
    color = 'black',
    linetype = 'dotdash',
    lwd = .1
  ) +
  theme(
    legend.position = 'right',
    plot.margin = unit(c(0.8, 0, .3, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial')
  ) +
  labs(x = NULL) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  stat_pvalue_manual(
    stat.test_histoscore4WPC,
    label = 'p adj. = {p.adj}',
    digits = 3,
    size = 3,
    tip.length = .005,
    vjust = -0.2,
    step.increase = 0.01
  )

count_df <- histoscore_data %>%
  count(Tissue, Treatment)

histoscore_plot +
  geom_text(
    data = count_df,
    aes(
      x = Tissue,
      y = 0,
      label = paste0('n = ', n),
      group = Treatment,
      fill = Treatment
    ),
    position = position_dodge(width = 0.8),
    vjust = 2,
    size = 3,
    family = 'Arial',
    color = 'black'
  )


ggsave(
  filename = 'path/to/data/HISTOSCORE4WPC_plot_FINAL.jpg',
  width = 3700,
  height = 2500,
  dpi = 300,
  units = 'px'
)

# histoscore stacked plots
histoscore_data <- readxl::read_xlsx('histoscore_data_updated.xlsx')

histoscore_data <- histoscore_data %>% 
  mutate(Treatment = factor(Treatment, levels = c('CONU', 'Plasmid', 'IV-LD', 'EOMES', 'GATA3'))) %>% 
  mutate(Sampling = factor( Sampling, levels = c('4WPC', '6WPC', '10WPC')))

levels(histoscore_data$Treatment)
kruskal.test(Histoscore ~ Treatment, data = dplyr::filter(histoscore_data, Tissue == 'Pancreas'))

histoscore_data %>%
  group_by(Treatment, Sampling) %>%
  summarise(max_value = max(Histoscore))  # checking max histoscores per treatment, per tissue

# Add a column to classify Histoscore
histoscore_data$Histoscore_Classification <- cut(
  histoscore_data$Histoscore,
  breaks = c(-Inf, 0.5, 1.5, 2.5, 3),
  labels = c("Low", "mediumLow", "mediumHigh", "High"),
  right = TRUE  # Intervals are inclusive on the right
)

histoscore_data_pancreas <- histoscore_data %>% dplyr::filter(., Tissue == 'Pancreas') %>% 
  group_by(Treatment, Sampling, Histoscore_Classification) %>% 
  summarise(
    Count = n(),
    .groups = 'drop'
  )

data_percentage_pancreas <- histoscore_data_pancreas %>%
  group_by(Treatment, Sampling) %>%
  mutate(Total = sum(Count),                          # Calculate the total Count per group
         Percentage = (Count / Total) * 100) %>%      # Calculate the percentage
  ungroup() 
data_percentage_pancreas %>% print(n = Inf)

# Colorblind-friendly traffic light-inspired palette
traffic_palette <- c("#009E60", "#F0E442", "#E69F00", "#D32F2F")

data_percentage_pancreas <- data_percentage_pancreas %>%
  mutate(Treatment = factor(
    Treatment,
    levels = c('CONU', 'Plasmid', 'IV-LD', 'EOMES', 'GATA3')
  ))

# Stacked barplot 4 WPC
pancreas_4wpc <- data_percentage_pancreas %>%
  filter(Sampling == '4WPC') %>%
  ggplot(aes(fill = Histoscore_Classification, y = Percentage, x = Treatment)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_cartesian(clip = 'off') +
  scale_fill_manual(
    name = 'Histoscore',
    labels = c('≤0.5', '0.5 - 1.5', '1.5 - 2.5', '2.5 - 3'),
    values = traffic_palette
  ) +
  # Add total sample size (n = X) on top of each stack
  geom_text(
    data = data_percentage_pancreas %>%
      filter(Sampling == '4WPC') %>%
      distinct(Treatment, Sampling, Total), # Use unique data for the labels
    aes(x = Treatment, y = 1, label = paste0("n = ", Total)), # Avoid inheriting unwanted aesthetics
    inherit.aes = F,
    vjust = -0.5, # Slightly above the stack
    size = 3,
    family = 'Arial'
  ) +
  theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = 'right',
    plot.margin = unit(c(.5, .5, .1, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial')
  ) +
  scale_y_continuous(labels = scales::percent) +
  ylab('Percentage of individuals') +
  labs(title = 'Pancreas \n 4 weeks post-challenge',
       tag = 'a') +
  theme(plot.title = element_text(hjust = .5),
        plot.tag.position = c(0.07, .99))

# Stacked barplot 6 WPC
pancreas_6wpc <- data_percentage_pancreas %>%
  filter(Sampling == '6WPC') %>%
  ggplot(aes(fill = Histoscore_Classification, y = Percentage, x = Treatment)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_cartesian(clip = 'off') +
  scale_fill_manual(
    name = 'Histoscore',
    labels = c('≤0.5', '0.5 - 1.5', '1.5 - 2.5', '2.5 - 3'),
    values = traffic_palette # Replace with your desired color palette
  ) +
  # Add total sample size (n = X) on top of each stack
  geom_text(
    data = data_percentage_pancreas %>%
      filter(Sampling == '6WPC') %>%
      distinct(Treatment, Sampling, Total), # Use unique data for the labels
    aes(x = Treatment, y = 1, label = paste0("n = ", Total)), # Avoid inheriting unwanted aesthetics
    inherit.aes = F,
    vjust = -0.5, # Slightly above the stack
    size = 3,
    family = 'Arial'
  ) +
  theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = 'right',
    plot.margin = unit(c(.8, .5, .1, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial')
  ) +
  scale_y_continuous(labels = scales::percent) +
  ylab('Percentage of individuals') +
  labs(title = 'Pancreas \n 6 weeks post-challenge',
       tag = 'c') +
  theme(plot.title = element_text(hjust = .5),
        plot.tag.position = c(0.07, .99))

# Stacked barplot 10 WPC
pancreas_10wpc <- data_percentage_pancreas %>%
  filter(Sampling == '10WPC') %>%
  ggplot(aes(fill = Histoscore_Classification, y = Percentage, x = Treatment)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_cartesian(clip = 'off') +
  scale_fill_manual(
    name = 'Histoscore',
    labels = c('≤0.5', '0.5 - 1.5', '1.5 - 2.5', '2.5 - 3'),
    values = traffic_palette # Replace with your desired color palette
  ) +
  # Add total sample size (n = X) on top of each stack
  geom_text(
    data = data_percentage_pancreas %>%
      filter(Sampling == '10WPC') %>%
      distinct(Treatment, Sampling, Total), # Use unique data for the labels
    aes(x = Treatment, y = 1, label = paste0("n = ", Total)), # Avoid inheriting unwanted aesthetics
    inherit.aes = F,
    vjust = -0.5, # Slightly above the stack
    size = 3,
    family = 'Arial'
  ) +
  theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = 'right',
    plot.margin = unit(c(.8, .5, .1, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial')
  ) +
  scale_y_continuous(labels = scales::percent) +
  ylab('Percentage of individuals') +
  labs(title = 'Pancreas \n 10 weeks post-challenge',
       tag = 'e') +
  theme(plot.title = element_text(hjust = .5),
        plot.tag.position = c(0.07, .99))

## Heart STACKED plots HISTOSCORE
histoscore_data_heart <- histoscore_data %>% dplyr::filter(., Tissue == 'Heart') %>% 
  group_by(Treatment, Sampling, Histoscore_Classification) %>% 
  summarise(
    Count = n(),
    .groups = 'drop'
  )

data_percentage_heart <- histoscore_data_heart %>%
  group_by(Treatment, Sampling) %>%
  mutate(Total = sum(Count),                          # Calculate the total Count per group
         Percentage = (Count / Total) * 100) %>%      # Calculate the percentage
  ungroup() 

data_percentage_heart %>% print(n = Inf)

data_percentage_heart <- data_percentage_heart %>%
  mutate(Treatment = factor(
    Treatment,
    levels = c('CONU', 'Plasmid', 'IV-LD', 'EOMES', 'GATA3')
  ))

# Stacked barplot 4 WPC
heart_4wpc <- data_percentage_heart %>%
  filter(Sampling == '4WPC') %>%
  ggplot(aes(fill = Histoscore_Classification, y = Percentage, x = Treatment)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_cartesian(clip = 'off') +
  scale_fill_manual(
    name = 'Histoscore',
    labels = c('≤0.5', '0.5 - 1.5', '1.5 - 2.5', '2.5 - 3'),
    values = traffic_palette
  ) +
  # Add total sample size (n = X) on top of each stack
  geom_text(
    data = data_percentage_heart %>%
      filter(Sampling == '4WPC') %>%
      distinct(Treatment, Sampling, Total),
    aes(x = Treatment, y = 1, label = paste0("n = ", Total)),
    inherit.aes = F,
    vjust = -0.5,
    size = 3,
    family = 'Arial'
  ) +
  theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = 'right',
    plot.margin = unit(c(.8, .5, .1, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial')
  ) +
  scale_y_continuous(labels = scales::percent) +
  ylab('Percentage of individuals') +
  labs(title = 'Heart \n 4 weeks post-challenge',
       tag = 'b') +
  theme(plot.title = element_text(hjust = .5),
        plot.tag.position = c(0.07, .99))

# Stacked barplot 6 WPC
heart_6wpc <- data_percentage_heart %>%
  filter(Sampling == '6WPC') %>%
  ggplot(aes(fill = Histoscore_Classification, y = Percentage, x = Treatment)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_cartesian(clip = 'off') +
  scale_fill_manual(
    name = 'Histoscore',
    labels = c('≤0.5', '0.5 - 1.5', '1.5 - 2.5', '2.5 - 3'),
    values = traffic_palette # Replace with your desired color palette
  ) +
  # Add total sample size (n = X) on top of each stack
  geom_text(
    data = data_percentage_heart %>%
      filter(Sampling == '6WPC') %>%
      distinct(Treatment, Sampling, Total), # Use unique data for the labels
    aes(x = Treatment, y = 1, label = paste0("n = ", Total)), # Avoid inheriting unwanted aesthetics
    inherit.aes = F,
    vjust = -0.5, # Slightly above the stack
    size = 3,
    family = 'Arial'
  ) +
  theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = 'right',
    plot.margin = unit(c(.8, .5, .1, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial'),
  ) +
  scale_y_continuous(labels = scales::percent) +
  ylab('Percentage of individuals') +
  labs(title = 'Heart \n 6 weeks post-challenge',
       tag = 'd') +
  theme(plot.title = element_text(hjust = .5),
        plot.tag.position = c(0.07, .99))

# Stacked barplot 10 WPC
heart_10wpc <- data_percentage_heart %>%
  filter(Sampling == '10WPC') %>%
  ggplot(aes(fill = Histoscore_Classification, y = Percentage, x = Treatment)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_cartesian(clip = 'off') +
  scale_fill_manual(
    name = 'Histoscore',
    labels = c('≤0.5', '0.5 - 1.5', '1.5 - 2.5', '2.5 - 3'),
    values = traffic_palette # Replace with your desired color palette
  ) +
  # Add total sample size (n = X) on top of each stack
  geom_text(
    data = data_percentage_heart %>%
      filter(Sampling == '10WPC') %>%
      distinct(Treatment, Sampling, Total), # Use unique data for the labels
    aes(x = Treatment, y = 1, label = paste0("n = ", Total)), # Avoid inheriting unwanted aesthetics
    inherit.aes = F,
    vjust = -0.5, # Slightly above the stack
    size = 3,
    family = 'Arial'
  ) +
  theme_pubr() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = 'right',
    plot.margin = unit(c(.8, .5, .1, .5), 'cm'),
    plot.title = element_text(
      hjust = .5,
      vjust = 6,
      size = 12
    ),
    text = element_text(size = 12, family = 'Arial')
  ) +
  scale_y_continuous(labels = scales::percent) +
  ylab('Percentage of individuals') +
  labs(title = 'Heart \n 10 weeks post-challenge',
       tag = 'f') +
  theme(plot.title = element_text(hjust = .5),
        plot.tag.position = c(0.07, .99))

# Suppress legends, leaving only one
heart_4wpc <- heart_4wpc + guides(fill = "none")
heart_6wpc <- heart_6wpc + guides(fill = "none")
heart_10wpc <- heart_10wpc + guides(fill = "none")
pancreas_4wpc <- pancreas_4wpc + guides(fill = "none")
pancreas_6wpc <- pancreas_6wpc + guides(fill = "legend")
pancreas_10wpc <- pancreas_10wpc + guides(fill = "none")

histoscore_stacked_plots <- (pancreas_4wpc / pancreas_6wpc / pancreas_10wpc | heart_4wpc / heart_6wpc / heart_10wpc) +
  plot_layout(guides = 'collect', axes = 'collect') &
  theme(
    legend.position = 'right',
    legend.justification = 'center',
    legend.box.just = 'center',
    plot.title = element_text(hjust = 0.5, size = 12),
    text = element_text(family = "Arial")
  )

ggsave('path/to/output/PhD/Papers/Paper II/data/definitive_plots/histoscore_stacked_plots.jpg', 
       width = 3700,
       height = 3700,
       dpi = 300,
       units = 'px'
)

# References ----
citation('readxl')    
citation('tidyverse')     
citation('rstatix')   
citation('dunn.test') 
citation('effectsize') 
citation('ggpubr')     
citation('patchwork')  
citation('scales')     
citation('MoMAColors') 
