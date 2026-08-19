
# Statistical analyses and figures
## Pharmaceutics, Figueiredo et al., 2026


# Packages ---------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(rstatix)
library(effectsize)
library(ggpubr)
library(patchwork)
library(scales)
library(MoMAColors)


# Paths ------------------------------------------------------------------------

data_dir <- "data"
figure_dir <- "figures"

dir.create(figure_dir, showWarnings = FALSE)


# Common settings --------------------------------------------------------------

treatment_colors <- moma.colors("Rattner")

treatment_labels <- c(
  "CONU",
  "Plasmid",
  "IV-LD",
  "EOMES",
  "GATA3"
)

theme_manuscript <- theme(
  legend.position = "right",
  plot.margin = unit(c(0.8, 0, 0.3, 0.5), "cm"),
  text = element_text(
    size = 12,
    family = "Times New Roman"
  )
)



# SGR ----
## Data

immunization_sgr <- read_excel(
  file.path(data_dir, "sgr_data.xlsx"),
  sheet = "SGR_immunization"
) %>%
  mutate(
    delta = as.numeric(difftime(sampling_date, start, units = "days")),
    sgr = (log(sampling_w) - log(starting_w)) / delta * 100
  )

challenge_sgr <- read_excel(
  file.path(data_dir, "sgr_data.xlsx"),
  sheet = "SGR_challenge"
) %>%
  mutate(
    delta = as.numeric(difftime(sampling_date, start, units = "days")),
    sgr = (log(sampling_w) - log(challenge_w)) / delta * 100
  )

df_sgr <- bind_rows(
  immunization_sgr %>%
    select(treatment, timepoint, assay, sgr, sampling_w),
  challenge_sgr %>%
    select(treatment, timepoint, assay, sgr, sampling_w)
) %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c("conu", "ptag", "ivld", "eomes", "gata")
    ),
    timepoint = factor(
      timepoint,
      levels = c("10WPI", "4WPC", "10WPC")
    )
  )

## Statistics
### Normality

shapiro_sgr <- df_sgr %>%
  group_by(timepoint, treatment) %>%
  group_modify(~ {
    test <- shapiro.test(.x$sgr)
    
    tibble(
      W = unname(test$statistic),
      p_value = test$p.value
    )
  }) %>%
  ungroup()

shapiro_sgr


### Homogeneity of variances

sgr_levene <- df_sgr %>%
  group_by(timepoint) %>%
  levene_test(sgr ~ treatment)

sgr_levene

### 10WPI: heterogeneous variances

sgr_welch_10WPI <- df_sgr %>%
  filter(timepoint == "10WPI") %>%
  welch_anova_test(sgr ~ treatment)

### 4WPC and 10WPC: homogeneous variances

sgr_anova_equal_var <- df_sgr %>%
  filter(timepoint %in% c("4WPC", "10WPC")) %>%
  group_by(timepoint) %>%
  anova_test(sgr ~ treatment)

sgr_welch_10WPI
sgr_anova_equal_var


### Post hoc comparisons

sgr_games_howell_10WPI <- df_sgr %>%
  filter(timepoint == "10WPI") %>%
  games_howell_test(sgr ~ treatment) %>%
  mutate(
    timepoint = factor(
      "10WPI",
      levels = levels(df_sgr$timepoint)
    )
  )

sgr_tukey_equal_var <- df_sgr %>%
  filter(timepoint %in% c("4WPC", "10WPC")) %>%
  group_by(timepoint) %>%
  tukey_hsd(sgr ~ treatment)

stat.test_sgr <- bind_rows(
  sgr_games_howell_10WPI,
  sgr_tukey_equal_var
) %>%
  filter(p.adj < 0.05) %>%
  mutate(
    p.label = if_else(
      p.adj < 0.001,
      "p < 0.001",
      sprintf("p = %.3f", p.adj)
    )
  ) %>%
  add_x_position(
    x = "timepoint",
    dodge = 0.8
  ) %>%
  mutate(
    y.position = case_when(
      group1 == "ptag" & group2 == "eomes" ~ 1.20,
      group1 == "ptag" & group2 == "gata"  ~ 1.30,
      TRUE ~ NA_real_
    ),
    xmin = case_when(
      group1 == "ptag" & group2 == "eomes" ~ 2.84,
      group1 == "ptag" & group2 == "gata"  ~ 2.84,
      TRUE ~ NA_real_
    ),
    xmax = case_when(
      group1 == "ptag" & group2 == "eomes" ~ 3.16,
      group1 == "ptag" & group2 == "gata"  ~ 3.32,
      TRUE ~ NA_real_
    )
  )

stat.test_sgr


### Plot

sgr_plot <- ggbarplot(
  data = df_sgr,
  x = "timepoint",
  y = "sgr",
  add = "mean_se",
  fill = "treatment",
  position = position_dodge(0.8),
  error.plot = "upper_pointrange",
  ylab = "Specific Growth Rate (%/day)"
) +
  scale_fill_manual(
    values = treatment_colors,
    name = "Treatment",
    labels = treatment_labels
  ) +
  geom_vline(
    xintercept = 1.5,
    color = "grey",
    linetype = "dashed",
    linewidth = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey10",
    linewidth = 0.1
  ) +
  stat_pvalue_manual(
    stat.test_sgr,
    label = "p.label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    size = 3,
    tip.length = 0.005,
    step.increase = 0,
    bracket.nudge.y = 0,
    family = "Times New Roman",
    inherit.aes = FALSE
  ) +
  scale_x_discrete(
    labels = c(
      "Immunization - 10WPI",
      "Challenge - 4WPC",
      "Challenge - 10WPC"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 2, by = 0.5)
  ) +
  labs(x = NULL) +
  theme_manuscript

print(sgr_plot)

ggsave(
  file.path(figure_dir, "SGR_plot.jpg"),
  plot = sgr_plot,
  width = 9,
  height = 5,
  dpi = 600
)


### Effect size

timepoints <- levels(df_sgr$timepoint)

effect_sizes <- map_dfr(timepoints, function(tp) {
  
  df_tp <- df_sgr %>%
    filter(timepoint == tp) %>%
    drop_na(sgr, treatment) %>%
    droplevels()
  
  treatment_pairs <- combn(
    levels(df_tp$treatment),
    2,
    simplify = FALSE
  )
  
  map_dfr(treatment_pairs, function(pair) {
    
    df_sub <- df_tp %>%
      filter(treatment %in% pair) %>%
      droplevels()
    
    d_result <- effectsize::cohens_d(
      sgr ~ treatment,
      data = df_sub,
      pooled_sd = FALSE,
      ci = 0.95
    )
    
    tibble(
      timepoint = tp,
      group1 = pair[1],
      group2 = pair[2],
      effect_size = d_result$Cohens_d,
      CI_low = d_result$CI_low,
      CI_high = d_result$CI_high
    )
  })
}) %>%
  mutate(
    interpretation = case_when(
      abs(effect_size) < 0.2 ~ "negligible",
      abs(effect_size) < 0.5 ~ "small",
      abs(effect_size) < 0.8 ~ "medium",
      TRUE ~ "large"
    ),
    comparison = paste(group1, "vs", group2),
    period = factor(
      timepoint,
      levels = c("10WPI", "4WPC", "10WPC")
    )
  )

effect_sizes %>%
  print(n = Inf)


# ELISA ----
## Data

elisa_data <- read.delim(
  file.path(data_dir, "ELISA_data.tsv"),
  check.names = FALSE
) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("CONU", "Plasmid", "IVLD", "EOMES", "GATA3")
    ),
    Sampling = factor(
      Sampling,
      levels = c("10WPI", "6WPC")
    )
  )


## Statistics
### Descriptive stats

elisa_summary <- elisa_data %>%
  group_by(Sampling, Treatment) %>%
  summarise(
    mean_absorbance = mean(Absorbance, na.rm = TRUE),
    sd_absorbance = sd(Absorbance, na.rm = TRUE),
    n = sum(!is.na(Absorbance)),
    sem_absorbance = sd_absorbance / sqrt(n),
    .groups = "drop"
  )

elisa_summary


### Normality

shapiro_elisa <- elisa_data %>%
  group_by(Sampling, Treatment) %>%
  group_modify(~ {
    test <- shapiro.test(.x$Absorbance)
    
    tibble(
      W = unname(test$statistic),
      p_value = test$p.value
    )
  }) %>%
  ungroup()

shapiro_elisa


### Treatment comparisons within sampling points

kruskal_elisa <- elisa_data %>%
  group_by(Sampling) %>%
  kruskal_test(Absorbance ~ Treatment)

kruskal_elisa

dunn_elisa_all <- elisa_data %>%
  group_by(Sampling) %>%
  dunn_test(
    Absorbance ~ Treatment,
    p.adjust.method = "BH",
    detailed = TRUE
  )

dunn_elisa_all


### Between-timepoint comparisons within treatments

elisa_between_timepoints <- elisa_data %>%
  group_by(Treatment) %>%
  wilcox_test(
    Absorbance ~ Sampling,
    detailed = TRUE
  ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

elisa_between_timepoints


### Plot annotations

stat.test.elisa <- dunn_elisa_all %>%
  filter(p.adj < 0.05) %>%
  mutate(
    p.label = if_else(
      p.adj < 0.001,
      "p < 0.001",
      sprintf("p = %.3f", p.adj)
    )
  ) %>%
  add_xy_position(
    x = "Sampling",
    dodge = 0.8,
    scales = "fixed"
  )

elisa_bracket_base <- elisa_data %>%
  group_by(Sampling, Treatment) %>%
  summarise(
    upper_error = mean(Absorbance, na.rm = TRUE) +
      sd(Absorbance, na.rm = TRUE) /
      sqrt(sum(!is.na(Absorbance))),
    .groups = "drop"
  ) %>%
  group_by(Sampling) %>%
  summarise(
    base_y = max(upper_error) + 0.025,
    .groups = "drop"
  )

stat.test.elisa <- stat.test.elisa %>%
  left_join(elisa_bracket_base, by = "Sampling") %>%
  group_by(Sampling) %>%
  arrange(xmax - xmin, .by_group = TRUE) %>%
  mutate(
    y.position = base_y + (row_number() - 1) * 0.045,
    y.position = if_else(
      Sampling == "6WPC",
      y.position - 0.08,
      y.position
    )
  ) %>%
  ungroup()

elisa_ymax <- max(
  stat.test.elisa$y.position,
  na.rm = TRUE
) + 0.035


### Plot

elisa_plot <- ggbarplot(
  data = elisa_data,
  x = "Sampling",
  y = "Absorbance",
  add = "mean_se",
  fill = "Treatment",
  position = position_dodge(0.8),
  error.plot = "upper_pointrange",
  ylab = "Absorbance at 492 nm"
) +
  scale_fill_manual(
    values = treatment_colors,
    name = "Treatment",
    labels = treatment_labels
  ) +
  geom_vline(
    xintercept = 1.5,
    color = "grey",
    linetype = "dashed",
    linewidth = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey10",
    linewidth = 0.1
  ) +
  stat_pvalue_manual(
    stat.test.elisa,
    label = "p.label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    size = 3,
    tip.length = 0.005,
    step.increase = 0,
    bracket.nudge.y = 0,
    family = "Times New Roman",
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      ceiling(elisa_ymax * 10) / 10,
      by = 0.1
    ),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(
    ylim = c(0, elisa_ymax)
  ) +
  labs(x = NULL) +
  theme_manuscript

print(elisa_plot)

ggsave(
  file.path(figure_dir, "ELISA_plot.jpg"),
  plot = elisa_plot,
  width = 9,
  height = 5,
  dpi = 600
)


# VNA ----
## Data

vna_data <- read.delim(
  file.path(data_dir, "VNA_data.tsv"),
  check.names = FALSE
) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Plasmid", "IVLD", "EOMES", "GATA3")
    ),
    Sampling = factor(
      Sampling,
      levels = c("10WPI", "6WPC")
    )
  )


## Statistics

vna_summary <- vna_data %>%
  group_by(Treatment, Sampling) %>%
  summarise(
    mean_dilution = mean(Dilution, na.rm = TRUE),
    n = sum(!is.na(Dilution)),
    sem = sd(Dilution, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

vna_summary


### Normality

shapiro_vna <- vna_data %>%
  group_by(Sampling, Treatment) %>%
  group_modify(~ {
    test <- shapiro.test(.x$Dilution)
    
    tibble(
      W = unname(test$statistic),
      p_value = test$p.value
    )
  }) %>%
  ungroup()

shapiro_vna


### Treatment comparisons within sampling points

kruskal_vna <- vna_data %>%
  group_by(Sampling) %>%
  kruskal_test(Dilution ~ Treatment)

kruskal_vna

dunn_vna_all <- vna_data %>%
  group_by(Sampling) %>%
  dunn_test(
    Dilution ~ Treatment,
    p.adjust.method = "BH",
    detailed = TRUE
  )

dunn_vna_all


### Overall comparison between sampling points

vna_between_sampling <- wilcox.test(
  Dilution ~ Sampling,
  data = vna_data
)

vna_between_sampling


### Plot annotations

stat.test_vna <- dunn_vna_all %>%
  filter(p.adj < 0.05) %>%
  mutate(
    p.label = if_else(
      p.adj < 0.001,
      "p < 0.001",
      sprintf("p = %.3f", p.adj)
    )
  ) %>%
  add_xy_position(
    x = "Sampling",
    dodge = 0.8,
    scales = "fixed"
  )

### Manually positioning significance
stat.test_vna$y.position <- 85

rattner_colors_4 <- treatment_colors[2:5]


### Plot

vna_plot <- ggbarplot(
  data = vna_data,
  x = "Sampling",
  y = "Dilution",
  add = "mean_se",
  fill = "Treatment",
  position = position_dodge(0.8),
  error.plot = "upper_pointrange",
  ylab = "Dilution 1:Y"
) +
  scale_fill_manual(
    values = rattner_colors_4,
    name = "Treatment",
    labels = c("Plasmid", "IV-LD", "EOMES", "GATA3")
  ) +
  geom_vline(
    xintercept = 1.5,
    color = "grey",
    linetype = "dashed",
    linewidth = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey10",
    linewidth = 0.1
  ) +
  stat_pvalue_manual(
    stat.test_vna,
    label = "p.label",
    size = 3,
    tip.length = 0.005,
    step.increase = 0,
    bracket.nudge.y = 0,
    family = "Times New Roman"
  ) +
  scale_y_continuous(
    breaks = seq(0, 250, by = 50)
  ) +
  labs(x = NULL) +
  theme_manuscript

print(vna_plot)

ggsave(
  file.path(figure_dir, "VNA_plot.jpg"),
  plot = vna_plot,
  width = 9,
  height = 5,
  dpi = 600
)


# Viral load ----
# Data

viral_load_data <- read_excel(
  file.path(data_dir, "viral_load_heart.xlsx"),
  sheet = "Sample sheet_shipped"
)

viral_load_df <- viral_load_data %>%
  transmute(
    treatment = factor(
      Treatment,
      levels = c("CONU", "p-TagRFP", "IV-LD", "Eomes", "GATA3")
    ),
    sampling = factor(
      Sampling,
      levels = c("10WPI", "4WPC")
    ),
    aqua_id = `Aqua-ID`,
    fam_avg = as.numeric(
      str_replace(as.character(`FAM-avg.`), ",", ".")
    ),
    rox_avg = as.numeric(
      str_replace(as.character(`ROX-avg.`), ",", ".")
    )
  ) %>%
  filter(sampling %in% c("10WPI", "4WPC"))


## Checking baseline

baseline_summary <- viral_load_df %>%
  filter(sampling == "10WPI") %>%
  group_by(treatment) %>%
  summarise(
    n = n(),
    mean_ct = mean(fam_avg, na.rm = TRUE),
    sd_ct = sd(fam_avg, na.rm = TRUE),
    n_at_cutoff = sum(fam_avg >= 40, na.rm = TRUE),
    n_below_cutoff = sum(fam_avg < 40, na.rm = TRUE),
    .groups = "drop"
  )

baseline_summary


## Treatment-specific baseline

baseline_ct <- viral_load_df %>%
  filter(sampling == "10WPI") %>%
  group_by(treatment) %>%
  summarise(
    baseline_ct = mean(fam_avg, na.rm = TRUE),
    .groups = "drop"
  )

viral_load_relative <- viral_load_df %>%
  left_join(
    baseline_ct,
    by = "treatment"
  ) %>%
  mutate(
    delta_delta_ct = fam_avg - baseline_ct,
    relative_load = 2^(-delta_delta_ct)
  )

viral_load_4wpc <- viral_load_relative %>%
  filter(sampling == "4WPC")


## Statistics
### Normality

shapiro_viral_load <- viral_load_4wpc %>%
  group_by(treatment) %>%
  group_modify(~ {
    test <- shapiro.test(.x$relative_load)
    
    tibble(
      W = unname(test$statistic),
      p_value = test$p.value
    )
  }) %>%
  ungroup()

shapiro_viral_load


### Kruskal-Wallis and Dunn

kruskal_viral_load <- viral_load_4wpc %>%
  kruskal_test(relative_load ~ treatment)

kruskal_viral_load

dunn_viral_load_all <- viral_load_4wpc %>%
  dunn_test(
    relative_load ~ treatment,
    p.adjust.method = "BH",
    detailed = TRUE
  )

dunn_viral_load_all


### Plot annotations

stat.test_viral_load <- dunn_viral_load_all %>%
  filter(p.adj < 0.05) %>%
  mutate(
    p.label = if_else(
      p.adj < 0.001,
      "p < 0.001",
      sprintf("p = %.3f", p.adj)
    ),
    y.position = c(7.1, 7.5, 6.4)
  ) %>%
  select(
    group1,
    group2,
    p.label,
    y.position
  )


### Plot

viral_load_plot <- ggplot(
  viral_load_4wpc,
  aes(
    x = treatment,
    y = relative_load,
    fill = treatment
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey"
  ) +
  geom_boxplot(
    width = 0.7,
    color = "black",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 1,
    alpha = 0.7
  ) +
  scale_y_log10(
    labels = trans_format(
      "log10",
      math_format(10^.x)
    ),
    breaks = log_breaks(n = 6)
  ) +
  scale_x_discrete(
    labels = c(
      "CONU" = "CONU",
      "p-TagRFP" = "Plasmid",
      "IV-LD" = "IV-LD",
      "Eomes" = "EOMES",
      "GATA3" = "GATA3"
    )
  ) +
  scale_fill_manual(
    values = treatment_colors
  ) +
  labs(
    x = NULL,
    y = expression(
      "Relative viral load (" *
        2^{-Delta*Delta*Ct} *
        ", log"[10] * " scale)"
    )
  ) +
  annotate(
    "text",
    x = -Inf,
    y = 1.12,
    label = "Treatment-specific 10WPI baseline",
    hjust = -0.02,
    vjust = 0,
    family = "Times New Roman",
    size = 4
  ) +
  stat_pvalue_manual(
    stat.test_viral_load,
    label = "p.label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    size = 3,
    tip.length = 0.005,
    vjust = -0.2,
    step.increase = 0,
    bracket.nudge.y = 0,
    hide.ns = TRUE,
    inherit.aes = FALSE,
    family = "Times New Roman"
  ) +
  theme_classic(
    base_family = "Times New Roman"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12)
  )

print(viral_load_plot)

ggsave(
  file.path(figure_dir, "viral_load_plot.jpg"),
  plot = viral_load_plot,
  width = 9,
  height = 5,
  dpi = 600
)


# Histoscore ----
# Data

histoscore_data <- read_xlsx(
  file.path(data_dir, "histoscore_data_updated.xlsx")
) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("CONU", "Plasmid", "IV-LD", "EOMES", "GATA3")
    ),
    Sampling = factor(
      Sampling,
      levels = c("4WPC", "6WPC", "10WPC")
    ),
    Tissue = factor(
      Tissue,
      levels = c("Pancreas", "Heart")
    )
  )

## Statistics
### Normality

shapiro_histoscore <- histoscore_data %>%
  group_by(Tissue, Sampling, Treatment) %>%
  group_modify(~ {
    test <- shapiro.test(.x$Histoscore)
    
    tibble(
      W = unname(test$statistic),
      p_value = test$p.value
    )
  }) %>%
  ungroup()

shapiro_histoscore


### 4WPC treatment comparisons

kruskal_histoscore <- histoscore_data %>%
  filter(Sampling == "4WPC") %>%
  group_by(Tissue) %>%
  kruskal_test(Histoscore ~ Treatment)

kruskal_histoscore

dunn_histoscore_all <- histoscore_data %>%
  filter(Sampling == "4WPC") %>%
  group_by(Tissue) %>%
  dunn_test(
    Histoscore ~ Treatment,
    p.adjust.method = "BH",
    detailed = TRUE
  )

dunn_histoscore_all


### Plot annotations

stat.test_histoscore4WPC <- dunn_histoscore_all %>%
  filter(p.adj < 0.05) %>%
  mutate(
    p.label = if_else(
      p.adj < 0.001,
      "p < 0.001",
      sprintf("p = %.3f", p.adj)
    )
  ) %>%
  add_xy_position(
    x = "Tissue",
    dodge = 0.8,
    scales = "fixed"
  )

histoscore_bracket_base <- histoscore_data %>%
  filter(Sampling == "4WPC") %>%
  group_by(Tissue, Treatment) %>%
  summarise(
    upper_error = mean(Histoscore) +
      sd(Histoscore) / sqrt(n()),
    .groups = "drop"
  ) %>%
  group_by(Tissue) %>%
  summarise(
    base_y = max(upper_error) + 0.10,
    .groups = "drop"
  )

stat.test_histoscore4WPC <- stat.test_histoscore4WPC %>%
  left_join(
    histoscore_bracket_base,
    by = "Tissue"
  ) %>%
  group_by(Tissue) %>%
  arrange(xmax - xmin, .by_group = TRUE) %>%
  mutate(
    y.position = base_y +
      (row_number() - 1) * 0.20
  ) %>%
  ungroup()

histoscore_ymax <- max(
  stat.test_histoscore4WPC$y.position
) + 0.15


### Plot

histoscore_plot <- ggbarplot(
  data = histoscore_data %>%
    filter(Sampling == "4WPC"),
  x = "Tissue",
  y = "Histoscore",
  add = "mean_se",
  fill = "Treatment",
  position = position_dodge(0.8),
  error.plot = "upper_pointrange",
  ylab = "Histoscore"
) +
  scale_fill_manual(
    values = treatment_colors,
    name = "Treatment",
    labels = treatment_labels
  ) +
  geom_vline(
    xintercept = 1.5,
    color = "grey",
    linetype = "dashed",
    linewidth = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey10",
    linewidth = 0.1
  ) +
  stat_pvalue_manual(
    stat.test_histoscore4WPC,
    label = "p.label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    size = 3,
    tip.length = 0.005,
    step.increase = 0,
    bracket.nudge.y = 0,
    family = "Times New Roman",
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      ceiling(histoscore_ymax * 2) / 2,
      by = 0.5
    ),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(
    ylim = c(0, histoscore_ymax)
  ) +
  labs(x = NULL) +
  theme_manuscript

print(histoscore_plot)

ggsave(
  file.path(figure_dir, "histoscore_plot.jpg"),
  plot = histoscore_plot,
  width = 9,
  height = 5,
  dpi = 600
)


# Histoscore stacked plots
### Classify histoscores

histoscore_stacked <- histoscore_data %>%
  mutate(
    Histoscore_Classification = cut(
      Histoscore,
      breaks = c(-Inf, 0.5, 1.5, 2.5, 3),
      labels = c(
        "Low",
        "mediumLow",
        "mediumHigh",
        "High"
      ),
      right = TRUE
    )
  ) %>%
  count(
    Tissue,
    Treatment,
    Sampling,
    Histoscore_Classification,
    name = "Count"
  ) %>%
  group_by(
    Tissue,
    Treatment,
    Sampling
  ) %>%
  mutate(
    Total = sum(Count),
    Percentage = Count / Total
  ) %>%
  ungroup()


### Color palette

traffic_palette <- c(
  "#009E60",
  "#F0E442",
  "#E69F00",
  "#D32F2F"
)


### Plot helper (GPT generated)

make_histoscore_stack <- function(data, tissue, sampling, title, tag) {
  
  plot_data <- data %>%
    filter(
      Tissue == tissue,
      Sampling == sampling
    )
  
  totals <- plot_data %>%
    distinct(
      Treatment,
      Total
    )
  
  ggplot(
    plot_data,
    aes(
      x = Treatment,
      y = Percentage,
      fill = Histoscore_Classification
    )
  ) +
    geom_col() +
    geom_text(
      data = totals,
      aes(
        x = Treatment,
        y = 1,
        label = paste0("n = ", Total)
      ),
      inherit.aes = FALSE,
      vjust = -0.5,
      size = 3,
      family = "Times New Roman"
    ) +
    coord_cartesian(
      ylim = c(0, 1),
      clip = "off"
    ) +
    scale_fill_manual(
      name = "Histoscore",
      labels = c(
        "≤0.5",
        "0.5 – 1.5",
        "1.5 – 2.5",
        "2.5 – 3"
      ),
      values = traffic_palette
    ) +
    scale_y_continuous(
      labels = percent
    ) +
    labs(
      x = NULL,
      y = "Percentage of individuals",
      title = title,
      tag = tag
    ) +
    theme_pubr() +
    theme(
      panel.grid.major = element_blank(),
      plot.margin = unit(
        c(0.8, 0.5, 0.1, 0.5),
        "cm"
      ),
      plot.title = element_text(
        hjust = 0.5,
        size = 12
      ),
      text = element_text(
        size = 12,
        family = "Times New Roman"
      ),
      plot.tag.position = c(0.07, 0.99)
    )
}


### Individual panels

pancreas_4wpc <- make_histoscore_stack(
  histoscore_stacked,
  "Pancreas",
  "4WPC",
  "Pancreas\n4 weeks post-challenge",
  "A"
)

heart_4wpc <- make_histoscore_stack(
  histoscore_stacked,
  "Heart",
  "4WPC",
  "Heart\n4 weeks post-challenge",
  "B"
)

pancreas_6wpc <- make_histoscore_stack(
  histoscore_stacked,
  "Pancreas",
  "6WPC",
  "Pancreas\n6 weeks post-challenge",
  "C"
)

heart_6wpc <- make_histoscore_stack(
  histoscore_stacked,
  "Heart",
  "6WPC",
  "Heart\n6 weeks post-challenge",
  "D"
)

pancreas_10wpc <- make_histoscore_stack(
  histoscore_stacked,
  "Pancreas",
  "10WPC",
  "Pancreas\n10 weeks post-challenge",
  "E"
)

heart_10wpc <- make_histoscore_stack(
  histoscore_stacked,
  "Heart",
  "10WPC",
  "Heart\n10 weeks post-challenge",
  "F"
)


### Combine panels

histoscore_stacked_plots <-
  (
    pancreas_4wpc /
      pancreas_6wpc /
      pancreas_10wpc |
      heart_4wpc /
      heart_6wpc /
      heart_10wpc
  ) +
  plot_layout(
    guides = "collect",
    axes = "collect"
  ) &
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.box.just = "center",
    text = element_text(
      family = "Times New Roman"
    )
  )

print(histoscore_stacked_plots)

ggsave(
  file.path(
    figure_dir,
    "histoscore_stacked_plots.jpg"
  ),
  plot = histoscore_stacked_plots,
  width = 3700,
  height = 3700,
  dpi = 300,
  units = "px"
)



# Package citations ----

citation("R")
citation("readxl")
citation("tidyverse")
citation("rstatix")
citation("effectsize")
citation("ggpubr")
citation("patchwork")
citation("scales")
citation("MoMAColors")