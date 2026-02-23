library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)

path <- "table_varyAsym_D2_coverage_from2025_vs_baseD1_33_allCoverage_levels_23_jan.xlsx"
sheet_names <- excel_sheets(path)

# helper: robust numeric conversion (handles "10%", "0.1", 10, etc.)
to_num <- function(x){
  x <- as.character(x)
  x <- str_replace_all(x, "%", "")
  as.numeric(x)
}

# Read + combine all sheets
df_all <- map_dfr(sheet_names, function(sh){
  
  dat <- read_excel(path, sheet = sh)
  
  # Asym from sheet name like "Asym_90pct" -> 0.9
  asym_val <- str_match(sh, "Asym_(\\d+)pct")[,2] %>% as.numeric() / 100
  
  dat %>%
    mutate(
      Asym = asym_val,
      cov_D2 = to_num(`durg.coverage.D2`),
      red_inc_r_1y  = to_num(red_inc_r_1y),
      red_inc_r_3y  = to_num(red_inc_r_3y),
      red_GR_inc_1y = to_num(red_GR_inc_1y),
      red_GR_inc_3y = to_num(red_GR_inc_3y)
    )
})

# Long format for plotting
df_long <- df_all %>%
  pivot_longer(
    cols = c(red_inc_r_1y, red_inc_r_3y, red_GR_inc_1y, red_GR_inc_3y),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    # order Asym: 0.9 -> 0.1 (top to bottom)
    Asym_f = factor(Asym, levels = sort(unique(Asym), decreasing = TRUE)),
    cov_f  = factor(cov_D2, levels = sort(unique(cov_D2)))
  )

# Heatmaps (facet by metric)
ggplot(df_long, aes(x = cov_f, y = Asym_f, fill = value)) +
  geom_tile() +
  facet_wrap(~ metric, ncol = 2, scales = "free") +
  labs(x = "cov_D2", y = "Asym (fraction)", fill = "value") +
  theme_minimal()

# single plot
metric_to_plot <- "red_inc_r_1y"

df_plot <- df_all %>%
  mutate(
    Asym_f = factor(Asym, levels = sort(unique(Asym), decreasing = TRUE)),
    cov_f  = factor(cov_D2, levels = sort(unique(cov_D2)))
  )

ggplot(df_plot,
       aes(x = cov_f,
           y = Asym_f,
           fill = .data[[metric_to_plot]])) +
  geom_tile() +
  labs(
    x = "cov_D2",
    y = "Asym (fraction)",
    fill = metric_to_plot,
    title = paste("Heatmap:", metric_to_plot)
  ) +
  theme_minimal()
