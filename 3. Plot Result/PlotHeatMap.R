library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

path <- "table_varyAsym_D2_coverage_from2025_vs_baseD1_33_allCoverage_levels_tanzania.xlsx"
sheet_names <- excel_sheets(path)

to_num <- function(x) as.numeric(str_remove(as.character(x), "%"))

y_vars <- c("red_inc_r_1y","red_inc_r_3y","red_GR_inc_1y","red_GR_inc_3y")

# ---- Read all sheets and bind ----
df_all <- map_dfr(sheet_names, function(sh){
  
  dat <- read_excel(path, sheet = sh)
  
  # from sheet name like "D2_80pct" -> 80
  covD2 <- as.numeric(str_match(sh, "D2_(\\d+)pct")[,2])
  
  dat %>%
    mutate(
      cov_D2 = covD2,
      Asym   = to_num(Asym),              # should be 90,80,... or 0.9,... depending on file
      across(any_of(y_vars), to_num)
    )
})

# ---- Choose ONE metric for the heatmap ----
metric_to_plot <- "red_GR_inc_3y"   # change to any of y_vars

# ---- Order axes ----
df_plot <- df_all %>%
  mutate(
    cov_f  = factor(cov_D2, levels = sort(unique(cov_D2))),                 # 0 -> 100
    Asym_f = factor(Asym,   levels = sort(unique(Asym), decreasing = TRUE)) # 90 -> 10 (top)
  )

# ---- Heatmap ----
ggplot(df_plot, aes(x = cov_f, y = Asym_f, fill = .data[[metric_to_plot]])) +
  geom_tile() +
  labs(
    x = "cov_D2 (%)",
    y = "Asym (%)",
    fill = metric_to_plot,
    title = paste("Heatmap:", metric_to_plot)
  ) +
  theme_minimal()
