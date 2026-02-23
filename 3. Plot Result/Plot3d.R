library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(plot3D)

path <- "table_varyAsym_D2_coverage_from2025_vs_baseD1_33_allCoverage_levels_tanzania.xlsx"
sheet_names <- excel_sheets(path)

to_num <- function(x) as.numeric(str_remove(as.character(x), "%"))

y_vars <- c("red_inc_r_1y","red_inc_r_3y","red_GR_inc_1y","red_GR_inc_3y")

# Read ALL sheets and bind rows
df_all <- map_dfr(sheet_names, function(sh){
  
  dat <- read_excel(path, sheet = sh)
  
  dat %>%
    mutate(
      cov_D2 = to_num(`D2`)*100,
      across(any_of(y_vars), to_num)
    )
})

df_all <- df_all %>% arrange(Asym, cov_D2)

# 3D scatter for each Y
for (y in y_vars) {
  scatter3D(
    x = df_all$cov_D2,
    y = df_all[[y]],
    z = df_all$Asym,
    xlab = "cov_D2",
    ylab = y,
    zlab = "Asym_XXpct",
    main = paste("3D:", y),
    pch = 19
  )
}



# ---- Line plot ----
library(plotly)

# Make sure points are ordered along the line
p <- plot_ly()

for (a in sort(unique(df_all$Asym), decreasing = TRUE)) {
  d <- df_all %>% filter(Asym == a) %>% arrange(cov_D2)
  
  p <- p %>% add_trace(
    data = d,
    x = ~cov_D2,
    y = ~.data[[metric_to_plot]],
    z = ~Asym,
    type = "scatter3d",
    mode = "lines+markers",
    name = paste0("Asym=", a)
  )
}

p %>% layout(
  scene = list(
    xaxis = list(title = "cov_D2"),
    yaxis = list(title = metric_to_plot),
    zaxis = list(title = "Asym")
  ),
  title = paste("3D line plot:", metric_to_plot)
)

df_long <- df_plot %>%
  pivot_longer(cols = all_of(y_vars), names_to = "metric", values_to = "value")

ggplot(df_long, aes(x = cov_f, y = Asym_f, fill = value)) +
  geom_tile() +
  facet_wrap(~metric, ncol = 2, scales = "free") +
  labs(x = "cov_D2 (%)", y = "Asym (%)", fill = "value") +
  theme_minimal()

# ---- 3D line plot: one line per Asym ----
p <- plot_ly()

for (a in sort(unique(df_all$Asym), decreasing = TRUE)) {
  d <- df_all %>% filter(Asym == a) %>% arrange(cov_D2)
  
  p <- p %>%
    add_trace(
      data = d,
      x = ~cov_D2,
      y = ~.data[[metric_to_plot]],
      z = ~Asym,
      type = "scatter3d",
      mode = "lines+markers",
      name = paste0("Asym=", a)
    )
}

p %>%
  layout(
    title = paste("3D line plot:", metric_to_plot),
    scene = list(
      xaxis = list(title = "cov_D2 (%)"),
      yaxis = list(title = metric_to_plot),
      zaxis = list(title = "Asym")
    )
  )
