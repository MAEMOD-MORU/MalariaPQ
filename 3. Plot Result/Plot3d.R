library(readxl)
library(dplyr)
library(stringr)

path <- "table_varyAsym_D2_coverage_from2025_vs_baseD1_33_allCoverage_levels_23_jan.xlsx"
sheet_names <- excel_sheets("table_varyAsym_D2_coverage_from2025_vs_baseD1_33_allCoverage_levels_23_jan.xlsx")
df <- read_excel(path) %>%
  mutate(
    cov_D2 = as.numeric(str_remove(as.character(durg.coverage.D2), "%")),
    Asym_XXpct = as.numeric(str_remove(as.character(Asym), "%"))
  )

y_vars <- c("red_inc_r_1y","red_inc_r_3y","red_GR_inc_1y","red_GR_inc_3y")

# y_vars <- data.frame(
#   red_inc_r_1y = as.numeric(str_remove(as.character(df$red_inc_r_1y), "%")),
#   red_inc_r_3y = as.numeric(str_remove(as.character(df$red_inc_r_3y), "%")),
#   red_GR_inc_1y = as.numeric(str_remove(as.character(df$red_GR_inc_1y), "%")),
#   red_GR_inc_3y = as.numeric(str_remove(as.character(df$red_GR_inc_3y), "%"))
# )

library(plot3D)

for (y in y_vars) {
  scatter3D(
    x = df$cov_D2,
    y = as.numeric(str_remove(as.character(df[[y]]), "%")),
    z = df$Asym_XXpct,
    xlab = "cov_D2",
    ylab = y,
    zlab = "Asym_XXpct",
    main = paste("3D:", y),
    pch = 19
  )
}
