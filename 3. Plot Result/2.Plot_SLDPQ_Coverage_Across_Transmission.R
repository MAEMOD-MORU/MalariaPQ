library(ggplot2)
library(dplyr)

# --- Filter to Asym Want only ----
asym_show <- c(20, 40, 60, 80)

df_plot <- df_results %>%
  filter(Asym %in% asym_show) %>%
  mutate(Asym_label = factor(paste0(Asym), levels = c("20","40","60","80")))

colors_manual <- c(
  "20" = "#4472C4",  # blue
  "40" = "#ED7D31",  # orange
  "60" = "#70AD47",  # green
  "80" = "#FF0000"   # red
)


df_results <- df_results %>%
  mutate(red_inc_r_3y = ifelse(D2_coverage == 0, 0, red_inc_r_3y))

# --- df_plot ---
df_plot <- df_results %>%
  filter(Asym %in% asym_show) %>%
  mutate(Asym_label = factor(paste0(Asym), levels = c("20","40","60","80"))) %>%
  arrange(Asym_label, D2_coverage)

# --- Plot ---
ggplot(df_plot, aes(x = D2_coverage, y = red_inc_r_3y,
                    color = Asym_label, group = Asym_label)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +  # เส้น y=0
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = colors_manual,
    name   = "Asymptomatic infections(%)"
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, 20),
    expand = c(0.01, 0),
    limits = c(0, 100)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 105),   # ขยาย limit ด้านล่างให้เห็นเส้นติดลบ
    expand = c(0, 0)
  ) +
  labs(
    title = "Effect of SLDPQ Coverage Across Transmission Settings",
    x     = "SLDPQ Coverage (%)",
    y     = "% change in resistant incidence\n(vs no-SLDPQ, same setting)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title        = element_text(face = "bold", size = 12),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 10),
    legend.key.width  = unit(1.5, "cm"),
    axis.title.y      = element_text(angle = 90, vjust = 0.5, size = 10),
    axis.text         = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5)))
