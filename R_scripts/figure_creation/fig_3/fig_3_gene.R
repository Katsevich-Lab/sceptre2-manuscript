shared_fig_script <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                            "sceptre2-manuscript/R_scripts/figure_creation/shared_figure_script.R")
source(shared_fig_script)
library(sn)
library(tidyverse)
conflict_prefer("filter", "dplyr")

xi <- 0
omega <- 1
alpha <- 3
n <- 10000
z_star <- 1.8

##########################
# panel c 
##########################

set.seed(1)
null_z <- rsn(n = n, xi = xi, omega = omega, alpha = alpha)

panel_c <- tibble(null_z) |>
  ggplot(aes(x = null_z)) +
  geom_histogram(aes(y = after_stat(density)),
    color = "black",
    fill = "grey85",
    bins = 20
  ) +
  stat_function(
    fun = function(x) (dsn(x, xi, omega, alpha)),
    color = "darkblue",
    linewidth = 0.75,
    xlim = c(-1, 4.2)
  ) +
  labs(x = expression(paste("Null z-scores ", tilde(z)[i]))) +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  my_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(filename = paste0(.get_config_path("LOCAL_CODE_DIR"),
                         "sceptre2-manuscript/R_scripts/figure_creation/fig_3/figure_3c.png"),
       plot = panel_c, device = "png", width = 3, height = 2.25)

##########################
# panel d
##########################
panel_d <- tibble(z = null_z) |>
  filter(z > z_star) |>
  ggplot(aes(x = z)) +
  stat_function(
    fun = function(z) (dsn(z, xi, omega, alpha)),
    color = "darkblue",
    linewidth = 0.75,
    xlim = c(-1, 4.2)
  ) +
  geom_ribbon(aes(ymax = dsn(z, xi, omega, alpha)), 
              ymin = 0, alpha = 0.5, fill = "darkblue") +
  geom_vline(xintercept = z_star, col = "purple", linewidth = 1) +
  labs(x = "Null z-score distribution") +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  my_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(filename = paste0(.get_config_path("LOCAL_CODE_DIR"),
                         "sceptre2-manuscript/R_scripts/figure_creation/fig_3/figure_3d.png"),
       plot = panel_d, device = "png", width = 3, height = 2.25)

