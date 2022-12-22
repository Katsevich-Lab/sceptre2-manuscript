z_grid <- seq(-4, 4, length.out = 1000)
density_df <- data.frame(density = dnorm(z_grid),
                         z_grid = z_grid)
histogram_df <- data.frame(z_null = correlated_res$resamp_dist$camp_null)
p3 <- ggplot() + geom_histogram(aes(x = z_null, y = after_stat(density)),
                                data = histogram_df,
                                boundary = 0,
                                fill = "grey85",
                                color = "black",
                                bins = 25) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  geom_line(aes(x = z_grid, y = density, col = "N(0,1) density"),
            data = density_df, linewidth = 0.7) +
  scale_color_manual(values = c("N(0,1) density" = "purple")) +
  xlab("z null") +
  ylab("") +
  theme(legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(),
        legend.margin=margin(t = -0.5, unit='cm')) +
  ggtitle("SCEPTRE null z-scores") +
  xlim(-4, 4)


histogram_df <- data.frame(z_null = uncorrelated_res$resamp_dist$camp_null)
p4 <- ggplot() + geom_histogram(aes(x = z_null, y = after_stat(density)),
                                data = histogram_df,
                                boundary = 0,
                                fill = "grey85",
                                color = "black",
                                bins = 25) +
  my_theme +
  scale_y_continuous(expand = expansion(mult = c(0.0, .01))) +
  geom_line(aes(x = z_grid, y = density),
            data = density_df, linewidth = 0.7, col = "purple") +
  xlab("z null") +
  ylab("") +
  ggtitle("SCEPTRE null z-scores") +
  xlim(-6, 6)