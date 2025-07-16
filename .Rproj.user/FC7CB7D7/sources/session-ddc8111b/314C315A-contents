### Visualize in forest plot


all_logp <- c(-log10(age_df$p_value), -log10(sex_df$p_value), -log10(both_df$p_value))
color_limits <- range(all_logp, na.rm = TRUE)

p_age <- age_df[1:12,] %>% 
  ggplot(aes(x=estimate, y=trait)) + 
  geom_effect(ggplot2::aes(xmin = CI_lower, xmax = CI_upper, colour = -log10(p_value), shape = 'fixed'), fatten = 6, position = ggstance::position_dodgev(height = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, colour = -log10(p_value)), height = 0.2) +
  theme_forest() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) + 
  scale_color_gradientn(colors = c("#33A02C", "#FB9A99", "#1F78B4", "#CAB2D6", "#FFF123") , limits = color_limits) +
  scale_shape_manual(values = c("fixed" = 18)) +
  xlim(-8,8) +
  labs(title = expression("age" %~% "glycosylation trait"), x = "std.estimates (95% CI)") +
  guides(shape = "none") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")

p_sex <- sex_df[c(1:6,15:20),] %>% 
  ggplot(aes(x=estimate, y=trait)) + 
  geom_effect(ggplot2::aes(xmin = CI_lower, xmax = CI_upper, colour = -log10(p_value), shape = 'fixed'), fatten = 6, position = ggstance::position_dodgev(height = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, colour = -log10(p_value)), height = 0.2) +
  theme_forest() + 
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +  
  scale_color_gradientn(colors = c("#33A02C", "#FB9A99", "#1F78B4", "#CAB2D6", "#FFF123") , limits = color_limits) +
  scale_shape_manual(values = c("fixed" = 18)) +
  xlim(0, 2) + labs(title = expression("sex" %~% "glycosylation trait"), x = "odds ratios (95% CI)") +
  guides(shape = "none") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")

p_both <- both_df[1:12,] %>% 
  ggplot(aes(x=estimate, y=trait)) + 
  geom_effect(ggplot2::aes(xmin = CI_lower, xmax = CI_upper, colour = -log10(p_value), shape = 'fixed'), fatten = 6, position = ggstance::position_dodgev(height = 0.5)) + 
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, colour = -log10(p_value)), height = 0.2) +
  theme_forest() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) + 
  scale_color_gradientn(colors = c("#33A02C", "#FB9A99", "#1F78B4", "#CAB2D6", "#FFF123") , limits = color_limits) +
  scale_shape_manual(values = c("fixed" = 18)) +
  xlim(-8, 8) + labs(title = expression("age" %~% "glycosylation trait * sex"), x = "std.estimates (95% CI)") +
  guides(shape = "none") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")

#scale_color_gradient2(low = "#33A02C", mid = "#FB9A99", high = "#1F78B4", midpoint = 40, limits = color_limits) 

p <- p_age + p_sex + p_both



## Merge plots

# Merge 4 forest plots stored in the list with patchwork
p <- (wrap_plots(c(list(p_age),list(p_sex),list(p_both)), ncol=3) + plot_annotation(tag_levels = 'A')) + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") 
#& guides(color = guide_legend(nrow = 1, override.aes = list(size = 3)), shape = guide_legend(nrow = 1, override.aes = list(size = 3)))

# Save the combined plot
ggsave("Results/Stat.png", plot = p, width = 12, height = 8 , units = "in", dpi = 500)
