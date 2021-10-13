# Belowground evolution Figure 4

#   mutate(year = as.numeric(year) + 2020) %>% 
#   mutate(color_code = case_when(iteration == 101 ~ 1,
#                                 iteration == 102 ~ 2,
#                                 iteration == 103 ~ 3,
#                                 iteration == 104 ~ 4,
#                                 T ~ 0),
#          size_code = case_when(iteration > 100 ~ "big",
#                                T ~ "small")) %>% 
#   ggplot(aes(x = year, y = value, group = iteration, color = factor(color_code), size = factor(size_code))) + 
#   geom_line(aes(alpha = size_code)) +
#   ylab("marsh elevation (cm NAVD88)") +
#   scale_color_manual(values = c("gray11", colors[1], colors[4], colors[1], colors[4])) +
#   scale_size_manual(values = c(1.5,0.8)) +
#   scale_alpha_manual(values = c(0.8, 0.2)) +
#   geom_point(aes(x = 2100, y = run_store_cohort[1,80]), color = colors[1], size = 3) +
#   geom_point(aes(x = 2100, y = run_store_cohort[2,80]), color = colors[4], size = 3) +
#   geom_point(aes(x = 2100, y = run_store_cohort[3,80]), color = colors[1], size = 3, shape = 17) +
#   geom_point(aes(x = 2100, y = run_store_cohort[4,80]), color = colors[4], size = 3, shape = 17) +
#   theme(legend.position = "none") +
#   ylim(22,33.5)-> Fig4_panelA
# 
# tibble(acc_rate = avg_accretion_rates1*10) %>% 
#   ggplot(aes(x = acc_rate)) +
#   geom_histogram(bins=10, color = "black", fill = "white") +
#   xlab(expression(paste("vertical accretion rate (mm ",yr^-1,")"))) +
#   geom_point(aes(x = cohort_summary$acc_v[1], y = 0), color = colors[1], size = 3) +
#   geom_point(aes(x = cohort_summary$acc_v[2], y = 0), color = colors[4], size = 3) +
#   geom_point(aes(x = cohort_summary$acc_v[3], y = 0), color = colors[1], size = 3, shape = 17)+
#   geom_point(aes(x = cohort_summary$acc_v[4], y = 0), color = colors[4], size = 3, shape = 17) -> Fig4_panelB
# 
# tibble(acc_rate = avg_C_accum_rate * 1e-6 / 1e-8) %>% 
#   ggplot(aes(x = acc_rate)) +
#   geom_histogram(bins = 10, fill = "white", color = "black") +
#   xlab(expression(paste("carbon accumulation rate (t C ", ha^-1, yr^-1,")"))) +
#   geom_point(aes(x = cohort_summary$acc_C[1], y = 0), color = colors[1], size = 3) +
#   geom_point(aes(x = cohort_summary$acc_C[2], y = 0), color = colors[4], size = 3) +
#   geom_point(aes(x = cohort_summary$acc_C[3], y = 0), color = colors[1], size = 3, shape = 17) +
#   geom_point(aes(x = cohort_summary$acc_C[4], y = 0), color = colors[4], size = 3, shape = 17) -> Fig4_panelC 
# 
# Fig4_panelsBC <- cowplot::plot_grid(Fig4_panelB, Fig4_panelC, nrow = 2, labels = c("b", "c"))
# Fig4_panelA_label <- cowplot::plot_grid(Fig4_panelA, labels = "a")


tibble(elevation_store_withcohorts2) %>% 
  mutate(iteration = 1:104) %>% 
  gather(key = year, value = value, `1`:`80`) %>% 
  mutate(year = as.numeric(year) + 2020) %>% 
  mutate(color_code = case_when(iteration == 101 ~ 1,
                                iteration == 102 ~ 2,
                                iteration == 103 ~ 3,
                                iteration == 104 ~ 4,
                                T ~ 0),
         size_code = case_when(iteration > 100 ~ "big",
                               T ~ "small")) %>% 
  ggplot(aes(x = year, y = value, group = iteration, color = factor(color_code), size = factor(size_code))) + 
  geom_line(aes(alpha = size_code)) +
  ylab("marsh elevation (cm NAVD88)") +
  scale_color_manual(values = c("gray67", colors[1], colors[4], colors[1], colors[4])) +
  scale_size_manual(values = c(1.5,0.8)) +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  geom_point(aes(x = 2100, y = run_store_cohort2[1,80]), color = colors[1], size = 3) +
  geom_point(aes(x = 2100, y = run_store_cohort2[2,80]), color = colors[4], size = 3,) +
  geom_point(aes(x = 2100, y = run_store_cohort2[3,80]), color = colors[1], size = 3, shape = 17) +
  geom_point(aes(x = 2100, y = run_store_cohort2[4,80]), color = colors[4], size = 3, shape = 17) +
  theme(legend.position = "none") +
  ylim(22,33.5)-> Fig4_panelD

## Last panel -- compare variation for agb + bgb and abg only models ####

quantile_agb_only <- quantile(run_store2[,80], c(0.025, 0.975))
quantile_agb_bgb <- quantile(run_store1[,80], c(0.025, 0.975))

tibble(y = c(run_store2[,80], run_store1[,80]),
       x = rep(c("agb only", "agb + bgb"), each = length(run_store1[,80]))) %>% 
  ggplot(aes(x = x, y = y, fill = x)) +
  geom_violin(draw_quantiles = c(0.025,0.5, 0.975), size = 0.5, alpha = 0.4) +
  scale_fill_manual(values = c("gray11", "gray67")) +
  theme(legend.position = "none") +
  ylab("elevation at t = 80 (cm NAVD88)") +
  xlab("scenario") -> Fig4_panelE

## Bring all plots together ####
Fig4_panelsDE <- cowplot::plot_grid(Fig4_panelD, Fig4_panelE, nrow = 1,
                                    labels = c("d", "e"), rel_widths = c(3,2))
Fig4_panelsABC <- cowplot::plot_grid(Fig4_panelA_label, Fig4_panelsBC, rel_widths = c(3,2))

Fig4_panelsABCDE <- cowplot::plot_grid(Fig4_panelsABC, Fig4_panelsDE, nrow = 2)

png("chp1/plots/MS_plots/Vahsen_Fig4.png", height = 6.8, width = 8, units = "in", res = 300)
Fig4_panelsABCDE
dev.off()
