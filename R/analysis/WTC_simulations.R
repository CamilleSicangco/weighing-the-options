# Model testing with the WTC4 dataset
# by Camille Sicangco
# Created 30 Sept 2024

# Load data frame
WTC4_data = read.csv("data/in/WTC4_data.csv")
WTC4_data$DateTime_hr <- as.POSIXct(WTC4_data$DateTime_hr,format="%Y-%m-%d %T",tz="GMT")

# Determine max VPDs during the heatwave
VPDs = WTC4_data %>% 
  dplyr::filter(DateTime_hr >= as.POSIXct("2016-11-01") & DateTime_hr <= as.POSIXct("2016-11-04")) %>% 
  group_by(chamber, HWtrt) %>% 
  summarise(maxVPD = max(VPD))
range(VPDs$maxVPD[VPDs$HWtrt == "C"])
range(VPDs$maxVPD[VPDs$HWtrt == "HW"])

# Fit models ###################################################################

## Control ---------------------------------------------------------------------

# Subset the model predictions for the ambient treatment  
control <- subset(WTC4_data,HWtrt=="C" & PPFD > 500)

## Prescribed E ----------------------------------------------------------------

# Test energy balance by prescribing E
pred_prescE.c = prescribedE_pred(df = filter(control, !is.na(gs)))

diff_long.c = pred_prescE.c %>% 
  select(ends_with(".diff")) %>%
  pivot_longer(cols = ends_with(".diff"), names_to = "var", values_to = "diff",
               names_pattern = "(.*).diff")
obs_long.c = pred_prescE.c %>% 
  select(ends_with(".obs")) %>%
  select(!E.obs) %>% 
  pivot_longer(cols = ends_with("obs"), names_to = "var", values_to = "obs",
               names_pattern = "(.*).obs")
prescE_long.c = cbind(diff_long.c, obs_long.c[2])

prescE_long.c %>% 
  ggplot(aes(x = obs, y = diff)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, color = "blue") +
  facet_wrap(vars(var), scales = "free") +
  theme_classic() +
  ggtitle("Control")

### Fitting --------------------------------------------------------------------

# GAMs
gam_E.c = gam(E ~ s(Tcan), data = control)
gam_A.c = gam(A ~ s(Tcan), data = control)

# Get gs model predictions
pred.c = get_predictions(df = control, b_USO = 0, Tcrit = 43.7, T50 = 48.6)

pred.c$Tdiff <- with(pred.c,Tleaf-Tair)

# Check if any model predicts negative gs values
sapply(unique(pred.c$Model), 
      function(model) table(pred.c$gs[pred.c$Model == model] >= 0)) # all predict positive gs

# Save control observations and predictions
save(control, pred.c, file = "data/out/control_runs.Rdata")
load("data/out/control_runs.Rdata")


# Heatwave  ####################################################################

# Subset the model predictions for the ambient treatment  
heatwave <- subset(WTC4_data,HWtrt=="HW" & PPFD > 500 & E >= 0)

## Prescribed E tests ---------------------------------------------------------- 
pred_prescE.hw = prescribedE_pred(df = filter(heatwave, !is.na(gs)))

diff_long.hw = pred_prescE.hw %>% 
  select(ends_with(".diff")) %>%
  pivot_longer(cols = ends_with(".diff"), names_to = "var", values_to = "diff",
               names_pattern = "(.*).diff")
obs_long.hw = pred_prescE.hw %>% 
  select(ends_with(".obs")) %>%
  select(!E.obs) %>% 
  pivot_longer(cols = ends_with("obs"), names_to = "var", values_to = "obs",
               names_pattern = "(.*).obs")
prescE_long.hw = cbind(diff_long.hw, obs_long.hw[2])

prescE_long.hw %>% 
  ggplot(aes(x = obs, y = diff)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, color = "blue") +
  facet_wrap(vars(var), scales = "free") +
  theme_classic() + ggtitle("Heatwave")

# Test leaf temperature predictions only
Tleaves.hw = try(calc_Tleaf(heatwave$Tair, VPD = heatwave$VPD, E = heatwave$E, 
                           PPFD = heatwave$PPFD, Wind = 4, Wleaf = 0.04, LeafAbs = 0.5))
Tleaves.hw[grep("Error", Tleaves.hw)] = NA
Tleaves.hw = as.numeric(Tleaves.hw)

# Tends to over-predict leaf temperature slightly
Tdiff.hw = Tleaves.hw - heatwave$Tcan
summary(Tdiff.hw)
hist(Tdiff.hw)
plot(x = heatwave$Tcan, y = Tdiff.hw, col = "deeppink")
abline(h = 0, col = "blue")

## Fitting ---------------------------------------------------------------------

# Fit GAMs
gam_E.hw = gam(E ~ s(Tcan), data = heatwave)
gam_A.hw = gam(A ~ s(Tcan), data = heatwave)

# Get gs model predictions
pred.hw = get_predictions(df = heatwave, b_USO = 0)

 
# Check if any model predicts negative gs values
sapply(unique(pred.hw$Model), 
       function(model) table(pred.hw$gs[pred.hw$Model == model] >= 0)) # all predict positive gs

# Save heatwave observations and predictions
save(heatwave, pred.hw, file = "data/out/heatwave_runs.Rdata")
load("data/out/heatwave_runs.Rdata")

# Plotting #####################################################################

## Prep outputs ---------------------

preds = list(list(control, pred.c), list(heatwave, pred.hw_DKparams))
out.l = lapply(preds, function(ls) {
  out = data.frame(datetime = c(ls[[1]]$DateTime_hr, rep(ls[[1]]$DateTime_hr, each = 4), 
                                ls[[1]]$DateTime_hr),
                   chamber = c(ls[[1]]$chamber, rep(ls[[1]]$chamber, each = 4), 
                                ls[[1]]$chamber),
                   Model = c(rep("observed", nrow(ls[[1]])), ls[[2]]$Model),
                   E = c(ls[[1]]$E, ls[[2]]$E),
                   A = c(ls[[1]]$A, ls[[2]]$A),
                   gs = c(ls[[1]]$gs, ls[[2]]$gs),
                   Dleaf = c(ls[[1]]$Dleaf, ls[[2]]$Dleaf),
                   Tleaf = c(ls[[1]]$Tcan, ls[[2]]$Tleaf),
                   Pleaf = c(rep(NA, nrow(ls[[1]])), ls[[2]]$P))
  
  out$Model = factor(out$Model, levels = c("observed", "Medlyn", "Sperry", 
                                               "Sperry + varkmax", "Sperry + CGnet", 
                                               "Sperry + CGnet + TC"))
  out = out %>% arrange(Model)
  return(out)
  
})

names(out.l) = c("control", "heatwave")

palette = c(Sperry = "#1E88E5", "Sperry + varkmax" = "#332288", "Sperry + CGnet" = "#e7d87d", "Sperry + CGnet + TC" = "orange", 
            Medlyn = "#D81B60", observed = "grey50")

## Figure 5: A, E vs Tleaf -------------------------

# Plot A and E versus canopy temperature for the control and heatwave treatments 
EvT.c = plot_AEvT_WTC(gam_E.c, out.l$control, "E")
AvT.c = plot_AEvT_WTC(gam_A.c, out.l$control, "A")
EvT.hw = plot_AEvT_WTC(gam_E.hw, out.l$heatwave, "E")
AvT.hw = plot_AEvT_WTC(gam_A.hw, out.l$heatwave, "A")

# Combine all plots (i.e. recreate Drake et al. Fig 5 with all models)
AEvT.plt = ggarrange(AvT.c + ylim(-2.5, 13) + 
                       labs(title = "Control",
                            subtitle = expression(bold("(a)"))) + 
                       theme(axis.title.x = element_blank(), 
                             plot.title = element_text(hjust = 0.5)), 
                     AvT.hw + ylim(-2.5, 13) +
                       labs(title = "Heatwave",
                            subtitle = expression(bold("(b)"))) + 
                       theme(axis.title.y = element_blank(), 
                             axis.title.x = element_blank(), 
                             plot.title = element_text(hjust = 0.5)), 
                     EvT.c + ylim(-1, 4) + ggtitle("(c)") + 
                       theme(axis.title.x = element_blank(), 
                             plot.title = element_text(face = "bold", hjust = 0, vjust = 5)), 
                     EvT.hw + ylim(-1,4) + ggtitle("(d)") + 
                       theme(axis.title.y = element_blank(), 
                             axis.title.x = element_blank(), 
                             plot.title = element_text(face = "bold", hjust = 0, vjust = 5)),
                     nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
AEvT.plt = 
  annotate_figure(AEvT.plt,
                  bottom = text_grob(expression("T"[leaf]*" (\u00B0C)"), hjust = 1))

ggsave(plot = AEvT.plt, filename = "figs/Fig5_AEvT_WTC.tiff", width = 8, height = 7, bg = "white")

## Figure 6: Tleaf predictions vs observations -------------------------

Tleaf_pred_obs.plt = 
  bind_rows(out.l, .id = "treatment") %>% 
  pivot_wider(names_from = Model, values_from = Tleaf, id_cols = c(chamber, datetime, treatment)) %>% 
  pivot_longer(cols = Medlyn:"Sperry + CGnet + TC", names_to = "Model", values_to = "Tleaf_pred") %>% 
  rename(Tleaf = observed) %>% 
  ggplot() +
  geom_point(aes(x = Tleaf, y = Tleaf_pred, color = Model), shape = 1) +
  theme_classic() +
  scale_color_manual(values = palette[-6],
                     labels = c("USO",
                                "ProfitMax",
                                expression("ProfitMax"[k[max](T)]),
                                expression("ProfitMax"[net]),
                                expression("ProfitMax"[TC])
                     )) +
  xlab(expression("observed T"[leaf]*" (\u00B0C)")) +
  ylab(expression("predicted T"[leaf]*" (\u00B0C)")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") + 
  geom_abline(slope = 1) +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2))) + 
  geom_vline(xintercept = 43.4, linetype = "dashed", colour = "darkorange") +
  geom_hline(yintercept = 43.4, linetype = "dashed", colour = "darkorange") +
  annotate("text", x=Tcrit, y = 10, label=expression("T"[crit]), hjust = -0.5, colour = "darkorange", size = 5)+ 
  geom_vline(xintercept = 49.6, linetype = "dashed", colour = "orangered3") +
  geom_hline(yintercept = 49.6, linetype = "dashed", colour = "orangered3") +
  annotate("text", x=T50, y = 10, label=expression("T"[50]), hjust = -0.5, colour = "orangered3", size = 5)+
  theme(plot.title = element_blank(),text = element_text(size = 14)) +  
  xlim(NA, 53)
ggsave("figs/Fig6_Tleaf_pred_vs_obs_WTC.tiff", Tleaf_pred_obs.plt, height = 7, width = 11, bg = "white")

## Figure 7: Pleaf vs Tleaf ----------------------------------------------------

Fig7_PleafvT = 
  bind_rows(out.l) %>% 
  filter(!is.na(Pleaf)) %>% 
  ggplot(aes(x = Tleaf, y = Pleaf, color = Model)) +
  geom_point(shape = 1) +
  scale_color_manual(values = palette[-c(6)],
                     labels = c("USO",
                                "ProfitMax",
                                expression("ProfitMax"[k[max](T)]),
                                expression("ProfitMax"[net]),
                                expression("ProfitMax"[TC])
                     )) +
  theme_classic() +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression(Psi[leaf]*" (-MPa)")) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2)),
         linetype = "none") +
  geom_hline(yintercept = 4.07, linetype = "dashed", colour = "darkorange") +
  annotate("text", x = 20, y = 4.07, label=expression("P"[50]), vjust = -0.5, colour = "darkorange", size = 5) +
  geom_hline(yintercept = 5.50, linetype = "dashed", colour = "orangered3") +
  annotate("text", x = 20, y = 5.50, label=expression("P"[88]), vjust = -0.5, colour = "orangered3", size = 5) +
  theme(text = element_text(size = 16)) +
  ylim(NA,5.7)
ggsave("figs/Fig7_Pleaf_vs_T_WTC.tiff", Fig7_PleafvT, height = 7, width = 10, bg = "white")

## E vs Tleaf -----------
out.l$heatwave %>% 
  filter(Model %in% c("observed", "Sperry", "Sperry + varkmax")) %>%
  ggplot(aes(x = Dleaf, y = E, color = Model)) +
geom_point(#size = .5, 
           alpha = 0.3) +
  theme_classic() +
  scale_color_manual(
    values = palette,
    labels = c("Observations",
               "ProfitMax",
               expression("ProfitMax"[k[max](T)])
    )) +
  xlab(expression("D"[leaf]*" (kPa)")) +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2, alpha = 1)),
         linetype = "none") +
  theme(plot.title = element_blank())

## gs vs Tleaf -----------

# GAMs
gam_gs.hw = gam(gs ~ s(Tcan), data = heatwave)
gam_gs.c = gam(gs ~ s(Tcan), data = control)

gs_vs_Tleaf.hw = plot_AEvT_WTC(gam_gs.hw, out.l$heatwave, "gs")
gs_vs_Tleaf.c = plot_AEvT_WTC(gam_gs.c, out.l$control, "gs")

gs_vs_Tleaf.plt = 
  ggarrange(gs_vs_Tleaf.c + #ylim(-2.5, 13) + 
              ggtitle(expression(bold("(a)")*" Control")) + 
              theme(axis.title.x = element_blank(), 
                    plot.title = element_text(hjust = 0)), 
            gs_vs_Tleaf.hw + #ylim(-2.5, 13) +
              ggtitle(expression(bold("(b)")*" Heatwave")) +
              theme(axis.title.y = element_blank(), 
                    axis.title.x = element_blank(), 
                    plot.title = element_text(hjust = 0)), 
            nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")

gs_vs_Tleaf.plt = annotate_figure(gs_vs_Tleaf.plt,
                bottom = text_grob(expression("T"[leaf]*" (\u00B0C)"), hjust = 1))

ggsave(plot = gs_vs_Tleaf.plt, filename = "figs/FigS7_gs_vs_Tleaf.tiff", width = 11, height = 6, bg = "white")

## Time series of A, E -----------

plt_timeseries = function(df,
                          heatwave = TRUE,
                          yvar = c("E", "A", "gs")) {
  palette = c(Sperry = "#1E88E5", "Sperry + varkmax" = "#332288", "Sperry + CGnet" = "#e7d87d", "Sperry + CGnet + TC" = "orange", 
              Medlyn = "#D81B60", observed = "grey50")
  
  ylabel = 
    if (yvar == "gs") {
      expression("g"[s]*" (mol m"^-2*"s"^-1*")")
    } else if (yvar == "P") {
      expression(psi[leaf]*" (-MPa)")
    } else if (yvar == "Tleaf") {
      expression("T"[leaf]*" (\u00B0C)")
    } else if (yvar == "A") {
      expression("A (" * mu * "mol m"^-2*"s"^-1*")")
    } else if (yvar == "Dleaf") {
      expression("VPD"[leaf]*" (kPa)")
    } else if (yvar == "E") {
      expression("E"*" (mmol m"^-2*"s"^-1*")")
    }
  
  plt = df %>% 
    ggplot() +
    geom_point(aes(x = datetime, y = !!sym(yvar), color = Model), alpha = 0.3, size = 0.5) +
    theme_classic() +
    scale_color_manual(
      values = palette,
      labels = c("Observations",
                 "USO",
                 "ProfitMax",
                 expression("ProfitMax"[k[max](T)]),
                 expression("ProfitMax"[net]),
                 expression("ProfitMax"[TC])
      )) +
    guides(color = guide_legend(override.aes = list(shape = 19, size = 2, alpha = 1)),
           linetype = "none") +
    ylab(ylabel) +
    xlab("Date") 
  
  if(isTRUE(heatwave)) {
    plt = plt +
      annotate("rect",
               xmin = as.POSIXct("2016-10-31", tz = "GMT"),
               xmax = as.POSIXct("2016-11-04", tz = "GMT"),
               ymin = -Inf,
               ymax = Inf,
               alpha = 0.2, fill = "red") 
  }
  
  return(plt)
}

yvars = c("A", "E", "gs")
timeseries.c = lapply(yvars, function(yvar) plt_timeseries(df =out.l$control, heatwave = FALSE, yvar))
timeseries.hw = lapply(yvars, function(yvar) plt_timeseries(df =out.l$heatwave, heatwave = TRUE, yvar))

# Create composite plot
legend = cowplot::get_legend(timeseries.c[[1]])
aligned_plots = cowplot::align_plots(timeseries.c[[1]] + ggpubr::rremove("xlab") +
                                       labs(title = "Control",
                                            subtitle = expression(bold("(a)"))) + 
                                       theme(axis.title.x = element_blank(), 
                                             plot.title = element_text(hjust = 0.5),
                                             legend.position = "none"), 
                                     timeseries.c[[2]] + ggpubr::rremove("xlab") + theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0, vjust = 5)) + ggtitle(bquote(bold("(b)"))),
                                     timeseries.c[[3]] + ggpubr::rremove("xlab") + theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0, vjust = 5)) + ggtitle(bquote(bold("(c)"))),
                                     timeseries.hw[[1]] + ggpubr::rremove("xlab") + ggpubr::rremove("ylab") +
                                       labs(title = "Heatwave",
                                            subtitle = expression(bold("(d)"))) + 
                                       theme(axis.title.x = element_blank(), 
                                             plot.title = element_text(hjust = 0.5),
                                             legend.position = "none"),
                                     timeseries.hw[[2]] + ggpubr::rremove("xlab") + ggpubr::rremove("ylab") + theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0, vjust = 5)) + ggtitle(bquote(bold("(e)"))),
                                     timeseries.hw[[3]] + ggpubr::rremove("xlab") + ggpubr::rremove("ylab") + theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0, vjust = 5)) + ggtitle(bquote(bold("(f)"))),
                                     align = "v", axis = "l")
plt_comp = cowplot::plot_grid(aligned_plots[[1]],
                              aligned_plots[[2]],
                              aligned_plots[[3]],
                              aligned_plots[[4]],
                              aligned_plots[[5]],
                              aligned_plots[[6]],
                              nrow = 3, ncol = 2, byrow = FALSE,
                              rel_heights = c(1, 0.75, 0.75, 1, 0.75, 0.75))
plt_comp = ggpubr::annotate_figure(plt_comp, 
                                   bottom = "Date")
aligned_plts2 = cowplot::align_plots(plt_comp, legend, align = "h", axis = "t")
timeseries.plt = cowplot::plot_grid(aligned_plts2[[1]], aligned_plts2[[2]],
                         rel_widths = c(1,0.2))
timeseries.plt

ggsave(plot = timeseries.plt, filename = "figs/FigS8_timeseries.tiff", width = 11, height = 7, bg = "white")


# Calculate TSM and HSM ########################################################

# TSM
TSM_summary = out.l$heatwave %>%
  filter(Model != "observed") %>% 
  group_by(Model) %>% 
  summarise(Tleaf_max = max(Tleaf)) %>% 
  mutate(TSM_Tcrit = 46.5 - Tleaf_max,
         TSM_T50 = 50.4 - Tleaf_max) %>% 
  arrange(factor(Model, levels = c("Medlyn", "Sperry", 
                                   "Sperry + varkmax", "Sperry + CGnet", 
                                   "Sperry + CGnet + TC"))) %>% 
  as.data.frame()
write.csv(TSM_summary, "data/out/TSM_summary.csv", row.names = FALSE)

# HSM
HSM_summary = out.l$heatwave %>%
  filter(Model != "observed") %>% 
  group_by(Model) %>% 
  summarise(Pleaf_min = max(Pleaf)) %>% 
  mutate(HSM_P50 = 4.07 - Pleaf_min,
         HSM_P88 = 5.50 - Pleaf_min) %>% 
  arrange(factor(Model, levels = c("Medlyn", "Sperry", 
                                   "Sperry + varkmax", "Sperry + CGnet", 
                                   "Sperry + CGnet + TC")))
write.csv(HSM_summary, "data/out/HSM_summary.csv", row.names = FALSE)


# Model evaluation #############################################################
eval_model = function(fit_group = c("low", "medium", "high", "none"),
                      treatment = c("control", "heatwave"),
                      model = c("Medlyn", "Sperry", "Sperry + varkmax", 
                                "Sperry + CGnet", "Sperry + CGnet + TC"),
                      variable = c("E", "A", "gs", "Dleaf", "Tleaf")
)
{
  df = if (treatment == "control") {out.l$control} else {out.l$heatwave}
  df = df %>% 
    select(!Pleaf) %>% 
    pivot_wider(names_from = Model, values_from = c(E, A, gs, Dleaf, Tleaf))
  
  
  df = if (fit_group == "low") {
    df %>% filter(Tleaf_observed < 25)
  } else if (fit_group == "medium" & treatment == "control") {
    df %>% filter(Tleaf_observed >= 25)
  } else if (fit_group == "medium" & treatment == "heatwave") {
    df %>% filter(Tleaf_observed >= 25 & Tleaf_observed < 35)
  } else if (fit_group == "high") {
    df %>% filter(Tleaf_observed >= 35)
  } else if (fit_group == "none") {
    df
  }
  
  # Create vector with observations
  obs = df[[paste0(variable, "_observed")]]
  
  # Create vector with predictions
  pred = df[[paste0(variable, "_", model)]]
  
  mae = mae_vec(obs, pred)
  r2 = rsq_vec(obs, pred)
  
  eval_metrics = c("MAE" = mae, "R2" = r2)
  
  return(eval_metrics)
}


eval_model(fit_group = "none",treatment = "heatwave", model = "Sperry", variable = "E")

## Data binned by temperature --------------------------------------------------

# Generate evaluation metrics
model_eval_df = data.frame(treatment = c(rep("control", 50), rep("heatwave", 75)),
                           fit_group = c(rep(c("low", "medium"), each = 25, times = 2), rep("high", 25)),
                           Model = rep(c("Medlyn", "Sperry", "Sperry + varkmax", 
                                         "Sperry + CGnet", "Sperry + CGnet + TC"), times = 25), 
                           variable = rep(c("E", "A", "gs", "Dleaf", "Tleaf"), each = 5, times = 5)
)
model_eval_df$Model = factor(model_eval_df$Model,
                             levels = c("observed", "Medlyn", "Sperry", 
                                        "Sperry + varkmax", "Sperry + CGnet", 
                                        "Sperry + CGnet + TC"))
model_evals = sapply(1:nrow(model_eval_df), 
                     function(i) eval_model(model_eval_df$fit_group[i],
                                            model_eval_df$treatment[i],
                                            model_eval_df$Model[i], 
                                            model_eval_df$variable[i]))
model_evals = as_data_frame(t(model_evals))
model_eval_df = cbind(model_eval_df, model_evals)

model_eval_df$fit_group = factor(model_eval_df$fit_group, 
                                 levels = c("low", "medium", "high"))

# Plot metrics
MAE_plt_binned.hw = plot_evals(model_eval_df, "heatwave", "MAE", binned = TRUE)
R2_plt_binned.hw = plot_evals(model_eval_df, "heatwave", "R2", binned = TRUE)
MAE_plt_binned.c = plot_evals(model_eval_df, "control", "MAE", binned = TRUE)
R2_plt_binned.c = plot_evals(model_eval_df, "control", "R2", binned = TRUE)

## Unbinned data----------------------------------------------------------------

# Generate evaluation metrics
model_eval_df_unbinned = data.frame(treatment = rep(c("control", "heatwave"), each = 25),
                                    fit_group = rep("none", 50),
                                    Model = rep(c("Medlyn", "Sperry", "Sperry + varkmax", 
                                                  "Sperry + CGnet", "Sperry + CGnet + TC"), 
                                                times = 10), 
                                    variable = rep(c("E", "A", "gs", "Dleaf", "Tleaf"), each = 5, times = 2)
)
model_evals_unbinned = sapply(1:nrow(model_eval_df_unbinned), 
                     function(i) eval_model(model_eval_df_unbinned$fit_group[i],
                                            model_eval_df_unbinned$treatment[i],
                                            model_eval_df_unbinned$Model[i], 
                                            model_eval_df_unbinned$variable[i]))
model_evals_unbinned = as_data_frame(t(model_evals_unbinned))
model_eval_df_unbinned = cbind(model_eval_df_unbinned, model_evals_unbinned)

# Plot metrics
MAE_plt.hw = plot_evals(model_eval_df_unbinned, "heatwave", "MAE", binned = FALSE)
R2_plt.hw = plot_evals(model_eval_df_unbinned, "heatwave", "R2", binned = FALSE)
MAE_plt.c = plot_evals(model_eval_df_unbinned, "control", "MAE", binned = FALSE)
R2_plt.c = plot_evals(model_eval_df_unbinned, "control", "R2", binned = FALSE)


ggarrange(MAE_plt.c + ggtitle("Control"), 
          MAE_plt.hw + ggtitle("Heatwave"), 
          R2_plt.c, 
          R2_plt.hw, 
          common.legend = TRUE)
