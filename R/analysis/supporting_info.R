# Supporting information figures
# by Camille Sicangco
# 19 June 2025

# TODO Fig S1: ProfitMax with CG via gross or net photosynthesis w/ constrained Ci ####

# Fig S2: Jmax-T response w/ and w/o TC ########################################
Jold = TJmax(Tair_vec,EaJ=33115,EdVJ=2e5,delsJ=635)
Jnew = TJmax_updated(Tair_vec,EaJ=33115,EdVJ=2e5,delsJ=635, Tcrit = 43.4, T50 = 49.6)
FigS2 = data.frame(Tair = rep(Tair_vec, 2), 
                   Jmax = c(Jold, Jnew),
                   Method = rep(c("old", "new"), each = length(Tair_vec))) %>%
  ggplot(aes(x = Tair, y = Jmax, linetype = Method)) +
  geom_line() +
  theme_classic() +
  ylab(expression("J"[max]*" / J"[max25])) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  theme(text = element_text(size = 14)) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c(
                          expression("(1 - "*italic("TC")*") "*italic("J")["peakedArr"]),
                          expression(italic("J")["peakedArr"])
                        )) +
  guides(linetype = guide_legend(title = expression("J"[max])))
ggsave("figs/FigS2_JvsT.tiff", FigS2, height = 6, width = 10)

# Fig S3: Predawn time series ##################################################
predawn_LWP =  
  predawn_df %>% 
  filter(tissue == "leaf") %>% 
  rename(Ps = LWP) %>% 
  mutate(date = as.POSIXct(Date)) %>% 
  as.data.frame()

FigS3 = ggplot(NULL, aes(x = date, y = Ps, color = chamber)) +
  geom_point(data = predawn_LWP) +
  geom_line(data = predawn_est) +
  theme_classic() +
  guides(color = guide_legend(title = "Chamber")) +
  xlab("Date") +
  ylab(expression("Predawn  "*psi[leaf]*" (MPa)")) +
  theme(text = element_text(size = 14))
ggsave("figs/FigS3_Pleaf_timeseries.tiff", FigS3, height = 6, width = 10)

# Fig S4: Cost, gain, and profit functions at Tair = 30 deg C ##################

cost_gain30 = calc_costgain_netorig(P, b, c, kmax_25 = kmax_25, 
                                    Tair = 30, PPFD = PPFD, VPD = VPD,
                                    Tcrit = Tcrit, T50 = T50, #constant_kmax = TRUE,
                                    Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
                                    Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                                    Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)
FigS4 = composite_plot(cost_gain30)
ggsave(filename = "figs/FigS4_30deg_inst_sim.tiff", FigS4, width = 12.25, height = 6.5)

# Fig S5: Predicted Ci vs Pleaf for Tair = 48 deg C, ProfitMax_net #############
# Set hydraulic parameters
Weibull = fit_Weibull(P50 = 4.07, P88 = 5.50)
b = Weibull[1,1]
c = Weibull[1,2]
Pcrit = calc_Pcrit(b, c)

# Environmental variables
Ps = 0.5
Tair = 48
VPD = 1.5
PPFD = 1500

# Calculate physiological variables
P = Ps_to_Pcrit(Ps, Pcrit)
E = trans_from_vc(P, kmax_25 = 0.5, Tair, b, c, constant_kmax = FALSE)
Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, 
                   E = E, Wind = 8, Wleaf = 0.02, 
                   LeafAbs = 0.5)
g_w = calc_gw(E, Tleaf = Tleaf, Tair = Tair, VPD = VPD, PPFD = PPFD,
              Wind = 8, Wleaf = 0.01)

# Calculate photosynthesis and associated variables
Photosyn_v = Vectorize(Photosyn_custom)
Photo_out = Photosyn_v(Ca = 420, GS = g_w, Tleaf = Tleaf, VPD = 1.5, PPFD = 1500,
                       Vcmax=36.4,EaV=62307,EdVC=2e5,delsC=639,
                       Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635,
                       Rd0 = 0.92
)
Ci_pred = unlist(Photo_out[1,])

# Plot Ci vs Pleaf
FigS5 = data.frame(P, Ci = Ci_pred) %>% 
  ggplot(aes(x = P, y = Ci)) +
  geom_point(shape = 1) +
  theme_classic() +
  xlab(expression(Psi[leaf]*" (-MPa)")) +
  ylab(expression("C"[i]*" (ppm)")) +
  theme(text = element_text(size = 14)) +
  geom_hline(yintercept = 420, linetype = "dashed")
ggsave("figs/FigS5_Ci_vs_Pleaf.tiff", FigS5, width = 9, height = 6)

# Fig S6: T50 sensitivity analysis #############################################

palette = scales::seq_gradient_pal("slategray1","darkslateblue", "Lab")(seq(0,1,length.out=4))
preds_varT50 = make_pred_Tthresholds(Tair_sim.df, 
                                     Wind = 8, Wleaf = 0.02, LeafAbs = 0.5, 
                                     kmax_25 = 0.5, constant_kmax = FALSE,
                                     hold_Tcrit = TRUE,
                                     Thold_val = 43.4,
                                     Tvar_vals = c(44.5, 45.5, 47.5, 49.6))

T50.plt = preds_varT50 %>% 
  ggplot(aes(x = Tleaf, y = gs, linetype = T50, color = T50)) + 
  geom_line(linewidth = 1) + 
  theme_classic() +
  ylab(expression("g"[s]*" (mol m"^-2*"s"^-1*")")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  guides(linetype = guide_legend(title = expression("T"[50]*" (\u00B0C)")),
         color = guide_legend(title = expression("T"[50]*" (\u00B0C)"))) +
  scale_colour_manual(values = palette) +
  theme(axis.title = element_text(size = 14))
ggsave("figs/FigS6_SA_T50.tiff", T50.plt, height = 7, width = 11)

# TODO Fig S7: kmax sensitivity analysis ############################################