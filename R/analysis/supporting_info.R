# Supporting information figures
# by Camille Sicangco
# 19 June 2025

# Fig S1: Jmax-T response w/ and w/o TC ########################################

Tair_vec = seq(20,60)

Jold = TJmax(Tair_vec,EaJ=33115,EdVJ=2e5,delsJ=635)
Jnew = TJmax_updated(Tair_vec,EaJ=33115,EdVJ=2e5,delsJ=635, Tcrit = 46.5, T50 = 50.4)

FigS1 = data.frame(Tair = rep(Tair_vec, 2), 
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
ggsave("figs/FigS1_JvsT.tiff", FigS1, height = 6, width = 10, bg = "white")

# Fig S2: Predawn time series ##################################################

# See R/processing/WTC4_data_processing.R

# Fig S3: Theoretical simulations with constant VPD ############################

# See R/analysis/T_range_testing.R

# Fig S4: Cost, gain, and profit functions at Tair = 30 deg C ##################

# See R/analysis/inst_sims.R

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
ggsave("figs/FigS5_Ci_vs_Pleaf.tiff", FigS5, width = 9, height = 6, bg = "white")

# Fig S6: T50 sensitivity analysis #############################################

# See R/analysis/Tthreshold_sensitivity.R
# Environment
Tair_vec = seq(20,60, by = 5)
PPFD = 1500
VPD = RHtoVPD(RH = 60, TdegC = Tair_vec)

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 0.5)
#out_Ps0.5_constRH = make_pred_compEB(Tair_sim.df, g1 = g1_alt)
out = bind_rows(lapply(1:nrow(Tair_sim.df), 
                       function(i) make_pred_compEB(
                         Tair = Tair_sim.df$Tair[i], Ps = Tair_sim.df$Ps[i], 
                         VPD = Tair_sim.df$VPD[i], PPFD = Tair_sim.df$PPFD[i],
                         )))
out$gs_feedback = rep(c("TRUE", "FALSE"), times = length(Tair_vec))
out %>% 
  #filter(gs_feedback == FALSE) %>% 
  ggplot(aes(x = Tair, y = E, color = gs_feedback)) +
  geom_point() +
  theme_classic()
