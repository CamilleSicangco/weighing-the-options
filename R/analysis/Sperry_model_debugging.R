# Single T run ##################
Tair = 52
PPFD = 1500
VPD = 1.5
test_25deg = make_pred(Tair = Tair, PPFD = PPFD, VPD = VPD, Ps = 0.5,
                       Vcmax = 92, Jmax = 165.53,
                       kmax_25 = 1.5)
# Hydraulics
kmax_25 = 3
Weibull = fit_Weibull(P50 = 4.07, P88 = 5.50)
b = Weibull[1,1]
c = Weibull[1,2]
Pcrit = calc_Pcrit(b, c)
P = Ps_to_Pcrit(Ps = 0.5, Pcrit, pts = 300)

# Photosynthesis
Vcmax = 92
Jmax = 165.53

# Calculate costs and gains
HC = hydraulic_cost(P, b, c, kmax_25, Tair, constant_kmax = TRUE)
CG_gross_52 = Cgain_alt_lowres(P, b, c, kmax_25 = kmax_25, 
                  Tair = Tair, PPFD = PPFD, VPD = VPD, 
                  Vcmax = Vcmax, Jmax = Jmax, 
                  constant_kmax = TRUE, net = FALSE, new_JT = FALSE)
CG_net = C_gain(P, b, c, kmax_25 = kmax_25, 
                Tair = Tair, PPFD = PPFD, VPD = VPD, 
                Vcmax = Vcmax, Jmax = Jmax,  
                constant_kmax = TRUE, net = TRUE, netOrig = TRUE, new_JT = FALSE)
E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, E = E)
g_w = calc_gw(E, Tleaf, Tair = Tair, VPD = VPD, PPFD = PPFD)
Dleaf = plantecophys::VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, 
                                   VPD = VPD, Pa = 101.325)
Photosyn_out = mapply(Photosyn_custom, 
                      VPD = Dleaf, PPFD = PPFD, Tleaf = Tleaf, 
                      GS = g_w, Jmax = Jmax, Vcmax = Vcmax)
Anet = as.numeric(Photosyn_out[2, ])
Rd = as.numeric(Photosyn_out[8, ])
#Agross = Anet + Rd
Ci = as.numeric(Photosyn_out[1, ])
out = data.frame(index = seq(0,1, length.out = 300),
                      P, E, Ci, Anet, Rd, CG_gross, CG_net, HC, 
                      profit_gross = CG_gross - HC, profit_net = CG_net - HC)
out_long =  out %>% 
  pivot_longer(P:profit_net)
out_long$name = factor(out_long$name, levels = c("P", "E", "Ci", "Agross", "Rd",
                                                 "CG_gross", "CG_net", "HC", 
                                                 "profit_gross", "profit_net"))
ggplot(out_long, aes(x = index, y = value)) +
  geom_point() +
  #theme_classic() +
  facet_wrap(vars(name), scales = "free")
ggsave("debugging_out/fig_Tair60_ProfitMax_constrCi.jpg", height = 8.5, width = 11)
write.csv(out,"debugging_out/outputs_Tair60_ProfitMax_constrCi.csv", row.names = FALSE)
# Multiple #################
# Environment
Tair_vec = seq(20,60, by = 8)
#VPD = RHtoVPD(RH = 60, TdegC = Tair_vec)


# Medlyn model: fit b
b_USO = log(((0.05 - 1e-5)/1.6*400/5.89 - 1) * sqrt(1.5) / 2.9)/(-0.5 + 0.2)
g1_alt = ((0.05 - 1e-5)/1.6*400/5.89 - 1)*sqrt(1.5) / exp(0.55*(-0.5 + 0.3))


# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 0.5)
test_0.5 = get_predictions(Tair_sim.df, g1 = g1_alt, 
                       Vcmax = 92, Jmax = 165.53,
                       kmax_25 = 3)
#test_0.5_constkmax$idx[test_0.5_constkmax$Model == "Sperry + CGnet"]
test_0.5 %>%
  #filter(Model %in% c("Sperry", "Sperry + CGnet")) %>% 
  ggplot(aes(x = Tleaf, y = E, color = Model)) +
  geom_point() +
  theme_classic()

test_0.5 %>% 
  filter(Model %in% c("Sperry", "Sperry + CGnet")) %>% 
  mutate(Pleaf = -P) %>% 
  select(Model, Tair, Pleaf, E, gs, A, Ci) %>% 
  pivot_longer(Pleaf:Ci) %>% 
  ggplot(aes(x = Tair, color = Model, y = value)) +
  geom_point() +
  facet_wrap(vars(name), scales = "free")
ggsave("debugging_out/multiple_Ts.jpg", height = 6, width = 9)

test_constrCi = get_predictions(Tair_sim.df, g1 = g1_alt, 
                       Vcmax = 92, Jmax = 165.53,
                       kmax_25 = 3, constr_Ci = TRUE)

# Ps = -2 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 2)
test_2 = get_predictions(Tair_sim.df, g1 = g1_alt, 
                         Vcmax = 92, Jmax = 165.53,
                         kmax_25 = 3)

# Ps = -4 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 4)
test_4 = get_predictions(Tair_sim.df, g1 = g1_alt, 
                         Vcmax = 92, Jmax = 165.53,
                         kmax_25 = 3)

#plot(test$Tleaf, test$A)

# Plot all
# Combine dataframes into one list
test_constVPD = list(test_0.5, test_2, test_4)

# Generate composite plots
plot_composite_Trange(test_constVPD, vars = "gs, A, Pleaf")
plot_composite_Trange(test_constVPD, vars = "E, Dleaf")


# Pleaf quite off compared to Manon's
test_constrCi %>%
  filter(Model != "Medlyn") %>% 
  #mutate(Ci = as.numeric(Ci)) %>% 
  #filter(Ci <1E4) %>% 
  #filter(Model == "Sperry + CGnet") %>% 
  ggplot(aes(x = Tleaf, y = Ci, color = Model)) +
  geom_point() +
  theme_classic() #+
  #geom_hline(yintercept = 400, color = "grey", linetype = "dashed")

test_USO = out.l$heatwave %>% filter(Model == "Sperry")
Photosyn_USO = mapply(Photosyn_custom, VPD = test_USO$Dleaf, 
                      Ca = 400, PPFD = heatwave$PPFD, Tleaf = test_USO$Tleaf, Patm = 101.325, 
                      GS = test_USO$gs, Vcmax = 100.52, Jmax = 165.53#, 
                      #g1 = g1, g0 = g0
                      )
plot(test_USO$A, as.numeric(Photosyn_USO[2, ]))

out_Ps0.5_constVPD_Dleaf_PM_simpleE %>% 
  mutate(partial_deriv = c(NaN, diff(A)/diff(E))) %>% 
  filter(!is.infinite(partial_deriv)& !is.nan(partial_deriv)) %>% 
  filter(Model != "Medlyn") %>% 
  ggplot(aes(x = E, y = partial_deriv, color = Model)) +
  geom_point()
calc_A()

Tair = 25
Ps = 0.5
P50 =4.07
P88 = 5.5
# Hydraulics
Weibull = fit_Weibull(P50, P88)
b = Weibull[1,1]
c = Weibull[1,2]
Pcrit = calc_Pcrit(b, c)
P = Ps_to_Pcrit(Ps, Pcrit, pts = 500)

E = trans_from_vc(P, kmax_25 = 0.5, Tair, b, c, constant_kmax = TRUE)

Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD,
                   E = E)
gs_in = calc_gw(E = E, Tleaf = Tleaf, Tair = Tair, VPD = VPD, PPFD = PPFD)
Dleaf = plantecophys::VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)
Photosyn_out = mapply(plantecophys::Photosyn, VPD = Dleaf,
                      Ca = 400, PPFD = PPFD, Tleaf = Tleaf, GS = gs_in)
A = as.numeric(Photosyn_out[2, ])
Ci = as.numeric(Photosyn_out[1, ])
gs_out1 = as.numeric(Photosyn_out[3, ])

Ci <- Ca - Am/GC
gs_out <- A/(400 - Ci)*1.57
plot(gs_in, col = "red")
points(gs_out)
plot(gs_in - gs_out)
