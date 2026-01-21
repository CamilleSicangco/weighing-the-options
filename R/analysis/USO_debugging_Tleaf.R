# Debugging of USO model Eleaf predictions

# Generate USO predictions with different LEB ---------------------------------
# Skip this section and proceed from line 34!

# Load data
WTC4_data = read.csv("data/in/WTC4_data.csv")
heatwave <- subset(WTC4_data,HWtrt=="HW" & PPFD > 500 & E >= 0)
df = heatwave

# Make Medlyn predictions
Medlyn_preds = plantecophys::PhotosynEB(Tair=df$Tair, VPD=df$VPD, PPFD = df$PPFD,
                                        g1 = 2.9, g0=1.e-5, b_USO = 0.55,
                                        Wind = 5, Wleaf = 0.025, LeafAbs = 0.5,
                                        Vcmax=100.52,EaV=62307,EdVC=2e5,delsC=639,
                                        Jmax = 165.53,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
                                        Ps = df$Ps
)
get_Pleaf_Medlyn_V = Vectorize(get_Pleaf_Medlyn, c("Ps", "Tair", "E"))
Medlyn_preds$P = get_Pleaf_Medlyn_V(
  Ps = df$Ps, Tair = df$Tair, E = Medlyn_preds$ELEAF, 
  P50 = 4.07, P88 = 5.50, kmax_25 = 1.5, constant_kmax = FALSE)

Medlyn_preds_sim_newLEB3 = Medlyn_preds %>%
  rename(E = ELEAF, Dleaf = VPDleaf, gs = GS, A = ALEAF) %>%
  mutate(Model = "Medlyn", .before = 1) %>%
  #mutate(P = NA) %>% 
  select(Model, Tair, E, Tleaf, Dleaf, gs, A, P, Ci
  )

save(Medlyn_preds_sim_oldLEB, Medlyn_preds_sim_newLEB, Medlyn_preds_sim_newLEB2,
     Medlyn_preds_sim_newLEB3, file = "data/out/USO_preds_LEB_test.Rdata")

# Compare outputs ----------------------------------------------
load("data/out/USO_preds_LEB_test.Rdata")

# VPD and PPFD responses are pretty much the same
plot(df$VPD, Medlyn_preds_sim_oldLEB$Tleaf)
points(df$VPD, Medlyn_preds_sim_newLEB$Tleaf, col = "red")

plot(df$PPFD, Medlyn_preds_sim_oldLEB$Tleaf)
points(df$PPFD, Medlyn_preds_sim_newLEB$Tleaf, col = "red")

# gs very similar - so why is E so much noisier????
## Answer: issue is in leaf energy balance calculation 
## - needed to remove old unit conversion for Cp in Penman-Monteith calc of E
plot(Medlyn_preds_sim_oldLEB$Tleaf, Medlyn_preds_sim_oldLEB$gs)
points(Medlyn_preds_sim_newLEB$Tleaf, Medlyn_preds_sim_newLEB$gs, col = "red")

# Compare E-Tleaf
plot(Medlyn_preds_sim_oldLEB$Tleaf, Medlyn_preds_sim_oldLEB$E)
points(Medlyn_preds_sim_newLEB$Tleaf, Medlyn_preds_sim_newLEB$E, col = "red")

# Add version with fixed units in P-M
## Fixed spread of E, but now magnitude is off?
## 3rd version is fixed! Issue was an extra AIRMA multiplier
points(Medlyn_preds_sim_newLEB$Tleaf, Medlyn_preds_sim_newLEB2$E, col = "blue")
points(Medlyn_preds_sim_newLEB$Tleaf, Medlyn_preds_sim_newLEB3$E, col = "green")

# Magnitude is off by a factor of around 3.75
hist(Medlyn_preds_sim_oldLEB$E/Medlyn_preds_sim_newLEB2$E)
summary(Medlyn_preds_sim_oldLEB$E/Medlyn_preds_sim_newLEB2$E)

# Not a problem for most other variables, e.g. A
plot(Medlyn_preds_sim_oldLEB$Tleaf, Medlyn_preds_sim_oldLEB$A)
points(Medlyn_preds_sim_newLEB$Tleaf, Medlyn_preds_sim_newLEB2$A, col = "blue")

# But is a problem for Pleaf, since calculated from E
plot(Medlyn_preds_sim_oldLEB$Tleaf, Medlyn_preds_sim_oldLEB$P)
points(Medlyn_preds_sim_newLEB$Tleaf, Medlyn_preds_sim_newLEB2$P, col = "blue")
