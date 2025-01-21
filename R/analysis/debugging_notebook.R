# Miscellaneous debugging

# Issues with Anet and CG calculation ##########################################
## Point test ##################################################################
Photosyn()
# Set parameters
P50 = 4.07
P88 = 5.50
Weibull = fit_Weibull(P50, P88)
b = Weibull[1,1]
c = Weibull[1,2]
Pcrit = calc_Pcrit(b, c)
Ps = 0.25
kmax_25 = 0.7
Tair = 36.5
VPD = 3
PPFD = 1750
Rd0 = 1.115
TrefR = 25

# Calculate physiological variables
P = Ps_to_Pcrit(Ps, Pcrit)
E_vec = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, 
                   E = E, Wind = 8, Wleaf = 0.01, 
                   LeafAbs = 0.5)
D_leaf = plantecophys::VPDairToLeaf(VPD = VPD, Tair = Tair, 
                                    Tleaf = Tleaf)
g_w = calc_gw(E, Tleaf = Tleaf, Tair = Tair, VPD = VPD, PPFD = PPFD,
              Wind = 8, Wleaf = 0.01)
Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2))

# Calculate photosynthesis

# Net photosynthesis
## "Jumps" when g_w becomes nonzero
## Notably at high temps (>40) the jump is still evident, but Anet < 0 for all gs
A_net = calc_A(Tair, VPD, PPFD, E = E, Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
           Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, 
           Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
           net = TRUE, netOrig = TRUE)
plot(g_w, A_net)

# Agross - Rd (with Agross calculated assuming Rd = 0)
## Gives continuous curve but considerably different from the "true" Anet
A_net2 = calc_A(Tair, VPD, PPFD, E = E, Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
               Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, 
               Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
               net = TRUE, netOrig = FALSE)
points(g_w, A_net2, col = "blue")

# Comparison of calc_A vs Photosyn for second g_w value (i.e. min g_w > 0)
calc_A(VPD, PPFD, 
       Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
       Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
       net = TRUE, netOrig = TRUE,
       Tleaf = Tleaf[2], g_w = g_w[2])

Photosyn_custom(VPD = VPD, PPFD = PPFD, 
         Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
         Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
         g1 = 2.9, g0 = 0.003,
         GS = g_w[2], Tleaf = Tleaf[2], Ca = 420, Patm = 101.325)

# Calculate costs and gains

## CGnet calculated based on Anet or Agross - Rd
CG_net = C_gain(P, b, c, Amax = NULL, kmax_25, Tair, VPD, PPFD, 
                Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
                Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, 
                Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
                constant_kmax = TRUE, 
                net = TRUE, netOrig = TRUE)
CG_net2 = C_gain(P, b, c, Amax = NULL, kmax_25, Tair, VPD, PPFD, 
                Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
                Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, 
                Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
                constant_kmax = TRUE, 
                net = TRUE, netOrig = FALSE)
CG_gross = C_gain(P, b, c, Amax = NULL, kmax_25, Tair, VPD, PPFD, 
                  Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
                  Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, 
                  Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
                  constant_kmax = TRUE, net = FALSE)
HC = hydraulic_cost(P, b, c, kmax_25, Tair, constant_kmax = TRUE)
TC = thermal_cost(P, b, c, kmax_25, Tair, VPD, PPFD, 
                  Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
                  Tcrit = 45.5, T50 = 47.5, constant_kmax = TRUE)

# Solve optimization
HC.m = marginal_GainCost(P, HC)
CC.m = marginal_GainCost(P, HC + TC)
CG_gross.m = marginal_GainCost(P, CG_gross)
CG_net.m = marginal_GainCost(P, CG_net)
CG_net2.m = marginal_GainCost(P, CG_net2)

i_Sperry = which.min(abs(CG_gross.m-HC.m)) # 144
i_Sicangco = which.min(abs(CG_net.m-CC.m)) # 146, near-equal to Sperry
i_Sicangco2 = which.min(abs(CG_net2.m-CC.m)) # 185, more open stomata than Sperry

# Plot gains minus costs
plot(abs(CG_gross.m-HC.m), ylim = c(0, 1.5)) # Sperry
points(abs(CG_net.m-CC.m), col = "orange") # Sicangco Anet
points(abs(CG_net2.m-CC.m), col = "blue") # Sicangco Agross - Rd

plot(CG_net)
points(CG_net2, col = "blue")

# Predicted E and A for Sperry and Sicangco2
# Sicangco2 predicts higher E, lower A
E[144] 
E[185]
A_net[144]
A_net2[185]

plot(pred2.hw$Tleaf, pred2.hw.test$Tleaf)
abline(a = 0, b = 1)

## Ac predictions ##############################################################
# Ac functions
.Rgas <- function()8.314
Tk <- function(x)x+273.15

# Arrhenius
arrh <- function(Tleaf, Ea){
  exp((Ea * (Tk(Tleaf) - 298.15)) / (298.15 * .Rgas() * Tk(Tleaf))) 
}

TGammaStar <- function(Tleaf, Patm,
                       Egamma=37830.0, 
                       value25=42.75){  
  
  value25*arrh(Tleaf,Egamma)*Patm/100
}

TKm <- function(Tleaf, Patm,
                Oi = 210,      # O2 concentration (mmol mol-1)
                Ec = 79430.0,  # activation energy for Kc 
                Eo = 36380.0,  # activation energy for Ko
                Kc25 = 404.9,  # Kc at 25C
                Ko25 = 278.4  # Ko at 25C
){
  
  Oi <- Oi * Patm / 100
  
  Ko <- Ko25*arrh(Tleaf, Eo)
  Kc <- Kc25*arrh(Tleaf, Ec)
  Km <- Kc * (1.0 + Oi / Ko)
  
  return(Km)
}

TVcmax <- function(Tleaf, EaV, delsC, EdVC){
  
  if(EdVC > 0){
    V1 <- 1+exp((delsC*(25 + 273.15)-EdVC)/(.Rgas()*(25 + 273.15)))
    V2 <- 1+exp((delsC*(Tleaf+273.15)-EdVC)/(.Rgas()*(Tk(Tleaf))))
    f <- V1/V2
  } else f <- 1
  
  exp((Tleaf-25)*EaV/(.Rgas()*Tk(Tleaf)*Tk(25))) * f
}

# Calculate Ac and Rd for a set temperature over a range of Cis
Vcmax = 34
Tleaf = 35
Vcmax = Vcmax * TVcmax(Tleaf,EaV, delsC, EdVC)
GammaStar = TGammaStar(Tleaf, Patm = 100)
Km = TKm(Tleaf, Patm = 100)
Cis = seq(0,800, by = 5)
Ac <- Vcmax*(Cis - GammaStar)/(Cis + Km)

Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2))

# Plot Ac vs Ci, plus Rd line
plot(Cis,Ac)
abline(h = Rd)
## Only for Ci < ~180 is Ac < Rd, i.e. Anet < 0 when Tleaf = 35
## If can get model to behave s.t. Ci stays low when gs is near 0, then would 
## give desired behaviour(?). But notably equations are defined s.t. the limit
## of Ci as gs approaches 0 is Ca!

## Photosynthetic predictions over a range, including Ci #######################
Photosyn_v = Vectorize(Photosyn_custom)
Photo_out = Photosyn_v(Tleaf = Tleaf, GS = g_w, VPD = VPD, PPFD = PPFD, 
           Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
           Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)
Cis_pred = unlist(Photo_out[1,]) # Calculated as Ci = Ca - Am/GC
CICs_pred = unlist(Photo_out[15,]) # CICs calculated with Leuning 1990 quadratic
As_pred = unlist(Photo_out[2,])

plot(g_w, Cis_pred)
plot(g_w, CICs_pred) # Near-constant, obviously not true! Needs gs dependency
plot(g_w, As_pred) # Same thing for both Ac and Aj

Photosyn_custom2(Tleaf = Tleaf[1], GS = g_w[1], VPD = VPD, PPFD = PPFD, 
           Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
           Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)

# Photosynthesis predictions at high temperatures ##############################

gs_vec = c(1, 1e-1, 1e-2, 1e-3, 1e-4)
A_v_gs_preds = Photosyn_v(GS = gs_vec, Tleaf=44.5,VPD=4.4,
                          PPFD=1800,g1=g1,g0=g0,
                          Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
                          Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633, Ca = 420)

# Lower gs gives higher Aleaf?!
A_preds = unlist(A_v_gs_preds[2,])
plot(log(gs_vec), A_preds, ylim = c(-2, 5))

# Note that this is only an issue with Ac, not Aj (which decreases as expected)
Ac_preds = unlist(A_v_gs_preds[5,])
points(log(gs_vec), Ac_preds, col = "blue")
Aj_preds = unlist(A_v_gs_preds[6,])
points(log(gs_vec), Aj_preds, col = "orange")
