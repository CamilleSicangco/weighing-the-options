# Miscellaneous debugging

# Point test ###################################################################
P50 = 4.31
P88 = 6.64
Weibull = fit_Weibull(P50, P88)
b = Weibull[1,1]
c = Weibull[1,2]
Pcrit = calc_Pcrit(b, c)

Ps = 0.25
kmax_25 = 0.7
Tair = 36.5
VPD = 3
PPFD = 1750
P = Ps_to_Pcrit(Ps, Pcrit)
E_vec = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)

E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, 
                   E = E, Wind = 8, Wleaf = 0.01, 
                   LeafAbs = 0.5)
D_leaf = plantecophys::VPDairToLeaf(VPD = VPD, Tair = Tair, 
                                    Tleaf = Tleaf)
g_w = calc_gw(E, D_leaf)

Rd0 = 1.115
TrefR = 25
# Tleaf = 38.04572
# Rd = 2.154639
Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2))
A_net = calc_A(Tair, VPD, PPFD, E = E, Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
           Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
           Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
           net = TRUE, netOrig = TRUE)

Photosyn(Tleaf = Tleaf[2], GS = g_w[2], VPD = VPD, PPFD = PPFD, 
         Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
         Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)
plot(A_net)

calc_A(Tair, VPD, PPFD, Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
       Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
       Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
       net = TRUE, netOrig = TRUE,
       Tleaf = Tleaf[2], g_w = g_w[2])

Photosyn_custom(VPD = VPD, PPFD = PPFD, 
         Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
         Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
         g1 = 2.9, g0 = 0.003,
         GS = g_w[2], Tleaf = Tleaf[2])


# Test predictions #############################################################
A_gross = calc_A(Tair, VPD, PPFD, E = E, Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
               Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
               Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
               net = FALSE)

CG_net = C_gain(P, b, c, Amax = NULL, kmax_25, Tair, VPD, PPFD, 
                Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
                Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
                Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
                constant_kmax = TRUE, 
                net = TRUE, netOrig = TRUE)

CG_gross = C_gain(P, b, c, Amax = NULL, kmax_25, Tair, VPD, PPFD, 
                  Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
                  Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
                  Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
                  constant_kmax = TRUE, net = FALSE)
HC = hydraulic_cost(P, b, c, kmax_25, Tair, constant_kmax = TRUE)
TC = thermal_cost(P, b, c, kmax_25, Tair, VPD, PPFD, 
                  Wind = 8, Wleaf = 0.01, LeafAbs = 0.5, 
                  Tcrit = 45.5, T50 = 47.5, constant_kmax = TRUE)

HC.m = marginal_GainCost(P, HC)
CC.m = marginal_GainCost(P, HC + TC)
CG_gross.m = marginal_GainCost(P, CG_gross)
CG_net.m = marginal_GainCost(P, CG_net)

which.min(abs(CG_gross.m-HC.m))
which.min(abs(CG_net.m-CC.m))

plot(abs(CG_gross.m-HC.m))
plot(abs(CG_net.m-CC.m), ylim = c(-.5, 1))
plot(CG_net, ylim = c(0,1))
plot(CG_gross)
E[144]
E[185]
A[144]
A[145]

plot(pred2.hw$Tleaf, pred2.hw.test$Tleaf)
abline(a = 0, b = 1)

# Ac predictions ###############################################################
# Ac function
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
Vcmax = 34
Tleaf = 35
GammaStar = TGammaStar(Tleaf, Patm = 100)
Km = TKm(Tleaf, Patm = 100)
Cis = seq(0,800, by = 5)
Ac <- Vcmax*(Cis - GammaStar)/(Cis + Km)

Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2))

plot(Cis,Ac)
abline(h = Rd)

# Predictions over a range, including Ci #######################################
Photosyn_v = Vectorize(Photosyn_custom2)
Photo_out = Photosyn_v(Tleaf = Tleaf, GS = g_w, VPD = VPD, PPFD = PPFD, 
           Vcmax=34,EaV=51780,EdVC=2e5,delsC=640, # add other VJ params
           Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)
Cis_pred = unlist(Photo_out[1,])
As_pred = unlist(Photo_out[2,])
plot(g_w, Cis_pred)
plot(g_w, As_pred) # Same thing for both Ac and Aj


# Model comparison #############################################################
# Medlyn (-0.319)
Photosyn(Tleaf=out.hw$Tleaf[1451],VPD=out.hw$VPD[1451],
           PPFD=out.hw$PPFD[1451],g1=g1,g0=g0,
           Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
           Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
         GS = out.hw$gs[1451], Ca = 420)
out.hw$A[1451]

# Sicangco
Photosyn(Tleaf=out.hw$Tleaf[2723],VPD=out.hw$VPD[2723],
         PPFD=out.hw$PPFD[2723],g1=g1,g0=g0,
         Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
         Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
         GS = out.hw$gs[2723], Ca = 420)
out.hw$A[2723]

# Sperry
Photosyn(Tleaf=out.hw$Tleaf[3995],VPD=out.hw$VPD[3995],
         PPFD=out.hw$PPFD[3995],g1=g1,g0=g0,
         Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
         Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
         GS = out.hw$gs[3995])
out.hw$A[3995]
