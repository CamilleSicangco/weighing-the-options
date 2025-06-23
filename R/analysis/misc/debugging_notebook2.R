# Calculate Photosyn for increasing Tleaf and constant gs ######################
Photosyn_v = Vectorize(Photosyn_custom)
Tleaves = seq(10,60, by = 1)
gs_decr = seq(0.1,0.001, length = length(Tleaves))

# If gs is unspecified, then A -> 0 as Tleaf increases
# If gs is set at a constant value, then has expected behaviour
# Likewise if gs decreases (with exception of increase at very end?)
Photo_out = Photosyn_v(Ca = 420, GS = 0.001, 
                       Tleaf = Tleaves, VPD = 3, PPFD = 1800,
                           Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                           Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635)
Vcmax <- 34 * TVcmax(Tleaf=Tleaves,EaV=62307, delsC=639, EdVC=2e5)
GammaStar = TGammaStar(Tleaf = Tleaves, Patm = 100)
Km = TKm(Tleaf = Tleaves, Patm = 100)
plot(Tleaves, Km)
points(Tleaves, GammaStar, col = "red")
#Cis = seq(0,800, by = 5)
Ac <- Vcmax*(Ci_pred - GammaStar)/(Ci_pred + Km)
As_pred = unlist(Photo_out[2,])
gs_pred = unlist(Photo_out[3,])
Ac_pred = unlist(Photo_out[5,])
Aj_pred = unlist(Photo_out[6,])
Rd_pred = unlist(Photo_out[8,])
Ci_pred = unlist(Photo_out[1,])

# If gs is unspecified, then A -> 0 as Tleaf increases
# If gs is set at a constant value, then has expected behaviour
# Likewise if gs decreases (with exception of increase at very end?)
plot(Tleaves, As_pred)
plot(Tleaves, Aj_pred, col = "blue")
points(Tleaves, Ac_pred, col = "orange")
points(Tleaves, Rd_pred, col = "grey")
plot(Tleaves,Ci_pred, ylim = c(300,500))
plot(Tleaves,Vcmax,ylim=c(0,55))
points(Tleaves,Ac_pred, col = "orange")
points(Tleaves,Ac, col = "darkgreen")
plot(Km)
# gs set to 0 whenever A <= 0 
plot(Tleaves, gs_pred)
plot(gs_pred, As_pred)
plot(Tleaves, gs_pred)

Photosyn(Ca = 420, Tleaf = 45, VPD = 3, PPFD = 1800,
         Vcmax=34,EaV=62307,EdVC=2e5,delsC=639, # add other VJ params
         Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635)

# Repeat with energy balance (increasing Tair but calculated Tleaf, gs) ########
# Note that can't provide gs directly to this function #########################
PhotosynEB_v = Vectorize(PhotosynEB)
Tair = Tleaves

PhotoEB_out = PhotosynEB_v(Ca = 420, Tair = Tair, VPD = 3, PPFD = 1800,
                           Wind=8,Wleaf=0.01,StomatalRatio=1,LeafAbs=0.5,
                           Vcmax=34,EaV=62307,EdVC=2e5,delsC=639, 
                           Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635,
                           g1 = 2.9, g0=0.003)
As_pred = unlist(PhotoEB_out[2,])
gs_pred = unlist(PhotoEB_out[3,])
Ac_pred = unlist(PhotoEB_out[5,])
Aj_pred = unlist(PhotoEB_out[6,])
Rd_pred = unlist(PhotoEB_out[8,])
Tleaf_pred = unlist(PhotoEB_out[10,])

plot(Tleaf_pred, As_pred)
plot(Tair, Tleaf_pred - Tair) # More decoupling of Tleaf-Tair as Tair increases

# Beyond intersection of Ac and Rd, gs = 0 and Ac ~ Rd s.t. A ~ 0
# Seems unrealistic, but matches observations around 40-45 deg?
plot(Tair, As_pred)
plot(Tair, Aj_pred, col = "blue")
points(Tair, Ac_pred, col = "orange")
points(Tair, Rd_pred, col = "grey")

plot(Tair, gs_pred)
plot(Tair, gs_pred)

Photosyn(Ca = 420, GS = 0.003, Tleaf = 55, VPD = 3, PPFD = 1800,
             #Wind=8,Wleaf=0.01,StomatalRatio=1,LeafAbs=0.5,
             Vcmax=34,EaV=62307,EdVC=2e5,delsC=639, # add other VJ params
             Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635)

PhotosynEB(Ca = 420, Tair = 60, VPD = 3, PPFD = 1800,
           Wind=8,Wleaf=0.01,StomatalRatio=1,LeafAbs=0.5,
           Vcmax=34,EaV=62307,EdVC=2e5,delsC=639, # add other VJ params
           Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635)


# Repeat for a range of values along transpiration supply stream ###############
# Set parameters
P50 = 4.07
P88 = 5.50
Weibull = fit_Weibull(P50, P88)
b = Weibull[1,1]
c = Weibull[1,2]
Pcrit = calc_Pcrit(b, c)
Ps = 0.25
kmax_25 = 0.7
Tair = 45
VPD = 3
PPFD = 1800
Rd0 = 1.115
TrefR = 25

# Calculate physiological variables
P = Ps_to_Pcrit(Ps, Pcrit)
E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, 
                   E = E, Wind = 8, Wleaf = 0.01, 
                   LeafAbs = 0.5)
D_leaf = plantecophys::VPDairToLeaf(VPD = VPD, Tair = Tair, 
                                    Tleaf = Tleaf)
g_w = calc_gw(E, Tleaf = Tleaf, Tair = Tair, VPD = VPD, PPFD = PPFD,
              Wind = 8, Wleaf = 0.01)
Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2))

# When calculating net A, gives bizarre T response
# Can be resolved by calculating gross and subtracting out Rd, but values are very different
Photo_out = Photosyn_v(Ca = 420, GS = g_w, Tleaf = Tleaf, VPD = 3, PPFD = 1800,
                       #Rd0 = 0,
                       Vcmax=34,EaV=62307,EdVC=2e5,delsC=639, # delsC causes weird Ac response
                       Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635
                       )
As_pred = unlist(Photo_out[2,])# - Rd
Ac_pred = unlist(Photo_out[5,])
Aj_pred = unlist(Photo_out[6,])
Rd_pred = unlist(Photo_out[8,])
Ci_pred = unlist(Photo_out[1,])


plot(E, As_pred)
plot(E, Aj_pred, col = "blue")
points(E, Ac_pred, col = "orange")
points(E, Rd_pred, col = "grey")

plot(E, g_w)
plot(E,Tleaf)
plot(E, Ci_pred)

# Limit of Ci as gs -> inf = Ca
# But concave down if A > 0 and concave up if A < 0 
A = -5
Ca = 420
gs = seq(0,0.01, length.out = 500)
Ci = Ca - A/gs
plot(gs, Ci)


Tleaf = 45
Vcmax = 34 * TVcmax(Tleaf = Tleaf,EaV=62307,EdVC=2e5,delsC=639)
GammaStar = TGammaStar(Tleaf = Tleaf, Patm = 100)
Km = TKm(Tleaf = Tleaf, Patm = 100)

plot(Vcmax*(Ci_pred - GammaStar)/(Ci_pred + Km))
plot(Rd_pred)
plot(Vcmax)
cis = seq(0,600, by = 1)
plot(cis, Vcmax[1]*(cis - GammaStar[1]) / (cis + Km[1]))

