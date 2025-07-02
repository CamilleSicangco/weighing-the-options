# Modified C gain functions

# Manon/Venturas implementation
C_gain_alt = function (P, b = -2.5, c = 2, Amax = NULL, kmax_25 = 4, Tair = 25, 
                       VPD = 1.5, PPFD = 1000, Patm = 101.325, Wind = 2, Wleaf = 0.01, 
                       LeafAbs = 0.5, Ca = 420, Jmax = 100, Vcmax = 50, constant_kmax = FALSE, 
                       Rd0 = 0.92, TrefR = 25, net = TRUE, ...) 
{
  # Calculate transpiration supply stream
  E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax)
  
  # Create 2D array of possible Ci's for each value of E
  Cis = array(seq(1, Ca, length.out = 500), dim = c(500, length(E)))
  
  # Get leaf temperature and stomatal conductance for each value of E
  Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, 
                     E = E, Wind = Wind, Patm = Patm, Wleaf = Wleaf, 
                     LeafAbs = LeafAbs)
  g_w = calc_gw(E, Tleaf, Patm, Tair, VPD, PPFD, Wind, 
                Wleaf)
  
  # Calculate supply A from Fick's law
  g_ws = t(array(g_w, dim = c(500, length(E))))
  A_supply = g_ws/1.6 * (Ca - Cis)
  
  Es = t(array(E, dim = c(500, length(E))))
  
  # Calculate demand A with Farqhuar model
  Tleaves = t(array(Tleaf, dim = c(500, length(E))))
  Photosyn_out = mapply(Photosyn_custom, VPD = VPD, 
                        Ca = Ca, PPFD = PPFD, Tleaf = Tleaves, 
                        Patm = Patm, 
                        Ci = Cis, Jmax = Jmax, Vcmax = Vcmax,
                        Rd0 = Rd0, TrefR = TrefR, new_JT = FALSE,
                        ...)
  
  if (isTRUE(net)) {
    A_demand = as.numeric(Photosyn_out[2,])
  } else {
    Anet = as.numeric(Photosyn_out)[2,]
    Rd = as.numeric(Photosyn_out)[8,]
    A_demand = Anet + Rd
  }
 dim(A_demand) = c(500, 500)
  
  # Find which values of supply and demand A match closest for each value of E
  idx = apply(abs(A_supply - A_demand), 2, FUN = which.min)
  
  # Select the corresponding values of Ci and A
  Ci = sapply(1:length(E), function(e) {Cis[idx[e], e]})
  A_P = sapply(1:length(E), function(e) {A_demand[idx[e], e]})
  
  # Calculate Amax
  Amax = if (is.null(Amax) & !(all(A_P <= 0))) {
    max(abs(A_P))
  } else if (is.null(Amax) & (all(A_P <= 0))) {
    max(abs(A_P))
  } else {
    Amax
  }
  
  # Normalise A to get CG
  gain = A_P/Amax
  return(gain)
}
