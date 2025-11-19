# Functions for analysis
# by Camille Sicangco
# Created 30 Sept 2024

# Make predictions with prescribed values of E
prescribedE_pred = function(df, Wind = 5, Wleaf = 0.025, LeafAbs=0.5, 
                            Rd0 = 0.92, TrefR = 25,...)
{
  Tleaf = try(calc_Tleaf(Tair = df$Tair, E = df$E, VPD = df$VPD, PPFD = df$PPFD, Wind = Wind, 
                     Wleaf = Wleaf, LeafAbs = LeafAbs))
  Tleaf[grep("Error", Tleaf)] = NA
  Tleaf = as.numeric(Tleaf)
  Dleaf = try(VPDairToLeaf(Tleaf = Tleaf, Tair = df$Tair, VPD = df$VPD))
  #gs = try(calc_gw(E = df$Trans, D_leaf = Dleaf))
  gs = try(calc_gw(E = df$E, Tleaf = Tleaf, Tair = df$Tair, VPD = df$VPD, 
                   PPFD = df$PPFD, Wind = Wind, Wleaf = Wleaf))
  Agr =  mapply(plantecophys::Photosyn, VPD = df$VPD, 
              Ca = 400, PPFD = df$PPFD, Tleaf = df$Tcan,  
              GS = gs,
              Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
              Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633, Rd0 = 0)
  Rd = Rd0 * exp(0.1012 * (Tleaf - TrefR) - 5e-04 * (Tleaf^2 - 
                                                       TrefR^2))
  A = as.numeric(Agr[2, ]) - Rd
  out = data.frame(Tleaf.pred = Tleaf, Dleaf.pred = Dleaf, gs.pred = gs, A.pred = A,
                   Tleaf.obs = df$Tcan, Dleaf.obs = df$Dleaf, gs.obs = df$gs, A.obs = df$A, E.obs = df$E)
  out = out %>% 
    mutate(Tleaf.diff = Tleaf.pred - Tleaf.obs, 
           Dleaf.diff = Dleaf.pred - Dleaf.obs,
           gs.diff = gs.pred - gs.obs,
           A.diff = A.pred - A.obs)
  return(out)
}


# Generate predictions for all models
make_pred = function(
    Tair, Ps, VPD, PPFD, 
    P50 = 4.07, P88 = 5.5,
    Wind = 5,
    Wleaf = 0.025,
    LeafAbs=0.5,
    Tcrit = 46.5,
    T50 = 59.6,
    kmax_25 = 0.5,
    Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
    Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635,
    Rd0 = 0.92,
    constr_Ci = FALSE,
    ...) {
  
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  P = Ps_to_Pcrit(Ps, Pcrit, pts = 500)
  
  
  # Calculate costs and gains
  cost_gain = calc_costgain(P, b, c, kmax_25 = kmax_25, Tair = Tair, VPD = VPD, 
                            PPFD = PPFD, Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                            Tcrit = Tcrit, T50 = T50, 
                            Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                            Jmax = Jmax,EaJ=EaJ,EdVJ=EdVC,delsJ=delsJ,
                            Rd0 = Rd0, constr_Ci = constr_Ci, ...)

  models = c("Sperry", "Sperry + varkmax", "Sperry + CGnet", "Sperry + CGnet + TC")

  
  # Generate model predictions
  preds = sapply(models, 
                 function(x) make_pred_fn(model = x, 
                                          cost_gain = cost_gain, 
                                          P = P, b = b, c = c,
                                          Tair = Tair, Ps = Ps, VPD = VPD, 
                                          PPFD = PPFD, kmax_25 = kmax_25, 
                                          Wind = Wind, Wleaf = Wleaf, 
                                          LeafAbs = LeafAbs,
                                          Tcrit = Tcrit, T50 = T50,
                                          ...)
  )
  preds = data.frame(t(preds))
  rownames(preds) = NULL
  preds = preds %>% mutate(across(Tair:A, as.numeric))
  return(preds)
}

# Make predictions with Sperry or Sicangco models for a single observation
make_pred_fn = function(P, b, c,
                         cost_gain,
                         Tair, 
                         Ps, 
                         VPD, 
                         PPFD, 
                         Wind, 
                         LeafAbs = 0.5,
                         Wleaf = 0.025,
                         Tcrit = 46.5,
                         T50 = 50.4,
                         P50 = 4.07,
                         P88 = 5.50,
                         kmax_25 = 0.5,
                         Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                         Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635,
                         Rd0 = 0.92,
                         model,
                         constr_Ci = FALSE,
                         #constant_kmax = TRUE,
                         ...) {
  
  
  # Solve optimization
  if (model == "Sperry") {
    i = which.max(cost_gain$CG_gross_constkmax - cost_gain$HC_constkmax)
  } else if (model == "Sperry + varkmax") {
    i = which.max(cost_gain$CG_gross_varkmax - cost_gain$HC_varkmax)
  } else if (model == "Sperry + CGnet") {
    i = which.max(cost_gain$CG_net - cost_gain$HC_varkmax)
  } else if (model == "Sperry + CGnet + TC") {
    i = which.max(cost_gain$CG_net_newJT - (cost_gain$HC_varkmax + cost_gain$TC))
  } else if (model == "Sperry + TC") {
    i = which.max(cost_gain$CG_gross_varkmax - (cost_gain$HC_varkmax + cost_gain$TC))
  } else {
    stop()
  }
  
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
    
    E_vec = if (model %in% c("Sperry", "Sperry + CGnet_constrCi")) {
      trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
    } else {
      trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = FALSE)
    }
    
    P = P[i]
    E = E_vec[i]
    Tleaf = calc_Tleaf(Tair = Tair, E = E, VPD = VPD, PPFD = PPFD, Wind = Wind, 
                       Wleaf = Wleaf, LeafAbs = LeafAbs)
    Dleaf = VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)
    gs = calc_gw(E = E, Tleaf = Tleaf, Tair = Tair, VPD = VPD, 
                 PPFD = PPFD, Wind = Wind, Wleaf = Wleaf)
    
    new_JT = if (model %in% c("Sperry + TC", "Sperry + CGnet + TC")) {
      TRUE
    } else {
      FALSE
    }
    
    A = calc_A(Tleaf = Tleaf, g_w = gs, VPD = VPD, net = TRUE, netOrig = TRUE,
               PPFD = PPFD, Wind = Wind, 
               Wleaf = Wleaf, LeafAbs = LeafAbs,
               Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
               Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, new_JT = new_JT,
               ...)
    
    out = c(model, Tair, Ps, P, E, Tleaf, Dleaf, gs, A)
    
    names(out) = c("Model", "Tair", "Ps", "P", "E", "Tleaf", "Dleaf", "gs", "A")
    return(out)
}

# Calculate Pleaf for the Medlyn model
calc_Pleaf0 = function(E, Ps, P50 = 4.07, P88 = 5.50) {
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  P = Ps_to_Pcrit(Ps, Pcrit, pts = 600)
  E_vec = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
  
  # Get corresponding Pleaf
  i = which.min(abs(E - E_vec))
  Pleaf = P[i]
  
  return(Pleaf)
}

# Custom version of Photosyn function using Heskel et al. 2017 equation for R(T)
Photosyn_custom <- function(VPD=1.5, 
                     Ca=400, 
                     PPFD=1500,
                     Tleaf=25,
                     Patm=101.325,
                     RH=NULL,
                     
                     gsmodel=c("BBOpti","BBLeuning","BallBerry","BBdefine"),
                     g1=4,
                     g0=0, 
                     gk=0.5,
                     vpdmin=0.5,
                     D0=5,
                     GS=NULL,
                     BBmult=NULL,
                     
                     alpha=0.24, 
                     theta=0.85, 
                     Jmax=100, 
                     Vcmax=50, 
                     gmeso=NULL,
                     TPU=1000,
                     alphag=0,
                     
                     Rd0 = 1.115,
                     Q10 = 1.92,
                     Rd=NULL,
                     TrefR = 25,
                     Rdayfrac = 1.0,
                     
                     EaV = 58550,
                     EdVC = 200000,
                     delsC = 629.26,
                     
                     EaJ = 29680,
                     EdVJ = 200000,
                     delsJ = 631.88,
                     
                     GammaStar = NULL,
                     Km = NULL,
                     
                     Ci = NULL,
                     Tcorrect=TRUE,  
                     returnParsOnly=FALSE,
                     whichA=c("Ah","Amin","Ac","Aj"),
                     
                     Tcrit = 46.5,
                     T50 = 50.4,
                     
                     new_JT = TRUE,
                     
                     b_USO = 0.55, # sensitivity of g1 to SWP
                     Ps = 0.5 # soil water potential, -MPa
                     ){
  
  
  whichA <- match.arg(whichA)
  gsmodel <- match.arg(gsmodel)
  if(gsmodel == "BBdefine" && is.null(BBmult)){
    Stop("When defining your own BB multiplier, set BBmult.")
  }
  inputCi <- !is.null(Ci)
  inputGS <- !is.null(GS)
  
  if(inputCi & inputGS)Stop("Cannot provide both Ci and GS.")
  
  if(is.null(TPU))TPU <- 1000
  
  if(is.null(VPD) && !is.null(RH)){
    VPD <- RHtoVPD(RH, Tleaf)
  } 
  if(is.null(VPD) && is.null(RH)){
    Stop("Need one of VPD, RH.")
  }
  
  #---- Constants; hard-wired parameters.
  Rgas <- .Rgas()
  GCtoGW <- 1.57     # conversion from conductance to CO2 to H2O
  
  
  #---- Do all calculations that can be vectorized
  
  # g1 and g0 are input ALWAYS IN UNITS OF H20
  # G0 must be converted to CO2 (but not G1, see below)
  g0 <- g0/GCtoGW
  
  # Leaf respiration
  if(is.null(Rd)){
    Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2))
  }
  
  # CO2 compensation point in absence of photorespiration
  if(is.null(GammaStar))GammaStar <- TGammaStar(Tleaf,Patm)
  
  # Michaelis-Menten coefficient
  if(is.null(Km))Km <- TKm(Tleaf,Patm)
  
  #-- Vcmax, Jmax T responses
  if(Tcorrect){
    Vcmax <- Vcmax * TVcmax(Tleaf,EaV, delsC, EdVC)
    Jmax <- 
      if(isTRUE(new_JT)) {
        Jmax * TJmax_updated(Tleaf, EaJ, delsJ, EdVJ, Tcrit, T50)
      } else {
        Jmax * TJmax(Tleaf, EaJ, delsJ, EdVJ)  
        }
  }
  
  # Electron transport rate
  J <- Jfun(PPFD, alpha, Jmax, theta)
  VJ <- J/4
  
  #--- Stop here if only the parameters are required
  if(returnParsOnly){
    return(list(Vcmax=Vcmax, Jmax=Jmax, Km=Km, GammaStar=GammaStar, VJ=VJ))
  }
  
  # Medlyn et al. 2011 model gs/A. NOTE: 1.6 not here because we need GCO2!
  if(gsmodel == "BBOpti"){
    vpduse <- VPD
    vpduse[vpduse < vpdmin] <- vpdmin
    
    g1 = g1 * exp(b_USO * (-Ps + 0.2))
    
    GSDIVA <- (1 + g1/(vpduse^(1-gk)))/Ca
  }
  
  # Leuning 1995 model, without gamma (CO2 compensation point)
  if(gsmodel == "BBLeuning"){
    GSDIVA <- g1 / Ca / (1 + VPD/D0)
    GSDIVA <- GSDIVA / GCtoGW   # convert to conductance to CO2
  }
  
  # Original Ball&Berry 1987 model.
  if(gsmodel == "BallBerry"){
    if(is.null(RH))RH <- VPDtoRH(VPD, Tleaf)
    RH <- RH / 100
    GSDIVA <- g1 * RH / Ca
    GSDIVA <- GSDIVA / GCtoGW   # convert to conductance to CO2
  } 
  
  # Multiplier is user-defined.
  if(gsmodel == "BBdefine"){
    GSDIVA <- BBmult / GCtoGW
  }
  
  if(inputGS){
    
    # Convert stomatal conductance to water vapor -> carbon
    Gsc <- GS / GCtoGW
    
    # Calculate boundary layer conductance
    Tair_k <- Tair + 273.15
    CMOLAR <- Patm * 1000/(8.314 * Tair_k)
    Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR
    GRASHOF <- 1.6e+08 * abs(Tleaf - Tair) * (Wleaf^3)
    Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25)/Wleaf * CMOLAR
    Gbh <- 2 * (Gbhfree + Gbhforced)
    Gbw = GbH / 0.93
    Gbc <- Gbw / 1.37 # Aphalo and Jarvis 1993
    
    GC = (Gsc * Gbc) / (Gsc + Gbc)
    
    if(GS > 0){
      
      # Solution when Rubisco activity is limiting
      A <- 1./GC
      B <- (Rd - Vcmax)/GC - Ca - Km
      C <- Vcmax * (Ca - GammaStar) - Rd * (Ca + Km)
      Ac <- QUADM(A,B,C)
      
      # Photosynthesis when electron transport is limiting
      B <- (Rd - VJ)/GC - Ca - 2*GammaStar
      C <- VJ * (Ca - GammaStar) - Rd * (Ca + 2*GammaStar)
      Aj <- QUADM(A,B,C)
      
      # NOTE: the solution above gives net photosynthesis, add Rd
      # to get gross rates (to be consistent with other solutions).
      Ac <- Ac + Rd
      Aj <- Aj + Rd
    } else {
      Ac <- Aj <- 0
    }
    
  } else {
    
    # If CI not provided, calculate from intersection between supply and demand
    if(!inputCi){
      
      #--- non-vectorized workhorse
      getCI <- function(VJ,GSDIVA,PPFD,VPD,Ca,Tleaf,vpdmin,g0,Rd,
                        Vcmax,Jmax,Km,GammaStar){
        
        if(identical(PPFD, 0) | identical(VJ, 0)){
          vec <- c(Ca,Ca)
          return(vec)
        }
        
        # Taken from MAESTRA.
        # Following calculations are used for both BB & BBL models.
        # Solution when Rubisco activity is limiting
        A <- g0 + GSDIVA * (Vcmax - Rd)
        B <- (1. - Ca*GSDIVA) * (Vcmax - Rd) + g0 * 
          (Km - Ca)- GSDIVA * (Vcmax*GammaStar + Km*Rd)
        C <- -(1. - Ca*GSDIVA) * (Vcmax*GammaStar + Km*Rd) - g0*Km*Ca
        
        CIC <- QUADP(A,B,C)
        
        # Solution when electron transport rate is limiting
        A <- g0 + GSDIVA * (VJ - Rd)
        B <- (1 - Ca*GSDIVA) * (VJ - Rd) + g0 * (2.*GammaStar - Ca)- 
          GSDIVA * (VJ*GammaStar + 2.*GammaStar*Rd)
        C <- -(1 - Ca*GSDIVA) * GammaStar * (VJ + 2*Rd) - 
          g0*2*GammaStar*Ca
        
        CIJ <- QUADP(A,B,C)
        return(c(CIJ,CIC))
      }
      
      # get Ci
      x <- mapply(getCI, 
                  VJ=VJ,
                  GSDIVA = GSDIVA,
                  PPFD=PPFD,
                  VPD=VPD,
                  Ca=Ca,
                  Tleaf=Tleaf,
                  vpdmin=vpdmin,
                  g0=g0,
                  Rd=Rd,
                  Vcmax=Vcmax,
                  Jmax=Jmax,
                  Km=Km,
                  GammaStar=GammaStar)
      
      CIJ <- x[1,]
      CIC <- x[2,]
    } else {
      
      # Rare case where one Ci is provided, and multiple Tleaf (Jena bug).
      if(length(Ci) == 1){
        Ci <- rep(Ci, length(Km))
      }
      
      # Ci provided (A-Ci function mode)
      CIJ <- Ci
      
      if(length(GammaStar) > 1){
        CIJ[CIJ <= GammaStar] <- GammaStar[CIJ < GammaStar]
      } else {
        CIJ[CIJ <= GammaStar] <- GammaStar
      }
      
      CIC <- Ci
      
    }
    
    # Photosynthetic rates, without or with mesophyll limitation
    if(is.null(gmeso) || gmeso < 0){
      # Get photosynthetic rate  
      Ac <- Vcmax*(CIC - GammaStar)/(CIC + Km)
      Aj <- VJ * (CIJ - GammaStar)/(CIJ + 2*GammaStar)
      
    } else {
      # Ethier and Livingston (2004) (Equation 10).
      A <- -1/gmeso
      BC <- (Vcmax - Rd)/gmeso + CIC + Km
      CC <- Rd*(CIC+Km)-Vcmax*(CIC-GammaStar)
      Ac <- mapply(QUADP, A=A,B=BC,C=CC)
      
      BJ <- (VJ - Rd)/gmeso + CIC + 2.0*GammaStar
      CJ <- Rd*(CIC+2.0*GammaStar) - VJ*(CIC - GammaStar)
      Aj <- mapply(QUADP, A=A,B=BJ,C=CJ)
      
      Ac <- Ac + Rd
      Aj <- Aj + Rd
      
    }
    
    
    # When below light-compensation points, assume Ci=Ca.
    if(!inputCi){
      lesslcp <- vector("logical", length(Aj))
      lesslcp <- Aj <= Rd + 1E-1
      
      if(length(Ca) == 1)Ca <- rep(Ca, length(CIJ))
      if(length(GammaStar) == 1)GammaStar <- rep(GammaStar, length(CIJ))
      if(length(VJ) == 1)VJ <- rep(VJ, length(CIJ))
      
      CIJ[lesslcp] <- Ca[lesslcp]
      Aj[lesslcp] <- VJ[lesslcp] * (CIJ[lesslcp] - GammaStar[lesslcp]) / 
        (CIJ[lesslcp] + 2*GammaStar[lesslcp])
      
      Ci <- ifelse(Aj < Ac, CIJ, CIC)
      
    }
  }
  
  # Limitation by triose-phosphate utilization
  if(!is.null(Ci)){
    Ap <- 3 * TPU * (Ci - GammaStar)/(Ci - (1 + 3*alphag)*GammaStar)
    Ap[Ci < 400] <- 1000  # avoid nonsense
  } else {
    Ap <- 1000  # This is when inputGS = TRUE; 
  }  
  
  
  # Hyperbolic minimum.
  Am <- -mapply(QUADP, A = 1 - 1E-04, B = Ac+Aj, C = Ac*Aj)
  
  # Another hyperbolic minimum with the transition to TPU
  tpulim <- any(Ap < Am)
  if(!is.na(tpulim) && tpulim){
    Am <- -mapply(QUADP, A = 1 - 1E-07, B = Am+Ap, C = Am*Ap)
  }
  
  # Net photosynthesis
  Am <- Am - Rd
  
  # Calculate conductance to CO2
  if(!inputCi && !inputGS){
    if(whichA == "Ah")GS <- g0 + GSDIVA*Am
    if(whichA == "Aj")GS <- g0 + GSDIVA*(Aj-Rd)
    if(whichA == "Ac")GS <- g0 + GSDIVA*(Ac-Rd)
  } 
  if(inputCi) {
    if(whichA == "Ah")GS <- Am/(Ca - Ci)
    if(whichA == "Aj")GS <- (Aj-Rd)/(Ca - Ci)
    if(whichA == "Ac")GS <- (Ac-Rd)/(Ca - Ci)
  }
  
  # Extra step here; GS can be negative
  GS[GS < g0] <- g0
  
  # Output conductance to H2O
  if(!inputGS){
    GS <- GS*GCtoGW
  }
  
  # Calculate Ci if GS was provided as input.
  if(inputGS){
    Ci <- Ca - Am/GC
    
    # For zero GC:
    Ci[!is.finite(Ci)] <- Ca
    # Stomata fully shut; Ci is not really Ca but that's 
    # how we like to think about it.
  }
  
  # Chloroplastic CO2 concentration
  if(!is.null(gmeso)){
    Cc <- Ci - Am/gmeso
  } else {
    Cc <- Ci
  }
  
  # Transpiration rate assuming perfect coupling.
  # Output units are mmol m-2 s-1
  E <- 1000*GS*VPD/Patm
  
  df <- data.frame( Ci=Ci,
                    ALEAF=Am,
                    GS=GS,
                    ELEAF=E,
                    Ac=Ac,
                    Aj=Aj,
                    Ap=Ap,
                    Rd=Rd,
                    VPD=VPD,
                    Tleaf=Tleaf,
                    Ca=Ca,
                    Cc=Cc,
                    PPFD=PPFD,
                    Patm=Patm)
  
  return(df)
}


# Function to fit temperature response of Vcmax 
# by Dushan Kumarathunge
# inputs; dat is a dataframe with Vcmax, Tleaf in Kelvin (named as "TsK")
# return; "Arr" for fitted Arrhenius model parameters, "Peak" for fitted Peak Arrhenius parameters
fitpeaked<-function(dat,return=c("Arr","Peak"),start=list(k25=100, Ea=60, delS = 0.64)){
  
  return <- match.arg(return)
  
  try(Vc<-nls(fVc, start=start, data=dat))
  Vc2<-summary(Vc)
  res1<-data.frame(cbind(Vc2$coefficients[1],Vc2$coefficients[2],Vc2$coefficients[3],Vc2$coefficients[4],Vc2$coefficients[5],
                         Vc2$coefficients[6]))
  
  names(res1)[1:6]<-c("Vcmax25","EaV","delsV","Vcmax25.se","EaV.se","delsV.se")
  
  try(Vr<-nls(fVc.a, start = list(k25=100, Ea=40), data = dat))
  Vr2<-summary(Vr)
  res<-data.frame(cbind(Vr2$coefficients[1],Vr2$coefficients[2],Vr2$coefficients[3],Vr2$coefficients[4]))
  
  
  names(res)[1:4]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se")
  
  
  an<-anova(Vr,Vc)
  AIC_Arr<-AIC(Vr)
  AIC_Peak<-AIC(Vc)
  prob_V<-an[[6]][[2]]
  
  
  r1<-cor(fitted(Vc),dat$Vcmax)
  R2_Peak<-r1*r1
  
  r2<-cor(fitted(Vr),dat$Vcmax)
  R2_Arr<-r2*r2
  
  #test for normality of residuals
  rest<-residuals(Vc)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  topt<-200/(res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  Topt<-topt-273.15
  
  Topt.se <- topt * 
    sqrt(res1[[6]]**2 + 
           (0.008314 * ((res1[[5]] / res1[[2]])**2 + (res1[[5]] / (200 - res1[[2]]))**2))**2) /
    (res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  
  
  param_peak<-cbind(res1,Topt,Topt.se,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_V)
  names(param_peak)[1:13]<-c("Vcmax25","EaV","delsV","Vcmax25.se","EaV.se","delsV.se","ToptV",
                             "ToptV.se","R2_Arr","R2_Peak","AIC_Arr","AIC_Peak","prob_V")
  
  param_arr<-cbind(res,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_V)
  names(param_arr)[1:9]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se","R2_Arr",
                           "R2_Peak","AIC_Arr","AIC_Peak","prob_V")
  
  
  if(return == "Peak")return(param_peak)
  if(return == "Arr")return(param_arr)
  
}

# Function to fit temperature response of Jmax 
# by Dushan Kumarathunge
# inputs; dat is a dataframe with Vcmax, Tleaf in Kelvin (named as "TsK")
# return; "Arr" for fitted Arrhenius model parameters, "Peak" for fitted Peak Arrhenius parameters

fitpeakedJ<-function(dat,return=c("Arr","Peak"),start=list(k25=100, Ea=60, delS = 0.64)){
  return <- match.arg(return)
  dat$TsK <- dat$Tleaf + 273.15
  
  try(Vj<-nls(fVJ, start = start, data = dat))
  Vj2<-summary(Vj)
  res1<-Vj2$coefficients[1:6]
  names(res1)[1:6]<-c("Jmax25","Ea","delsJ","Jmax25.se","EaJ.se","delsJ.se")
  
  try(Vr<-nls(fVJ.a, start = list(k25=100, Ea=40), data = dat))
  Vr2<-summary(Vr)
  res<-Vr2$coefficients[1:4]
  names(res)[1:4]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se")
  
  
  an<-anova(Vr,Vj)
  AIC_Arr<-AIC(Vr)
  AIC_Peak<-AIC(Vj)
  prob_J<-an[[6]][[2]]
  
  
  r1<-cor(fitted(Vj),dat$Jmax)
  R2_Peak<-r1*r1
  
  r2<-cor(fitted(Vr),dat$Jmax)
  R2_Arr<-r2*r2
  
  
  
  #test for normality of residuals
  rest<-residuals(Vj)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  topt<-200/(res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  Topt<-topt-273.15
  Topt.se <- topt * 
    sqrt(res1[[6]]**2 + 
           (0.008314 * ((res1[[5]] / res1[[2]])**2 + (res1[[5]] / (200 - res1[[2]]))**2))**2) /
    (res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  
  param_peak<-c(res1,Topt,Topt.se,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_J)
  names(param_peak)[1:13]<-c("Jmax25","EaJ","delsJ","Jmax25.se","EaJ.se","delsJ.se","ToptJ","ToptJ.se","R2_Arr","R2_Peak",
                             "AIC_Arr","AIC_Peak","prob_J")
  
  param_arr<-c(res,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_J)
  names(param_arr)[1:9]<-c("Jmax25","EaJ","Jmax25.se","EaJ.se","R2_Arr",
                           "R2_Peak","AIC_Arr","AIC_Peak","prob_J")
  
  
  if(return == "Peak")return(param_peak)
  if(return == "Arr")return(param_arr)
  
  #return(c(res,r2,s,pvalue))
  #return(list(res))
}

# define peak Arrhenius model for Vcmax
fVc <- as.formula(Vcmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * 
                    (1+exp((298.15*delS - 200)/(298.15*0.008314))) / 
                    (1+exp((TsK*delS-200)/(TsK*0.008314))))

# define peak Arrhenius model for Jmax
fVJ <- as.formula(Jmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * 
                    (1+exp((298.15*delS - 200)/(298.15*0.008314))) / 
                    (1+exp((TsK*delS-200)/(TsK*0.008314))))

# define standard Arrhenius model for Vcmax

fVc.a<-as.formula(Vcmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)))


# define standard Arrhenius model for Vcmax

fVJ.a<-as.formula(Jmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)))

# Generate predictions with all models
get_predictions = 
  function(
    df, Tcrit = 46.5, T50 = 50.4, P50 = 4.07, P88 = 5.50,
    Wind = 5, Wleaf = 0.025, LeafAbs = 0.5,
    Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
    Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
    kmax_25 = 0.5, #net = TRUE, netOrig = TRUE,
    g1 = 2.9, g0=1.e-5, b_USO = 0.55,
    constr_Ci = FALSE,
    ...
  ) {
    out = bind_rows(lapply(1:nrow(df), 
                 function(i) make_pred(
                   Tair = df$Tair[i], Ps = df$Ps[i], VPD = df$VPD[i], PPFD = df$PPFD[i],
                   Tcrit = Tcrit, T50 = T50, P50 = P50, P88 = P88,
                   Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                   Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC, 
                   Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, 
                   Rd0 = Rd0, kmax_25 = kmax_25, constr_Ci = constr_Ci,
                   ...)))
    
    # Make Medlyn predictions
    Medlyn_preds = plantecophys::PhotosynEB(Tair=df$Tair, VPD=df$VPD, PPFD = df$PPFD,
                                            Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                                            g1 = g1,g0=g0,
                                            Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC, 
                                            Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, Rd0 = Rd0,
                                            b_USO = b_USO, Ps = df$Ps
                                            )
    get_Pleaf_Medlyn_V = Vectorize(get_Pleaf_Medlyn, c("Ps", "Tair", "E"))
    Medlyn_preds$P = get_Pleaf_Medlyn_V(
      Ps = df$Ps, Tair = df$Tair, E = Medlyn_preds$ELEAF, 
      P50 = P50, P88 = P88, kmax_25 = kmax_25, constant_kmax = FALSE)
    
    Medlyn_preds_sim = Medlyn_preds %>%
      rename(E = ELEAF, Dleaf = VPDleaf, gs = GS, A = ALEAF) %>%
      mutate(Model = "Medlyn", .before = 1) %>%
      #mutate(P = NA) %>% 
      select(Model, Tair, E, Tleaf, Dleaf, gs, A, P)
    
    # Combine all predictions
    out_all = bind_rows(out, Medlyn_preds_sim)
    return(out_all)
  }

# Estimate Pleaf for the USO model
get_Pleaf_Medlyn = function(Ps, E, Tair, P50, P88, kmax_25, constant_kmax) {
  
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  
  P = Ps_to_Pcrit(Ps, Pcrit, pts = 600)
  E_vec = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax)
  j = which.min(abs(E - E_vec))
  Pleaf = P[j]
  
  return(Pleaf)
}


# Make predictions with different Tcrit/T50 values
make_pred_Tthresholds = function(df, Wind = 5, Wleaf = 0.025, LeafAbs=0.5,
                                 hold_Tcrit = FALSE,
                                 Thold_val = 50.4,
                                 Tvar_vals = c(46.5, 47.5, 48.5, 49.5),
                                 constr_Ci = FALSE,
                                 ...) {
  if(isTRUE(hold_Tcrit)) {
    Tcrit = Thold_val
    T50s = Tvar_vals
    
    sims.l = lapply(T50s, function(T50,...) {
      preds = get_predictions(df, Tcrit = Tcrit, T50 = T50) %>% 
        filter(Model == "Sperry + CGnet + TC")
    }
    )
    names(sims.l) = paste0("T50_", Tvar_vals)
    preds = bind_rows(sims.l, .id = "ID") %>% 
      mutate(Tcrit = Thold_val,
             T50 = str_sub(ID, 5)) %>% 
      select(!ID)
  } else {
    T50 = Thold_val
    Tcrits = Tvar_vals
    
    sims.l = lapply(Tcrits, function(Tcrit,...) {
      preds = get_predictions(df, Tcrit = Tcrit, T50 = T50) %>% 
        filter(Model == "Sperry + CGnet + TC")
    }
    )
    names(sims.l) = paste0("T50_", Tvar_vals)
    preds = bind_rows(sims.l, .id = "ID") %>% 
      mutate(T50 = Thold_val,
             Tcrit = str_sub(ID, 5)) %>% 
      select(!ID)
  }
  preds = preds %>% mutate(across(Tair:A, as.numeric),
                           across(T50:Tcrit, as.factor))
  return(preds)
}

# From plantecophys
# Jmax temperature response (Arrhenius)
TJmax <- function(Tleaf, EaJ, delsJ, EdVJ){
  J1 <- 1+exp((298.15*delsJ-EdVJ)/8.314/298.15)
  J2 <- 1+exp(((Tleaf + 273.15)*delsJ-EdVJ)/8.314/(Tleaf + 273.15))
  exp(EaJ/8.314*(1/298.15 - 1/(Tleaf + 273.15)))*J1/J2
}
