# FORGAP.MOD
# FORest GAP model
# (c) O. van Dam, Utrecht University, Tropenbos-Guyana
# Version 8.1, jan 2001
# combination of GAPRAD.MOD and POTEVAP.MOD and WATBAL.MOD
 
# -----------------------------------------------------
#                 Model for calculation
# incoming potential light energy, potential and actual
#     evapotranspiration and the soil water balance
#  in a Tropical Rain Forest Gap and surrounding edges
# -----------------------------------------------------
 
# INPUT:  - DEM study area
#         -location attributes
# OUTPUT: - hour reports on the sample locations
#         - daily and monthly sums of:
#           - potential global vegetation radiation [MJ/m2]
#           - potential global soil radiation [MJ/m2]
#           - potential and actual vegetation interception evaporation,
#             litter interception evaporation, soil evaporation and
#             transpiration by the vegetation
#           - sums of soil water suction levels per time step
# NOTE:
# all evaporation calculations in POTEVAP.MOD in mm
# all soil water calculations in WATBAL.MOD in cm
# all reports of actual evaporations in WATBAL.MOD in mm
 
binding
#___________________________________________________________________________
# All models: basic maps
 Gap  = e:\pcrp\maps\g5cal.map;  # Boolean map gap: 1 forest, 0 gap
 Loc  = e:\pcrp\maps\g5fdr2.loc;    # Nominal map with sample locations
 Alti = scalar(50);                # Altitude study area [m]
 Tree = scalar(28);                # Tree height forest [m]
# Time series
 Meteo  = e:\pcrp\meteo\m98d147.txt; # Timeseries with climate data
 StartDay = 136;        # 96:292, 97:1, 98:1
 Year   = 1998;       # 1996, 1997, 1998 or 9198
 LastTS = 96;       # 96:1800, 97:8760, 98:7565, 99:1920 calib g5:3193,1162 fdr 504
 
# Calculated maps from previous years
## Initial soil moisture
 ThetaI1 = .184847;
 ThetaI2 = .300000;                  # use in 96
 ThetaI3 = .400000;
 ThetaI4 =  ThetaI3; ThetaI5 =  ThetaI3;  # use in 96
# ThetaI1 = Theta1.97; ThetaI2 = Theta2.97; # use in 97,98,99
# ThetaI3 = Theta3.97; ThetaI4 = Theta4.97; # use in 97,98,99
# ThetaI5 = Theta5.97;                         # use in 97,98,99
## Net radiation sum  [MJ/m2]
# GrowI = 0;               # use in 96
 GrowI = Grow.98;         # use in 97,98,99
## Litter mass
# LMini   = 0;             # use in 96
 LMini   = LM.98;         # use in 97,98,99
 
# GAPRAD.MOD
# ========================================================================= Lat      = 5.03;            # latitude study area
 SlopMap  = 0.01;            # slope area when no slope calculated: default
 pi       = 3.1415;          # pi
 Sc       = 1367.0;          # Solar constant (Duffie & Beckmen 1991) [W/m2]
 Trans    = 0.6;             # Transmissivity tau (Gates, 1980)
 Cext     = 0.295;            # Gap edge extincting factor
 Kext     = 0.29;            # radiation extinction factor
 
# POTEVAP.MOD
# =========================================================================# Penman-Monteith
 Cp       = 1013;             # spec. heat moist air at constant P [J/kg/K]
 Epsilon  = 0.622;            # ratio molecul weight water vap / dry air
 Roa      = 1.2047;           # air density [kg/m3]
 StefBo   = 5.67E-8;          # Stefan Boltzman constant [W/m2/K4]
# Aerodynamic resistance, wind, canopy resistance
 Hwind    = 10.0;             # height windspeed measurements
 Ka       = 0.41;             # von Karman constant [-]
 RaMaxVeg = 80.0;             # maximum aerodynamic resistance
 RaMaxSap = 150.0;
# Rainfall interception
 CanStor  = 0.1505;           # Can stors 0.89[mm]/LAI 5.913[m2/m2]
 StemStor = 0.0373;           # Stem storage [mm] (Jetten 1994)
 StemCoef = 0.0118;
# Litter mass
 LMforest = 228.1;            # avg litter mass [g/m2]:
# Vegetation parameters
# MMF sigmoid model for sapling growth par.a
 MMFa     = 0.2569; MMFb = 71403490; MMFc = 28.3784; MMFd = 2.2682;
# parameters in edge litter fall logistic func.
 LFa      = 5.6133; LFb = -0.0439; LFc = 0.2037;
 
# WATBAL.MOD
# =========================================================================# Constants all layers
 R        = 8.3144;         # Universal gas constant [J/mol/K]
 M        = 0.018;          # Molecule weight of water [kg/mol]
 G        = 9.8;            # Gravitation acceleration [m/s2]
 H50      = 3500.0;         # Press head with 50% reduction transpi [cm]
 RedT     = 1.5;            # Transpiration reduction function steepness [-]
 Hmax     = 100000.0;       # Max press head: H > Hmax -> no ET [cm]
 TeMin    = 0.001;          # min theta eff.
# Constants per layer (4 layers)
# Top layer
 ThetaR  = 0.01;            # Residual soil moisture [cm3/cm3] L2
 ThetaS1g = 0.361;            # Saturated soil moisture [cm3/cm3] L1
 ThetaS1f = 0.450;            # Saturated soil moisture [cm3/cm3] L1
 Ks1g     = 8.342;          # Sat hydr cond [cm/h] L1
 Ks1f     = 15.000;          # Sat hydr cond [cm/h] L1
 n1g      = 1.500;           # Mualem n [-] L1
 n1f      = 1.500;           # Mualem n [-] L1
 a1g      = 0.100;           # Mualem aplha [-] L1
 a1f      = 0.070;           # Mualem aplha [-] L1
 D1       = 15.0;            # Depth first layer [cm]
 R1g      = 0.79;           # Fraction of all roots in gap in L1
 R1f      = 0.78;        # Fraction of all roots in forest in L1
# Other layers
 ThetaS2g = 0.501;  ThetaS3g = 0.460;   ThetaS4g = ThetaS3g;
 ThetaS2f = 0.512;  ThetaS3f = 0.437;   ThetaS4f = ThetaS3f;
 Ks2g     = 9.390;  Ks3g     = 7.830;   Ks4g     = Ks3g;
 Ks2f     = 13.00;  Ks3f     = 10.988;  Ks4f     = Ks3f;
 n2g      = 1.300;  n3g      = 1.300;   n4g      = n3g;
 n2f      = 1.300;  n3f      = 1.300;   n4f      = n3f;
 a2g      = 0.051;  a3g      = 0.070;   a4g      = a3g;
 a2f      = 0.020;  a3f      = 0.070;   a4f      = a3f;
 D2       = 25.0;     D3       = 40.0;      D4       = 60.0;
 R2g      = 0.11;  R3g      = 0.05;   R4g      = 0.05;
 R2f      = 0.10;  R3f      = 0.07;   R4f      = 0.05;
 
areamap
#___________________________________________________________________________
 Gap;
 
timer
 1 LastTS 1;
# day = 24+24..endtime;                                             # Day sum
# d30 = 720+720..endtime;
# grow = 888, 1800;              # 1996
# grow = 2232, 5616, 8760;      # 1997
# grow = 3192;
# grow = 960, 4368, 7560, 8760; # 1998
# grow = 1920;                  # 1999
# mon = 744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760; # Month sum
# yr  = endtime;                                                    # Year sum
 
initial
# Time calculation loops & DEM
#---------------------------------
 Day = StartDay;
 Hour = scalar(0); Month = scalar(0);
 GapScal = scalar(Gap);                      # scalar gap: 1 forest, 0 gap
 ForH    = Tree * GapScal;                   # forest height [m]
 ScalInv = if(Gap eq 0,scalar(1),scalar(0)); # scalar gap: 0 forest, 1 gap
 BoolInv = boolean(ScalInv);                 # boolean gap: 0 forest, 1 gap
 test    = GapScal+ScalInv;
 
# GAPRAD.MOD
# =========================================================================# initial maps radiation
 DEdge    = spread(Gap,0,1);                # distance from edge into gap
 DForest  = spread(BoolInv,0,1);            # distance from edge into forest
 DFE      = DEdge - DForest;
 Edge     = if(DEdge gt 3, scalar(0), scalar(1));
 Edge     = Edge * ScalInv;                   # edge gap which lies in gap
 EdgeBool = boolean(Edge);
 ForGap   = GapScal + Edge;                   # forest + edge=1, gap=0
 
# LAI forest & can
 LAIfor   = test*(0.2076 * Tree + 0.1);              # LAI dense forest
 LAIcan   = 0.2076 * (Tree - 0.5) + 0.1;      # LAI canopy above 0.5m sap
 
# LDD for soil radiation
 Ldd1 = ldd(1); Ldd1 = lddrepair(Ldd1); Ldd2 = ldd(2); Ldd2 = lddrepair(Ldd2);
 Ldd3 = ldd(3); Ldd3 = lddrepair(Ldd3); Ldd4 = ldd(4); Ldd4 = lddrepair(Ldd4);
 Ldd6 = ldd(6); Ldd6 = lddrepair(Ldd6); Ldd7 = ldd(7); Ldd7 = lddrepair(Ldd7);
 Ldd8 = ldd(8); Ldd8 = lddrepair(Ldd8); Ldd9 = ldd(9); Ldd9 = lddrepair(Ldd9);
# begin maps
 Sold    = scalar(0);  SsapOld = scalar(0);  SsoilO   = scalar(0);
 VegRpY  = scalar(0);  VegRpM  = scalar(0);  VegRpD   = scalar(0);
 SapRpY  = scalar(0);  SapRpM  = scalar(0);  SapRpD   = scalar(0);
 SoilRpY = scalar(0);  SoilRpM = scalar(0);  SoilRpD  = scalar(0);
 
# POTEVAP.MOD
# =========================================================================# Sapling growth with net rad function
 Grow    = GrowI;
# Forest floor litter LM LF LO
 LMgap   = 183.83*exp(-0.1002*DFE);                # LM in gap R2=0.70
 LM      = if(LMgap ge LMforest, LMforest, LMgap); # 1996 only
 LM      = if(Year eq 1996, LM, LMini);
 Kfactor = LFa*100 / LMforest;                     # k-factor for decomp
 
# begin maps
 ETpOld = scalar(0);
 PonStemOver = scalar(0); ESTD = scalar (0);
 ESTM = scalar(0); ESTY = scalar(0); RainY = scalar(0);
 OverCan = scalar(0); OverLM = scalar(0);
 RHo = scalar(0.99);  TempO = scalar(22.0);
 RsVegY = scalar(0);  RlVegY = scalar(0);
 RnetVeg = scalar(0); RnetSoil = scalar(0); EWy = scalar(0);
 RnVegY = scalar(0);  RnVegM = scalar(0);  RnVegD = scalar(0);
 RnSapY = scalar(0);  RnSapM = scalar(0);  RnSapD = scalar(0);
 RnSoilY = scalar(0); RnSoilM = scalar(0); RnSoilD = scalar(0);
 EIY = scalar(0);     EIM = scalar(0);     EID = scalar(0);
 ELY = scalar(0);     ELM = scalar(0);     ELD = scalar(0);
 ETpY = scalar(0);    ETpM = scalar(0);    ETpD = scalar(0);
 ETsapY = scalar(0);  ETsapM = scalar(0);  ETsapD = scalar(0);
 ESpY = scalar(0);    ESpM = scalar(0);    ESpD = scalar(0);
 PSY = scalar(0);     PSM = scalar(0);     PSD = scalar(0);
 EPY = scalar(0);     EPM = scalar(0);     EPD = scalar(0);
 STY  = scalar(0);    STM = scalar(0);     STD = scalar(0);
 THY  = scalar(0);    THM = scalar(0);     THD = scalar(0);
 ETaL1Y = scalar(0);  ETaL2Y = scalar(0);  ETaL3Y = scalar(0);
 ETaL1M = scalar(0);  ETaL2M = scalar(0);  ETaL3M = scalar(0);
 ETaL1D = scalar(0);  ETaL2D = scalar(0);  ETaL3D = scalar(0);
 ETaL4D = scalar(0);  ETaL4M = scalar(0);  ETaL4Y = scalar(0);
 ETaY = scalar(0);    ETaM = scalar(0);    ETaD = scalar(0);
 RSap = scalar(0);    RSoil = scalar(0);   DECY = scalar(0);
 LFY = scalar(0);     LFedgY = scalar(0);  LFgapY = scalar(0);
 
# WATBAL.MOD
# ========================================================================= R1 = GapScal*R1f + ScalInv*R1g;  R2 = GapScal*R2f + ScalInv*R2g;
 R3 = GapScal*R3f + ScalInv*R3g;  R4 = GapScal*R4f + ScalInv*R4g;
 ThetaS1 = GapScal*ThetaS1f + ScalInv*ThetaS1g;
 ThetaS2 = GapScal*ThetaS2f + ScalInv*ThetaS2g;
 ThetaS3 = GapScal*ThetaS3f + ScalInv*ThetaS3g;
 ThetaS4 = GapScal*ThetaS4f + ScalInv*ThetaS4g;
 Ks1 = GapScal*Ks1f + ScalInv*Ks1g;  Ks2 = GapScal*Ks2f + ScalInv*Ks2g;
 Ks3 = GapScal*Ks3f + ScalInv*Ks3g;  Ks4 = GapScal*Ks4f + ScalInv*Ks4g;
 n1 = GapScal*n1f + ScalInv*n1g;  n2 = GapScal*n2f + ScalInv*n2g;
 n3 = GapScal*n3f + ScalInv*n3g;  n4 = GapScal*n4f + ScalInv*n4g;
 a1 = GapScal*a1f + ScalInv*a1g;  a2 = GapScal*a2f + ScalInv*a2g;
 a3 = GapScal*a3f + ScalInv*a3g;  a4 = GapScal*a4f + ScalInv*a4g;
 Theta1 = ThetaI1;  Theta2 = ThetaI2;  Theta3 = ThetaI3;
 Theta4 = ThetaI4;  Theta5 = ThetaI5;
 m1 = 1-1/n1;  m2 = 1-1/n2;  m3 = 1-1/n3;  m4 = 1-1/n4;  # Mualem m [-]
 SMI = Theta1*D1 + Theta2*D2 + Theta3*D3 + Theta4*D4;
 ThetaSR1 = (ThetaS1 - ThetaR);         ThetaSR2 = (ThetaS2 - ThetaR);
 ThetaSR3 = (ThetaS3 - ThetaR);         ThetaSR4 = (ThetaS4 - ThetaR);
 ThetaE1 = (Theta1-ThetaR)/(ThetaSR1);  ThetaE2 = (Theta2-ThetaR)/(ThetaSR2);
 ThetaE3 = (Theta3-ThetaR)/(ThetaSR3);  ThetaE4 = (Theta4-ThetaR)/(ThetaSR4);
 Psi4     = (((1/min(1,ThetaE4))**(1/m4)-1)**(1/n4))/a4;
 Kpsi4    = Ks4*(1-(a4*Psi4)**(n4-1)*(1+(a4*Psi4)**n4)**-m4)**2
            /(1+(a4*Psi4)**n4)**(m4/2);
 PF02L1 = scalar(0); PF01L1 = scalar(0); PF1L1  = scalar(0);
 PF2L1  = scalar(0); PF3L1  = scalar(0); PF4L1  = scalar(0);
 PF02L2 = scalar(0); PF01L2 = scalar(0); PF1L2  = scalar(0);
 PF2L2  = scalar(0); PF3L2  = scalar(0); PF4L2  = scalar(0);
 PF02L3 = scalar(0); PF01L3 = scalar(0); PF1L3  = scalar(0);
 PF2L3  = scalar(0); PF3L3  = scalar(0); PF4L3  = scalar(0);
 PF02L4 = scalar(0); PF01L4 = scalar(0); PF1L4  = scalar(0);
 PF2L4  = scalar(0); PF3L4  = scalar(0); PF4L4  = scalar(0);
 ETa1Y  = scalar(0); ETa1M  = scalar(0); ETa1D  = scalar(0);
 ETa2Y  = scalar(0); ETa2M  = scalar(0); ETa2D  = scalar(0);
 ETa3Y  = scalar(0); ETa3M  = scalar(0); ETa3D  = scalar(0);
 ETa4Y  = scalar(0); ETa4M  = scalar(0); ETa4D  = scalar(0);
 ETaY   = scalar(0); ETaM   = scalar(0); ETaD   = scalar(0);
 ESaY   = scalar(0); ESaM   = scalar(0); ESaD   = scalar(0);
 LoY    = scalar(0); LoM    = scalar(0); LoD    = scalar(0);
 Pond   = scalar(0);
 PondD  = scalar(0); PondM  = scalar(0); PondY  = scalar(0);
 RunOD  = scalar(0); RunOM  = scalar(0); RunOY  = scalar(0);
 EpondD = scalar(0); EpondM = scalar(0); EpondY = scalar(0);
 Theta1old = Theta1; Theta2old = Theta2;
 Theta3old = Theta3; Theta4old = Theta4;
 
dynamic
# Daily calculation loop all models
# ---------------------------
 Day    = if(Hour ge 24, Day+1, Day);
 Hour   = if(Hour ge 24, scalar(1), Hour+1);
 Month  = if(Month ge 720, scalar(1), Month+1) ;
 
# Sapling height growth as function of RnetVeg
# -----------------------------------------
 SapH     = (MMFa*MMFb+MMFc*Grow**MMFd)/(MMFb+Grow**MMFd);
 VegH     = if(ForH+SapH*ScalInv<=0,0.001,ForH+SapH*ScalInv);  # vegetation height [m]
 DEM      = VegH + Alti;                       # DEM with vegetation [m]
 LAI      = 0.2076 * VegH + 0.1;               # Leaf Area Index all vegetat
 LAIsap   = 0.2076 * SapH + 0.1;               # Leaf Area Index saplings
 LAIsapHi = LAI-LAIsap;              # LAI above saplings (gap=0)
#report(day) SapHD = SapH;    #report(day) VegHD = VegH;
#report(mon) SapHM= SapH;     #report(mon) VegHM = VegH;
#report(yr)  SapHY = SapH;     #report(yr)  VegHY = VegH;
 
# GAPRAD.MOD
# =========================================================================# When there is a slope in the research area, better to make a manual slope
# map, because edge effects due to the trees will occur:
# SlopMap = scalar(atan(slope(DEM)));
# SlopMap = if(SlopMap eq 0 then scalar(0.01) else SlopMap);
 AtmPcor = ((288-0.0065*DEM)/288)**5.256; # atm pressure corr [-]
 AspMap  = scalar(aspect(DEM));           # aspect [deg]
 AspMap  = if(AspMap le 0 then scalar(0.01) else AspMap);
 
# Solar geometry
# ----------------------------
# SolDec  :declination sun per day  between +23 & -23 [deg]
# HourAng :hour angle [-] of sun
# SolAlt  :solar altitude [deg], height of sun above horizon
 SolDec  = -23.4*cos(360*(Day+10)/365);
 HourAng = 15*(Hour-12.01);
 SolAlt  = scalar(asin(scalar(sin(Lat)*sin(SolDec)+cos(Lat)*
           cos(SolDec)*cos(HourAng))));
 
# Solar azimuth
# ----------------------------
# SolAzi  :angle solar beams to N-S axes earth [deg]
 SolAzi = scalar(acos((sin(SolDec)*cos(Lat)-cos(SolDec)*
          sin(Lat)*cos(HourAng))/cos(SolAlt)));
 SolAzi = if(Hour le 12 then SolAzi else 360 - SolAzi);
# Additonal extra correction by R.Sluiter, Aug '99
 SolAzi = if(SolAzi gt 89.994 and SolAzi lt 90, 90, SolAzi);
 SolAzi = if(SolAzi gt 269.994 and SolAzi lt 270, 270, SolAzi);
 
# Surface azimuth
# ----------------------------
# cosIncident :cosine of angle of incident; angle solar beams to angle Surface
 cosIncident = sin(SolAlt)*cos(SlopMap)+cos(SolAlt)*sin(SlopMap)
               *cos(SolAzi-AspMap);
 
# Critical angle sun
# ----------------------------
# HoriAng  :tan maximum angle over DEM in direction sun, 0 if neg
# CritSun  :tan of maximum angle in direction solar beams
# Shade    :cell in sun 1, in shade 0
 HoriAng   = horizontan(DEM,directional(SolAzi));
 HoriAng   = if(HoriAng lt 0 then scalar(0) else HoriAng);
 CritSun   = if(SolAlt gt 90 then scalar(0) else scalar(atan(HoriAng)));
 Illum     = if(SolAlt gt CritSun then scalar(1) else scalar(0));
 Shade     = if(SolAlt gt CritSun then scalar(0) else scalar(1));
 
# Radiation outer atmosphere
# ----------------------------
 OpCorr = Trans**((sqrt(1229+(614*sin(SolAlt))**2)-614*sin(SolAlt))
          *AtmPcor);                      # correction for air masses [-]
 Sout   = Sc*(1+0.034*cos(360*Day/365)); # radiation outer atmosphere [W/m2]
 Snor   = Sout*OpCorr;                    # rad on Surface normal to the beam [W/m2]
 
# Radiation on vegetation
# ----------------------------
# Sdir   :direct sunlight on a horizontal Surface [W/m2] if no shade
# Sdiff  :diffuse light [W/m2] for shade and no shade
 Sdir    = if(Snor*cosIncident*Illum<0,0.0,Snor*cosIncident*Illum);
 Sdiff   = if(Sout*(0.271-0.294*OpCorr)*sin(SolAlt)<0, 0.0,
           Sout*(0.271-0.294*OpCorr)*sin(SolAlt));
 Stot    = cover(Sdir+Sdiff,mapmaximum(Sdir+Sdiff));                   # Tot rad with dir rad
 CANex   = if(sin(SolAlt)<=0, 0.00001, exp(-Kext*LAIcan/sin(SolAlt))); # LAI rad ext dense forest
 SmapMax = mapmaximum(Stot);
 Rcan    = SmapMax * CANex;                  # LAI rad ext canopy
 Sveg    = cover(Illum*Stot + Shade*(Rcan+Sdiff),
           SmapMax);                         # Radiation on vegetat
 VegRad  = (Sold + Sveg)/2;                           # Radiation hour intval
 Sold    = Sveg;
#report(yr) VegRpY = VegRpY + VegRad*0.0036;        # year rad [MJ/d]
#report(mon) VegRpM = if(Month eq 1,VegRad*0.0036,VegRpM+VegRad*0.0036);
#report(day) VegRpD = if(Hour eq 1,VegRad*0.0036,VegRpD+VegRad*0.0036);
 
# Sapling and Soil radiation
# ----------------------------
# First: Assign ldd to direction solar beams direction: azimuth
 SolBeam = if(SolAzi > 337.5, Ldd8, if(SolAzi > 292.5, Ldd7,
           if(SolAzi > 247.5, Ldd4, if(SolAzi > 202.5, Ldd1,
           if(SolAzi > 157.5, Ldd2, if(SolAzi > 112.5, Ldd3,
           if(SolAzi > 67.5, Ldd6, if(SolAzi > 22.5, Ldd9, Ldd8))))))));
 BeamID  = catchment(SolBeam,pit(SolBeam));
 
# Second: distance from gap edge with illumination
 SolEdge  = cover(if(Sdir*Edge < 0.0001, 0.0, Sveg*Edge),0.0); #rad at edge
 SolEdgeB = boolean(if(SolEdge eq 0,0,1));
 LDDEdge  = cover(ldddist(SolBeam,SolEdgeB,1),0);      # dist from edge
 DistEdge = if(abs(HourAng)<1,LDDEdge*ForGap,LDDEdge); # -neg effects 12h
 SolRarea = if(DistEdge > 0, scalar(1), scalar(0));    # area with Sdir influ
 SolRinv  = if(SolRarea eq 0, scalar(1), scalar(0));
 MaxRadE  = areamaximum(SolEdge,BeamID) * SolRarea;
 
# Third: radiation extinction coefficients
 LAIex    = if(sin(SolAlt)<=0, 0.00001, exp(-Kext*LAI/sin(SolAlt)));
 SAPex    = if(sin(SolAlt)<=0, 0.00001, exp(-Kext*LAIsapHi/sin(SolAlt)));
 EdgeEx   = if(sin(SolAlt) <= 0, 0.00001, if(DistEdge eq 0, 0.00001,
            exp(-Cext*DistEdge/sin(SolAlt))));
 
# Fourth: radiation on saplings
 Rsap     = SolRinv * Sveg * SAPex;       # radiation on saplings not in edge
 EdgeSap  = SolRarea * if(SolAlt<0, 0.0, (MaxRadE-Rsap)*EdgeEx+Rsap);
 Ssap    = Rsap + EdgeSap;
 SapRad  = (SsapOld + Ssap)/2;
 SsapOld = Ssap;
#report(yr) SapRpY = SapRpY + SapRad*0.0036; # year rad [MJ/d]
#report(mon) SapRpM = if(Month eq 1,SapRad*0.0036,SapRpM+SapRad*0.0036);
#report(day) SapRpD = if(Hour eq 1,SapRad*0.0036,SapRpD+SapRad*0.0036);
 
# Fith: radiation on soil
 Rsoil    = SolRinv * Sveg * LAIex;       # radiation on soil not in edge
 EdgeSoil = SolRarea * if(SolAlt<0, 0.0, (MaxRadE-Rsoil)*EdgeEx+Rsoil);
 Ssoil    = Rsoil + EdgeSoil;
 SoilRad  = (Ssoil + SsoilO)/2;                       # hour interval
 SsoilO   = Ssoil;
#report(yr) SoilRpY = SoilRpY + SoilRad*0.0036; # year rad [MJ/d]
#report(mon) SoilRpM = if(Month eq 1,SoilRad*0.0036,SoilRpM+SoilRad*0.0036);
#report(day) SoilRpD = if(Hour eq 1,SoilRad*0.0036,SoilRpD+SoilRad*0.0036);
 
# POTEVAP.MOD
# =========================================================================# INPUT METEO DATA
# Humidity and Temperature averaged over the hour
 RHn   = timeinputscalar(Meteo,1)/100;      # Rel Hum [% to fraction]
 RHn   = if(RHn ge 1, 0.999, RHn);
 RH    = (RHo + RHn)/2;
 RHo   = RHn;
 TempN = timeinputscalar(Meteo,2);            # Temperature [oC]
 Temp  = (TempN + TempO)/2;
 TempO = TempN;
# Wind and cloud are measured as average of last hour
 Wind  = timeinputscalar(Meteo,3);            # Windspeed [m/s]
 Wind  = if(Wind eq 0, 0.0001, Wind);
 Cloud = timeinputscalar(Meteo,5);            # cloud factor [-]
# Air pressure [mbar to Pa] not avg over hour because almost no fluctuation
 Press = timeinputscalar(Meteo,4);            # in mbar
 Rain  = timeinputscalar(Meteo,6);            # Rainfall sum [mm] of last hour
#report(yr) RainY = (Rain + RainY)*(GapScal+ScalInv);
# Meteo for understory
 RHS   = if(Gap eq 1, 0.0965*RH+0.89892, RH);                   # RH und [-]
 TempS = if(Gap eq 1, 0.5243*Temp+11.206, Temp);                # T und [deg C]
 WindS = if(Gap eq 1, 0.1889*Wind**2+0.0523*Wind+0.0013, Wind); # Wind und [m/s]
 
# Litter Mass
# -----------------------------------------
 LFedg  = (LFa/(1+exp(LFb+LFc*DFE)))*100/8760;    # edge litterfall
 LFgap  = (0.7616 * SapH + 0.1)*100/8760;         # gap litterfall [g/m2/h]
 LF     = if(LFedg+LFgap ge 560/8760, 560/8760, LFedg+LFgap);
 DEC    = LM * Kfactor / 8760;                      # decomposition [g/m2/h]
 LM     = LM - DEC + LF;                            # litter mass [g/m2]
 WHCLM  = 0.065*(1-exp(-0.0033*LM));             # Stor cap. LM [mm] R2=0.47
 LO     = exp(-0.0143*LM);                       # Litter open [-]:R2=0.6814
#report(yr) LFedgY  = LFedg + LFedgY;
#report(yr) LFgapY  = LFgap + LFgapY;
#report(yr) LFY     = LF + LFY;
#report(yr) DECY    = DEC + DECY;
#report(yr) LMY     = LM;
 
# Vapour pressure ao
# -----------------------------------------
# Forest & Gap
 Lv    = 2.501E6 - 2361 * Temp;                # lat heat vap of water [MJ/kg]
 Es    = 6.107*exp((17.27*Temp)/(Temp+237.3)); # sat vap pressure [mbar]
 Ea    = RH * Es;                              # actual vap pressure [mbar]
 Delta = 409.8 * Es / (Temp+237.3)**2;        # slope sat vap pres [mbar/K]
 Gamma = (Cp*Press)/(Epsilon*Lv);                  # psychrometric cons. [mbar/K]
# Understory, soil
 LvS    = 2.501E6 - 2361 * TempS;
 EsS    = 6.107*exp((17.27*TempS)/(TempS+237.3));
 EaS    = RHS * EsS;
 DeltaS = 409.8 * EsS / (TempS+237.3)**2;
 GammaS = (Cp*Press)/(Epsilon*LvS);
 
# Wind correction
# -----------------------------------------
 Hmet  = if(Year eq 1996, 0.0019 * Day - 0.3653,
         if(Year eq 1997, 0.0048 * Day + 0.3301,
         if(Year eq 1998, 0.0037 * Day + 2.0821,
         0.0037 * Day + 4.9791)));            # SapH veg under meteotower
 Hmet  = if(Hmet le 0.2, scalar(0.2), Hmet);
 Dveg  = 0.7 * VegH;                          # zero-plane displacement
 Dsap  = 0.7 * SapH;
 Z0veg = 0.1 * VegH;                          # roughness length vegetation
 Z0sap = 0.1 * SapH;
 Z0met = 0.1 * Hmet;                          # roughness length at meteo
 WcorF = (ln((VegH-Dveg)/Z0veg)*ln(60/Z0met))/
         (ln(60/Z0veg)*ln(Hwind/Z0met));      # wind cor factor forest
 WcorG = (ln(VegH/Z0veg)*ln(60/Z0met))/
         (ln(60/Z0veg)*ln(Hwind/Z0met));      # wind cor factor gap
 Wcor  = if(VegH > 20.0, WcorF, WcorG);
 
# CO & storage capacities
# -----------------------------------------
 WHCCan = CanStor * LAI;                      # Stor capacity canopy [mm]
 CO = 1 - 0.0259 * VegH;                      # Cell openness [-]
 
# Surface roughness parameters
# -----------------------------------------
 RaVeg = (ln((VegH-Dveg)/Z0veg))**2/(Ka**2*Wind*Wcor);
 RaVeg = if(RaVeg > RaMaxVeg, RaMaxVeg, RaVeg);            # Aerodynamic resistance [s/m]
 RaSap = if(Gap eq 1, (ln((SapH)/Z0sap))**2/(Ka**2*WindS),
         (ln((SapH)/Z0sap))**2/(Ka**2*Wind*Wcor));
 RaSap = if(RaSap > RaMaxSap, RaMaxSap, RaSap);            # Aerodynamic resistance [s/m]
 Ra    = if(Gap eq 0, RaSap, RaVeg);
 
# Radiation
# -----------------------------------------
# RsBal: shortwave rad [W/m2], RlBal: longwave rad [W/m2]
# Rnet: net rad (Rs + Rl) [MJ/m2]
# Veg: on vegetation,
 Albedo  = min(0.1509-0.00136*SolAlt+0.0000123*SolAlt**2,0.15); # Shuttleworth '94
 RsVeg   = (1-Albedo)*VegRad*Cloud;               # short wave on veg [W/m2]
 RlVeg   = if(SolAlt < 0, 0, StefBo*((Temp+273.15)**4)*(0.56-0.079*
           sqrt(Ea))*(0.1+0.9*Cloud)); # long wave on veg [W/m2]
 RnetVeg = if(RsVeg - RlVeg < 0, scalar(0), RsVeg - RlVeg);# net rad [W/m2]
#report(yr) RnVegY  = RnetVeg*0.0036 + RnVegY;       # [MJ/m2]
#report(mon) RnVegM = if(Month eq 1,RnetVeg*0.0036,RnVegM+RnetVeg*0.0036);
#report(day) RnVegD = if(Hour eq 1,RnetVeg*0.0036,RnVegD+RnetVeg*0.0036);
#report(yr) RsVegY = 0.0036*RsVeg + RsVegY;
#report(yr) RlVegY = 0.0036*RlVeg*test + RlVegY;
 
# Sap: on saplings
 RnetSap = SolRarea*RnetVeg*EdgeEx + SolRinv*RnetVeg*SAPex;     # W/m2
#report(yr)  RnSapY = RnetSap*0.0036 + RnSapY;          # MJ/h
#report(mon) RnSapM = if(Month eq 1,RnetSap*0.0036,RnSapM+RnetSap*0.0036);
#report(day) RnSapD = if(Hour eq 1,RnetSap*0.0036,RnSapD+RnetSap*0.0036);
 
# Soil: on soil
 RnetSoil = SolRarea*RnetVeg*EdgeEx + SolRinv*RnetVeg*LAIex;     # W/m2
#report(yr) RnSoilY = RnetSoil*0.0036 + RnSoilY;
#report(mon) RnSoilM = if(Month eq 1,RnetSoil*0.0036,RnSoilM+RnetSoil*0.0036);
#report(day) RnSoilD = if(Hour eq 1,RnetSoil*0.0036,RnSoilD+RnetSoil*0.0036;
 
# No evapotranspiration during rainfall
# -----------------------------------------
 NoRain   = 1-(3.6542*Rain/60);                  # Time no rain Ei&Et [min/h]
 NoRain   = if(NoRain < 0, scalar(0), NoRain);
 
# Evaporation from a wet vegetated Surface
# -----------------------------------------
 InbetVeg = 3600*NoRain*(Delta*RnetVeg+Roa*Cp*(Es-Ea)/Ra); # inbetween
 InbetSap = 3600*NoRain*(DeltaS*RnetSap+Roa*Cp*(EsS-EaS)/RaSap);
 EW       = InbetVeg/(Lv*(Delta+Gamma));           # Evap wet Surface [mm/h]
 EW       = if(EW < 0, scalar(0), EW);
#report(yr) EWy = EW + EWy;
 
# Stemflow & Stem evaporation
# -----------------------------------------
 PonStem     = if(Rain + PonStemOver eq 0, 0, if(Rain + PonStemOver > StemStor,
               StemStor, Rain + PonStemOver));           #stem storage [mm]
 StemFlow    = max((1-CO)*Rain*StemCoef-StemStor,0);    # Stemflow [mm/h]
 Estem       = NoRain * min(EW, PonStem);
 PonStemOver = max(PonStem - Estem, 0);
#report(yr) STY = StemFlow + STY;
#report(mon) STM = if(Month eq 1,StemFlow,STM+StemFlow);
#report(day) STD = if(Hour eq 1,StemFlow,STD+StemFlow);
#report(yr) ESTY = Estem + ESTY;
#report(mon) ESTM = if(Month eq 1,Estem,ESTM+Estem);
#report(day) ESTD = if(Hour eq 1,Estem,ESTD+Estem);
 
#Evaporation of intercepted water on canopy
# -----------------------------------------
 PonCan   = (1-CO) * Rain - StemFlow + OverCan;      # Total water [mm/h]
 IncepCan = min(PonCan,WHCCan);                      # Water storred [mm/h]
 DripCan  = if(PonCan > WHCCan, PonCan - WHCCan, 0); # Canopy drip [mm/h]
 Through  = CO * Rain + DripCan;                     # Throughflow [mm/h]
 EI       = min(IncepCan,EW);                        # Intercep evap [mm/h]
 OverCan  = if(EI < IncepCan, IncepCan - EI, 0);     # left over water [mm/h]
#report(yr) THY = Through + THY;
#report(mon) THM = if(Month eq 1,Through,THM+Through);
#report(day) THD = if(Hour eq 1,Through,THD+Through);
#report(yr) EIY = EI + Estem + EIY;
#report(mon) EIM = if(Month eq 1,EI,EI+EIM);
#report(day) EID = if(Hour eq 1,EI,EI+EID);
 
# Potential Transpiration of dry vegetation
# -----------------------------------------
 RcVeg = exp(0.6251-0.0019*VegRad+0.00155*(Es-Ea)+
         0.02247*(Temp)-0.5699*Wcor*Wind)*100;     # adjusted Schellekens 2000
 ETp   = InbetVeg/(Lv*(Delta+Gamma*(1+RcVeg/Ra))); # Pot Trans dry veg [mm/h]
 ETp   = if(ETp < 0, scalar(0), ETp);
#report(yr) ETpY = ETp + ETpY;
#report(mon) ETpM = if(Month eq 1,ETp,ETp+ETpM);
#report(day) ETpD = if(Hour eq 1,ETp,ETp+ETpD);
 
# Potential Transpiration of dry seedlings
# -----------------------------------------
# WcorSap = (ln(SapH/Z0sap)*ln(60/Z0met))/
#            (ln(60/Z0sap)*ln(Hwind/Z0met));
# WcorSap = if(VegH > 20.0, WcorSap, WcorG);
# RaSap   = if (Gap eq 0 then (ln((SapH-Dsap)/Z0sap))**2/(Ka**2*Wind*WcorSap)
#           else (ln((SapH-Dsap)/Z0sap))**2/(Ka**2*WindS*WcorSap));
# RaSap   = if(RaSap > RaMax , RaMax, RaSap);
 RcSap   = exp(0.6251-0.0019*SapRad+0.00155*(EsS-EaS)+
           0.02247*(TempS)-0.5699*WindS)*100;  # adjusted Schellekens 2000
 ETsap   = NoRain * if(Gap eq 0, ETp,
           InbetSap/(LvS*(DeltaS+GammaS*(1+RcSap/RaSap))));
 ETsap   = if(ETsap < 0, scalar(0), ETsap);
#report(yr) ETsapY = ETsap + ETsapY;
#report(mon) ETsapM = if(Month eq 1,ETsap,ETsap+ETsapM);
#report(day) ETsapD = if(Hour eq 1,ETsap,ETsap+ETsapD);
 
# Potential Soil evaporation
# -----------------------------------------
 ESp = 3600*NoRain * if(GapScal eq 0, (Delta*RnetSoil)/(Lv*(Delta+Gamma)),
       (DeltaS*RnetSoil)/(LvS*(DeltaS+GammaS)));
 ESp = if(ESp < 0, scalar(0), ESp);
#report(yr) ESpY = ESp + ESpY;
#report(mon) ESpM = if(Month eq 1,ESp,ESp+ESpM);
#report(day) ESpD = if(Hour eq 1,ESp,ESp+ESpD);
 
# Drainage canopy & Litter Mass Evaporation
# -----------------------------------------
 PonLM   = (1-LO) * Through + OverLM;            # Water on littermass [mm/h]
 IncepLM = if(PonLM > WHCLM, WHCLM, PonLM);      # Storrage on LM [mm/h]
 DripLM  = if(PonLM > WHCLM, PonLM - WHCLM, 0);  # Drip from LM [mm/h]
 EL      = if(IncepLM > ESp, ESp, IncepLM);      # Litter Mass evap [mm/h]
 OverLM  = if(EL < IncepLM, IncepLM - EL, 0);    # left over water [mm/h]
 Psoil   = LO * Through + StemFlow + DripLM;     # tot drainage to soil [mm/h]
#report(yr) PSY = Psoil + PSY;
#report(mon) PSM = if(Month eq 1,Psoil,Psoil+PSM);
#report(day) PSD = if(Hour eq 1,Psoil,Psoil+PSD);
#report(yr) ELY = EL + ELY;
#report(mon) ELM = if(Month eq 1,EL,EL+ELM);
#report(day) ELD = if(Hour eq 1,EL,EL+ELD);
 
# Potential Evapotranspiration
# -----------------------------------------
 EP = EI + EL + ETp + ESp;          # Tot Pot evapotranspi [mm/h]
#report(yr) EPY = EP + EPY;
#report(mon) EPM = if(Month eq 1,EP,EP+EPM);
#report(day) EPD = if(Hour eq 1,EP,EP+EPD);
 
# WATBAL.MOD
# =========================================================================# Soil matrix suction 'Psi' is directed positive downwards, all in cm
# Top boundary condition equals atmospheric suction
# Bottom boundary condition equals condition layer 3 previous time step
 
# Calculation of psi (suction) positive (cm) Kutílek & Nielsen, 1994, eq5.50
 Psi5    = Psi4;
 Psi1    = min((((1/min(0.999,ThetaE1))**(1/m1)-1)**(1/n1))/a1, Hmax);
 Psi2    = min((((1/min(0.999,ThetaE2))**(1/m2)-1)**(1/n2))/a2, Hmax);
 Psi3    = min((((1/min(0.999,ThetaE3))**(1/m3)-1)**(1/n3))/a3, Hmax);
 Psi4    = min((((1/min(0.999,ThetaE4))**(1/m4)-1)**(1/n4))/a4, Hmax);
 
# Calculation of Kpsi (cm/hour): Mualem (1976) Kutílek & Nielsen eq.5.52
 Kpsi5    = Kpsi4;
 Kpsi1    = Ks1*(1-(a1*Psi1)**(n1-1)*(1+(a1*Psi1)**n1)**-m1)**2
            /(1+(a1*Psi1)**n1)**(m1/2);
 Kpsi2    = Ks2*(1-(a2*Psi2)**(n2-1)*(1+(a2*Psi2)**n2)**-m2)**2
            /(1+(a2*Psi2)**n2)**(m2/2);
 Kpsi3    = Ks3*(1-(a3*Psi3)**(n3-1)*(1+(a3*Psi3)**n3)**-m3)**2
            /(1+(a3*Psi3)**n3)**(m3/2);
 Kpsi4    = Ks4*(1-(a4*Psi4)**(n4-1)*(1+(a4*Psi4)**n4)**-m4)**2
            /(1+(a4*Psi4)**n4)**(m4/2);
 
# Infiltration
# ------------------------
# Potential & actual infiltration rate, cannot be less than zero (cm/hr)
 Ipot = max(-sqrt(Ks1*Kpsi1)*((-Psi1)/(D1/2)-1), 0);
 Iact = max(min(Ipot, Psoil/10+Pond, (ThetaS1-Theta1)*D1), 0);
 Pond = max(Psoil/10 + Pond - Iact,0);
 RunOff = if(SlopMap > 0.5, Pond, 0);
 Pond = Pond - RunOff;
#report(yr) RunOY = RunOff*10 + RunOY;
#report(mon) RunOM = if(Month eq 1,RunOff*10,RunOff*10+RunOM);
#report(day) RunOD = if(Hour eq 1,RunOff*10,RunOff*10+RunOD);
 
# Evaporation of ponded water
# ----------------------
 Epond = min(ESp,Pond);
#report(yr) EpondY = Epond*10 + EpondY;
#report(mon) EpondM = if(Month eq 1,Epond*10,Epond*10+EpondM);
#report(day) EpondD = if(Hour eq 1,Epond*10,Epond*10+EpondD);
 Pond  = Pond - Epond;
#report(yr) PondY = Pond*10 + PondY;
#report(mon) PondM = if(Month eq 1,Pond*10,Pond*10+PondM);
#report(day) PondD = if(Hour eq 1,Pond*10,Pond*10+PondD);
 
# Percolation & evapotranspiration
# ------------------------
# Calculation of Perculation between two layers (mm)
 Perc12 = min(-sqrt(Kpsi1*Kpsi2)*((Psi1-Psi2)/((D1+D2)/2)-1),
             (ThetaS2-Theta2)*D2, (Theta1-ThetaR)*D1);
 Perc12 = if(Pond > 0 and Perc12 < 0, 0, Perc12);
 Perc12 = if(Iact-Perc12 > (ThetaS1-Theta1)*D1, 0, Perc12);
 Perc23 = min(-sqrt(Kpsi2*Kpsi3)*((Psi2-Psi3)/((D2+D3)/2)-1),
             (ThetaS3-Theta3)*D3, (Theta2-ThetaR)*D2);
 Perc34 = min(-sqrt(Kpsi3*Kpsi4)*((Psi3-Psi4)/((D3+D4)/2)-1),
             (ThetaS4-Theta4)*D4, (Theta3-ThetaR)*D3);
 Perc45 = min(-sqrt(Kpsi4*Kpsi5)*((Psi4-Psi5)/(D4)-1),
             (ThetaS4-Theta5)*D4, (Theta5-ThetaR)*D4);
#              (ThetaSR4-ThetaE5*ThetaSR4)*D4);
# Loss below third layer [mm]
#report(yr) LoY = Perc45*10 + LoY;
#report(mon) LoM = if(Month eq 1,Perc45*10,Perc45*10+LoM);
#report(day) LoD = if(Hour eq 1,Perc45*10,Perc45*10+LoD);
 
# Soil evaporation
 PsiEvap  = -((R*TempS+273.15)/(M*G))*ln(RHS);   # Press head atm is neg upwards
 KpsiEvap = Ks1*(1-(a1*PsiEvap)**(n1-1)*(1+(a1*PsiEvap)**n1)**-m1)**2
            /(1+(a1*PsiEvap)**n1)**(m1/2);
 Esoil    = max(min(-sqrt(KpsiEvap*Kpsi1)*((Psi1-PsiEvap)/(D1/2)-1),
                   (Theta1-ThetaR)*D1+Iact-Perc12), 0);
 Esoil    = if(Pond > 0, 0, Esoil);
 ESa      = min(ESp/10, Esoil);                          # Act soil evap [cm]
#report(yr) ESaY = ESa*10 + ESaY;
#report(mon) ESaM = if(Month eq 1,Esoil,Esoil+ESaM);
#report(day) ESaD = if(Hour eq 1,Esoil,Esoil+ESaD);
 
# Transpiration
# Layer 1
 Red11    = 0.101 * Psi1 - 0.0101;                     # Trans red H<10
 Red21    = 1 / (1 + (Psi1/H50)**RedT);                # Trans red 10<H<Hmax
 Reduc1   = if(Psi1 < 0.1, 0.0, if(Psi1 < 1.0, Red11,
            if(Psi1 < Hmax, Red21, 0.0)));
 ETa1     = R1 * Reduc1 * ETp/10;                # Actual transpi [cm]
 ETa1     = max(min(ETa1, (Theta1-ThetaR)*D1+Iact-Perc12-ESa), 0);
#report(yr) ETa1Y = ETa1*10 + ETa1Y;
#report(mon) ETa1M = if(Month eq 1,ETa1,ETa1+ETa1M);
#report(day) ETa1D = if(Hour eq 1,ETa1,ETa1+ETa1D);
 
# Layer 2
 Red12    = 0.101 * Psi2 - 0.0101;                     # Trans red H<10
 Red22    = 1 / (1 + (Psi2/H50)**RedT);                # Trans red 10<H<Hmax
 Reduc2   = if(Psi2 < 0.1, 0.0, if(Psi2 < 1.0, Red12,
            if(Psi2 < Hmax, Red22, 0.0)));
 ETa2     = R2 * Reduc2 * ETp/10;                # Actual transpi [cm]
#report(yr) ETa2Y = ETa2*10 + ETa2Y;
#report(mon) ETa2M = if(Month eq 1,ETa2,ETa2+ETa2M);
#report(day) ETa2D = if(Hour eq 1,ETa2,ETa2+ETa2D);
 
# Layer 3
 Red13    = 0.101 * Psi3 - 0.0101;                     # Trans red H<10
 Red23    = 1 / (1 + (Psi3/H50)**RedT);                # Trans red 10<H<Hmax
 Reduc3   = if(Psi3 < 0.1, 0.0, if(Psi3 < 1.0, Red13,
             if(Psi3 < Hmax, Red23, 0.0)));
 ETa3     = R3 * Reduc3 * ETp/10;                # Actual transpi [cm]
#report(yr) ETa3Y = ETa3*10 + ETa3Y;
#report(mon) ETa3M = if(Month eq 1,ETa3,ETa3+ETa3M);
#report(day) ETa3D = if(Hour eq 1,ETa3,ETa3+ETa3D);
 
# Layer 4
 Red14    = 0.101 * Psi4 - 0.0101;                     # Trans red H<10
 Red24    = 1 / (1 + (Psi4/H50)**RedT);                # Trans red 10<H<Hmax
 Reduc4   = if(Psi4 < 0.1, 0.0, if(Psi4 < 1.0, Red14,
             if(Psi4 < Hmax, Red24, 0.0)));
 ETa4     = R4 * Reduc4 * ETp/10;                # Actual transpi [cm]
#report(yr) ETa4Y = ETa4*10 + ETa4Y;
#report(mon) ETa4M = if(Month eq 1,ETa4,ETa4+ETa4M);
#report(day) ETa4D = if(Hour eq 1,ETa4,ETa4+ETa4D);
 
# Total actual transpiration [mm]
 ETa = (ETa1 + ETa2 + ETa3 + ETa4)*10;
#report(yr) ETaY = ETa + ETaY;
#report(mon) ETaM = if(Month eq 1,ETa,ETa+ETaM);
#report(day) ETaD = if(Hour eq 1,ETa,ETa+ETaD);
# Act. trans. in gap for growth MMF model and pot. SAP trans. in forest
Grow = GrowI + ScalInv*ETaY + GapScal*ETsapY;
 
# New theta and counting of pF occurance
# ------------------------
# New moisture content in layers
 Theta5 = Theta4;
 Theta1 = max(min(ThetaS1, Theta1+(Iact-Perc12-ESa-ETa1)/D1), ThetaR);
 Theta2 = max(min(ThetaS2, Theta2+(Perc12-Perc23-ETa2)/D2), ThetaR);
 Theta3 = max(min(ThetaS3, Theta3+(Perc23-Perc34-ETa3)/D3), ThetaR);
 Theta4 = max(min(ThetaS4, Theta4+(Perc34-Perc45-ETa4)/D4), ThetaR);
#report(day) SMD = (SMI-(Theta1*D1+Theta2*D2+Theta3*D3+Theta4*D4))*10;
#report(mon) SMM = (SMI-(Theta1*D1+Theta2*D2+Theta3*D3+Theta4*D4))*10;
#report(yr)  SMY = (SMI-(Theta1*D1+Theta2*D2+Theta3*D3+Theta4*D4))*10;
# initial soil moisture for next year
#report(yr) Theta1Y = Theta1;  report(yr) Theta2Y = Theta2;
#report(yr) Theta3Y = Theta3;  report(yr) Theta4Y = Theta4;
#report(yr) Theta5Y = Theta5;
 
# Theta eff, psi(theta), k(psi) relations
# ------------------------
# Effective moisture content in layers
 ThetaE1 = max((Theta1-ThetaR)/(ThetaSR1), TeMin);
 ThetaE2 = max((Theta2-ThetaR)/(ThetaSR2), TeMin);
 ThetaE3 = max((Theta3-ThetaR)/(ThetaSR3), TeMin);
 ThetaE4 = max((Theta4-ThetaR)/(ThetaSR4), TeMin);
 
# Time steps with drought or water excess
 PF02L1 = if(Psi1 le 1, PF02L1 +1, PF02L1);
 PF01L1 = if(Psi1 gt 1 and Psi1 le 10, PF01L1 +1, PF01L1);
 PF1L1  = if(Psi1 gt 10 and Psi1 le 100, PF1L1 +1, PF1L1);
 PF2L1  = if(Psi1 gt 100 and Psi1 le 1000, PF2L1 +1, PF2L1);
 PF3L1  = if(Psi1 gt 1000 and Psi1 le 10000, PF3L1 +1, PF3L1);
 PF4L1  = if(Psi1 gt 10000, PF01L1 +1, PF4L1);
 
 PF02L2 = if(Psi2 le 1, PF02L2 +1, PF02L2);
 PF01L2 = if(Psi2 gt 1 and Psi2 le 10, PF01L2 +1, PF01L2);
 PF1L2  = if(Psi2 gt 10 and Psi2 le 100, PF1L2 +1, PF1L2);
 PF2L2  = if(Psi2 gt 100 and Psi2 le 1000, PF2L2 +1, PF2L2);
 PF3L2  = if(Psi2 gt 1000 and Psi2 le 10000, PF3L2 +1, PF3L2);
 PF4L2  = if(Psi2 gt 10000, PF01L2 +1, PF4L2);
 
#report(yr) PF02 = (PF02L1 + PF02L2)/2;
#report(yr) PF01 = (PF01L1 + PF01L2)/2;
#report(yr) PF1  = (PF1L1  + PF1L2)/2;
#report(yr) PF2  = (PF2L1  + PF2L2)/2;
#report(yr) PF3  = (PF3L1  + PF3L2)/2;
#report(yr) PF4  = (PF4L1  + PF4L2)/2;
 
#report(yr) PF02L3 = if(Psi3 le 1, PF02L3 +1, PF02L3);
#report(yr) PF01L3 = if(Psi3 gt 1 and Psi3 le 10, PF01L3 +1, PF01L3);
#report(yr) PF1L3  = if(Psi3 gt 10 and Psi3 le 100, PF1L3 +1, PF1L3);
#report(yr) PF2L3  = if(Psi3 gt 100 and Psi3 le 1000, PF2L3 +1, PF2L3);
#report(yr) PF3L3  = if(Psi3 gt 1000 and Psi3 le 10000, PF3L3 +1, PF3L3);
#report(yr) PF4L3  = if(Psi3 gt 10000, PF01L3 +1, PF4L3);
 
#report(yr) PF02L4 = if(Psi4 le 1, PF02L4 +1, PF02L4);
#report(yr) PF01L4 = if(Psi4 gt 1 and Psi4 le 10, PF01L4 +1, PF01L4);
#report(yr) PF1L4  = if(Psi4 gt 10 and Psi4 le 100, PF1L4 +1, PF1L4);
#report(yr) PF2L4  = if(Psi4 gt 100 and Psi4 le 1000, PF2L4 +1, PF2L4);
#report(yr) PF3L4  = if(Psi4 gt 1000 and Psi4 le 10000, PF3L4 +1, PF3L4);
#report(yr) PF4L4  = if(Psi4 gt 10000, PF01L4 +1, PF4L4);
 
# Water balance error
# ----------------------
 WBerr1    = Iact-Perc12-ESa-ETa1-(Theta1-Theta1old)*D1;
 WBerr2    = Perc12-Perc23-ETa2-(Theta2-Theta2old)*D2;
 WBerr3    = Perc23-Perc34-ETa3-(Theta3-Theta3old)*D3;
 WBerr4    = Perc34-Perc45-ETa4-(Theta4-Theta4old)*D4;
 WBErrtot  = WBerr1 + WBerr2 + WBerr3 + WBerr4;
 WBerr     = Iact-ESa-ETa/10-Perc45-(Theta1-Theta1old)*D1-
             (Theta2-Theta2old)*D2-(Theta3-Theta3old)*D3-
             (Theta4-Theta4old)*D4;
#report(day) WBerrD = WBerr;
#report(mon)  WBerrM = WBerr;
#report(yr)   WBerrY = WBerr;
 Theta1old = Theta1; Theta2old = Theta2;
 Theta3old = Theta3; Theta4old = Theta4;
 
# Timeoutput of sample locations
# ====================================#report saph.tim     = timeoutput(Loc,SapH);
#report vegh.tim     = timeoutput(Loc,VegH);
 
# POTRAD
#report soldec.tim  = timeoutput(Loc,SolDec);
#report hourang.tim = timeoutput(Loc,HourAng);
#report solalt.tim  = timeoutput(Loc,SolAlt);
#report solazi.tim  = timeoutput(Loc,SolAzi);
#report cosi.tim    = timeoutput(Loc,cosIncident);
 
#report crit.tim    = timeoutput(Loc,CritSun);
#report illum.tim   = timeoutput(Loc,Illum);
#report shade.tim   = timeoutput(Loc,Shade);
#report snor.tim    = timeoutput(Loc,Snor);
 
#report vegrad.tim  = timeoutput(Loc,VegRad);
#report sdir.tim    = timeoutput(Loc,Sdir);
#report sdif.tim    = timeoutput(Loc,Sdiff);
 
#report sveg.tim    = timeoutput(Loc,Sveg);
#report ssap.tim    = timeoutput(Loc,Ssap);
#report soilrad.tim = timeoutput(Loc,Ssoil);
 
# POTEVAP
#report ra.tim = timeoutput(Loc,Ra);
#report rc.tim = timeoutput(Loc,Rc);
 
#report albedo.tim = timeoutput(Loc,Albedo);
#report rnveg.tim  = timeoutput(Loc,RnetVeg);
#report rlveg.tim  = timeoutput(Loc,RlVeg);
#report rsveg.tim  = timeoutput(Loc,RsVeg);
 
#report rnsap.tim  = timeoutput(Loc,RnetSap);
#report rssap.tim  = timeoutput(Loc,RsSap);
#report grow.tim  = timeoutput(Loc,Grow);
#report rnsoil.tim = timeoutput(Loc,RnetSoil);
#report rssoil.tim  = timeoutput(Loc,RsSoil);
 
#report LFedg.tim = timeoutput(Loc,LFedg);
#report LFgap.tim = timeoutput(Loc,LFgap);
#report LF.tim = timeoutput(Loc,LF);
#report DEC.tim = timeoutput(Loc,DEC);
#report LM.tim = timeoutput(Loc,LM);
#report WHC.tim = timeoutput(Loc,WHCLM);
#report LO.tim = timeoutput(Loc,LO);
 
#report ew.tim       = timeoutput(Loc,EW);
#report ei.tim       = timeoutput(Loc,EI);
#report el.tim       = timeoutput(Loc,EL);
#report etp.tim      = timeoutput(Loc,ETp);
#report esp.tim      = timeoutput(Loc,ESp);
#report ep.tim       = timeoutput(Loc,EP);
#report esap.tim     = timeoutput(Loc,ETsap);
 
#report stem.tim     = timeoutput(Loc, StemFlow);
#report through.tim  = timeoutput(Loc, Through);
#report psoil.tim    = timeoutput(Loc,Psoil);
 
# WATBAL
#report ipot.tim  = timeoutput(Loc,Ipot*10);
#report iact.tim  = timeoutput(Loc,Iact*10);
#report epond.tim = timeoutput(Loc,Epond*10);
#report pond.tim  = timeoutput(Loc,Pond*10);
 
report t1.tim = timeoutput(Loc,Theta1*100);
#report t2.tim = timeoutput(Loc,Theta2*100);
#report t3.tim = timeoutput(Loc,Theta3*100);
#report t4.tim = timeoutput(Loc,Theta4*100);
 
#report te1.tim = timeoutput(Loc,ThetaE1);
#report te2.tim = timeoutput(Loc,ThetaE2);
#report te3.tim = timeoutput(Loc,ThetaE3);
#report te4.tim = timeoutput(Loc,ThetaE4);
 
#report psi1.tim = timeoutput(Loc,Psi1);
#report psi2.tim = timeoutput(Loc,Psi2);
#report psi3.tim = timeoutput(Loc,Psi3);
#report psi4.tim = timeoutput(Loc,Psi4);
 
#report kpsi1.tim = timeoutput(Loc,Kpsi1);
#report kpsi2.tim = timeoutput(Loc,Kpsi2);
#report kpsi3.tim = timeoutput(Loc,Kpsi3);
#report kpsi4.tim = timeoutput(Loc,Kpsi4);
 
#report perc12.tim = timeoutput(Loc,Perc12*10);
#report perc23.tim = timeoutput(Loc,Perc23*10);
#report perc34.tim = timeoutput(Loc,Perc34*10);
#report perc45.tim = timeoutput(Loc,Perc45*10);
 
#report esoil.tim  = timeoutput(Loc,Esoil*10);
#report esa.tim    = timeoutput(Loc,ESa*10);
#report eta1.tim   = timeoutput(Loc,ETa1*10);
#report eta2.tim   = timeoutput(Loc,ETa2*10);
#report eta3.tim   = timeoutput(Loc,ETa3*10);
#report eta4.tim   = timeoutput(Loc,ETa4*10);
#report eta.tim    = timeoutput(Loc,ETa*10);
 
#report wberr1.tim   = timeoutput(Loc,WBerr1);
#report wberr2.tim   = timeoutput(Loc,WBerr2);
#report wberr3.tim   = timeoutput(Loc,WBerr3);
#report wberr4.tim   = timeoutput(Loc,WBerr4);
#report wberrtot.tim = timeoutput(Loc,WBErrtot);
#report wberr.tim    = timeoutput(Loc,WBerr);
 
