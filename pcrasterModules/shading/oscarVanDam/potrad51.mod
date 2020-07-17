# POTRAD51.MOD
# Potential Radiation Equator model
# (c) O. van Dam, UU, Tropenbos-Guyana
# Version 5.1, June 2000
# HORIZONTAN per 5 deg
# NOTE: Copyright: This program is free to use provided 
#       you refer to the manualfor citation.
#       Do not distribute without prior approval of the author.
#       Manual and additional info: O.vanDam@geog.uu.nl

# -----------------------------------------------------
#                 Model for calculation
#             incoming potential light energy
# -----------------------------------------------------
# INPUT:  - DEM study area
#         - location attributes
# OUTPUT: - hour reports on the sample locations
#         - potential global radiation per hour [W/m2] (optional)
#         - potential global radiation per day, month & year [MJ/m2]

binding
# INPUT BY USER
 DEM  = hilldem.map;         # DEM research area
# Loc  = hillloc.map;         # nominal map with sample locations
 Lat  = 52.0;                # latitude study area
 NoTSteps = 8760;             # No of times steps
 HourStep =  1.0;            # New Hour = Old Hour + HourStep
 DayStep  =  1.0;            # New Day = Old Day + DayStep   
 StartDay =  1.0;            # Julian dayno startday

# OUTPUT BY USER
 Slope  = SlopMap;         # Slope of DEM
 Aspect = AspMap;          # Aspect of DEM

# constants
 pi       = 3.1415;          # pi
 Sc       = 1367.0;          # Solar constant (Gates, 1980) [W/m2]
 Trans    = 0.6;             # Transmissivity tau (Gates, 1980)

areamap
 DEM;

timer
 1 NoTSteps 1;                   # see manual for nr of timesteps in a year
 day = 24+24..endtime;                                        # Day totals
 mon = 744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016; # Month totals
 d30 = 72+72..endtime;   # FILL IN: 720/HourStep*DayStep for 30 days total
 yr  = endtime;                                              # Tstep total
# NOTE: PradM only valid for HourStep & DayStep = 1

initial
# Time calculation loops 
# ---------------------------
 Day = StartDay;    Hour = 0;

 Sold = 0;  PRadY = 0; PRadM = 0; PRadMD = 0; PRadD = 0;
 SlopMap = scalar(atan(slope(DEM)));
 report SlopMap = if(SlopMap eq 0 then scalar(0.001) else SlopMap);
 AspMap  = scalar(aspect(DEM));                      # aspect [deg]
 report AspMap  = if(AspMap le 0 then scalar(0.001) else AspMap);
 AtmPcor = ((288-0.0065*DEM)/288)**5.256;            # atm pressure corr [-]

# Calculation of Horizonts for DEM per 5ø
# NOTE: for a changing DEM in time DO NOT use the following statements
#       but put a # for these statements

 Hori01=horizontan(DEM,directional(3));
 Hori02=horizontan(DEM,directional(8));
 Hori03=horizontan(DEM,directional(13));
 Hori04=horizontan(DEM,directional(18));
 Hori05=horizontan(DEM,directional(23));
 Hori06=horizontan(DEM,directional(28));
 Hori07=horizontan(DEM,directional(33));
 Hori08=horizontan(DEM,directional(38));
 Hori09=horizontan(DEM,directional(43));
 Hori10=horizontan(DEM,directional(48));
 Hori11=horizontan(DEM,directional(53));
 Hori12=horizontan(DEM,directional(58));
 Hori13=horizontan(DEM,directional(63));
 Hori14=horizontan(DEM,directional(68));
 Hori15=horizontan(DEM,directional(73));
 Hori16=horizontan(DEM,directional(78));
 Hori17=horizontan(DEM,directional(83));
 Hori18=horizontan(DEM,directional(88));
 Hori19=horizontan(DEM,directional(93));
 Hori20=horizontan(DEM,directional(98));
 Hori21=horizontan(DEM,directional(103));
 Hori22=horizontan(DEM,directional(108));
 Hori23=horizontan(DEM,directional(113));
 Hori24=horizontan(DEM,directional(118));
 Hori25=horizontan(DEM,directional(123));
 Hori26=horizontan(DEM,directional(128));
 Hori27=horizontan(DEM,directional(133));
 Hori28=horizontan(DEM,directional(138));
 Hori29=horizontan(DEM,directional(143));
 Hori30=horizontan(DEM,directional(148));
 Hori31=horizontan(DEM,directional(153));
 Hori32=horizontan(DEM,directional(158));
 Hori33=horizontan(DEM,directional(163));
 Hori34=horizontan(DEM,directional(168));
 Hori35=horizontan(DEM,directional(173));
 Hori36=horizontan(DEM,directional(178));
 Hori37=horizontan(DEM,directional(183));
 Hori38=horizontan(DEM,directional(188));
 Hori39=horizontan(DEM,directional(193));
 Hori40=horizontan(DEM,directional(198));
 Hori41=horizontan(DEM,directional(203));
 Hori42=horizontan(DEM,directional(208));
 Hori43=horizontan(DEM,directional(213));
 Hori44=horizontan(DEM,directional(218));
 Hori45=horizontan(DEM,directional(223));
 Hori46=horizontan(DEM,directional(228));
 Hori47=horizontan(DEM,directional(233));
 Hori48=horizontan(DEM,directional(238));
 Hori49=horizontan(DEM,directional(243));
 Hori50=horizontan(DEM,directional(248));
 Hori51=horizontan(DEM,directional(253));
 Hori52=horizontan(DEM,directional(258));
 Hori53=horizontan(DEM,directional(263));
 Hori54=horizontan(DEM,directional(268));
 Hori55=horizontan(DEM,directional(273));
 Hori56=horizontan(DEM,directional(278));
 Hori57=horizontan(DEM,directional(283));
 Hori58=horizontan(DEM,directional(288));
 Hori59=horizontan(DEM,directional(293));
 Hori60=horizontan(DEM,directional(298));
 Hori61=horizontan(DEM,directional(303));
 Hori62=horizontan(DEM,directional(308));
 Hori63=horizontan(DEM,directional(313));
 Hori64=horizontan(DEM,directional(318));
 Hori65=horizontan(DEM,directional(323));
 Hori66=horizontan(DEM,directional(328));
 Hori67=horizontan(DEM,directional(333));
 Hori68=horizontan(DEM,directional(338));
 Hori69=horizontan(DEM,directional(343));
 Hori70=horizontan(DEM,directional(348));
 Hori71=horizontan(DEM,directional(353));
 Hori72=horizontan(DEM,directional(358));

 Hori01=if(Hori01 lt 0, scalar(0), Hori01);
 Hori02=if(Hori02 lt 0, scalar(0), Hori02);
 Hori03=if(Hori03 lt 0, scalar(0), Hori03);
 Hori04=if(Hori04 lt 0, scalar(0), Hori04);
 Hori05=if(Hori05 lt 0, scalar(0), Hori05);
 Hori06=if(Hori06 lt 0, scalar(0), Hori06);
 Hori07=if(Hori07 lt 0, scalar(0), Hori07);
 Hori08=if(Hori08 lt 0, scalar(0), Hori08);
 Hori09=if(Hori09 lt 0, scalar(0), Hori09); 
 Hori10=if(Hori10 lt 0, scalar(0), Hori10);
 Hori11=if(Hori11 lt 0, scalar(0), Hori11);
 Hori12=if(Hori12 lt 0, scalar(0), Hori12);
 Hori13=if(Hori13 lt 0, scalar(0), Hori13);
 Hori14=if(Hori14 lt 0, scalar(0), Hori14);
 Hori15=if(Hori15 lt 0, scalar(0), Hori15);
 Hori16=if(Hori16 lt 0, scalar(0), Hori16);
 Hori17=if(Hori17 lt 0, scalar(0), Hori17);
 Hori18=if(Hori18 lt 0, scalar(0), Hori18);
 Hori19=if(Hori19 lt 0, scalar(0), Hori19);
 Hori20=if(Hori20 lt 0, scalar(0), Hori20);
 Hori21=if(Hori21 lt 0, scalar(0), Hori21);
 Hori22=if(Hori22 lt 0, scalar(0), Hori22);
 Hori23=if(Hori23 lt 0, scalar(0), Hori23);
 Hori24=if(Hori24 lt 0, scalar(0), Hori24);
 Hori25=if(Hori25 lt 0, scalar(0), Hori25);
 Hori26=if(Hori26 lt 0, scalar(0), Hori26);
 Hori27=if(Hori27 lt 0, scalar(0), Hori27);
 Hori28=if(Hori28 lt 0, scalar(0), Hori28);
 Hori29=if(Hori29 lt 0, scalar(0), Hori29);
 Hori30=if(Hori30 lt 0, scalar(0), Hori30);
 Hori31=if(Hori31 lt 0, scalar(0), Hori31);
 Hori32=if(Hori32 lt 0, scalar(0), Hori32);
 Hori33=if(Hori33 lt 0, scalar(0), Hori33);
 Hori34=if(Hori34 lt 0, scalar(0), Hori34);
 Hori35=if(Hori35 lt 0, scalar(0), Hori35);
 Hori36=if(Hori36 lt 0, scalar(0), Hori36);
 Hori37=if(Hori37 lt 0, scalar(0), Hori37);
 Hori38=if(Hori38 lt 0, scalar(0), Hori38);
 Hori39=if(Hori39 lt 0, scalar(0), Hori39);
 Hori40=if(Hori40 lt 0, scalar(0), Hori40);
 Hori41=if(Hori41 lt 0, scalar(0), Hori41);
 Hori42=if(Hori42 lt 0, scalar(0), Hori42);
 Hori43=if(Hori43 lt 0, scalar(0), Hori43);
 Hori44=if(Hori44 lt 0, scalar(0), Hori44);
 Hori45=if(Hori45 lt 0, scalar(0), Hori45);
 Hori46=if(Hori46 lt 0, scalar(0), Hori46);
 Hori47=if(Hori47 lt 0, scalar(0), Hori47);
 Hori48=if(Hori48 lt 0, scalar(0), Hori48);
 Hori49=if(Hori49 lt 0, scalar(0), Hori49);
 Hori50=if(Hori50 lt 0, scalar(0), Hori50);
 Hori51=if(Hori51 lt 0, scalar(0), Hori51);
 Hori52=if(Hori52 lt 0, scalar(0), Hori52);
 Hori53=if(Hori53 lt 0, scalar(0), Hori53);
 Hori54=if(Hori54 lt 0, scalar(0), Hori54);
 Hori55=if(Hori55 lt 0, scalar(0), Hori55);
 Hori56=if(Hori56 lt 0, scalar(0), Hori56);
 Hori57=if(Hori57 lt 0, scalar(0), Hori57);
 Hori58=if(Hori58 lt 0, scalar(0), Hori58);
 Hori59=if(Hori59 lt 0, scalar(0), Hori59);
 Hori60=if(Hori60 lt 0, scalar(0), Hori60);
 Hori61=if(Hori61 lt 0, scalar(0), Hori61);
 Hori62=if(Hori62 lt 0, scalar(0), Hori62);
 Hori63=if(Hori63 lt 0, scalar(0), Hori63);
 Hori64=if(Hori64 lt 0, scalar(0), Hori64);
 Hori65=if(Hori65 lt 0, scalar(0), Hori65);
 Hori66=if(Hori66 lt 0, scalar(0), Hori66);
 Hori67=if(Hori67 lt 0, scalar(0), Hori67);
 Hori68=if(Hori68 lt 0, scalar(0), Hori68); 
 Hori69=if(Hori69 lt 0, scalar(0), Hori69);
 Hori70=if(Hori70 lt 0, scalar(0), Hori70);
 Hori71=if(Hori71 lt 0, scalar(0), Hori71);
 Hori72=if(Hori72 lt 0, scalar(0), Hori72);

dynamic                               
# Daily calculation loop all models
# ---------------------------
 Day    = if(Hour ge 24, Day + DayStep, Day);
 Hour   = if(Hour ge 24, scalar(1), Hour + HourStep);

# Solar geometry
# ----------------------------
# SolDec  :declination sun per day  between +23 & -23 [deg]
# HourAng :hour angle [-] of sun during day
# SolAlt  :solar altitude [deg], height of sun above horizon
 SolDec  = -23.4*cos(360*(Day+10)/365);
 HourAng = 15*(Hour-12.01);
 SolAlt  = scalar(asin(scalar(sin(Lat)*sin(SolDec)+cos(Lat)*
           cos(SolDec)*cos(HourAng))));
# report soldec.tss  = timeoutput(Loc,SolDec);
# report hourang.tss = timeoutput(Loc,HourAng);
# report solalt.tss  = timeoutput(Loc,SolAlt);

# Solar azimuth 
# ----------------------------
# SolAzi  :angle solar beams to N-S axes earth [deg]
 SolAzi = scalar(acos((sin(SolDec)*cos(Lat)-cos(SolDec)*
          sin(Lat)*cos(HourAng))/cos(SolAlt)));
 SolAzi = if(Hour le 12 then SolAzi else 360 - SolAzi);
# Additonal extra correction by R.Sluiter, Aug '99
 SolAzi = if(SolAzi gt 89.994 and SolAzi lt 90, 90, SolAzi);
 SolAzi = if(SolAzi gt 269.994 and SolAzi lt 270, 270, SolAzi);
# report solazi.tss = timeoutput(Loc,SolAzi);

# Surface azimuth
# ----------------------------
# cosIncident :cosine of angle of incident; angle solar beams to angle surface
 cosIncident = sin(SolAlt)*cos(SlopMap)+cos(SolAlt)*sin(SlopMap)
               *cos(SolAzi-AspMap);
# report cosi.tss = timeoutput(Loc,cosIncident);

# Critical angle sun
# ----------------------------
# HoriAng  :tan maximum angle over DEM in direction sun, 0 if neg 
# CritSun  :tan of maximum angle in direction solar beams
# Shade    :cell in sun 1, in shade 0
# NOTE: for a changing DEM in time use following 3 statements and put a #
#       for the 4th CritSun statement
# HoriAng   = horizontan(DEM,directional(SolAzi));
# HoriAng   = if(HoriAng lt 0 then scalar(0) else HoriAng);
# CritSun   = if(SolAlt gt 90 then scalar(0) else scalar(atan(HoriAng)));
 CritSun = if(SolAlt gt 90.0 then scalar(0) else
           if(SolAzi le 5.0 then scalar(atan(Hori01)) else
           if(SolAzi le 10.0 then scalar(atan(Hori02)) else
           if(SolAzi le 15.0 then scalar(atan(Hori03)) else
           if(SolAzi le 20.0 then scalar(atan(Hori04)) else
           if(SolAzi le 25.0 then scalar(atan(Hori05)) else
           if(SolAzi le 30.0 then scalar(atan(Hori06)) else
           if(SolAzi le 35.0 then scalar(atan(Hori07)) else
           if(SolAzi le 40.0 then scalar(atan(Hori08)) else
           if(SolAzi le 45.0 then scalar(atan(Hori09)) else
           if(SolAzi le 50.0 then scalar(atan(Hori10)) else
           if(SolAzi le 55.0 then scalar(atan(Hori11)) else
           if(SolAzi le 60.0 then scalar(atan(Hori12)) else
           if(SolAzi le 65.0 then scalar(atan(Hori13)) else
           if(SolAzi le 70.0 then scalar(atan(Hori14)) else
           if(SolAzi le 75.0 then scalar(atan(Hori15)) else
           if(SolAzi le 80.0 then scalar(atan(Hori16)) else
           if(SolAzi le 85.0 then scalar(atan(Hori17)) else
           if(SolAzi le 90.0 then scalar(atan(Hori18)) else
           if(SolAzi le 95.0 then scalar(atan(Hori19)) else
           if(SolAzi le 100.0 then scalar(atan(Hori20)) else
           if(SolAzi le 105.0 then scalar(atan(Hori21)) else
           if(SolAzi le 110.0 then scalar(atan(Hori22)) else
           if(SolAzi le 115.0 then scalar(atan(Hori23)) else
           if(SolAzi le 120.0 then scalar(atan(Hori24)) else
           if(SolAzi le 125.0 then scalar(atan(Hori25)) else
           if(SolAzi le 130.0 then scalar(atan(Hori26)) else
           if(SolAzi le 135.0 then scalar(atan(Hori27)) else
           if(SolAzi le 140.0 then scalar(atan(Hori28)) else
           if(SolAzi le 145.0 then scalar(atan(Hori29)) else
           if(SolAzi le 150.0 then scalar(atan(Hori30)) else
           if(SolAzi le 155.0 then scalar(atan(Hori31)) else
           if(SolAzi le 160.0 then scalar(atan(Hori32)) else
           if(SolAzi le 165.0 then scalar(atan(Hori33)) else
           if(SolAzi le 170.0 then scalar(atan(Hori34)) else
           if(SolAzi le 175.0 then scalar(atan(Hori35)) else
           if(SolAzi le 180.0 then scalar(atan(Hori36)) else
           if(SolAzi le 185.0 then scalar(atan(Hori37)) else
           if(SolAzi le 190.0 then scalar(atan(Hori38)) else
           if(SolAzi le 195.0 then scalar(atan(Hori39)) else
           if(SolAzi le 200.0 then scalar(atan(Hori40)) else
           if(SolAzi le 205.0 then scalar(atan(Hori41)) else
           if(SolAzi le 210.0 then scalar(atan(Hori42)) else
           if(SolAzi le 215.0 then scalar(atan(Hori43)) else
           if(SolAzi le 220.0 then scalar(atan(Hori44)) else
           if(SolAzi le 225.0 then scalar(atan(Hori45)) else
           if(SolAzi le 230.0 then scalar(atan(Hori46)) else
           if(SolAzi le 235.0 then scalar(atan(Hori47)) else
           if(SolAzi le 240.0 then scalar(atan(Hori48)) else
           if(SolAzi le 245.0 then scalar(atan(Hori49)) else
           if(SolAzi le 250.0 then scalar(atan(Hori50)) else
           if(SolAzi le 255.0 then scalar(atan(Hori51)) else
           if(SolAzi le 260.0 then scalar(atan(Hori52)) else
           if(SolAzi le 265.0 then scalar(atan(Hori53)) else
           if(SolAzi le 270.0 then scalar(atan(Hori54)) else
           if(SolAzi le 275.0 then scalar(atan(Hori55)) else
           if(SolAzi le 280.0 then scalar(atan(Hori56)) else
           if(SolAzi le 285.0 then scalar(atan(Hori57)) else
           if(SolAzi le 290.0 then scalar(atan(Hori58)) else
           if(SolAzi le 295.0 then scalar(atan(Hori59)) else
           if(SolAzi le 300.0 then scalar(atan(Hori60)) else
           if(SolAzi le 305.0 then scalar(atan(Hori61)) else
           if(SolAzi le 310.0 then scalar(atan(Hori62)) else
           if(SolAzi le 315.0 then scalar(atan(Hori63)) else
           if(SolAzi le 320.0 then scalar(atan(Hori64)) else
           if(SolAzi le 325.0 then scalar(atan(Hori65)) else
           if(SolAzi le 330.0 then scalar(atan(Hori66)) else
           if(SolAzi le 335.0 then scalar(atan(Hori67)) else
           if(SolAzi le 340.0 then scalar(atan(Hori68)) else
           if(SolAzi le 345.0 then scalar(atan(Hori69)) else
           if(SolAzi le 350.0 then scalar(atan(Hori70)) else
           if(SolAzi le 355.0 then scalar(atan(Hori71)) else
              scalar(atan(Hori72))
 ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))));  
 Shade   = if(SolAlt gt CritSun then scalar(1) else scalar(0));
# Horiang.tss  = timeoutput(Loc,HoriAng);
# report crit.tss  = timeoutput(Loc,CritSun);
# report shade.tss = timeoutput(Loc,Shade);

# Radiation outer atmosphere
# ----------------------------
 OpCorr = Trans**((sqrt(1229+(614*sin(SolAlt))**2)
          -614*sin(SolAlt))*AtmPcor);     # correction for air masses [-] 
 Sout   = Sc*(1+0.034*cos(360*Day/365)); # radiation outer atmosphere [W/m2]
 Snor   = Sout*OpCorr;                    # rad on surface normal to the beam [W/m2]
# report snor.tss = timeoutput(Loc,Snor);

# Radiation at DEM
# ----------------------------
# Sdir   :direct sunlight on a horizontal surface [W/m2] if no shade
# Sdiff  :diffuse light [W/m2] for shade and no shade
# Stot   :total incomming light Sdir+Sdiff [W/m2] at Hour
# PotRad :avg of Stot(Hour) and Stot(Hour-HourStep)
# NOTE: PradM only valid for HourStep & DayStep = 1
 Sdir   = if(Snor*cosIncident*Shade<0,0.0,Snor*cosIncident*Shade);
 Sdiff  = if(Sout*(0.271-0.294*OpCorr)*sin(SolAlt)<0, 0.0,
          Sout*(0.271-0.294*OpCorr)*sin(SolAlt));
# Stot   = cover(Sdir+Sdiff,windowaverage(Sdir+Sdiff,3));        # Rad [W/m2]
 Stot   = Sdir + Sdiff;                                         # Rad [W/m2]
 PotRad = (Sold + Stot)/2;                                # Rad interval
#report(day) PRadD  = PRadD + (PotRad*0.0036*HourStep);          # day rad
report(mon) PRadM  = PRadM + (PotRad*0.0036*HourStep*DayStep);  # month rad
report(d30) PRadMD = PRadMD + (PotRad*0.0036*HourStep*DayStep); # 30days rad
report(yr)  PRadY  = PRadY + (PotRad*0.0036*HourStep*DayStep);  # year rad
#report potrad.tss = timeoutput(Loc,PotRad);
#report sdir.tss=timeoutput(Loc,Sdir);
#report Sdif.tss=timeoutput(Loc,Sdiff);
#report Stotal.tss=timeoutput(Loc,Stot);

# Day & month loops & reassign maps for next timestep
 PRadD  = if(Hour ge 24 then scalar(0) else PRadD);
 PRadM  = if(time() eq 744 or time() eq 1416 or time() eq 2160 or
          time() eq 2880 or time() eq 3624 or time() eq 4344 or
          time() eq 5088 or time() eq 5832 or time() eq 6552 or
          time() eq 7296 or time() eq 8016, scalar(0), PRadM);
 PRadMD = if(time() idiv (720/(HourStep*DayStep)) eq
            time()/(720/(HourStep*DayStep)), scalar(0), PRadMD);
 Sold = Stot;

