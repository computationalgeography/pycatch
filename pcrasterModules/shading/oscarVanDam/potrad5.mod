# POTRAD5.MOD
# Potential Radiation Equator model
# (c) O. van Dam, UU, Tropenbos-Guyana
# Version 5, June 2000
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
 HoriAng   = horizontan(DEM,directional(SolAzi));
 HoriAng   = if(HoriAng lt 0 then scalar(0) else HoriAng);
 CritSun   = if(SolAlt gt 90 then scalar(0) else scalar(atan(HoriAng)));
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

