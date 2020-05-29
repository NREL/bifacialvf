# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 11:16:40 2018

@author: Silvana
"""

import math
import csv
import pvlib
import os
import pytz

from vf import getBackSurfaceIrradiances, getFrontSurfaceIrradiances, getGroundShadeFactors
from vf import getSkyConfigurationFactors, trackingBFvaluescalculator, rowSpacing
from sun import hrSolarPos, perezComp, solarPos, sunIncident
import pandas as pd
import datetime
from readepw import readepw

def simulateNSRDB(NSRDBdata, writefiletitle,  beta = 0, sazm = 180, C = 0.5, D = None,
             rowType = 'interior', transFactor = 0.01, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = 0.2,  
             tracking = False, backtrack = True, rtr = None, Cv = None, offset = 0, 
             lat = 39.742, lng = -105.179, tz = -7):


        dataInterval = 60
        
        if tracking == True:
            axis_tilt = 0  # algorithm only allows for zero north-south tilt with SAT
            max_angle = 90  # maximum tracker rotation 
            axis_azimuth=sazm    # axis_azimuth is degrees east of North
            beta = 0            # start with tracker tilt = 0
            hub_height = C      # Ground clearance at tilt = 0.  C >= 0.5
            if hub_height < 0.5:
                print('Warning: tracker hub height C < 0.5 may result in ground clearance errors')
        
        if (D == None) & (rtr != None):
            D = rtr - math.cos(beta / 180.0 * math.pi)
        elif (rtr == None) & (D != None):
            rtr = D + math.cos(beta / 180.0 * math.pi)
        elif (D == None) & (rtr == None):
            raise Exception('No row distance specified in either D or rtr') 
        else:
            print('Warning: Gap D and rtr passed in. Using ' + ('rtr' if tracking else 'D') )
        
        noRows, noCols = NSRDBdata.shape
        
        ## Distance between rows for no shading on Dec 21 at 9 am
        print " "
        print "********* "
        print "Running Simulation for SRRL DataFile: "
        print "Location:  ", name
        print "Lat: ", lat, " Long: ", lng, " Tz ", tz
        print "Parameters: beta: ", beta, "  Sazm: ", sazm, "  Height: ", C, "  rtr separation: ", rtr, "  Row type: ", rowType, "  Albedo: ", albedo
        print "Saving into", writefiletitle
        print " "
        print " "
        
        DD = rowSpacing(beta, sazm, lat, lng, tz, 9, 0.0);          ## Distance between rows for no shading on Dec 21 at 9 am
        print "Distance between rows for no shading on Dec 21 at 9 am solar time = ", DD
        print "Actual distance between rows = ", D  
        print " "
    
        if tracking==False:        
            ## Sky configuration factors are the same for all times, only based on geometry and row type
            [rearSkyConfigFactors, frontSkyConfigFactors, ffConfigFactors] = getSkyConfigurationFactors(rowType, beta, C, D);       ## Sky configuration factors are the same for all times, only based on geometry and row type
    
     
        ## Create WriteFile and write labels at this time
        
        #check that the save directory exists
        if not os.path.exists(os.path.dirname(writefiletitle)):
            os.makedirs(os.path.dirname(writefiletitle))
        
        with open (writefiletitle,'wb') as csvfile:
            sw = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            # Write Simulation Parameters (from setup file)
            
            outputheader=['Latitude(deg)','Longitude(deg)', 'Time Zone','Tilt(deg)', 
                         'PV Azimuth(deg)','GroundClearance(panel slope lengths)', 'Row-to-Row-Distance rtr', 'RowType(first interior last single)',
                         'TransmissionFactor(open area fraction)','CellRows(# hor rows in panel)', 
                         'PVfrontSurface(glass or AR glass)', 'PVbackSurface(glass or AR glass)',
                         'CellOffsetFromBack(panel slope lengths)','Albedo',  'Tracking']
            outputheadervars=[lat, lng, tz, beta, sazm, C, rtr, rowType, transFactor, cellRows, PVfrontSurface,
                             PVbackSurface, offset, albedo, tracking]
            
            
            if tracking==True:
                outputheader+=['Backtracking']
                outputheadervars.append(backtrack)
                
            sw.writerow(outputheader)
            sw.writerow(outputheadervars)
            
            # Write Results names"
            allrowfronts=[]
            allrowbacks=[]
            for k in range(0, cellRows):
                allrowfronts.append("No_"+str(k+1)+"_RowFrontGTI")
                allrowbacks.append("No_"+str(k+1)+"_RowBackGTI")      
            outputtitles=['Year', 'Month', 'Day', 'Hour', 'Minute', 'DNI', 'DHI', 
                         'decHRs', 'ghi', 'inc', 'zen', 'azm', 'pvFrontSH', 
                         'aveFrontGroundGHI', 'GTIfrontBroadBand', 'pvBackSH', 
                         'aveBackGroundGHI', 'GTIbackBroadBand', 'maxShadow', 'Tamb', 'Vwind']
            outputtitles+=allrowfronts
            outputtitles+=allrowbacks
            if tracking == True:
                print " ***** IMPORTANT --> THIS SIMULATION Has Tracking Activated"
                print "Backtracking Option is set to: ", backtrack
                outputtitles+=['beta']
                outputtitles+=['sazm']
                outputtitles+=['height']
                outputtitles+=['D']
                    
            sw.writerow(outputtitles)
            
            ## Loop through file
            rl = 0
            
            while (rl < noRows):
          
                myTimestamp=NSRDBdata.index[rl]
                year = myTimestamp.year
                month = myTimestamp.month
                day = myTimestamp.day
                hour = myTimestamp.hour
                minute = myTimestamp.minute
                dni = NSRDBdata.DNI[rl]#get_value(rl,5,"False")
                dhi = NSRDBdata.DHI[rl]#get_value(rl,8,"False")
                Tamb=NSRDBdata.Temp[rl]#get_value(rl,29,"False")
                Vwind=NSRDBdata.Wspd[rl]#get_value(rl,44,"False")
           #     
                rl = rl+1   # increasing while count
                            
                azm = 9999.0; zen = 9999.0; elv = 9999.0;
                if (dataInterval == 60):
                    azm, zen, elv, dec, sunrise, sunset, Eo, tst, suntime = hrSolarPos(year, month, day, hour, lat, lng, tz)
                    
                elif (dataInterval == 1 or dataInterval == 5 or dataInterval == 15):
                    azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos(year, month, day, hour, minute - 0.5 * dataInterval, lat, lng, tz) 
                else :  
                    print("ERROR: data interval not 1, 5, 15, or 60 minutes.");
            
                #123 check abouve this for reading / printing functions
            
                if (zen < 0.5 * math.pi):    # If daylight hours
                
                    # a. CALCULATE THE IRRADIANCE DISTRIBUTION ON THE GROUND *********************************************************************************************
                    #double[] rearGroundGHI = new double[100], frontGroundGHI = new double[100]; ;   # For global horizontal irradiance for each of 100 ground segments, to the rear and front of front of row edge         
                    # Determine where on the ground the direct beam is shaded for a sun elevation and azimuth
                    #int[] rearGroundSH = new int[100], frontGroundSH = new int[100]; # Front and rear row-to-row spacing divided into 100 segments, (later becomes 1 if direct beam is shaded, 0 if not shaded)
                    #double pvFrontSH = 0.0, pvBackSH = 0.0, maxShadow;     # Initialize fraction of PV module front and back surfaces that are shaded to zero (not shaded), and maximum shadow projected from front of row.
                    
                    # TRACKING ROUTINE CALULATING GETSKYCONFIGURATION FACTORS
                    if tracking == True:        
                        
                        #solpos = pvlib.solarposition.get_solarposition(myTimestamp, lat, lng)
                        #aazi= solpos['azimuth']
                        #azen= solpos['zenith']
                        aazi = pd.Series([azm*180.0/math.pi], index =[myTimestamp])                        
                        azen = pd.Series([zen*180.0/math.pi], index =[myTimestamp])

                        
                        gcr=1/rtr  
                        trackingdata = pvlib.tracking.singleaxis(azen, aazi, axis_tilt, axis_azimuth, max_angle, backtrack, gcr)
                                 ## Sky configuration factors are not the same for all times, since the geometry is changing with the tracking.
                        beta=trackingdata['surface_tilt'][0] # Trackingdata tracker_theta
                        sazm = trackingdata['surface_azimuth'][0]
                        if math.isnan(beta):
                            beta=90
    
                        # Rotate system if past sun's zenith ~ #123 Check if system breaks withot doing this.
                        if beta<0:
                            #sazm = sazm+180    # Rotate detectors
                            beta = -beta;
                            
                        [C, D] = trackingBFvaluescalculator(beta, hub_height, rtr)
                        [rearSkyConfigFactors, frontSkyConfigFactors] = getSkyConfigurationFactors(rowType, beta, C, D);       ## Sky configuration factors are the same for all times, only based on geometry and row type
    
    
    
                    rearGroundGHI=[]
                    frontGroundGHI=[]
                    pvFrontSH, pvBackSH, maxShadow, rearGroundSH, frontGroundSH = getGroundShadeFactors (rowType, beta, C, D, elv, azm, sazm)
            
                    # Sum the irradiance components for each of the ground segments, to the front and rear of the front of the PV row
                    #double iso_dif = 0.0, circ_dif = 0.0, horiz_dif = 0.0, grd_dif = 0.0, beam = 0.0;   # For calling PerezComp to break diffuse into components for zero tilt (horizontal)                           
                    ghi, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, zen, 0.0, zen)
                    
                    
                    for k in range (0, 100):
                    
                        rearGroundGHI.append(iso_dif * rearSkyConfigFactors[k]);       # Add diffuse sky component viewed by ground
                        if (rearGroundSH[k] == 0):
                            rearGroundGHI[k] += beam + circ_dif;                    # Add beam and circumsolar component if not shaded
                        else:
                            rearGroundGHI[k] += (beam + circ_dif) * transFactor;    # Add beam and circumsolar component transmitted thru module spacing if shaded
            
                        frontGroundGHI.append(iso_dif * frontSkyConfigFactors[k]);     # Add diffuse sky component viewed by ground
                        if (frontGroundSH[k] == 0):
                            frontGroundGHI[k] += beam + circ_dif;                   # Add beam and circumsolar component if not shaded 
                        else:
                            frontGroundGHI[k] += (beam + circ_dif) * transFactor;   # Add beam and circumsolar component transmitted thru module spacing if shaded
                    
            
                    # b. CALCULATE THE AOI CORRECTED IRRADIANCE ON THE FRONT OF THE PV MODULE, AND IRRADIANCE REFLECTED FROM FRONT OF PV MODULE ***************************
                    #double[] frontGTI = new double[cellRows], frontReflected = new double[cellRows];
                    #double aveGroundGHI = 0.0;          # Average GHI on ground under PV array
                    aveGroundGHI, frontGTI, frontReflected = getFrontSurfaceIrradiances(rowType, maxShadow, PVfrontSurface, beta, sazm, dni, dhi, C, D, albedo, zen, azm, cellRows, pvFrontSH, frontGroundGHI)
    
                    #double inc, tiltr, sazmr;
                    inc, tiltr, sazmr = sunIncident(0, beta, sazm, 45.0, zen, azm)	    # For calling PerezComp to break diffuse into components for 
                    save_inc=inc
                    gtiAllpc, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, inc, tiltr, zen)   # Call to get components for the tilt
                    save_gtiAllpc=gtiAllpc
                    #sw.Write(strLine);
                    #sw.Write(",{0,6:0.00}", hour - 0.5 * dataInterval / 60.0 + minute / 60.0);
                    #sw.Write(",{0,6:0.0},{1,5:0.0},{2,5:0.0},{3,5:0.0},{4,4:0.00},{5,6:0.0},{6,6:0.0}",
                        #dni * Math.Cos(zen) + dhi, inc * 180.0 / Math.PI, zen * 180.0 / Math.PI, azm * 180.0 / Math.PI, pvFrontSH, aveGroundGHI, gtiAllpc);
            
                    # CALCULATE THE AOI CORRECTED IRRADIANCE ON THE BACK OF THE PV MODULE,
                    #double[] backGTI = new double[cellRows];
                    backGTI, aveGroundGHI = getBackSurfaceIrradiances(rowType, maxShadow, PVbackSurface, beta, sazm, dni, dhi, C, D, albedo, zen, azm, cellRows, pvBackSH, rearGroundGHI, frontGroundGHI, frontReflected, offset=0)
               
                    inc, tiltr, sazmr = sunIncident(0, 180.0-beta, sazm-180.0, 45.0, zen, azm)       # For calling PerezComp to break diffuse into components for 
                    gtiAllpc, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, inc, tiltr, zen)   # Call to get components for the tilt
                    
                    
                    ## Write output
                    decHRs = hour - 0.5 * dataInterval / 60.0 + minute / 60.0
                    ghi_calc = dni * math.cos(zen) + dhi 
                    incd = save_inc * 180.0 / math.pi
                    zend = zen * 180.0 / math.pi
                    azmd = azm * 180.0 / math.pi
                    outputvalues=[year, month, day, hour, minute, dni, dhi, decHRs, 
                                  ghi_calc, incd, zend, azmd, pvFrontSH, aveGroundGHI, 
                                  save_gtiAllpc, pvBackSH, aveGroundGHI, 
                                  gtiAllpc, maxShadow, Tamb, Vwind]
                    frontGTIrow=[]
                    backGTIrow=[]
                    for k in range(0, cellRows):
                        frontGTIrow.append(frontGTI[k])
                        backGTIrow.append(backGTI[k])      
                    outputvalues+=frontGTIrow
                    outputvalues+=backGTIrow
                    
                    
                    if tracking==True:
                        outputvalues.append(beta)
                        outputvalues.append(sazm)
                        outputvalues.append(C)
                        outputvalues.append(D)
                                    
                    sw.writerow(outputvalues)
    
        	# End of daylight if loop 
    
        #strLine = sr.ReadLine();        # Read next line of data
       # End of while strLine != null loop
       
     
        print( "Finished")
        
        return;



#NSRDBfile='C:\Users\Silvana\Desktop\SOLTEC\\557152_31.21_-102.18_tmy.csv'
def readNSRDB(NSRDB_Filetitle, decimate= False, decimatestyle = "Mean" , decimatevalue = 5 ):      
    # decimatestyle = "Mean" or "Skip"
    day_all = []; month_all = []; year_all = []; hour_all = []; minute_all= []
    global_all=[]; direct_all=[]; diffuse_all=[]; 
    albedo_all=[]; temp_all=[]; pressure_all=[]; wind_all=[]
    
    timestamps_all = [];
    
    
    headeracquired= 0
    headererror = 0
    
    yearloc=0
    monthloc=1
    dayloc=2
    hourloc=3
    minuteloc=4
    directloc = 5
    diffuseloc = 6
    globalloc= 7
    albedoloc=8
    windspeedloc=9
    temperatureloc=10
    pressureloc=11
    
    
    with open(NSRDBfile, "r") as filestream:
    
        print "Reading NSRDBfile File: ", NSRDBfile
        print "Decimation options: ", decimate
        if decimate == True:
            print "Decimate Style: ", decimatestyle, " to interval ", decimatevalue, " minutes."
        headerline = filestream.readline() # read headers labels, read values
        head1 = headerline.strip().split(',')
        
        meta = dict(zip(head1, filestream.readline().rstrip('\n').split(",")))
            
        for line in filestream:
            if headeracquired == 0:
                header = line.split(",")
                        
                if header[yearloc] != 'Year': print "Issue reading" + header [yearloc] ; headererror = 1
                if header[directloc] != 'DHI': print "Issue reading" + header [directloc] ; headererror = 1
                if header[diffuseloc] != 'DNI': print "Issue reading" + header [diffuseloc] ; headererror = 1
                if header[globalloc] != 'GHI': print "Issue reading" + header [globalloc] ; headererror = 1
            
                headeracquired = 1
                
                if headererror == 1:
                    print "STOPPING File Read because of headers issue (expected data might not be where we think it is! Stop roll and check!"
                    continue
                
            else:
                
                if headererror == 1:
                    continue
            
                currentline=line.split(",")                
                year=int(currentline[yearloc])
                month=int(currentline[monthloc])
                day=int(currentline[dayloc])
                hour=int(currentline[hourloc])
                minute=int(currentline[minuteloc])
                
                if decimate == True:
                    if decimatestyle == "Skip":
                        if minute%decimatevalue != 0:
                              continue
    
                year_all.append(int(currentline[yearloc]))
                month_all.append(int(currentline[monthloc]))
                day_all.append(int(currentline[dayloc]))
                hour_all.append(int(currentline[hourloc]))
                minute_all.append(int(currentline[minuteloc]))
                direct_all.append(float(currentline[directloc]))
                diffuse_all.append(float(currentline[diffuseloc]))
                global_all.append(float(currentline[globalloc]))
                albedo_all.append(float(currentline[albedoloc]))
                wind_all.append(float(currentline[windspeedloc]))
                temp_all.append(float(currentline[temperatureloc]))
                pressure_all.append(float(currentline[pressureloc]))
                myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)                
                timestamps_all.append(myTimestamp)
                
    NSRDB = ({'Month': month_all, 'Day': day_all, 'Year': year_all, 'Hour': hour_all, 'Minute': minute_all,
                 'GHI': global_all, 'DNI': direct_all, 'DHI': diffuse_all, 
                 'Temp': temp_all, 'Presure': pressure_all, 'Wspd': wind_all, 'Albedo': albedo_all})
    
    NSRDB = pd.DataFrame.from_records(NSRDB, index=timestamps_all)
    
    
    if decimate == True:
        if decimatestyle == "Mean":
            if decimatevalue == 5:
                    NSRDB=NSRDB.resample('5Min', base=0).mean()
                    print "Data decimated to 5 min Interval by Average"
    
            if decimatevalue == 10:
                    NSRDB=NSRDB.resample('10Min').mean()
                    print "Data decimated to 10 min Interval by Average"
    
            if decimatevalue == 15:
                    NSRDB=NSRDB.resample('15Min').mean()
                    print "Data decimated to 15 min Interval by Average"
    
            if decimatevalue == 60:
                    NSRDB=NSRDB.resample('60Min').mean()
                    print "Data decimated to 1 Hr Interval by Average"
    
    return NSRDB, meta;                 



    
if __name__ == "__main__":

    NSRDBfile='C:\Users\Silvana\Desktop\SOLTEC\\557152_31.21_-102.18_tmy.csv'
    NSRDB, metdata = readNSRDB(NSRDBfile, decimate= False, decimatestyle = "Mean" , decimatevalue = 5 ) 
    TMYtoread='C:\Users\Silvana\Desktop\SOLTEC\\USA_TX_Midland-Odessa.722650_TMY2.epw'
    TMYtoread='C:\Users\Silvana\Desktop\SOLTEC\\Midland_722650TYA.csv'
    useTMY = True  # either use Midland - Odessa EPW (False) or NSRDB TMY2 file (True)
    
        ## Read TMY3 data and start loop ~  
    if TMYtoread.endswith('.csv') :  
        (myTMY3,meta)=pvlib.tmy.readtmy3(TMYtoread)
    elif TMYtoread.endswith('.epw') : 
        (myTMY3,meta) = readepw(TMYtoread)
        # rename different field parameters to match DNI, DHI, DryBulb, Wspd
        myTMY3.rename(columns={'Direct normal radiation in Wh/m2':'DNI','Diffuse horizontal radiation in Wh/m2':'DHI','Dry bulb temperature in C':'DryBulb','Wind speed in m/s':'Wspd', 'Global horizontal radiation in Wh/m2':'GHI'}, inplace=True)
    else:
        raise Exception('Incorrect extension for TMYtoread. Either .csv (TMY3) .epw or None')

   
   # NSRDB.index=NSRDB.index.tz_localize('Etc/GMT+6')
    #myAxisTitles=myTMY3.axes
    noRows, noCols = myTMY3.shape
    lat = meta['latitude']; lng = meta['longitude']; tz = meta['TZ']
    name = meta['Name']

    if useTMY == True:  # if true, substitute NSRDB TMY2 data for EPW file.
        NSRDB.DHI = list(myTMY3.DHI)
        NSRDB.DNI = list(myTMY3.DNI)
        NSRDB.GHI = list(myTMY3.GHI)
        NSRDB.DNI = NSRDB.DNI.shift(1)
        NSRDB.DHI = NSRDB.DHI.shift(1)
        NSRDB.GHI = NSRDB.GHI.shift(1)
        NSRDB.DNI[0] = 0.0; NSRDB.DHI[0] = 0.0; NSRDB.GHI[0] = 0.0
    
    # SIMULATION 1:
    # 1. 2-up portrait, no torque tube, standard simulation. There is a 15 cm gap between the modules. Clearance is 2.35m, GCR 0.33
    module_height=(1.996*2)+0.15
    gcr=0.33
    rtr = 1/gcr
    hubheight= 2.35
    C = hubheight/module_height
    albedo = 0.28
    lat = 31.21
    lng = -102.18
    tz = -6
    writefiletitle ='C:\Users\Silvana\Desktop\SOLTEC\\Simul1_VF_Results_UseTMYCsvTrue_downshift.csv'
        
    simulateNSRDB(NSRDB, writefiletitle,  beta = 0, sazm = 180, C = C, D = None,
             rowType = 'interior', transFactor = 0.01, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = albedo,  
             tracking = True, backtrack = True, rtr = rtr, Cv = None, offset = 0, 
             lat = lat, lng = lng, tz = tz)


    # SIMULATION 1:
    # 2. 1-up portrait, no torque tube, standard simulation. 1.35m clearance height, 0.33 GCR
    module_height=1.996
    gcr=0.33
    rtr = 1/gcr
    hubheight= 1.35
    C = hubheight/module_height
    albedo = 0.28
    lat = 31.21
    lng = -102.18
    tz = -6
    writefiletitle ='C:\Users\Silvana\Desktop\SOLTEC\\Simul2_VF_Results_UseTMYCsvTrue_downshift.csv'
        
    simulateNSRDB(NSRDB, writefiletitle,  beta = 0, sazm = 180, C = C, D = None,
             rowType = 'interior', transFactor = 0.01, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = albedo,  
             tracking = True, backtrack = True, rtr = rtr, Cv = None, offset = 0, 
             lat = lat, lng = lng, tz = tz)


   