# -*- coding: utf-8 -*-
"""
Created on Tue Nov 03 12:11:07 2017

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
        
def readSRRL(SRRL_Filetitle, decimate= True, decimatestyle = "Mean" , decimatevalue = 5 ):      
    # decimatestyle = "Mean" or "Skip"
    date_all=[]; time_all=[]; global_all=[]; direct_all=[]
    diffuse_all=[]; temp_all=[]; relhum_all=[]; wind_all=[]
    day_all = []; month_all = []; year_all = []; hour_all = []; minute_all= []
    timestamps_all = [];
    
    headeracquired= 0
    headererror = 0
    
    dateloc = 0
    timeloc = 1
    globalloc= 2
    directloc = 4
    diffuseloc = 9
    temploc = 13
    relhumloc = 14
    windloc = 15
    
    
    with open(SRRL_Filetitle, "r") as filestream:
    
        print "Reading SRRL File: ", SRRL_Filetitle
        print "Decimation options: ", decimate
        if decimate == True:
            print "Decimate Style: ", decimatestyle, " to interval ", decimatevalue, " minutes."
            
        for line in filestream:
            if headeracquired == 0:
                header = line.split(",")
                        
                if header[dateloc] != 'DATE (MM/DD/YYYY)': print "Issue reading" + header [dateloc] ; headererror = 1
                if header[timeloc] != 'MST': print "Issue reading" + header [timeloc] ; headererror = 1
                if header[globalloc] != 'Global CMP22 (vent/cor) [W/m^2]': print "Issue reading" + header [globalloc] ; headererror = 1
                if header[directloc] != 'Direct CHP1-1 [W/m^2]': print "Issue reading" + header [directloc] ; headererror = 1
                if header[diffuseloc] != 'Diffuse 8-48 (vent) [W/m^2]': print "Issue reading" + header [diffuseloc] ; headererror = 1
                if header[temploc] != 'Tower Dry Bulb Temp [deg C]': print "Issue reading" + header [temploc] ; headererror = 1 
                if header[relhumloc] != 'Tower RH [%]': print "Issue reading" + header [relhumloc] ; headererror = 1
                if header[windloc] != 'Avg Wind Speed @ 6ft [m/s]\n': print "Issue reading" + header [windloc] ; headererror = 1
            
                headeracquired = 1
                
                if headererror == 1:
                    print "STOPPING File Read because of headers issue"
                    continue
                
            else:
                
                if headererror == 1:
                    continue
            
                currentline=line.split(",")                
                month, day, year = currentline[dateloc].split("/")
                month = int(month); day = int(day); year = int(year)
                hour, minute = currentline[timeloc].split(":")
                hour = int(hour); minute = int(minute)
                wind, trash = currentline[windloc].split("\n")

                if decimate == True:
                    if decimatestyle == "Skip":
                        if minute%decimatevalue != 0:
                              continue

                date_all.append(currentline[dateloc])
                time_all.append(currentline[timeloc])
                global_all.append(float(currentline[globalloc]))
                direct_all.append(float(currentline[directloc]))
                diffuse_all.append(float(currentline[diffuseloc]))
                temp_all.append(float(currentline[temploc]))
                relhum_all.append(float(currentline[relhumloc]))
                day_all.append(day)
                month_all.append(month)
                year_all.append(year)
                hour_all.append(hour)
                minute_all.append(minute)
                wind_all.append(float(wind))
                myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)                
                timestamps_all.append(myTimestamp)
                
    SRRL = ({'Month': month_all, 'Day': day_all, 'Year': year_all, 'Hour': hour_all, 'Minute': minute_all,
                 'Global': global_all, 'Direct': direct_all, 'Diffuse': diffuse_all, 
                 'Temp': temp_all, 'Humidity': relhum_all, 'Wind': wind_all})

    SRRL = pd.DataFrame.from_records(SRRL, index=timestamps_all)


    if decimate == True:
        if decimatestyle == "Mean":
            if decimatevalue == 5:
                    SRRL=SRRL.resample('5Min', base=0).mean()
                    print "Data decimated to 5 min Interval by Average"
    
            if decimatevalue == 10:
                    SRRL=SRRL.resample('10Min').mean()
                    print "Data decimated to 10 min Interval by Average"
    
            if decimatevalue == 15:
                    SRRL=SRRL.resample('15Min').mean()
                    print "Data decimated to 15 min Interval by Average"
    
            if decimatevalue == 60:
                    SRRL=SRRL.resample('60Min').mean()
                    print "Data decimated to 1 Hr Interval by Average"

    return SRRL;                 

def readSRRL(SRRL_Filetitle, decimate= True, decimatestyle = "Mean" , decimatevalue = 5 ):      
    # decimatestyle = "Mean" or "Skip"
    date_all=[]; time_all=[]; global_all=[]; direct_all=[]
    diffuse_all=[]; temp_all=[]; relhum_all=[]; wind_all=[]
    day_all = []; month_all = []; year_all = []; hour_all = []; minute_all= []
    timestamps_all = [];
    
    headeracquired= 0
    headererror = 0
    
    dateloc = 0
    timeloc = 1
    globalloc= 2
    directloc = 4
    diffuseloc = 9
    temploc = 13
    relhumloc = 14
    windloc = 15
    
    
    with open(SRRL_Filetitle, "r") as filestream:
    
        print "Reading SRRL File: ", SRRL_Filetitle
        print "Decimation options: ", decimate
        if decimate == True:
            print "Decimate Style: ", decimatestyle, " to interval ", decimatevalue, " minutes."
            
        for line in filestream:
            if headeracquired == 0:
                header = line.split(",")
                        
                if header[dateloc] != 'DATE (MM/DD/YYYY)': print "Issue reading" + header [dateloc] ; headererror = 1
                if header[timeloc] != 'MST': print "Issue reading" + header [timeloc] ; headererror = 1
                if header[globalloc] != 'Global CMP22 (vent/cor) [W/m^2]': print "Issue reading" + header [globalloc] ; headererror = 1
                if header[directloc] != 'Direct CHP1-1 [W/m^2]': print "Issue reading" + header [directloc] ; headererror = 1
                if header[diffuseloc] != 'Diffuse 8-48 (vent) [W/m^2]': print "Issue reading" + header [diffuseloc] ; headererror = 1
                if header[temploc] != 'Tower Dry Bulb Temp [deg C]': print "Issue reading" + header [temploc] ; headererror = 1 
                if header[relhumloc] != 'Tower RH [%]': print "Issue reading" + header [relhumloc] ; headererror = 1
                if header[windloc] != 'Avg Wind Speed @ 6ft [m/s]\n': print "Issue reading" + header [windloc] ; headererror = 1
            
                headeracquired = 1
                
                if headererror == 1:
                    print "STOPPING File Read because of headers issue"
                    continue
                
            else:
                
                if headererror == 1:
                    continue
            
                currentline=line.split(",")                
                month, day, year = currentline[dateloc].split("/")
                month = int(month); day = int(day); year = int(year)
                hour, minute = currentline[timeloc].split(":")
                hour = int(hour); minute = int(minute)
                wind, trash = currentline[windloc].split("\n")

                if decimate == True:
                    if decimatestyle == "Skip":
                        if minute%decimatevalue != 0:
                              continue

                date_all.append(currentline[dateloc])
                time_all.append(currentline[timeloc])
                global_all.append(float(currentline[globalloc]))
                direct_all.append(float(currentline[directloc]))
                diffuse_all.append(float(currentline[diffuseloc]))
                temp_all.append(float(currentline[temploc]))
                relhum_all.append(float(currentline[relhumloc]))
                day_all.append(day)
                month_all.append(month)
                year_all.append(year)
                hour_all.append(hour)
                minute_all.append(minute)
                wind_all.append(float(wind))
                myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)                
                timestamps_all.append(myTimestamp)
                
    SRRL = ({'Month': month_all, 'Day': day_all, 'Year': year_all, 'Hour': hour_all, 'Minute': minute_all,
                 'Global': global_all, 'Direct': direct_all, 'Diffuse': diffuse_all, 
                 'Temp': temp_all, 'Humidity': relhum_all, 'Wind': wind_all})

    SRRL = pd.DataFrame.from_records(SRRL, index=timestamps_all)


    if decimate == True:
        if decimatestyle == "Mean":
            if decimatevalue == 5:
                    SRRL=SRRL.resample('5Min', base=0).mean()
                    print "Data decimated to 5 min Interval by Average"
    
            if decimatevalue == 10:
                    SRRL=SRRL.resample('10Min').mean()
                    print "Data decimated to 10 min Interval by Average"
    
            if decimatevalue == 15:
                    SRRL=SRRL.resample('15Min').mean()
                    print "Data decimated to 15 min Interval by Average"
    
            if decimatevalue == 60:
                    SRRL=SRRL.resample('60Min').mean()
                    print "Data decimated to 1 Hr Interval by Average"

    return SRRL;                 

def readpvdaqSNL(pvdaq_SNL_Filetitle, decimate= False, decimatestyle = "Mean" , decimatevalue = 5 ):      

    date_all=[]; time_all=[]
    temp_all=[]; relhum_all=[]; wind_all=[]
    day_all = []; month_all = []; year_all = []; hour_all = []; minute_all= []
    poa_irradiation_all = []; system_id_all = []; reference_yield_all = [];
    timestamps_all = [];
    
    headeracquired= 0
    headererror = 0
    
    measdatetime_loc = 0
    availability_loc = 5
    poa_irradiation_loc = 14
    reference_yield_loc = 15
    system_id_loc = 16
    
    with open(pvdaq_SNL_Filetitle, "r") as filestream:
        
        print "Reading pvDAQ SNL File: ", pvdaq_SNL_Filetitle
        print "Decimation options: ", decimate
        if decimate == True:
            print "Decimate Style: ", decimatestyle, " to interval ", decimatevalue, " minutes."

        for line in filestream:
            if headeracquired == 0:
                header = line.split(",")
                        
                if header[measdatetime_loc] != 'measdatetime': print "Issue reading" + header [measdatetime_loc] ; headererror = 1
      #          if header[timeloc] != 'MST': print "Issue reading" + header [timeloc] ; headererror = 1
                if header[availability_loc] != 'availability': print "Issue reading" + header [availability_loc] ; headererror = 1
                if header[poa_irradiation_loc] != 'poa_irradiation': print "Issue reading" + header [poa_irradiation_loc] ; headererror = 1
                if header[reference_yield_loc] !=  'reference_yield': print "Issue reading" + header [reference_yield_loc] ; headererror = 1
                if header[system_id_loc] != 'system_id': print "Issue reading" + header [system_id_loc] ; headererror = 1 
                
                headeracquired = 1
    
                if headererror == 1:
                    print "STOPPING File Read because of headers issue"
                    continue
            else:
    
                if headererror == 1:
                    continue
    
                currentline=line.split(",")                
    
                # Deleting datapoints where availability was not 1.
                # It's only like 6 points for 2017~
                if currentline[availability_loc] != '1':
                    continue
                
                # Ignoring Empty Values now
                if currentline[poa_irradiation_loc] == '':
                    continue
                else:                
                    measdatetime = currentline[measdatetime_loc]
                    date, mytime, timezone = measdatetime.split(" ")
                    year, month, day = date.split("-")
                    month = int(month); day = int(day); year = int(year)
                    hour, minute, second = mytime.split(":")
                    hour = int(hour); minute = int(minute); second = int(second)
    
                    if decimate == True:
                        if decimatestyle == "Skip":
                            if minute%decimatevalue != 0:
                                  continue
        
                    date_all.append(date)
                    time_all.append(mytime)
                    poa_irradiation_all.append(currentline[poa_irradiation_loc])
                    system_id_all.append(currentline[system_id_loc])
                    reference_yield_all.append(currentline[reference_yield_loc])
                    day_all.append(day)
                    month_all.append(month)
                    year_all.append(year)
                    hour_all.append(hour)
                    minute_all.append(minute)
                    wind_all.append(0)
                    temp_all.append(0)
                    relhum_all.append(0)
                    myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)                
                    timestamps_all.append(myTimestamp)
                
    
    pvDaq_SNL = ({'Month': month_all, 'Day': day_all, 'Year': year_all, 'Hour': hour_all, 'Minute': minute_all,
                 'System_Id': system_id_all, 'POA_Irradiation': poa_irradiation_all, 'Reference_Yield': reference_yield_all, 
                 'Temp': temp_all, 'Humidity': relhum_all, 'Wind': wind_all})
    
    pvDaq_SNL = pd.DataFrame.from_records(pvDaq_SNL, index=timestamps_all)
    
    # Decimation by Average
    if decimate == True:
        if decimatestyle == "Mean":
            if decimatevalue == 5:
                    pvDaq_SNL=pvDaq_SNL.resample('5Min', base=0).mean()
                    print "Data decimated to 5 min Interval by Average"
        
            if decimatevalue == 10:
                    pvDaq_SNL=pvDaq_SNL.resample('10Min').mean()
                    print "Data decimated to 10 min Interval by Average"
        
            if decimatevalue == 15:
                    pvDaq_SNL=pvDaq_SNL.resample('15Min').mean()
                    print "Data decimated to 15 min Interval by Average"
        
            if decimatevalue == 60:
                    pvDaq_SNL=pvDaq_SNL.resample('60Min').mean()
                    print "Data decimated to 1 Hr Interval by Average"

    pvDaq_SNL.dropna()

    return pvDaq_SNL;


def readpvdaqSNL_CORRECT(pvdaq_SNL_Filetitle, decimate= False, decimatestyle = "Mean" , decimatevalue = 5 ):      

    date_all=[]; time_all=[]
    direct_all=[]; diffuse_all=[]
    temp_all=[]; relhum_all=[]; wind_all=[]
    day_all = []; month_all = []; year_all = []; hour_all = []; minute_all= []
    timestamps_all = [];
    
    
    headeracquired= 0
    headererror = 0
    
    measdatetime_loc = 1
    directloc = 5
    diffuseloc = 4
    temploc = 8
    relhumloc = 10
    windloc = 17
    
    with open(pvdaq_SNL_Filetitle, "r") as filestream:
        
        print "Reading pvDAQ SNL File: ", pvdaq_SNL_Filetitle
        print "Decimation options: ", decimate
        if decimate == True:
            print "Decimate Style: ", decimatestyle, " to interval ", decimatevalue, " minutes."
    
        rl = 0
        
        for line in filestream:
            rl = rl+1
            if headeracquired == 0:
                header = line.split(",")
                        
                if header[measdatetime_loc] != 'Date-Time': print "1 Issue reading" + header [measdatetime_loc] ; headererror = 1
      #          if header[timeloc] != 'MST': print "Issue reading" + header [timeloc] ; headererror = 1
                if header[directloc] != 'Direct_Wm2': print "2 Issue reading" + header [directloc] ; headererror = 1
                if header[diffuseloc] != 'Diffuse_Wm2': print "3 Issue reading" + header [diffuseloc] ; headererror = 1
                if header[temploc] != 'Panel_Temp_C': print "4 Issue reading" + header [temploc] ; headererror = 1
                if header[relhumloc] != 'RH_pct': print "5 Issue reading" + header [relhumloc] ; headererror = 1
                if header[windloc] != 'WS_ms_Mean': print "6 Issue reading" + header [windloc] ; headererror = 1
    
                headeracquired = 1
    
                if headererror == 1:
                    print "STOPPING File Read because of headers issue"
                    continue
            else:
    
                if headererror == 1:
                    continue
    
                currentline=line.split(",")                
    
                # Deleting datapoints where availability was not 1.
                # It's only like 6 points for 2017~
                
                # Ignoring Empty Values now
                measdatetime = currentline[measdatetime_loc]
                date, mytime = measdatetime.split(" ")
                year, month, day = date.split("-")
                month = int(month); day = int(day); year = int(year)
                hour, minute, second = mytime.split(":")
                hour = int(hour); minute = int(minute); second = int(second)
    
                if decimate == True:
                    if decimatestyle == "Skip":
                        if minute%decimatevalue != 0:
                              continue
    
                date_all.append(date)
                time_all.append(mytime)

                direct = 0.0
                diffuse = 0.0
                wind = 0.0
                relhum = 0.0
                temp = 0.0
                
                if (currentline[directloc]) == '':
                    print "Empty Direct at line ", rl
                else:                    
                    if float(currentline[directloc]) >= 0:
                        direct = float(currentline[directloc])
                
                if (currentline[diffuseloc]) == '':
                    print "Empty Diffuse at line ", rl

                else:
                    if float(currentline[diffuseloc]) >= 0:
                        diffuse = float(currentline[diffuseloc])

                if (currentline[relhumloc]) == '':
                    print "Empty Rel. Humidity at line ", rl
                else:                    
                    relhum = float(currentline[relhumloc])

                if (currentline[windloc]) == '':
                    print "Empty Wind at line ", rl
                else:                    
                    wind = float(currentline[windloc])
                
                if (currentline[temploc]) == '':
                    print "Empty Temperature at line ", rl
                else:                    
                    temp = float(currentline[temploc])
                                

                direct_all.append(direct)
                diffuse_all.append(diffuse)
                temp_all.append(temp)
                relhum_all.append(relhum)
                wind_all.append(wind)
                day_all.append(day)
                month_all.append(month)
                year_all.append(year)
                hour_all.append(hour)
                minute_all.append(minute)
                myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)                
                timestamps_all.append(myTimestamp)
                
    
    pvDaq_SNL = ({'Month': month_all, 'Day': day_all, 'Year': year_all, 'Hour': hour_all, 'Minute': minute_all,
                 'Direct': direct_all, 'Diffuse': diffuse_all, 
                 'Temp': temp_all, 'Humidity': relhum_all, 'Wind': wind_all})
    
    pvDaq_SNL = pd.DataFrame.from_records(pvDaq_SNL, index=timestamps_all)
    
    # Decimation by Average
    if decimate == True:
        if decimatestyle == "Mean":
            if decimatevalue == 5:
                    pvDaq_SNL=pvDaq_SNL.resample('5Min', base=0).mean()
                    print "Data decimated to 5 min Interval by Average"
        
            if decimatevalue == 10:
                    pvDaq_SNL=pvDaq_SNL.resample('10Min').mean()
                    print "Data decimated to 10 min Interval by Average"
        
            if decimatevalue == 15:
                    pvDaq_SNL=pvDaq_SNL.resample('15Min').mean()
                    print "Data decimated to 15 min Interval by Average"
        
            if decimatevalue == 60:
                    pvDaq_SNL=pvDaq_SNL.resample('60Min').mean()
                    print "Data decimated to 1 Hr Interval by Average"
    
    pvDaq_SNL.dropna()

    return pvDaq_SNL;


def IrradfromPOA(timestamps_all, POA_Measured_all, poa_sazm=180.0, poa_tilt=35.0, lat= 35.0549, lng=-106.5433, 
                 elevation = 1658, tzstr='US/Mountain', name=  "[1429] RTC,SNL,Baseline", albedo = 0.4, transmittance = 0.75):
    
    scalingfactor_all=[]               
    dni_scaled_all=[]
    dhi_scaled_all=[]
    timestamp_all=[]
    poa2_all=[]
    
    MYLOC=pvlib.location.Location(lat, lng, tzstr, elevation, name)
    
    rowsno = len(timestamps_all)
    
    for rl in range(0, rowsno):   
        if float(POA_Measured_all[rl]) <= 0:
            scalingfactor_all.append(0)
            dni_scaled_all.append(0)
            dhi_scaled_all.append(0)
            timestamp_all.append(0)
            poa2_all.append(0)
        else:
            times_loc= timestamps_all[rl].replace(tzinfo=pytz.timezone('US/Mountain'))
            solpos = pvlib.solarposition.get_solarposition(times_loc, lat, lng)
            AM=MYLOC.get_airmass(times_loc, solpos, model='kastenyoung1989')        
            aazir=math.radians(float(solpos['azimuth'])); azenr=math.radians(float(solpos['apparent_zenith']))
            irradiance=pvlib.irradiance.liujordan(float(solpos['apparent_zenith']), transmittance, float(AM['airmass_relative']), pressure=101325.0, dni_extra=1367.0)
            dni=float(irradiance['dni']); dhi=float(irradiance['dhi']); ghi=float(irradiance['ghi'])
    
            # pvDAQ sensor 1429 is located in SNL at a tilt of 35 deg.:        
            inc, tiltr, sazmr = sunIncident(0, poa_tilt, poa_sazm, 180.0, azenr, aazir)	    # For calling PerezComp to break diffuse into components for 
            poa, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, inc, tiltr, azenr)   # Call to get components for the tilt
            
            # Calculate scaling factor
            scalingfactor = 0
            if poa > 0:
                scalingfactor = float(POA_Measured_all[rl])/((poa)*0.001)               
            if scalingfactor > 1.3:
                scalingfactor = 1.3
                
            #Sanity Check
            dni_scaled=dni*scalingfactor; dhi_scaled=dhi*scalingfactor
            poa2, iso_dif2, circ_dif2, horiz_dif2, grd_dif2, beam2 = perezComp(dni_scaled, dhi_scaled, albedo, inc, tiltr, azenr)   # Call to get components for the tilt
            
            # Save params
            scalingfactor_all.append(scalingfactor)
            dni_scaled_all.append(dni_scaled)
            dhi_scaled_all.append(dhi_scaled)
            poa2_all.append(poa2)
            timestamp_all.append(times_loc)
            
    #Append to pvDaq.row
    
    IRRAD = ({'Direct': dni_scaled_all, 'Diffuse': dhi_scaled_all, 'POA_Calculated': poa2_all, 'ScalingFactor': scalingfactor_all})
    IRRAD = pd.DataFrame.from_records(IRRAD, index=timestamps_all)
    IRRAD = IRRAD.fillna(0)
    
    return IRRAD;

def simulateSRRL(SRRLtoread, writefiletitle,  beta = 0, sazm = 180, C = 0.5, D = None,
             rowType = 'interior', transFactor = 0.01, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = 0.2,  
             tracking = False, backtrack = True, rtr = None, Cv = None, offset = 0, 
             decimate = True, decimatestyle = "Mean" , decimatevalue = 5,      
             lat = 39.742, lng = -105.179, tz = -7, name = 'NRELGolden', dataInterval= 1,
             datestart = '2000-11-2 00:00:00', dateend = '3000-11-24 23:59:59'):


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
        
        ## Read SRRL data and start loop ~  
        SRRL = readSRRL(SRRLtoread, decimate= decimate, decimatestyle = decimatestyle , decimatevalue = decimatevalue)
#        in_range_SRRL = SRRL[SRRL.isin(pd.date_range(datestart, dateend))]
        mask = (SRRL.index > datestart) & (SRRL.index <= dateend )
        SRRL=SRRL.loc[mask]
        SRRL=SRRL.dropna()

        noRows, noCols = SRRL.shape
        
        ## Distance between rows for no shading on Dec 21 at 9 am
        print " "
        print "********* "
        print "Running Simulation for SRRL DataFile: ", SRRLtoread
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
          
                dni = SRRL['Direct'][rl]
                dhi = SRRL['Diffuse'][rl]
                Tamb = SRRL['Temp'][rl]
                Vwind = SRRL['Wind'][rl]
                #tzinfo='pytz.FixedOffset(-300)'
                #myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)e, 0, 0)                
                myTimestamp= SRRL.index[rl]
                hour = myTimestamp.hour
                minute = myTimestamp.minute
                year = myTimestamp.year
                day = myTimestamp.day
                month = myTimestamp.month
                
                rl = rl+1   # increasing while count

                # decimation                
                            
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
                        [rearSkyConfigFactors, frontSkyConfigFactors, ffConfigFactors] = getSkyConfigurationFactors(rowType, beta, C, D);       ## Sky configuration factors are the same for all times, only based on geometry and row type
    
    
    
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
                    backGTI, aveGroundGHI = getBackSurfaceIrradiances(rowType, maxShadow, PVbackSurface, beta, sazm, dni, dhi, C, D, albedo, zen, azm, cellRows, pvBackSH, rearGroundGHI, frontGroundGHI, frontReflected, offset)
               
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
       
     
        print "Finished"
        
        return;
      
def simulatepvDAQ(pvDAQ, writefiletitle,  beta = 0, sazm = 180, C = 0.5, D = None,
             rowType = 'interior', transFactor = 0.01, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = 0.2,  
             tracking = True, backtrack = True, rtr = None, Cv = None, offset = 0, 
             lat = 35.0549, lng = -106.5433, tz = -7, name = "[1429] RTC,SNL,Baseline", dataInterval= 1):

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
        
        ## Read SRRL data and start loop ~  
        noRows, noCols = pvDAQ.shape
        
        
        ## Distance between rows for no shading on Dec 21 at 9 am
        print " "
        print "********* "
        print "Running Simulation for pvDAQ DataFile"
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
          
                dni = pvDAQ['Direct'][rl]
                dhi = pvDAQ['Diffuse'][rl]
                Tamb = pvDAQ['Temp'][rl]
                Vwind = pvDAQ['Wind'][rl]
                #tzinfo='pytz.FixedOffset(-300)'
                #myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)e, 0, 0)                
                myTimestamp= pvDAQ.index[rl]
                hour = myTimestamp.hour
                minute = myTimestamp.minute
                year = myTimestamp.year
                day = myTimestamp.day
                month = myTimestamp.month
                
                rl = rl+1   # increasing while count

                # decimation                
                            
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


                        
                        gcr=1.0/rtr  
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
                        [rearSkyConfigFactors, frontSkyConfigFactors, ffConfigFactors] = getSkyConfigurationFactors(rowType, beta, C, D);       ## Sky configuration factors are the same for all times, only based on geometry and row type
    
    
    
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
                    backGTI, aveGroundGHI = getBackSurfaceIrradiances(rowType, maxShadow, PVbackSurface, beta, sazm, dni, dhi, C, D, albedo, zen, azm, cellRows, pvBackSH, rearGroundGHI, frontGroundGHI, frontReflected, offset)
               
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
       
     
        print "Finished"
        
        return;
    
def readOTFMeasuredFile(OTF_Filetitle, decimate= True, decimatestyle = "Mean" , decimatevalue = 5 ):      
    
#    OTF_Filetitle= "data\COMPLETE_OTF_Roof_Irrad.csv"
#    decimate= False
#    decimatestyle = "Mean"
#    decimatevalue = 5
    
    # decimatestyle = "Mean" or "Skip"
    day_all = []; month_all = []; year_all = []; hour_all = []; minute_all= []
    timestamps_all = [];
    CM11_1all = []; CM11_2all = []
    IMT1_all = []; IMT2_all = []; IMT3_all = []; IMT4_all = []; IMT5_all = []; IMT6_all = []

    
    headeracquired= 0
    headererror = 0

    timestamploc = 0    
    CM11_1loc = 2
    CM11_2loc = 3
    IMT1loc = 4
    IMT2loc = 5
    IMT3loc = 6
    IMT4loc = 7
    IMT5loc = 8
    IMT6loc= 9
    
    with open(OTF_Filetitle, "r") as filestream:
    
        print "Reading OTF File: ", OTF_Filetitle
        print "Decimation options: ", decimate
        if decimate == True:
            print "Decimate Style: ", decimatestyle, " to interval ", decimatevalue, " minutes."
            
        for line in filestream:
            if headeracquired == 0:
                header = line.split(",")
                        
                if header[timestamploc] != 'TIMESTAMP': print "Issue reading" + header [timestamploc] ; headererror = 1
                if header[CM11_1loc] != 'CM11_1': print "Issue reading" + header [CM11_1loc] ; headererror = 1
                if header[IMT1loc] != 'IMT1': print "Issue reading" + header [IMT1loc] ; headererror = 1

                headeracquired = 1
                
                if headererror == 1:
                    print "STOPPING File Read because of headers issue"
                    continue
                
            else:
                
                if headererror == 1:
                    continue
            
                currentline=line.split(",")                
                date, hour = currentline[timestamploc].split(" ")
                month, day, year = date.split("/")
                month = int(month); day = int(day); year = int(year)
                hour, minute = hour.split(":")
                hour = int(hour); minute = int(minute)
                
                if decimate == True:
                    if decimatestyle == "Skip":
                        if minute%decimatevalue != 0:
                              continue

                
                CM11_1all.append(float(currentline[CM11_1loc]))
                CM11_2all.append(float(currentline[CM11_2loc]))
                IMT1_all.append(float(currentline[IMT1loc]))
                IMT2_all.append(float(currentline[IMT2loc]))
                IMT3_all.append(float(currentline[IMT3loc]))
                IMT4_all.append(float(currentline[IMT4loc]))
                IMT5_all.append(float(currentline[IMT5loc]))
                IMT6_all.append(float(currentline[IMT6loc]))
                day_all.append(day)
                month_all.append(month)
                year_all.append(year)
                hour_all.append(hour)
                minute_all.append(minute)
                myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)                
                timestamps_all.append(myTimestamp)
                
    OTF = ({'Month': month_all, 'Day': day_all, 'Year': year_all, 'Hour': hour_all, 'Minute': minute_all,
                 'IMT1': IMT1_all, 'IMT2': IMT2_all, 'IMT3': IMT3_all, 'IMT4': IMT4_all, 'IMT5': IMT5_all, 'IMT6': IMT6_all,
                 'CM11_1': CM11_1all, 'CM11_2': CM11_2all})

    OTF = pd.DataFrame.from_records(OTF, index=timestamps_all)


    if decimate == True:
        if decimatestyle == "Mean":
            if decimatevalue == 5:
                    OTF=OTF.resample('5Min', base=0).mean()
                    print "Data decimated to 5 min Interval by Average"
    
            if decimatevalue == 10:
                    OTF=OTF.resample('10Min').mean()
                    print "Data decimated to 10 min Interval by Average"
    
            if decimatevalue == 15:
                    OTF=OTF.resample('15Min').mean()
                    print "Data decimated to 15 min Interval by Average"
    
            if decimatevalue == 60:
                    OTF=OTF.resample('60Min').mean()
                    print "Data decimated to 1 Hr Interval by Average"

    return OTF;                 


def readSetupResult(SetupResult_Filetitle, decimate= False, decimatestyle = "Mean" , decimatevalue = 5 ):      

    
    # decimatestyle = "Mean" or "Skip"
    day_all = []; month_all = []; year_all = []; hour_all = []; minute_all= []
    myTimestamp_all = [];
    IMT1_all = []; IMT2_all = []; IMT3_all = []; IMT4_all = []; IMT5_all = []; IMT6_all = []
    GTIFrontavg_all=[]; GTIBackavg_all=[]
    
    headeracquired= 0
    headererror = 0
    
    Yearloc = 1
    Monthloc = 2
    Dayloc = 3
    Hourloc = 4
    Minuteloc = 5
    DNIloc = 6
    DHIloc = 7
    ghiloc = 9
    IMT1loc = 22
    IMT2loc = 23
    IMT3loc = 24
    IMT4loc = 25
    IMT5loc = 26
    IMT6loc = 27
    GTIFrontavgloc = 28
    GTIBackavgloc = 29
        
    with open(SetupResult_Filetitle, "r") as filestream:
    
        print "Reading VF Setup Result File: ", SetupResult_Filetitle
        print "Decimation options: ", decimate
        if decimate == True:
            print "Decimate Style: ", decimatestyle, " to interval ", decimatevalue, " minutes."
            
        for line in filestream:
            if headeracquired == 0:
                header = line.split(",")
                        
                if header[Yearloc] != 'Year': print "Issue reading" + header [Yearloc] ; headererror = 1
                if header[IMT1loc] != 'IMT1': print "Issue reading" + header [IMT1loc] ; headererror = 1
    
                headeracquired = 1
                
                if headererror == 1:
                    print "STOPPING File Read because of headers issue"
                    continue
                
            else:
                
                if headererror == 1:
                    continue
            
                currentline=line.split(",")                
                
                year = int(currentline[Yearloc])
                month = int(currentline[Monthloc])
                day = int(currentline[Dayloc])
                hour = int(currentline[Hourloc])
                minute = int(currentline[Minuteloc])
    
                if decimate == True:
                    if decimatestyle == "Skip":
                        if minute%decimatevalue != 0:
                              continue
                
                IMT1_all.append(float(currentline[IMT1loc]))
                IMT2_all.append(float(currentline[IMT2loc]))
                IMT3_all.append(float(currentline[IMT3loc]))
                IMT4_all.append(float(currentline[IMT4loc]))
                IMT5_all.append(float(currentline[IMT5loc]))
                IMT6_all.append(float(currentline[IMT6loc]))
                GTIFrontavg_all.append(float(currentline[GTIFrontavgloc]))
                GTIBackavg_all.append(float(currentline[GTIBackavgloc]))
                day_all.append(day)
                month_all.append(month)
                year_all.append(year)
                hour_all.append(hour)
                minute_all.append(minute)
                myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)                
                myTimestamp_all.append(myTimestamp)
                
    SetupResults = ({'Month_M': month_all, 'Day_M': day_all, 'Year_M': year_all, 'Hour_M': hour_all, 'Minute_M': minute_all,
                 'IMT1_M': IMT1_all, 'IMT2_M': IMT2_all, 'IMT3_M': IMT3_all, 'IMT4_M': IMT4_all, 'IMT5_M': IMT5_all, 'IMT6_M': IMT6_all,
                 'GTIFrontavg_M': GTIFrontavg_all, 'GTIBackavg_M': GTIBackavg_all})
    
    SetupResults = pd.DataFrame.from_records(SetupResults, index=myTimestamp_all)
    
    
    if decimate == True:
        if decimatestyle == "Mean":
            if decimatevalue == 5:
                    SetupResults=SetupResults.resample('5Min', base=0).mean()
                    print "Data decimated to 5 min Interval by Average"
    
            if decimatevalue == 10:
                    SetupResults=SetupResults.resample('10Min').mean()
                    print "Data decimated to 10 min Interval by Average"
    
            if decimatevalue == 15:
                    SetupResults=SetupResults.resample('15Min').mean()
                    print "Data decimated to 15 min Interval by Average"
    
            if decimatevalue == 60:
                    SetupResults=SetupResults.resample('60Min').mean()
                    print "Data decimated to 1 Hr Interval by Average"

    return SetupResults;   

def readRadSetupResult(RadSetupResult_Filetitle, decimate= False, decimatestyle = "Mean" , decimatevalue = 5 ):      

    
    # decimatestyle = "Mean" or "Skip"
    day_all = []; month_all = []; year_all = []; hour_all = []; minute_all= []
    myTimestamp_all = [];
    IMT1_all = []; IMT2_all = []; IMT3_all = []; IMT4_all = []; IMT5_all = []; IMT6_all = []
   
    headeracquired= 0
    headererror = 0
    
    Yearloc = 0
    Monthloc = 1
    Dayloc = 2
    Hourloc = 3
    Minuteloc = 4
    DNIloc = 5
    DHIloc = 6
    Tambloc= 7
    VWindloc= 8
    IMT1loc = 9
    IMT2loc = 10
    IMT3loc = 11
    IMT4loc = 12
    IMT5loc = 13
    IMT6loc = 14
        
    with open(RadSetupResult_Filetitle, "r") as filestream:
    
        print "Reading Rad Setup Result File: ", RadSetupResult_Filetitle
        print "Decimation options: ", decimate
        if decimate == True:
            print "Decimate Style: ", decimatestyle, " to interval ", decimatevalue, " minutes."
            
        for line in filestream:
            if headeracquired == 0:
                header = line.split(",")
                        
                if header[Yearloc] != 'Year': print "Issue reading" + header [Yearloc] ; headererror = 1
                if header[IMT1loc] != 'IMT1_Rad': print "Issue reading" + header [IMT1loc] ; headererror = 1
    
                headeracquired = 1
                
                if headererror == 1:
                    print "STOPPING File Read because of headers issue"
                    continue
                
            else:
                
                if headererror == 1:
                    continue
            
                currentline=line.split(",")                
                
                year = int(currentline[Yearloc])
                month = int(currentline[Monthloc])
                day = int(currentline[Dayloc])
                hour = int(currentline[Hourloc])
                minute = int(currentline[Minuteloc])-2
    
                if decimate == True:
                    if decimatestyle == "Skip":
                        if minute%decimatevalue != 0:
                              continue
                
                IMT1_all.append(float(currentline[IMT1loc]))
                IMT2_all.append(float(currentline[IMT2loc]))
                IMT3_all.append(float(currentline[IMT3loc]))
                IMT4_all.append(float(currentline[IMT4loc]))
                IMT5_all.append(float(currentline[IMT5loc]))
                IMT6_all.append(float(currentline[IMT6loc]))
                day_all.append(day)
                month_all.append(month)
                year_all.append(year)
                hour_all.append(hour)
                minute_all.append(minute)
                myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)                
                myTimestamp_all.append(myTimestamp)
                
    SetupResults = ({'Month_R': month_all, 'Day_R': day_all, 'Year_R': year_all, 'Hour_R': hour_all, 'Minute_R': minute_all,
                 'IMT1_Rad': IMT1_all, 'IMT2_Rad': IMT2_all, 'IMT3_Rad': IMT3_all, 'IMT4_Rad': IMT4_all, 'IMT5_Rad': IMT5_all, 'IMT6_Rad': IMT6_all})
    
    SetupResults = pd.DataFrame.from_records(SetupResults, index=myTimestamp_all)
    
    
    if decimate == True:
        if decimatestyle == "Mean":
            if decimatevalue == 5:
                    SetupResults=SetupResults.resample('5Min', base=0).mean()
                    print "Data decimated to 5 min Interval by Average"
    
            if decimatevalue == 10:
                    SetupResults=SetupResults.resample('10Min').mean()
                    print "Data decimated to 10 min Interval by Average"
    
            if decimatevalue == 15:
                    SetupResults=SetupResults.resample('15Min').mean()
                    print "Data decimated to 15 min Interval by Average"
    
            if decimatevalue == 60:
                    SetupResults=SetupResults.resample('60Min').mean()
                    print "Data decimated to 1 Hr Interval by Average"

    return SetupResults;   