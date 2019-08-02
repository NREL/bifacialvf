#!/usr/bin/env python2
# -*- coding: utf-8 -*-
        #          This program calculates irradiances on the front and back surfaces of bifacial PV modules.
        #          Key dimensions and nomenclature:
        #          beta = PV module tilt angle from horizontal, in degrees
        #          sazm = PV module surface azimuth from north, in degrees
        #          1.0 = normalized PV module/panel slant height
        #          C = ground clearance of PV module, in PV module/panel slant heights
        #          D = distance between rows, from rear of module to front of module in next row, in PV module/panel slant heights
        #          h = sin(beta), vertical PV module dimension, in PV module/panel slant heights
        #          x1 = cos(beta), horizontal PV module dimension, in PV module/panel slant heights
        #          rtr = x1 + D, row-to-row distance, from front of module to front of module in next row, in PV module/panel slant heights
        #          cellRows = number of horzontal rows of cells in a PV module/panel
        #          PVfrontSurface = PV module front surface material type, either "glass" or "ARglass"
        #          PVbackSurface = PV module back surfac ematerial type, either "glass" or "ARglass"
        #        
        #         Program flow consists of:
        #          a. Calculate irradiance distribution on ground
        #          b. Calculate AOI corrected irradiance on front of PV module, and irradiance reflected from front of PV module
        #          c. Calculate irradiance on back of PV module

# ensure python3 compatible division and printing
from __future__ import division, print_function, absolute_import
 
import math
import csv
import pvlib
import os
import sys

from .vf import getBackSurfaceIrradiances, getFrontSurfaceIrradiances, getGroundShadeFactors
from .vf import getSkyConfigurationFactors, trackingBFvaluescalculator, rowSpacing
from .sun import hrSolarPos, perezComp, solarPos, sunIncident
import pandas as pd
from readepw import readepw

# Electrical Mismatch Calculation 
import scipy.io as sio
sys.path.insert(0, 'BF_BifacialIrradiances')
from PortraitSingleHour import PortraitSingleHour    # For calculateBilInterpol
from LandscapeSingleHour import LandscapeSingleHour # For calculateBilInterpol
from pvmismatch import *  # this imports everything we need
import numpy as np

def simulate(TMYtoread=None, writefiletitle=None, beta=0, sazm=180, C=0.5, D=None,
             rowType='interior', transFactor=0.01, cellRows=6, 
             PVfrontSurface='glass', PVbackSurface='glass', albedo=0.2,  
             tracking=False, backtrack=True, rtr=None, max_angle=45,
             calculatePVMismatch=False, portraitorlandscape='landscape',
             calculateBilInterpol=False, interpolA=0.005, IVArray=None, 
             beta_voc_all=None, m_all=None, bee_all=None):
                       
        '''
      
        Description
        -----------
        Main function to run the bifacialvf routines 
    
        Parameters
        ----------
        TMYtoread: TMY3 .csv weather file, which can be downloaded at http://rredc.nrel.gov/solar/old_data/nsrdb/1991-2005/tmy3/by_state_and_city.html
                   Also .epw weather files, which can be downloaded here: https://energyplus.net/weather and here: http://re.jrc.ec.europa.eu/pvg_tools/en/tools.html#TMY
                 
        writefiletitle:  name of output file
        beta:    tilt angle in degrees.  Not used for tracking
        sazm:    surface azimuth orientation in degrees east of north. For tracking this is the tracker axis orientation
        C:       normalized ground clearance.  For trackers, this is the module height at zero tilt
        D:       normalized gap between PV module rows.  For trackers use rtr
        rtr:     row-to-row normalized distance.  = 1/GCR
        transFactor:   PV module transmission fraction.  Default 1% (0.01)
        cellRows:      Number of points along the module chord to return irradiance values.  Default 6 (1-up landscape module)
        max_angle:     1-axis tracking maximum limits of rotation
        tracking, backtrack:  boolean to enable 1-axis tracking and pvlib backtracking algorithm, respectively
        
        
        Returns
        -------
        none
        '''    
        
        if tracking == True:
            axis_tilt = 0  # algorithm only allows for zero north-south tilt with SAT
            #max_angle = 45  # maximum tracker rotation 
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
        if writefiletitle == None:
            writefiletitle = "data/Output/TEST.csv"
        
        ## Read TMY3 data and start loop ~  
        if TMYtoread is None: # if no file passed in, the readtmy3 graphical file picker will open.
            (myTMY3,meta)=pvlib.tmy.readtmy3(TMYtoread)        
        elif TMYtoread.lower().endswith('.csv') :  
            (myTMY3,meta)=pvlib.tmy.readtmy3(TMYtoread)
        elif TMYtoread.lower().endswith('.epw') : 
            (myTMY3,meta) = readepw(TMYtoread)
            # rename different field parameters to match DNI, DHI, DryBulb, Wspd
            myTMY3.rename(columns={'Direct normal radiation in Wh/m2':'DNI','Diffuse horizontal radiation in Wh/m2':'DHI','Dry bulb temperature in C':'DryBulb','Wind speed in m/s':'Wspd'}, inplace=True)
        else:
            raise Exception('Incorrect extension for TMYtoread. Either .csv (TMY3) .epw or None')
            
        #myAxisTitles=myTMY3.axes
        noRows, noCols = myTMY3.shape
        lat = meta['latitude']; lng = meta['longitude']; tz = meta['TZ']
        name = meta['Name']
        
        ## infer the data frequency in minutes
        dataInterval = (myTMY3.index[1]-myTMY3.index[0]).total_seconds()/60
    
        ## Distance between rows for no shading on Dec 21 at 9 am
        print( " ")
        print( "********* ")
        print( "Running Simulation for TMY3: ", TMYtoread)
        print( "Location:  ", name)
        print( "Lat: ", lat, " Long: ", lng, " Tz ", tz)
        print( "Parameters: beta: ", beta, "  Sazm: ", sazm, "  Height: ", C, "  rtr separation: ", rtr, "  Row type: ", rowType, "  Albedo: ", albedo)
        print( "Saving into", writefiletitle)
        print( " ")
        print( " ")
        
        DD = rowSpacing(beta, sazm, lat, lng, tz, 9, 0.0);          ## Distance between rows for no shading on Dec 21 at 9 am
        print( "Distance between rows for no shading on Dec 21 at 9 am solar time = ", DD)
        print( "Actual distance between rows = ", D  )
        print( " ")
    
        if tracking==False:        
            ## Sky configuration factors are the same for all times, only based on geometry and row type
            [rearSkyConfigFactors, frontSkyConfigFactors] = getSkyConfigurationFactors(rowType, beta, C, D);       ## Sky configuration factors are the same for all times, only based on geometry and row type
    
     
        ## Create WriteFile and write labels at this time
        
        #check that the save directory exists, unless it's in root
        savedirectory = os.path.dirname(writefiletitle)
        if ( (not os.path.exists(savedirectory)) and (savedirectory is not '')):
            os.makedirs(savedirectory)
        
        with open (writefiletitle,'w') as csvfile:
            sw = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
            # Write Simulation Parameters (from setup file)
            
            outputheader=['Latitude(deg)','Longitude(deg)', 'Time Zone','Tilt(deg)', 
                         'PV Azimuth(deg)','GroundClearance(panel slope lengths)', 'Row-to-Row-Distance rtr', 'RowType(first interior last single)',
                         'TransmissionFactor(open area fraction)','CellRows(# hor rows in panel)', 
                         'PVfrontSurface(glass or ARglass)', 'PVbackSurface(glass or ARglass)',
                         'Albedo',  'Tracking', 'CalculatePVOutput (Bilinear Interpol)','CalculatePVOutput (PVMismatch)']
            outputheadervars=[lat, lng, tz, beta, sazm, C, rtr, rowType, transFactor, cellRows, PVfrontSurface,
                             PVbackSurface, albedo, tracking, calculateBilInterpol, calculatePVMismatch]
            
            
            if tracking==True:
                outputheader+=['Backtracking']
                outputheadervars.append(backtrack)

            if calculateBilInterpol==True:  
                outputheader+=['interpolA']
                outputheadervars.append(interpolA)               
                outputheader+=['PortraitorLandscape']
                outputheadervars.append(portraitorlandscape)       
            
                cellCenterBI=[]
                if portraitorlandscape == 'portrait':
                    if cellRows != 6:
                        for i in range (0, 6):
                            cellCenterBI.append((i*cellRows/6.0+(i+1)*cellRows/6.0)/2)
                            
                if portraitorlandscape == 'landscape':
                    if cellRows != 6:
                        for i in range (0, 6):
                            cellCenterBI.append((i*cellRows/6.0+(i+1)*cellRows/6.0)/2)
                            
                            
            if calculatePVMismatch == True:
                outputheader+=['PortraitorLandscape']
                outputheadervars.append(portraitorlandscape)  
                
                stdpl=np.array([[0,	23,	24,	47,	48,	71,	72,	95],
                [1,	22,	25,	46,	49,	70,	73,	94],
                [2,	21,	26,	45,	50,	69,	74,	93],
                [3,	20,	27,	44,	51,	68,	75,	92],
                [4,	19,	28,	43,	52,	67,	76,	91],
                [5,	18,	29,	42,	53,	66,	77,	90],
                [6,	17,	30,	41,	54,	65,	78,	89],
                [7,	16,	31,	40,	55,	64,	79,	88],
                [8,	15,	32,	39,	56,	63,	80,	87],
                [9,	14,	33,	38,	57,	62,	81,	86],
                [10,	13,	34,	37,	58,	61,	82,	85],
                [11,	12,	35,	36,	59,	60,	83,	84]])

                cellCenterPVM=[]
                
                if portraitorlandscape == 'portrait':
                    if cellRows != 12:
                        for i in range (0, 12):
                            cellCenterPVM.append((i*cellRows/12.0+(i+1)*cellRows/12.0)/2)
                            
                if portraitorlandscape == 'landscape':
                    stdpl = stdpl.transpose()
                    if cellRows != 8:
                        for i in range (0, 8):
                            cellCenterPVM.append((i*cellRows/8.0+(i+1)*cellRows/8.0)/2)
                
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
                         'aveBackGroundGHI', 'GTIbackBroadBand', 'maxShadow', 'Tamb', 'VWind']
            outputtitles+=allrowfronts
            outputtitles+=allrowbacks
            if tracking == True:
                print( " ***** IMPORTANT --> THIS SIMULATION Has Tracking Activated")
                print( "Backtracking Option is set to: ", backtrack)
                outputtitles+=['beta']
                outputtitles+=['sazm']
                outputtitles+=['height']
                outputtitles+=['D']
                    
            if calculateBilInterpol==True:
                outputtitles+=['BilInterpol FRONT + BACK (Averaged) PmaxIdeal [W]']
                outputtitles+=['BilInterpol FRONT + BACK (Detailed) PmaxUnmatched [W]']

            if calculatePVMismatch == True:
                outputtitles+=['PVMismatch FRONT + BACK (Averaged) PmaxIdeal [W]']
                outputtitles+=['PVMismatch FRONT + BACK (Detailed) PmaxUnmatched [W]']
                
            sw.writerow(outputtitles)
            
            ## Loop through file.  TODO: replace this with for loop.
            rl = 0
            
            while (rl < noRows):
                
                myTimestamp=myTMY3.index[rl]
                year = myTimestamp.year
                month = myTimestamp.month
                day = myTimestamp.day
                hour = myTimestamp.hour
                minute = myTimestamp.minute
                dni = myTMY3.DNI[rl]#get_value(rl,5,"False")
                dhi = myTMY3.DHI[rl]#get_value(rl,8,"False")
                Tamb=myTMY3.DryBulb[rl]#get_value(rl,29,"False")
                VWind=myTMY3.Wspd[rl]#get_value(rl,44,"False")
                
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
                                  gtiAllpc, maxShadow, Tamb, VWind]
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

                    if calculateBilInterpol==True:
                        
                        if portraitorlandscape=='portrait':
                        
                            if cellRows != 6:                        
                                cellCenterValFront= np.interp(cellCenterBI, list(range(0,cellRows)), frontGTIrow)
                                cellCenterValBack= np.interp(cellCenterBI, list(range(0,cellRows)), backGTIrow)
                            else:
                                cellCenterValFront = frontGTIrow
                                cellCenterValBack = backGTIrow
                            
                            [PmaxIdeal, PmaxUnmatched, PmaxAvg] = PortraitSingleHour(cellCenterValFront, cellCenterValBack, Tamb, VWind, 6, interpolA,IVArray,beta_voc_all,m_all,bee_all)
                            outputvalues.append(PmaxIdeal)
                            outputvalues.append(PmaxUnmatched)

                        if portraitorlandscape=='landscape':
                            
                            if cellRows != 6:                        
                                cellCenterValFront= np.interp(cellCenterBI, list(range(0,cellRows)), frontGTIrow)
                                cellCenterValBack= np.interp(cellCenterBI, list(range(0,cellRows)), backGTIrow)

                            else:
                                cellCenterValFront = frontGTIrow
                                cellCenterValBack = backGTIrow
                                
                            [PmaxIdeal, PmaxUnmatched, PmaxAvg] = LandscapeSingleHour(cellCenterValFront, cellCenterValBack, Tamb, VWind, 6, interpolA,IVArray,beta_voc_all,m_all,bee_all)
                            outputvalues.append(PmaxIdeal)
                            outputvalues.append(PmaxUnmatched)


                    if calculatePVMismatch==True:

                        pvsys = pvsystem.PVsystem(numberStrs=1, numberMods=1)  # makes the system  # 1 module, in portrait mode. 
                        pmp_ideal=pvsys.Pmp   # Panel ideal. Monofacial.                            

                        if portraitorlandscape == 'portrait':                    
                        
                            if cellRows != 12:                        
                                cellCenterValFront= np.interp(cellCenterPVM, list(range(0,cellRows)), frontGTIrow)
                                cellCenterValBack= np.interp(cellCenterPVM, list(range(0,cellRows)), backGTIrow)
                            else:
                                cellCenterValFront = frontGTIrow
                                cellCenterValBack = backGTIrow
                                
                            sunmatDetailed=[]
                            sunmatAveraged=[]
                            
                            cellCenterValues_FrontPlusBack = cellCenterValFront+cellCenterValBack
                            AveFront=cellCenterValFront.mean()                
                            AveBack=cellCenterValBack.mean()
                                     
                            # Repeat to create a matrix to pass matrix.
                            for j in range (0, len(cellCenterValues_FrontPlusBack)):
                                sunmatDetailed.append([cellCenterValues_FrontPlusBack[j]/1000]*8)
                                
                            for j in range (0, len(cellCenterValFront)):
                                sunmatAveraged.append([(AveFront+AveBack)/1000]*8)
                            
                            # ACtually do calculations
                            pvsys.setSuns({0: {0: [sunmatAveraged, stdpl]}})
                            PowerAveraged=pvsys.Pmp
                            
                            pvsys.setSuns({0: {0: [sunmatDetailed, stdpl]}})
                            PowerDetailed=pvsys.Pmp

                            # Append Values
                            outputvalues.append(PowerAveraged)
                            outputvalues.append(PowerDetailed)
                            
                        if portraitorlandscape == 'landscape':  
                            
                            if cellRows != 8:                        
                                cellCenterValFront= np.interp(cellCenterPVM, list(range(0,cellRows)), frontGTIrow)
                                cellCenterValBack= np.interp(cellCenterPVM, list(range(0,cellRows)), backGTIrow)
                            else:
                                cellCenterValFront = frontGTIrow
                                cellCenterValBack = backGTIrow
   
                            
                            sunmatDetailed=[]
                            sunmatAveraged=[]
                            
                            cellCenterValues_FrontPlusBack = cellCenterValFront+cellCenterValBack
                            AveFront=cellCenterValFront.mean()                
                            AveBack=cellCenterValBack.mean()
                                     
                            # Repeat to create a matrix to pass matrix.
                            for j in range (0, len(cellCenterValues_FrontPlusBack)):
                                sunmatDetailed.append([cellCenterValues_FrontPlusBack[j]/1000]*12)
                                
                            for j in range (0, len(cellCenterValFront)):
                                sunmatAveraged.append([(AveFront+AveBack)/1000]*12)
                            
                            # ACtually do calculations
                            pvsys.setSuns({0: {0: [sunmatAveraged, stdpl]}})
                            PowerAveraged=pvsys.Pmp
                            
                            pvsys.setSuns({0: {0: [sunmatDetailed, stdpl]}})
                            PowerDetailed=pvsys.Pmp
                            

                            
                            # Append Values
                            outputvalues.append(PowerAveraged)     
                            outputvalues.append(PowerDetailed)
                                         
                    sw.writerow(outputvalues)
    
        	# End of daylight if loop 
    
        #strLine = sr.ReadLine();        # Read next line of data
       # End of while strLine != null loop
       
     
        print( "Finished")
        
        return;
        
if __name__ == "__main__":    
    #import time
    #start_time = time.time()


    beta = 10                   # PV tilt (deg)
    sazm = 180                  # PV Azimuth(deg) or tracker axis direction
    C = 1.0                      # GroundClearance(panel slope lengths). For tracking this is tilt = 0 hub height 
    D = 0.51519                 # DistanceBetweenRows(panel slope lengths)
    rowType = "interior"        # RowType(first interior last single)
    transFactor = 0.013         # TransmissionFactor(open area fraction)
    cellRows = 6                # CellRows(# hor rows in panel)   <--> THIS ASSUMES LANDSCAPE ORIENTATION 
    PVfrontSurface = "glass"    # PVfrontSurface(glass or ARglass)
    PVbackSurface = "glass"     # PVbackSurface(glass or ARglass)
    albedo = 0.62               # ground albedo
    calculateDetailedMismatch = False

     #BILINEAR INTERPOLATION VALUES    
    calculateBilInterpol = True
    calculatePVMismatch = True
    # PORTRAIT BILINEAR INTERPOLATION NOT WORKING ATM!!!!! 
    portraitorlandscape='landscape'   # portrait or landscape
    mat_contents = sio.loadmat(r'BF_BifacialIrradiances\BilinearInterpParams\IVArrayYingli.mat')
    IVArray=mat_contents['IVArray']        
    mat_contents = sio.loadmat(r'BF_BifacialIrradiances\BilinearInterpParams\newBilinearParamsYingLi.mat')
    beta_voc_all=mat_contents['beta_voc_all']        
    m_all=mat_contents['m_all']
    bee_all=mat_contents['bee_all']
    interpolA = 0.005  # More accurate interpolation. Do 0.01 as an option.

    
    # Tracking instructions
    tracking=False
    backtrack=True
    rtr = 1.5                   # row to row spacing in normalized panel lengths. 
    max_angle = 60

    TMYtoread="data/724010TYA.csv"   # VA Richmond
    writefiletitle="data/Output/RICHMOND_1.0.csv"

    simulate(TMYtoread, writefiletitle, beta, sazm, C, D, 
                     rowType, transFactor, cellRows, 
                     PVfrontSurface, PVbackSurface, albedo, 
                     tracking, backtrack, rtr, max_angle,
                     calculatePVMismatch, portraitorlandscape, 
                     calculateBilInterpol, interpolA, IVArray, beta_voc_all,
                     m_all, bee_all)
                
    #Load the results from the resultfile
    from loadVFresults import loadVFresults
    (data, metadata) = loadVFresults(writefiletitle)
    #print data.keys()
    # calculate average front and back global tilted irradiance across the module chord
    data['GTIFrontavg'] = data[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
    data['GTIBackavg'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)
    
    # Print the annual bifacial ratio.
    frontIrrSum = data['GTIFrontavg'].sum()
    backIrrSum = data['GTIBackavg'].sum()
    print('The bifacial ratio for ground clearance {} and rtr spacing {} is: {:.1f}%'.format(C,rtr,backIrrSum/frontIrrSum*100))
    #print("--- %s seconds ---" % (time.time() - start_time))
