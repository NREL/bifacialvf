#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 07:59:37 2017

@author: sayala
Modified 07.28 by Sara MacAlpine
Calls matlab routines from here rather than in the single hour code
Pulls in front and rear irradiance directly rather than using ratios as I did in matlab code
Interpolation just 0.005A for current (half as many steps) - I can now get ~1min run times on my laptop
Has numsens sensors, passed to the hourly code
Also passes Tamb and Vwind to the hourly calcs to get cell temperature
"""

import scipy.io as sio
import numpy as np
#Sara added this line because like to see entire lists when debugging
np.set_printoptions(threshold=np.inf)
#from pvmismatch import *  # this imports everything we need
import csv  
import sys 
import math
sys.path.insert(0, '../BF_BifacialIrradiances')
from loadVFresults import loadVFresults
from PortraitSingleHour import PortraitSingleHour
from LandscapeSingleHour import LandscapeSingleHour
from matplotlib import pyplot as plt  # now lets make some plots
#plt.ion()  # this turns on interactive plotting
import time

mat_contents = sio.loadmat('BilinearInterpParams//IVArrayYingli.mat')
IVArray=mat_contents['IVArray']
    
mat_contents = sio.loadmat('BilinearInterpParams//newBilinearParamsYingLi.mat')
#BilinearParams=mat_contents['BilinearParams']
#print mat_contents.keys()
    
#mat_contents = sio.loadmat('BilinearInterpParams//beta_voc_Yingli.mat')
beta_voc_all=mat_contents['beta_voc_all']
    
#mat_contents = sio.loadmat('BilinearInterpParams//m_Yingli.mat')
m_all=mat_contents['m_all']

#mat_contents = sio.loadmat('BilinearInterpParams//bee_Yingli.mat')
bee_all=mat_contents['bee_all']

#myfolder = '//Users//sayala//Desktop//Output//SolarWorld'   # Folder Path where I have the results
myfolder = '/Users/smacalpi/Documents/Python Scripts/Silvana_BifacialIrradiances/Silvana_BifacialIrradiances'   # Folder Path where I have the results

myf=(['landscape_25.csv'])#,
#'SolarWorldResults_Height_0.5_InteriorRow.csv',
#'SolarWorldResults_Height_1_InteriorRow.csv']) 


#for mainloop in range (0, len(myf)):
for mainloop in range (0, 1):
#for mainloop in range (18, len(myf)):
    myfall=myfolder+'/'+myf[mainloop]
    #data, meta = loadVFresults(myfall)   
    writefiletitle=myfolder+"testtest.csv"
    print "Working on ", myf[mainloop]
    
    with open(myfall,'r') as csvinput:
        with open(writefiletitle, 'w') as csvoutput:
            writer = csv.writer(csvoutput, lineterminator='\n')
            reader = csv.reader(csvinput)
            
            tic = time.clock()       

            all=[]
            row = next(reader)
            all.append(row)
            row = next(reader)
            all.append(row)
            row = next(reader)
            row.append('PmaxIdeal [W]')
            row.append('PmaxUnmatched [W]')
            row.append('PmaxAvg [W]')
            
            all.append(row)
            writer.writerows(all)
           
            ct=0;
            for row in reader:
                ct=ct+1;
                #print "ct...",ct
                #if ct<111:
                #    continue
                #if ct>111:
                #    break
                
                Tamb=float(row[19])
                Vwind=float(row[20])
                # frontGTI = data.loc[:,"No_1_RowFrontGTI":'No_9_RowFrontGTI']
                # backGTI = data.loc[:,"Noc_1_RowBackGTI":'No_9_RowBackGTI']
                # cellCenterValFront= np.interp(cellCenter, [0, 1, 2, 3, 4, 5, 6, 7, 8], frontGTI.loc[i,"No_1_RowFrontGTI":'No_9_RowFrontGTI'])
                # cellCenterValBack= np.interp(cellCenter, [0, 1, 2, 3, 4, 5, 6, 7, 8], backGTI.loc[i,"No_1_RowBackGTI":'No_9_RowBackGTI'])
           
                #FrontIrrad=row[19:27]
                #RearIrrad=row[28:36]      #row[28]:row[36] --> Back ["No_1_RowBackGTI":'No_9_RowBackGTI']
                 
                #  FrontIrrad=([float(row[19])/1000, float(row[20])/1000, float(row[21])/1000, float(row[22])/1000, 
                # float(row[23])/1000, float(row[24])/1000, float(row[25])/1000, float(row[26])/1000, float(row[27])/1000])
                
                #For Portrait
                #FrontIrrad=(float(row[19]), float(row[20]), float(row[21]), float(row[22]), float(row[23]), 
                #             float(row[24]), float(row[25]) , float(row[26]), float(row[27]))  # Dividing by 1 sun (1000 W/m2)
                #RearIrrad=(float(row[28]), float(row[29]), float(row[30]), float(row[31]), 
                #            float(row[32]), float(row[33]), float(row[34]), float(row[35]), float(row[36]))
                #For Landscape
                FrontIrrad=(float(row[21]), float(row[22]), float(row[23]), float(row[24]), float(row[25]), 
                             float(row[26]))
                RearIrrad=(float(row[27]), float(row[28]), float(row[29]), float(row[30]), 
                            float(row[31]), float(row[32]))
                #print "Front Rad", FrontIrrad
                #print "Rear Rad", RearIrrad
                interpolA = 0.005  # More accurate interpolation. Do 0.01 as an option.
                #numsens=9 #For Portrait
                numsens=6 #For Landscape
                
                [PmaxIdeal, PmaxUnmatched, PmaxAvg] = LandscapeSingleHour(FrontIrrad, RearIrrad, Tamb, Vwind, numsens, interpolA,IVArray,beta_voc_all,m_all,bee_all)
                               
               
                row.append(PmaxIdeal)
                row.append(PmaxUnmatched)
                row.append(PmaxAvg)
                writer.writerow(row)

        #writer.writerows(all)
        toc = time.clock()
        print "Results Ready in ", (toc-tic)/60.0 , " minutes"
        

