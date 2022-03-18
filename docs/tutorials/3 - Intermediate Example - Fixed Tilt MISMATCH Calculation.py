#!/usr/bin/env python
# coding: utf-8

# # 1 - Introductory Example: Fixed-Tilt simple setup
# 
# This jupyter journal will walk us through the creation of the most basic fixed-tilt simulation possible with bifacialvf.
# We will simulate a 1-up landscape system over a white rooftop.
# 
# Steps include:
# 
# <ol>
#     <li> <a href='#step1'> Create a folder for your simulation, and Load bifacialvf </a></li> 
#     <li> <a href='#step2'> Download Weather File and Read </a></li> 
#     <li> <a href='#step3'> Set all your Parameters </a></li> 
#     <li> <a href='#step4'> Run simulation </a></li> 
#     <li> <a href='#step5'> Read Results </a></li> 
#     <li> <a href='#step6'> Plot Results </a></li> 
# </ol>
# 

# <a id='step1'></a>

# In[1]:


from pathlib import Path
import os
import bifacialvf

# IO Files
testfolder = Path().resolve().parent.parent / 'bifacialvf' / 'TEMP' / 'Tutorial_03'
if not os.path.exists(testfolder):
    os.makedirs(testfolder)


# In[2]:


# Download and Read input
TMYtoread=bifacialvf.getEPW(lat=37.5407,lon=-77.4360, path = testfolder)
myTMY3, meta = bifacialvf.readInputTMY(TMYtoread)
deltastyle = 'TMY3'  # 


# In[3]:


# Variables
tilt = 10                   # PV tilt (deg)
sazm = 180                  # PV Azimuth(deg) or tracker axis direction
albedo = 0.62               # ground albedo
clearance_height=0.4
pitch = 1.5                   # row to row spacing in normalized panel lengths. 
rowType = "interior"        # RowType(first interior last single)
transFactor = 0.013         # TransmissionFactor(open area fraction)
sensorsy = 6                # sensorsy(# hor rows in panel)   <--> THIS ASSUMES LANDSCAPE ORIENTATION 
PVfrontSurface = "glass"    # PVfrontSurface(glass or ARglass)
PVbackSurface = "glass"     # PVbackSurface(glass or ARglass)

# Calculate PV Output Through Various Methods    
# This variables are advanced and explored in other tutorials.
calculateBilInterpol = True         # Only works with landscape at the moment.
calculatePVMismatch = True
portraitorlandscape='landscape'   # portrait or landscape
cellsnum = 72
bififactor = 1.0
agriPV = False                       # Returns ground irradiance values

# Tracking instructions
tracking=False
backtrack=False
limit_angle = 60

writefiletitle = os.path.join(testfolder, 'Tutorial3_Results.csv')
myTMY3 = myTMY3.iloc[0:24].copy()  # Simulate just the first 24 hours of the data file for speed on this example
bifacialvf.simulate(myTMY3, meta, writefiletitle=writefiletitle, 
         tilt=tilt, sazm=sazm, pitch=pitch, clearance_height=clearance_height, 
         rowType=rowType, transFactor=transFactor, sensorsy=sensorsy, 
         PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
         albedo=albedo, tracking=tracking, backtrack=backtrack, 
         limit_angle=limit_angle, deltastyle=deltastyle,
        calculateBilInterpol=calculateBilInterpol, calculatePVMismatch=calculatePVMismatch,
        cellsnum = cellsnum, bififactor=bififactor, agriPV=agriPV)


# # 5. Load the results from the resultfile
# 

# In[4]:


from bifacialvf import loadVFresults
mismatchResultstitle = os.path.join(testfolder, 'Tutorial3_Results_PVMismatch.csv')
(data, metadata) = loadVFresults(mismatchResultstitle)


# In[5]:


data.keys()

