#!/usr/bin/env python
# coding: utf-8

# # PinPV Fig 1
# 
# This journal shows how to use the Analysis function inside bifacialVF to call for PV Mismatch IV curve calculation and plotting for different irradiance values.
# 
# Estimating and parameterizing mismatch power loss in bifacial photovoltaic systems
# Chris Deline,Silvana Ayala Pelaez,Sara MacAlpine,Carlos Olalla
# First published: 06 March 2020 https://doi.org/10.1002/pip.3259Citations: 16
# 

# In[1]:


import bifacialvf
import pandas as pd


# In[2]:


# SETUP:
City='Cairo'
WeatherFile = 'EGY_Cairo.Intl.Airport.623660_ETMY'
Latitude = 30.13
Longitude= 31.4
Timezone = 2
Tilt = 10
Azimuth = 180
GroundClearance = 0.15
RtR = 1.5
RowType = 'interior'
TransmissionFactor = 0.013
CellRows = 6
PVfrontSurface = 'glass'
PVbackSurface = 'glass'
Albedo = 0.62
Tracking = False
calculatePVMismatch = True
PortraitorLandscape = 'landscape'


# In[3]:


# PVMismatch calculation
# BifacialVF front and rear irradiance values for June 21st. at 2 PM
frontGTIrow = [927.6982631, 928.4082993, 928.8050623, 931.0767397, 933.0536854, 933.8064969]
backGTIrow = [131.7201933, 47.90228722, 38.49283059, 54.58803934, 100.6857589, 194.5133033]

portraitorlandscape='landscape'
sensorsy=len(frontGTIrow)
debug=True
numcells=72
plotflag = True

stdpl, cellsx, cellsy = bifacialvf.analysis.setupforPVMismatch(portraitorlandscape=portraitorlandscape, sensorsy=sensorsy, numcells=numcells)

PowerAveraged, PowerDetailed, sunmatAveraged, sunmatDetailed = bifacialvf.analysis.calculateVFPVMismatch(stdpl=stdpl, cellsx=cellsx, cellsy=cellsy, sensorsy=sensorsy, frontGTIrow=frontGTIrow, backGTIrow=backGTIrow, bififactor=1.0, debug=debug, plotflag = plotflag)

print("Results: ", round(PowerAveraged,3), round(PowerDetailed,2))

print(portraitorlandscape, "Cellsx", cellsx, "Cellsy", cellsy)
print(pd.DataFrame(stdpl))
print(pd.DataFrame(sunmatAveraged).round(3))
pd.DataFrame(sunmatDetailed).round(3)


# In[ ]:




