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
# </ol>
# 

# <a id='step1'></a>

# 
# ## 1. Create a folder for your simulation, and load bifacial_radiance 
# 
# First let's set the folder where the simulation will be saved. By default, this is the TEMP folder in the bifacialvf distribution.
# 
# The lines below find the location of the folder relative to this Jupyter Journal. You can alternatively point to an empty directory (it will open a load GUI Visual Interface) or specify any other directory in your computer, for example:
# 
# #### testfolder = r'C:\Users\sayala\Documents\bifacialVFResults'
# 
# 

# In[1]:


from pathlib import Path
import os

# IO Files
testfolder = Path().resolve().parent.parent / 'bifacialvf' / 'TEMP' / 'Tutorial_01'
if not os.path.exists(testfolder):
    os.makedirs(testfolder)


# This will load bifacial_radiance and other libraries from python that will be useful for this Jupyter Journal:

# In[2]:


import bifacialvf


# <a id='step2'></a>

# ## 2. Download and Load Weather Files
# 
# There are various options provided in bifacialVF to load weatherfiles. getEPW is useful because you just set the latitude and longitude of the location and it donwloads the meteorologicla data for any location. 

# In[3]:


# Download and Read input
TMYtoread=bifacialvf.getEPW(lat=30,lon=20, path = testfolder)
myTMY3, meta = bifacialvf.readInputTMY(TMYtoread)


# We can also specify if data is right labeled (like TMY3 files), or left labeled (like SAM input files) to calculate teh sun position. By default, 'TMY3' format is selected 

# In[4]:


deltastyle = 'TMY3'  # 


# <a id='step3'></a>

# # 3. Set Parameters

# In[5]:


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
#calculateBilInterpol = True         # Only works with landscape at the moment.
#calculatePVMismatch = True
#portraitorlandscape='landscape'   # portrait or landscape
#cellsnum = 72
#bififactor = 1.0
#agriPV = True                       # Returns ground irradiance values

# Tracking instructions
tracking=False
backtrack=False
limit_angle = 60


# <a id='step4'></a>

# # 4. Run Simulation

# In[6]:


writefiletitle = os.path.join(testfolder, 'Tutorial1_Results.csv')
myTMY3_2 = myTMY3.iloc[0:24].copy()  # Simulate just the first 24 hours of the data file for speed on this example
bifacialvf.simulate(myTMY3_2, meta, writefiletitle=writefiletitle, 
         tilt=tilt, sazm=sazm, pitch=pitch, clearance_height=clearance_height, 
         rowType=rowType, transFactor=transFactor, sensorsy=sensorsy, 
         PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
         albedo=albedo, tracking=tracking, backtrack=backtrack, 
         limit_angle=limit_angle, deltastyle=deltastyle)


# <a id='step5'></a>

# # 5. Load the results from the resultfile
# 

# In[7]:


from bifacialvf import loadVFresults
(data, metadata) = loadVFresults(writefiletitle)


# In[8]:


data

