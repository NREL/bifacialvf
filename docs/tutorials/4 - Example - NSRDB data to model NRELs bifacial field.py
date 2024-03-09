#!/usr/bin/env python
# coding: utf-8

# # 4 - Example - NSRDB data to model NRELs bifacial field
# 

# <a id='step1'></a>

# In[1]:


from pathlib import Path
import os

testfolder = Path().resolve().parent.parent / 'bifacialvf' / 'TEMP' / 'BARN'

if not os.path.exists(testfolder):
    os.makedirs(testfolder)


# In[3]:


import pandas as pd
import bifacialvf


# In[4]:


# This information helps with debugging and getting support :)
import sys, platform
print("Working on a ", platform.system(), platform.release())
print("Python version ", sys.version)
print("Pandas version ", pd.__version__)
print("bifacialVF version ", bifacialvf.__version__)


# In[5]:


NREL_API_KEY = None  # <-- please set your NREL API key here
# note you must use "quotes" around your key as it is a string.

if NREL_API_KEY is None:
       NREL_API_KEY = 'DEMO_KEY'  # OK for this demo, but better to get your own key


# In[6]:


import pvlib

# NREL Coords: 39.7407° N, 105.1686° W

metdata, metadata = pvlib.iotools.get_psm3(
    latitude=39.7407, longitude=-105.1686,
    api_key=NREL_API_KEY,
    email='silvana.ovaitt@nrel.gov',  # <-- any email works here fine
    names='2022', map_variables=True, leap_day=False)


# # 3. Set Parameters

# In[7]:


# Variables
sazm = 180                  # PV Azimuth(deg) or tracker axis direction
cw = 2 # m . 1-up portrait
pitch=5.7
hub_height=1.5 # m
norm_hub_height = 1.5/cw
norm_pitch = 5.7/cw
rowType = "interior"        
sensorsy = 12                # sampling areas in the module. edges would be 1 and 12.
PVfrontSurface = "glass"    # options: glass or ARglass
PVbackSurface = "glass"     # options: glass or ARglass

# Tracking instructions
tracking=True
backtrack=True
limit_angle = 50

deltastyle = 'exact'  # NSRDB downloads data at center hour 11:30, 12:30, etc... 


# In[8]:


# We have outdated names that do not match NSRDB. renaming.
#name in NSRDB : name in bifacialVF
metdata.rename(columns={"dni": "DNI", "dhi": "DHI", 'temp_air':'DryBulb', 'wind_speed':'VWind'}, inplace=True)
metadata['TZ'] = metadata['Time Zone']
metadata['city'] = metadata['City']


# In[9]:


writefiletitle = os.path.join(testfolder, 'Dirk_2022_Results.csv')
bifacialvf.simulate(metdata, metadata, writefiletitle=writefiletitle, 
            sazm=sazm, pitch=norm_pitch, hub_height=norm_hub_height, 
         rowType=rowType, sensorsy=sensorsy, 
         PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
          tracking=tracking, backtrack=backtrack, 
         limit_angle=limit_angle, deltastyle=deltastyle)


