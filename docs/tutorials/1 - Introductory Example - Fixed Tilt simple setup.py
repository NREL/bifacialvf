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
TMYtoread=bifacialvf.getEPW(lat=37.5407,lon=-77.4360, path = testfolder)
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

# In[17]:


writefiletitle = os.path.join(testfolder, 'Tutorial1_Results.csv')
#myTMY3 = myTMY3.iloc[0:24].copy()  # Simulate just the first 24 hours of the data file for speed on this example
bifacialvf.simulate(myTMY3, meta, writefiletitle=writefiletitle, 
         tilt=tilt, sazm=sazm, pitch=pitch, clearance_height=clearance_height, 
         rowType=rowType, transFactor=transFactor, sensorsy=sensorsy, 
         PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
         albedo=albedo, tracking=tracking, backtrack=backtrack, 
         limit_angle=limit_angle, deltastyle=deltastyle)


# <a id='step5'></a>

# # 5. Load the results from the resultfile
# 

# In[23]:


from bifacialvf import loadVFresults
(data, metadata) = loadVFresults(writefiletitle)

# calculate average front and back global tilted irradiance across the module chord
data['GTIFrontavg'] = data[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
data['GTIBackavg'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)


# Print the annual bifacial ratio
frontIrrSum = data['GTIFrontavg'].sum()
backIrrSum = data['GTIBackavg'].sum()
print('The bifacial ratio for ground clearance {} and rtr spacing {} is: {:.1f}%'.format(clearance_height,pitch,backIrrSum/frontIrrSum*100))


# In[24]:


data


# # 6. Plot some results 

# In[41]:


# plot the rear irradiance distribution for a single point in time. 1999-07-06
import matplotlib.pyplot as plt
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')

data['GTIBackstd'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].std(axis=1)
data.set_index(pd.to_datetime(data['date']), inplace=True, drop=True)
data.index = data.index.map(lambda t: t.replace(year=2021))   # Chagning to be the same year
singleday = (data.index > '2021-07-09') & (data.index<'2021-07-10')
singleday2 = (data.index > '2021-07-15') & (data.index<'2021-07-16')

fig, ax = plt.subplots()
ax1 = ax.twinx()
ax1.plot(data.index[singleday],data['GTIFrontavg'][singleday],'k')
ax1.set_ylabel('Front Irradiance (Wm-2)')
ax.set_ylabel('Rear Irradiance (Wm-2)')
ax.plot(data.index[singleday], data['No_1_RowBackGTI'][singleday],'r' , alpha =0.5)
ax.plot(data.index[singleday], data['No_2_RowBackGTI'][singleday], 'b', alpha = 0.5)
ax.plot(data.index[singleday], data['No_6_RowBackGTI'][singleday], 'g', alpha = 0.5)
ax.set_title('Sunny day')
fig.autofmt_xdate()
fig.tight_layout()


fig2, ax2 = plt.subplots()
ax3 = ax2.twinx()
ax3.plot(data.index[singleday2],data['GTIFrontavg'][singleday2],'k')
ax3.set_ylabel('Front Irradiance (Wm-2)')
ax2.set_ylabel('Rear Irradiance (Wm-2)')
ax2.plot(data.index[singleday2], data['No_1_RowBackGTI'][singleday2],'r' , alpha =0.5)
ax2.plot(data.index[singleday2], data['No_2_RowBackGTI'][singleday2], 'b', alpha = 0.5)
ax2.plot(data.index[singleday2], data['No_6_RowBackGTI'][singleday2], 'g', alpha = 0.5)
ax2.set_title('Cloudy day')
fig2.autofmt_xdate()
fig2.tight_layout()

