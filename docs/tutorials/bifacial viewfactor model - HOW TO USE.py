#!/usr/bin/env python
# coding: utf-8

# In[1]:


# STEP 1: Install and import the program
#
# clone the repo to your local directory
# navigate to the \bifacialvf directory which contains setup
# run pip install -e .

import bifacialvf    
import os
# change directory to \bifacialvf\ root
os.chdir(os.path.dirname(bifacialvf.__file__))
get_ipython().run_line_magic('pwd', '')


# In[2]:


#2. Set the Values of your test
# Remember all values are normalized to your panel length (slope, which will equal 1).
# If your slope is different than 1 m, desired C and rtr (or D) will need to be 
# divided by the slope length.
# i.e.: panel 1mx1.59m, in portrait mode means slope = 1.59. For a height C of 1m, C = 1/1.59. 
#         For a rtr of 1.5m, D=0.51519/1.59 if beta tilt angle = 10 

# Set mandatory variables
beta = 10                   # PV tilt (deg)
sazm = 180                  # PV Azimuth(deg)
clearance_height = 1                      # GroundClearance(panel slope lengths)
pitch = 1.5              # row to row spacing in panel lengths. 
GCR = 1.0/pitch             # Ground Cover Ratio

TMYtoread = "data/724010TYA.csv"   # VA Richmond
writefiletitle = "data/Output/test.csv"

# Set optional variables.  These are the default values
rowType = "interior"        # RowType(first interior last single)
transFactor = 0.013         # TransmissionFactor(open area fraction)
sensorsy = 6                # CellRows(# hor rows in panel)   This is the number of irradiance values returned along module chord
PVfrontSurface = "glass"    # PVfrontSurface(glass or ARglass)
PVbackSurface = "glass"     # PVbackSurface(glass or ARglass)
albedo = 0.62               # ground albedo beneath system
# 1-axis tracking instructions (optional)
tracking=False
backtrack=True       # backtracking optimization as defined in pvlib


deltastyle='TMY3'


# Notice that fixed systems use the parameter D, which is the space between the end of one panel and the beginning of the end, so it relates to the inclination angle beta. To go from a rtr value to D, you can use
# 
# from rtr_and_D_calculation import rtr2D
# D = rtr2D(beta, rtr)

# In[3]:


myTMY3, meta = bifacialvf.bifacialvf.readInputTMY(TMYtoread)


# In[4]:


#3. Call the function.


bifacialvf.simulate(myTMY3, meta, writefiletitle, beta, sazm, 
                clearance_height=clearance_height, pitch=pitch, rowType=rowType, transFactor=transFactor, sensorsy=sensorsy,
                PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, albedo=albedo, 
                tracking=tracking, backtrack=backtrack, deltastyle=deltastyle)
    


# In[5]:


#4. Load the results from the resultfile
(data, metadata) = bifacialvf.loadVFresults(writefiletitle)
#print data.keys()
# calculate average front and back global tilted irradiance across the module chord
data['GTIFrontavg'] = data[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
data['GTIBackavg'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)

# Print the annual bifacial ratio
frontIrrSum = data['GTIFrontavg'].sum()
backIrrSum = data['GTIBackavg'].sum()
print('The bifacial ratio for ground clearance {} and rtr spacing {} is: {:.1f}%'.format(clearance_height,pitch,backIrrSum/frontIrrSum*100))


# In[ ]:


metdata


# In[ ]:


# plot the rear irradiance distribution for a single point in time. 1999-07-06
import matplotlib.pyplot as plt
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'inline')

data['GTIBackstd'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].std(axis=1)
data['measdatetime'] = pd.to_datetime(data[['Year', 'Month', 'Day','Hour','Minute']])
singleday = (data['measdatetime'] > '1999-07-06') & (data['measdatetime']<'1999-07-07')
singleday2 = (data['measdatetime'] > '1999-07-12') & (data['measdatetime']<'1999-07-13')

fig, ax = plt.subplots()
ax1 = ax.twinx()
ax1.plot(data['measdatetime'][singleday],data['GTIFrontavg'][singleday],'k')
ax1.set_ylabel('Front Irradiance (Wm-2)')
ax.set_ylabel('Rear Irradiance (Wm-2)')
ax.plot(data['measdatetime'][singleday], data['No_1_RowBackGTI'][singleday],'r' , alpha =0.5)
ax.plot(data['measdatetime'][singleday], data['No_2_RowBackGTI'][singleday], 'b', alpha = 0.5)
ax.plot(data['measdatetime'][singleday], data['No_6_RowBackGTI'][singleday], 'g', alpha = 0.5)
ax.set_title('Sunny day')
fig.autofmt_xdate()
fig.tight_layout()


fig2, ax2 = plt.subplots()
ax3 = ax2.twinx()
ax3.plot(data['measdatetime'][singleday2],data['GTIFrontavg'][singleday2],'k')
ax3.set_ylabel('Front Irradiance (Wm-2)')
ax2.set_ylabel('Rear Irradiance (Wm-2)')
ax2.plot(data['measdatetime'][singleday2], data['No_1_RowBackGTI'][singleday2],'r' , alpha =0.5)
ax2.plot(data['measdatetime'][singleday2], data['No_2_RowBackGTI'][singleday2], 'b', alpha = 0.5)
ax2.plot(data['measdatetime'][singleday2], data['No_6_RowBackGTI'][singleday2], 'g', alpha = 0.5)
ax2.set_title('Cloudy day')
fig2.autofmt_xdate()
fig2.tight_layout()



# # Example of 1-axis tracking
# 

# In[11]:


# tracker geometry options:
module_height = 4  # module portrait dimension in meters (2-up portrait 'module')
gcr = 0.35   # ground cover ratio,  = module_height / pitch
albedo = 0.3     # ground albedo
hub_height_meters = 2.5   # tracker height at 0 tilt in meters (hub height)
limit_angle = 45 # tracker rotation limit angle

# Set mandatory variables
beta = 0                   # Tracking doesn't allow tilted trackers at the moment.
sazm = 180                  # PV Azimuth(deg) - tracker axis orientation for E-W tracked system.
hub_height = hub_height_meters / module_height                      # GroundClearance(panel slope lengths)
pitch = 1 / gcr              # row to row spacing in panel lengths. 


TMYtoread = "data/724010TYA.csv"   # VA Richmond
writefiletitle = "data/Output/1AxisTest.csv"

# Set optional variables.  These are the default values
rowType = "interior"        # RowType(first interior last single)
transFactor = 0.013         # TransmissionFactor(open area fraction)
sensorsy = 6                # CellRows(# hor rows in panel)   This is the number of irradiance values returned along module chord
PVfrontSurface = "ARglass"    # PVfrontSurface(glass or ARglass)
PVbackSurface = "glass"     # PVbackSurface(glass or ARglass)
# 1-axis tracking instructions (optional)
tracking=True
backtrack=True       # backtracking optimization as defined in pvlib
deltastyle='TMY3'


# In[14]:


#3. Call the function.
myTMY3, meta = bifacialvf.bifacialvf.readInputTMY(TMYtoread)


bifacialvf.simulate(myTMY3, meta, writefiletitle, beta, sazm, 
                hub_height=hub_height, rowType=rowType, transFactor=transFactor, sensorsy=sensorsy,
                PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, albedo=albedo, 
                tracking=True, backtrack=backtrack, pitch=pitch, limit_angle = limit_angle, deltastyle=deltastyle)
  


# In[ ]:


#4. Load the results from the resultfile
(data, metadata) = bifacialvf.loadVFresults(writefiletitle)
#print data.keys()
# calculate average front and back global tilted irradiance across the module chord in 6 points
data['GTIFrontavg'] = data[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
data['GTIBackavg'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)

# Print the annual bifacial ratio
frontIrrSum = data['GTIFrontavg'].sum()
backIrrSum = data['GTIBackavg'].sum()
print('The 1-axis tracked bifacial ratio for ground clearance {} and row to row {} is: {:.1f}%\n\n *****'.format(C,rtr,backIrrSum/frontIrrSum*100))


# 
# 
