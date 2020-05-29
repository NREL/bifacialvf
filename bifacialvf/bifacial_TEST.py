# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 11:27:08 2018

@author: Silvana
"""

import datetime
import pvlib
import math
import matplotlib.pyplot as plt
import pytz

# Tracker Parameters
rtr=1.5
gcr=1.0/rtr  
max_angle=45

# Date Parameters
year=2018
month=6
day=21
lat=32.133	
lng=-110.95
axis_tilt=0.0
axis_azimuth=0.0 # North-south

# Arrays to save data
betaall=[]
betaall2=[]
timestampall=[]
sazmall=[]

# Loop to calculate tracking angles for algorithm WITH backtracking and WITHOUT backtracking.
for i in range (4, 20):
    for j in range (0, 59):
        hour=i;         minute=j
        myTimestamp=datetime.datetime(year, month, day, hour, minute, 0, 0)                
        
        # Tucson AZ Timezone
        mytz = pytz.timezone('US/Mountain')
        # Assign timezone for Daylight-Savings, March 11 to Nov. 4,
        # Arizona doesn't follow Daylight-Savings so on this times
        # It becomes "Pacific" time on pytz zones...
        if month == 3:
            if day >= 11:
                mytz = pytz.timezone('US/Pacific')
        
        if 4 <= month < 11:
            mytz = pytz.timezone('US/Pacific')
            
        if month == 11:
            if day < 4:
                mytz = pytz.timezone('US/Pacific')    
        
        
        myTimestamp = myTimestamp.replace(tzinfo=mytz)
        
        solpos = pvlib.solarposition.get_solarposition(myTimestamp, lat, lng)
        aazi= solpos['azimuth']
        azen= solpos['zenith']

        trackingdata = pvlib.tracking.singleaxis(azen, aazi, axis_tilt, axis_azimuth, max_angle, backtrack=True, gcr=gcr)
        beta=trackingdata['surface_tilt'][0] # Trackingdata tracker_theta
        sazm = trackingdata['surface_azimuth'][0]
        if math.isnan(beta):
            beta=0

        # Rotate system if past sun's zenith
        if beta<0:
            beta = -beta;

        trackingdata2 = pvlib.tracking.singleaxis(azen, aazi, axis_tilt, axis_azimuth, max_angle, backtrack=False, gcr=gcr)
        beta2=trackingdata2['surface_tilt'][0] # Trackingdata tracker_theta
        sazm2 = trackingdata2['surface_azimuth'][0]
        if math.isnan(beta2):
            beta2=0

        if beta2<0:
            beta2 = -beta2;

        betaall.append(beta)
        betaall2.append(beta2)
        timestampall.append(myTimestamp)
        sazmall.append(sazm)
        
import matplotlib.dates as mdates
import matplotlib

# Plots Comparison of Backtrack and NO Backtrack
matplotlib.rcParams.update({'font.size': 22})
fig = plt.figure(figsize=(10,5))        
fig, ax = plt.subplots(figsize=(10,5))
ax.fmt_xdata = mdates.DateFormatter('%H')
plt.plot(timestampall, betaall,linewidth=7.0)
plt.plot(timestampall, betaall2, linewidth=3.0)
#matplotlib.rcParams['timezone'] = 'US/Pacific'
#mdates.DateFormatter('%H:%M', tz=tz.gettz('Europe/Berlin'))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H', tz=mytz))
plt.xlabel('Hour')
plt.ylabel('Tracking Angle [$^\circ$]')
#fig.autofmt_xdate() # makes ticks sideays
plt.show()
