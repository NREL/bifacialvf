# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 13:02:33 2018

@author: cdeline


Bonanza data analysis - look at raw data in \Bonanza Trina raw data and 
\Bonanza bifacial raw data

Site raw data - coming from egauge Keeton site
Bifacial data: http://egauge34875.egaug.es/58C32/
Trina data: http://egauge34871.egaug.es/58C32/     

Version1:  just look at site data.  Filter from 8/16 onward.
 - Note.  Inverter D & F are swapped starting 2/16/18
V.2:  Import TWC weather data for PR and modeling. Save in \Bonanza weather data
  - merge weather and modeled data with inverter B and E data.  Save in directory\Bonanza_Merge.pickle
  - start version control with github


Version3: Teague data analysis from C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon\Teague bifacial
5/8/18: Additional site bifacial data:
    http://egauge34873.egaug.es/58C32/   Teague #1  100 kw  Silfab bifacial
    http://egauge34874.egaug.es/58C32/    Teague #2 100 kw Silfab bifacial
    
    Similar output to Keeton bifacial
    
    Strategy:  Use the Teague site as the bifacial reference, since it doesn't have nearly as many outages
    and changes in ground cover.  Compare with Keeton inverter E for the enhanced albedo study.  Use the 
    Keeton Trina modules as the monofacial reference.
    
    Also - Keeton Silfab modules has white fabric deployed beneath them starting February.
    
    Variable decoder:  
        df:  raw 15-minute measured site data from Teague, read from TeagueBifacial.csv
        
        psm: Weather.com hourly weather data, along with modeled Trina and Silfab.  read from Bonanza_TWC.csv
             psm['silfab_meas'] merged from  df['bifacial'] which is measured Teague bifacial data.
             Combined file saved as 'Teague_Merge.pickle'
             
        data:  hourly Keeton measured and modeled data read from 'Keeton_Merge.pickle'.  Silfab_meas is from InvE.  Trina_meas is from InvB.
            bi1, bi2 and bi3 are inverters D,E and F respectively.
        
    
"""
from __future__ import division
import pandas as pd
import numpy as np
import os
import glob
#import matplotlib.pyplot as plt, mpld3
import matplotlib.pyplot as plt

#from matplotlib import style
#style.use('fivethirtyeight')
import pvlib

# Set default font for plots
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Times']})
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']
plt.rcParams['axes.linewidth'] = 0.2 #set the value globally


# In[16]:
# Save combined merged document with psm model and Teague measured.
directory = r'C:\Users\Silvana\Documents\GitHub\bifacialvf\bifacialvf\bonanza'

psm_filtered = pd.read_pickle(os.path.join(directory,'Bonanza_Merge.pickle'))


TWC = ({'Month': psm_filtered.index.month, 
        'Day': psm_filtered.index.day, 
        'Year': psm_filtered.index.year, 
        'Hour': psm_filtered.index.hour, 
        'Minute': psm_filtered.index.minute,
        'Global': psm_filtered.ghi, 
        'Direct': psm_filtered.dni, 
        'Diffuse': psm_filtered.dhi, 
        'Temp': psm_filtered['SurfaceTemperatureCelsius'], 
        'Humidity': psm_filtered['RelativeHumidityPercent'], 
        'Wind': psm_filtered['WindSpeedKph']*1000/3600}) 
TWC = pd.DataFrame.from_records(TWC, index=psm_filtered.index)
# set up bifacial simulation
writefiletitle = 'Silfab_model_1axis_albedo0p2.csv'       
gcr = 0.35;  rtr = 1.0/gcr
C = 0.75  # clearance height: 2m table, 1.5m hub height

# NOTE: this uses simulateTWC below which needs to be read into memory.
#os.chdir(r'C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon')
os.chdir(r'C:\Users\Silvana\Documents\GitHub\bifacialvf\bifacialvf\bonanza')

# In[]:
    
os.chdir(r'C:\Users\Silvana\Documents\GitHub\bifacialvf\bifacialvf')
writefiletitle = r'C:\Users\Silvana\Documents\GitHub\\bifacialvf\\bifacialvf\Silfab_model_1axis_albedo0p2.csv'       


from modelPaperPVSC import simulateTWC

# Albedo range: 0.15 - 0.2 ?
simulateTWC(TWC, writefiletitle,  beta = 0, sazm = 180, C = C, D = None,
             rowType = 'interior', transFactor = 0.0, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = 0.2,  
             tracking = True, backtrack = True, rtr = rtr, Cv = None, offset = 0, max_angle = 45,
             lat = 42.07823, lng = -121.27869, tz = -7, name = 'BONANZA', dataInterval= 60)
# In[]:
    
os.chdir(r'C:\Users\Silvana\Documents\GitHub\bifacialvf\bifacialvf')
writefiletitle = r'C:\Users\Silvana\Documents\GitHub\\bifacialvf\\bifacialvf\Silfab_model_1axis_albedo0p25.csv'       


from modelPaperPVSC import simulateTWC

# Albedo range: 0.15 - 0.2 ?
simulateTWC(TWC, writefiletitle,  beta = 0, sazm = 180, C = C, D = None,
             rowType = 'interior', transFactor = 0.0, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = 0.25,  
             tracking = True, backtrack = True, rtr = rtr, Cv = None, offset = 0, max_angle = 45,
             lat = 42.07823, lng = -121.27869, tz = -7, name = 'BONANZA', dataInterval= 60)


# In[]:
    
os.chdir(r'C:\Users\Silvana\Documents\GitHub\bifacialvf\bifacialvf')
writefiletitle = r'C:\Users\Silvana\Documents\GitHub\\bifacialvf\\bifacialvf\Silfab_model_1axis_albedo0p44.csv'       


from modelPaperPVSC import simulateTWC

# Albedo range: 0.15 - 0.2 ?
simulateTWC(TWC, writefiletitle,  beta = 0, sazm = 180, C = C, D = None,
             rowType = 'interior', transFactor = 0.0, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = 0.44,  
             tracking = True, backtrack = True, rtr = rtr, Cv = None, offset = 0, max_angle = 45,
             lat = 42.07823, lng = -121.27869, tz = -7, name = 'BONANZA', dataInterval= 60)

# In[]:
    
os.chdir(r'C:\Users\Silvana\Documents\GitHub\bifacialvf\bifacialvf')
writefiletitle = r'C:\Users\Silvana\Documents\GitHub\\bifacialvf\\bifacialvf\Silfab_model_1axis_albedo0p15.csv'       


from modelPaperPVSC import simulateTWC

# Albedo range: 0.15 - 0.2 ?
simulateTWC(TWC, writefiletitle,  beta = 0, sazm = 180, C = C, D = None,
             rowType = 'interior', transFactor = 0.0, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = 0.15,  
             tracking = True, backtrack = True, rtr = rtr, Cv = None, offset = 0, max_angle = 45,
             lat = 42.07823, lng = -121.27869, tz = -7, name = 'BONANZA', dataInterval= 60)


# In[]:
writefiletitle = r'C:\Users\Silvana\Documents\GitHub\\bifacialvf\\bifacialvf\Silfab_model_1axis_albedo0p2.csv'       

from loadVFresults import loadVFresults
(model, metadata) = loadVFresults(writefiletitle)
# calculate average front and back global tilted irradiance across the module chord
model['GTIFrontavg'] = model[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
model['GTIBackavg'] = model[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)
# index by measdatetime
model.index = pd.to_datetime({'Year':model.Year,'Month':model.Month,'Day':model.Day,'Hour':model.Hour})

# monthly cumulative production for each string
mindate = psm_filtered.index[0]
minmonth = psm_filtered.index.month[0] + psm_filtered.index.year[0] * 12
minday = psm_filtered.index.dayofyear[0] + psm_filtered.index.year[0] * 365
model['sequential_month'] = (model.index.month + model.index.year * 12-minmonth)
model['month'] = model.index.month
model['sequential_day'] = (model.index.dayofyear + model.index.year * 365-minday)

# In[]:
psmmonth = psm_filtered.groupby(psm_filtered.sequential_month).sum()

def add_month_values(dfmonth,mindate):
    '''
    monthly cumulative total.  Import: dfmonth monthly aggregated dataframe 
    e.g. df_filter.groupby(df_filter.sequential_month).sum()  in kWh (need to divide by 4 for 15-min data)
    
    '''    

    dfmonth['sequential_month'] = dfmonth.index
    # save the month and year in mm-yy format
    from dateutil.relativedelta import relativedelta
    dfmonth['mmyy'] = dfmonth['sequential_month'].apply(lambda x: (mindate + relativedelta(months=+x)).strftime("%m-%y"))

    # also save the month name
    dfmonth['month'] = dfmonth['sequential_month'].apply(lambda x: (mindate + relativedelta(months=+x)).strftime("%b"))
    return dfmonth



psmmonth = add_month_values(psmmonth,mindate = psm_filtered.index[0])  #defined above

# performance ratio - define as separate function. Apply to PSM data either hourly, daily monthly or cumulative.
def addPR(df):
    df['Yr'] = df.poa_irradiance / 1000
    df['YfSilfab'] = df.silfab_meas / 285 
    df['YfSilfabModel'] = df.silfab_model / 285
    df['SilfabPR'] = df.YfSilfab / df.Yr
    df['SilfabModelPR'] = df.YfSilfabModel / df.Yr
    try:
        df['YfSilfabKeeton'] = df.silfab_Keeton_meas / 285
        df['SilfabKeetonPR'] = df.YfSilfabKeeton / df.Yr
        df['YfTrina'] = df.trina_meas / 300
        df['YfTrinaModel'] = df.trina_model   / 300
        df['BifiRatio'] = df.YfSilfab  / df.YfTrina
        df['TrinaPR'] = df.YfTrina / df.Yr
        df['TrinaModelPR'] = df.YfTrinaModel / df.Yr
    except:
        pass
    return df


psmmonth = addPR(psmmonth)

bifimonth = model.groupby(model.sequential_month).sum()
# add modeled bifi gain and measured PR / model PR difference
bifimonth['BGE_model'] = bifimonth.GTIBackavg / bifimonth.GTIFrontavg * 0.95  # Silfab bifaciality: 0.95?
bifimonth['PR_meas_diff'] = (psmmonth.SilfabPR / psmmonth.SilfabModelPR * psmmonth.TrinaModelPR / psmmonth.TrinaPR)-1
bifimonth = add_month_values(bifimonth,mindate)

# In[]:

psmyear = psm_filtered.sum()
psmyear = addPR(psmyear)  # Teague hourly data + modeled

bifiyear = model.sum()
bifiyear['BGE_model'] = bifiyear.GTIBackavg / bifiyear.GTIFrontavg * 0.95  # Silfab bifaciality: 0.95?
bifiyear['PR_meas_diff'] = (psmyear.SilfabPR / psmyear.SilfabModelPR * psmyear.TrinaModelPR / psmyear.TrinaPR)-1

# plot the modeled bifacial gain vs. monthly measured PR difference between Teague (psm bifacial) and Keeton (data monofacial)


# In[]:  
    # ALL OF THIS CELL IS JUST TO OBTAIN THE MONTHS NAMES FOR THE NEXT PLOT......... SRSLY
    
def analyzeMonth(df_filter,plotflag = True):
    '''
    analyzeMonth:  aggregate monthly and plot for filtered Teague or Keeton data
    input: df_filter:  filtered dataframe including bi1 through bi6 or bi1-bi3 and mono1-mono3
    return: dfmonth:  monthly sum in kWh
    '''
    #sum monthly values
    dfmonth = df_filter.groupby(df_filter.sequential_month).sum() / 4  #15 minute data from August 18 onward.
    mindate = df_filter.index[0]
    dfmonth = add_month_values(dfmonth,mindate)  #defined above

    
    if plotflag:
        #  compare monthly production of 5 inverters
        fig4 = plt.figure()
        width = .1
        #plt.rcParams['figure.figsize'] = (6, 3.5)
        plt.bar(dfmonth.index, dfmonth.bi1, width = width, label = 'bi1')
        plt.bar(dfmonth.index+width, dfmonth.bi2, width = width, label = 'bi2')
        plt.bar(dfmonth.index+width*2, dfmonth.bi3, width = width, label = 'bi3')
        try:
            plt.bar(dfmonth.index+width*3, dfmonth.bi5, width = width, label = 'bi5')
            plt.bar(dfmonth.index+width*4, dfmonth.bi6, width = width, label = 'bi6')
        except:
            pass
        try:
            plt.bar(dfmonth.index+width*3, dfmonth.mono1, width = width, label = 'mono1')
            plt.bar(dfmonth.index+width*4, dfmonth.mono2, width = width, label = 'mono2')
        except:
            pass
        plt.legend()
        plt.title('Monthly kWh/kWp production  ', fontsize=16)
        plt.xticks(dfmonth.index, dfmonth.mmyy, rotation='horizontal')
        fig4.autofmt_xdate()
        plt.ylim([0,300])
        plt.ylabel('Monthly kWh / kWp ')
        plt.show()
    
    return dfmonth

def add_sequential_counters(df):
    #
    minmonth = df.index.month[0] + df.index.year[0] * 12
    minday = df.index.dayofyear[0] + df.index.year[0] * 365
    df['sequential_month'] = (df.index.month + df.index.year * 12-minmonth)
    df['month'] = df.index.month
    df['sequential_day'] = (df.index.dayofyear + df.index.year * 365-minday)
    return df


def analyzeKeeton(df,plotflag = True):
    '''
    analyzeKeeton:  import dataframe df and append normalized bifacial and monofacial series.  return df_filter (filtered data)
    
    Input:  df: datafraome including 'Date & Time', Inverter A [kW] through Inverter F [kW]
    Output: df_filter:  data filtered
    
    Adapted from Bonanza_data_analysis.py
    
    Updated notes:  
        - mono3 has a stuck tracker in Aug.2018 (possibly before too)
        - bi3 has AC clipping issues from April onward
    
    '''
    
    #  weird switch around 2/16. Inverter F is now 37.62 kW, inverter D is 31.35kW
    dateOfChange = '2018-02-17'

    df['mono1'] = df['INVERTER A [kW]'] / 33.3
    df['mono2'] = df['INVERTER B [kW]'] / 33.3
    df['mono3'] = df['INVERTER C [kW]'] / 33.3
    df['bi1'] = df['INVERTER D [kW]'] / 37.62
    df['bi2'] = df['INVERTER E [kW]'] / 31.35
    df['bi3'] = df['INVERTER F [kW]'] / 31.35
    
    df.loc[df.index < dateOfChange, 'bi3'] = df['INVERTER D [kW]'] /37.62
    df.loc[df.index > dateOfChange, 'bi3'] = df['INVERTER F [kW]'] /37.62
    df.loc[df.index < dateOfChange, 'bi1'] = df['INVERTER F [kW]'] /31.35
    df.loc[df.index > dateOfChange, 'bi1'] = df['INVERTER D [kW]'] /31.35
    
    df['monofacial'] = df['mono2']
    df['bifacial'] = df['bi2']
    
    df = add_sequential_counters(df) #add ['sequential_month'], ['sequential_day'], ['month']

    df_filter = filterdf(df) 
    df_filter.index = df_filter.index.tz_localize('Etc/GMT+7') # set to pacific standard time
    
    if plotflag:
        fig = plt.figure()
        plt.rcParams['figure.figsize'] = (6, 4)
        plt.plot(df.index,df.mono1,label = 'InverterA', alpha = 0.5)
        plt.plot(df.index,df.mono2,label = 'InverterB', alpha = 0.4)
        plt.plot(df.index,df.mono3,label = 'InverterC', alpha = 0.5)
        #plt.xlim(['2017-06-10','2017-06-17'])
        #plt.xlim(['2018-02-15','2018-02-25'])
        fig.autofmt_xdate()
        plt.legend(loc='best')
        plt.title('Trina inverters')
        #plt.show()
        mpld3.show()
        
                
        
        # Bifacial inverters 
        fig = plt.figure()
        plt.rcParams['figure.figsize'] = (6, 4)
        plt.plot(df.index,df.bi1 ,label = 'Interior bifacial', alpha = 0.5)#/ 37.62
        plt.plot(df.index,df.bi2 ,label = 'Middle bifacial', alpha = 0.4)#/ 31.35
        plt.plot(df.index,df.bi3 ,label = 'Edge bifacial', alpha = 0.5)#/ 31.35
        #plt.xlim(['2018-02-15','2018-02-25'])
        fig.autofmt_xdate()
        plt.legend(loc='best')
        plt.title('Silfab inverters')
        #plt.show()
        mpld3.show()
        
               

    return df_filter

def filterdf(df):
    '''
    apply uniform filtering to all experiment data to allow everything to line up
        
    '''
    # Startup problems before 6/8/18.  Match to Keeton site starting 8/16/17
    mask1 = (df.index > '2017-08-17')  #((df.index < '2017-08-03') | (df.index > '2017-08-18'))
    
    #Keeton bifacial outage on 12/16
    mask6 = ((df.index < '2017-12-16' ) | ( df.index > '2017-12-17') )
    
    #Keeton  missing bifacial data from 3/22 - 4/2
    mask7 = ((df.index < '2018-3-22' ) | ( df.index > '2018-4-2') )
    
    # filter data with startup issues, missing partial times
    df_filter = df[mask1 & mask6 & mask7]#.dropna()
    
    return df_filter

data = pd.read_pickle(r'C:\Users\Silvana\Documents\GitHub\bifacialvf\bifacialvf\bonanza\CombinedData.pickle')

data_filtered = analyzeKeeton(data,plotflag = False)  # set plotflag = False if you don't want to plot the data increment


datamonth = analyzeMonth(data_filtered)



# In[]:
fig = plt.figure()
fig.set_size_inches(4, 2.5)
ax = fig.add_axes((0.15,0.15,0.75,0.7))
#plt.rc('font', family='sans-serif')
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
plt.rc('axes',labelsize=8)
width = 0.3
ax.bar(bifimonth.index-width/2, bifimonth.PR_meas_diff*100, width, label=r"$BG_E$,meas", color = 'steelblue')
ax.bar(bifimonth.index+width/2, bifimonth.BGE_model*100, width, label="$BG_E$,model", alpha = 0.5, color = 'darkorange')
plt.ylabel(r'$BG_E$,bifacial [%]')
plt.ylim([0, 15])
plt.title('Monthly Bifacial gain - Measured vs Modeled ',fontsize=9)
plt.legend(loc='upper left',fontsize=8)
plt.xticks(datamonth.index, datamonth.mmyy, rotation='horizontal')
fig.autofmt_xdate()
plt.show()

print('Cumulative bifacial model: +{:.2}%'.format(bifiyear['BGE_model']*100))
print('Total measured PR gain: +{:.2f}% .   Bifacial gain: +{:.2f}%'.format((psmyear.SilfabPR / psmyear.TrinaPR*100)-100,
      (psmyear.SilfabPR / psmyear.SilfabModelPR / psmyear.TrinaPR * psmyear.TrinaModelPR*100)-100))
# In[18]:
# compared with natural ground cover at Teague site, what is the benefit of the white ground cover at Keeton?


# print each row production relative to Teague for February and cumulative (Feb - Sept)
print('White fabric gain: +{:.2}% in Feb, +{:.2}% cumulative  '.format((psmmonth.SilfabKeetonPR[10] /  psmmonth.SilfabPR[10]*100)-100,(psmmonth.SilfabKeetonPR[10:].sum()/  psmmonth.SilfabPR[10:].sum()*100)-100)) 
