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
import matplotlib.pyplot as plt, mpld3
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

def _interactive_load(title = None):
    # Tkinter file picker
    import Tkinter
    from tkFileDialog import askopenfilename
    root = Tkinter.Tk()
    root.withdraw() #Start interactive file input
    root.attributes("-topmost", True) #Bring window into foreground
    return askopenfilename(parent = root, title = title) #initialdir = data_dir

def add_sequential_counters(df):
    #
    minmonth = df.index.month[0] + df.index.year[0] * 12
    minday = df.index.dayofyear[0] + df.index.year[0] * 365
    df['sequential_month'] = (df.index.month + df.index.year * 12-minmonth)
    df['month'] = df.index.month
    df['sequential_day'] = (df.index.dayofyear + df.index.year * 365-minday)
    return df

# group by month values - define as separate function. 
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


def readKeeton():
    '''
    readKeeton():   routine to read raw data from Keeton system and save to the root directory as CombinedData.pickle
    '''
    
    def readData(directory):
        # subroutine to read different .csv files in the Bonanza (Keeton) raw data directory and compile into a
        # single dataframe
        os.chdir(directory)
        filenames = glob.glob('*.csv')
        data = pd.concat((pd.read_csv(os.path.join(directory,f),index_col='Date & Time',skipfooter = 1, engine='python') for f in filenames))
        
        #re-order the data in increasing chronological order
        #data.index = data['Date & Time']
        data.index = pd.to_datetime(data.index)
        data = data.sort_index()
        data = data.drop_duplicates()
        data = data[~data.index.duplicated(keep='first')] #de-dupe where you have multiple instances
        data = data.resample('15T').mean()
        
        return data
    
    directory1 = r'C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon\Bonanza bifacial raw data'
    data = readData(directory1)
    data.to_csv('BonanzaBifacial.csv')
    
    
    directory2 = r'C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon\Bonanza Trina raw data'
    data2 = readData(directory2)
    data2.to_csv('BonanzaTrina.csv')
    
    directory = r'C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon'
    data3 = pd.concat([data,data2], axis=1, join='inner')
    data3.to_pickle(os.path.join(directory,'CombinedData.pickle'))
    #data3.to_csv(os.path.join(directory,'CombinedData.csv'))
    
    return data3

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


def readTeague(directory):
    # subroutine to read different .csv files in the Teague raw directory and compile into a
    # single dataframe
    os.chdir(directory)
    filenames = glob.glob('*T1.csv')
    data = pd.concat((pd.read_csv(os.path.join(directory,f),index_col='Date & Time',skipfooter = 1, engine='python' ) for f in filenames))
    
    filenames2 = glob.glob('*T2.csv')
    data2 = pd.concat((pd.read_csv(os.path.join(directory,f),index_col='Date & Time',skipfooter = 1, engine='python' ) for f in filenames2))
    
        
    #re-order the data in increasing chronological order
    #data.index = data['Date & Time']
    data.index = pd.to_datetime(data.index)
    data = data.sort_index()
    data = data.drop_duplicates()
    data = data[~data.index.duplicated(keep='first')] #de-dupe where you have multiple instances
    data = data.resample('15T').mean()
    
    data2.index = pd.to_datetime(data2.index)
    data2 = data2.sort_index()
    data2 = data2.drop_duplicates()
    data2 = data2[~data2.index.duplicated(keep='first')] #de-dupe where you have multiple instances
    data2 = data2.resample('15T').mean()    
    
    data3 = pd.concat([data,data2], axis=1)
    
    return data3

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


    
def analyzeTeague(df,plotflag = True):
    '''
    analyzeTeague:  import dataframe df and append normalized bifacial series.  Return df_filter (filtered data) 
    
    Input: df:  dataframe including 'Date & Time', 'Inverter A [kW]' thru 'Inverter F [kW]'
    Output: df_filter: data filtered for date > 8-16-2017
    
    '''
    
    df.index = pd.to_datetime(df['Date & Time'])
    # inverter A, B, C: Sunpreme Bifacial.  
    # inverter D, E, F: Sunpreme Bifacial
    # ignore data when (B-A)/B > 10%
    
    '''
    system details, Teague system:
       Silfab system: Inverter A: 6 strings of 22 (37.62) Inverter B: 5 strings of 22 (31.35). Inverter C: 108 modules (30.78)
      Inverter D: 6 strings of 22 (37.62) Inverter E: 5 strings of 22 (31.35). Inverter F: 108 modules (30.78)
    
      NOTES: Inverter bi4 underperforming relative to others (String outage?).  Neglect it from combined analysis
    '''
    
    df['bi1'] = df['INVERTER A [kW]'] / 37.62
    df['bi2'] = df['INVERTER B [kW]'] / 31.35
    df['bi3'] = df['INVERTER C [kW]'] / 30.78
    
    df['bi4'] = df['INVERTER D [kW]'] / 37.62
    df['bi5'] = df['INVERTER E [kW]'] / 31.35
    df['bi6'] = df['INVERTER F [kW]'] / 30.78
    
    df['bifacial'] = (df['INVERTER A [kW]'] + df['INVERTER B [kW]'] + df['INVERTER C [kW]']  
        + df['INVERTER E [kW]'] + df['INVERTER F [kW]']) / (37.62 + 2*31.35 + 2*30.78)
    
         
    fig = plt.figure()
    plt.plot(df.index,df.bi1,alpha = 0.4);plt.plot(df.index,df.bi4,alpha = 0.4);
    plt.legend(); plt.xlabel('measdatetime'); plt.ylabel('Power (normalized)')
    
    
    # monthly cumulative production for each string
    df = add_sequential_counters(df)  #add ['sequential_month']
    
    df_filter = filterdf(df)
    
    df_filter.index = df_filter.index.tz_localize('Etc/GMT+7') # set to pacific standard time
    
    # PLOTTING SECTION
    if plotflag == True:

            
        # plot 5 bifacial inverters, normalized, filtered
        fig = plt.figure()
        plt.rcParams['figure.figsize'] = (6, 4)
        plt.plot(df_filter.index,df_filter.bi1,label = 'InverterA', alpha = 0.5)
        plt.plot(df_filter.index,df_filter.bi2,label = 'InverterB', alpha = 0.4)
        plt.plot(df_filter.index,df_filter.bi3,label = 'InverterC', alpha = 0.5)
        plt.plot(df_filter.index,df_filter.bi5,label = 'InverterE', alpha = 0.4)
        plt.plot(df_filter.index,df_filter.bi6,label = 'InverterF', alpha = 0.5)
        #plt.xlim(['2017-06-10','2017-06-17'])
        plt.xlim(['2018-02-15','2018-02-25'])
        fig.autofmt_xdate()
        plt.legend(loc='best')
        plt.title('Silfab inverters')
        plt.show()
        
    return df_filter

def analyzeDay(df_filter):
    '''
    analyzeDay:  aggregate daily and return for filtered Teague or Keeton data
    
    '''
    mindate = df_filter.index[0]
    dfday = df_filter.groupby(df_filter.sequential_day).sum() / 4  #15 minute data from August 18 onward.
    dfday['sequential_day'] = dfday.index
    from dateutil.relativedelta import relativedelta
    dfday['mmdd'] = dfday['sequential_day'].apply(lambda x: (mindate + relativedelta(days=+x)).strftime("%m-%d"))

    return dfday


 
    
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

# In[13]:

def modelTrackedIrradiance(bonanza_TWC_data_file):
    '''
    From analysis_v2 - look at modeled PR and satellite weather data. Load Bonanza_TWC from local csv
    Hourly modeled data saved as psm (from NSRDB physical solar model)
    
    input: bonanza_TWC_data_file:  directory filename where weather channel data .csv data exists.  
        relevant fields of csv file from TheWeatherChannel: DateHrGmt, Latitude, Longitude, DiffuseHorizontalRadiationWsqm, DirectNormalIrradianceWsqm, DownwardSolarRadiationWsqm
    
    return: psm with hourly irradiance 'poa_irradiance'
    '''
       
    # work with weather data
    # PVSC-specific analysis:
    #psm = pd.read_csv(os.path.join(r'C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon\Teague bifacial','Bonanza_TWC_PVSC.csv')  )
    if bonanza_TWC_data_file is None:
        bonanza_TWC_data_file =  _interactive_load(title = 'Select Bonanza weather data .csv')       
    
    psm = pd.read_csv( bonanza_TWC_data_file  )
    
    psm.index = pd.to_datetime(psm.DateHrGmt)
    psm.index = psm.index.tz_localize('UTC')
    # calculate solar position and tracker tilt/azimuth based on DateHrGmt to get POA
    #pvlib.solarposition.get_solarposition(time, latitude, longitude)
    sun = pvlib.solarposition.get_solarposition(time = psm.index,latitude = psm.Latitude[0], longitude = psm.Longitude[0] )
    #singleaxis(apparent_zenith, apparent_azimuth,axis_tilt=0, axis_azimuth=0, max_angle=90, backtrack=True, gcr=2.0/7.0):
    tracker = pvlib.tracking.singleaxis(sun.zenith, sun.azimuth, max_angle = 60, gcr = 0.35, backtrack = True)
    # fill na's with zero
    tracker.fillna(value = 0, inplace = True)
    # calculate extraterrestrial  for the hay-davies model
    #dni_extra = pvlib.irradiance.extraradiation(psm.index.dayofyear)
    dni_extra = pvlib.irradiance.extraradiation(psm.index)
    dni_extra = pd.Series(dni_extra, index=psm.index)
    #haydavies(surface_tilt, surface_azimuth, dhi, dni, dni_extra, solar_zenith=None, solar_azimuth=None, projection_ratio=None):
    #perez(surface_tilt, surface_azimuth, dhi, dni, dni_extra, solar_zenith, solar_azimuth, airmass,  model='allsitescomposite1990', return_components=False):
    psm['sun_zenith'] = sun.zenith
    psm['sun_azimuth'] = sun.azimuth
    psm['surface_tilt'] = tracker.surface_tilt
    psm['surface_azimuth'] = tracker.surface_azimuth
    psm['dhi'] = psm['DiffuseHorizontalRadiationWsqm']
    psm['dni'] = psm['DirectNormalIrradianceWsqm']
    psm['ghi'] = psm['DownwardSolarRadiationWsqm']
    
    # calculate tilted diffuse
    #psm['poa_irradiance'] = pvlib.irradiance.total_irrad(psm.surface_tilt, psm.surface_azimuth,
    #                psm.sun_zenith, psm.sun_azimuth,
    #                psm.dni, psm.ghi, psm.dhi, model='haydavies', dni_extra = dni_extra  )
    
    total_irradiance = pvlib.irradiance.total_irrad(
        psm.surface_tilt,
        psm.surface_azimuth,
        psm['sun_zenith'],
        psm['sun_azimuth'],
        psm['dni'],
        psm['ghi'],
        psm['dhi'],
        model='haydavies',
        dni_extra=dni_extra)
    psm['poa_irradiance'] = total_irradiance['poa_global']
    return psm

# In[14]
def modelSilfabTrina(psm):
    '''
    model bifacial (silfab) and monofacial (trina) modules based on psm dataframe.  
    
    input parameter:  'psm' dataframe with relevant fields: poa_irradiance, WindSpeedKph, SurfaceTemperatureCelcius
       (poa_irradiance calculated from modelTrackedIrradiance)
    
    Returns:  psm dataframe of hourly satellite data with appended fields 'silfab_model' and 'trina_model'
    
    '''    
    # create modeled values based on CEC 5 parameter model parameters for Trina TSM-300D14.002 x33.3 kW and Silfab SLAx285 x 31.35 or 37.62 kW
    def findCEC(moduletext):
        if moduletext == None:
            moduletext = ''
        #return a list of CEC module names in the SAM database
        #moduletext: module text to search for.  cast as lowercase
        cec_modules = pvlib.pvsystem.retrieve_sam(path = r"C:\SAM\2017.9.5\libraries\CEC Modules.csv")
        #cec_modules = pvlib.pvsystem.retrieve_sam(path = r"C:\SAM\2017.1.17\libraries\CEC Modules.csv")
        modules = [col for col in cec_modules.columns if moduletext.lower() in col.lower()]
        print modules
        return cec_modules
    
    
    # create modeled values based on CEC 5 parameter model parameters for  Silfab SLAx285 x 31.35 or 37.62 kW
    #cec_modules = findCEC('trina')  # example of search by module type
    cec_modules = pvlib.pvsystem.retrieve_sam(path = r"C:\SAM\2017.9.5\libraries\CEC Modules.csv")
    
    silfab = cec_modules['Silfab_SLAX_285_Clear']  
    trina = cec_modules['Trina_Solar_TSM_300PD14_0x2']
    system_modules = pd.concat([trina,silfab],axis = 1)
    
    
    for col in system_modules:
        # calculate 5-parameter values for each point in time
        module_params = system_modules[col]
        #model_out = pd.DataFrame(index = data.measdatetime)
        # Calculate temp_cell from POA,  Windspeed and Tamb
        temp_cell = pvlib.pvsystem.sapm_celltemp(psm['poa_irradiance'],
                                                 psm['WindSpeedKph']*1000/3600,
                                                 psm['SurfaceTemperatureCelsius'])
        
        (photocurrent, saturation_current, resistance_series, 
        resistance_shunt, nNsVth) = pvlib.pvsystem.calcparams_desoto(psm['poa_irradiance'],
                          
                                         temp_cell = temp_cell['temp_cell'],               
                                         alpha_isc=module_params['alpha_sc'],
                                         module_parameters=module_params,
                                         EgRef=1.121,
                                         dEgdT=-0.0002677) 
    
        dc = pvlib.pvsystem.singlediode(photocurrent, saturation_current,
                                                      resistance_series,
                                                      resistance_shunt,
                                                      nNsVth)
        # apply DC losses.  SAM: 5% soiling, 4.4% DC loss, 1% AC loss
        # Total loss: 100% *(1- product of (1- loss/100%))
        effLoss = 0.14
        # inverter size: Solectria 36.7kW relative to 33kWdc string, scaled to 300W module power
        ac = pvlib.pvsystem.pvwatts_ac(pdc = dc['p_mp']*(1-effLoss),pdc0 = 330, eta_inv_nom = 0.96  ) 
    
        # join the power back to original dataframe
        ac.fillna(value = 0, inplace = True)
        psm = psm.join(ac)
    
        psm.rename(columns = {'p_mp':col+'_modeled'}, inplace = True)
    
    
    
    # remove modeled data if you need to re-run above
    #psm.drop(['Trina_Solar_TSM_300PD14_0x2_modeled','Silfab_SLAX_285_Clear_modeled'],inplace=True,axis=1)
            
    # Switch back to local wall time for index.  Merge inverter strings B and E (averaged over 1 hour)
    psm.index = psm.index.tz_convert('Etc/GMT+7')
    psm.rename(columns={'Silfab_SLAX_285_Clear_modeled':'silfab_model'}, inplace=True)
    psm.rename(columns={'Trina_Solar_TSM_300PD14_0x2_modeled':'trina_model'}, inplace=True)
  
    return psm




# In[1]:  START ANALYSIS
#  RAW 15 MINUTE DATA ANALYSIS OF TEAGUE DATA - saved as df

'''
re-load raw Teague measured data (if necessary). otherwise just load TeagueBifacial.csv
'''
directory = r'C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon'
#df = readTeague(os.path.join(directory,'raw'))  #  Only necessary with new data increment.

# save and re-load the raw csv joined file:
#df.to_csv(os.path.join(directory,'TeagueBifacial.csv'))
#df = pd.read_csv(os.path.join(directory,'TeagueBifacial_PVSC.csv'))
df = pd.read_csv(os.path.join(directory,'Teague bifacial','TeagueBifacial.csv'))

# Append normalized bifacial results to 'bifacial' and 'bi1' through 'bi6' keys.
# plot results (optional)
df_filter = analyzeTeague(df,plotflag = False) # use plotflag = True if you want to plot the data increment

# In[12]:  Monthly and daily aggregate analysis of Teague
dfmonth = analyzeMonth(df_filter)
dfday = analyzeDay(df_filter)



# In[]:   Raw 15 minute data analysis of Keeton
'''
re-load raw Keeton data using readKeeton().  otherwise just load CombinedData.pickle
'''
#data = readKeeton()
data = pd.read_pickle(r'C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon\CombinedData.pickle')
data_filtered = analyzeKeeton(data,plotflag = False)  # set plotflag = False if you don't want to plot the data increment
datamonth = analyzeMonth(data_filtered)
dataday = analyzeDay(data_filtered)
# In[]

# save and re-load weather data for the site.  Run performance model
bonanza_TWC_data_file = os.path.join(r'C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon\Teague bifacial','Bonanza_TWC.csv')
psm = modelTrackedIrradiance(bonanza_TWC_data_file)
psm = modelSilfabTrina(psm)  # return psm with additional keys: trina_model, silfab_model

# Startup problems before 6/8/18. Match to 8/16 start date of Keeton site.
psm_mask = (psm.index > '2017-8-16') 
psm_filtered = psm[psm_mask]#.dropna()
psm_filtered = filterdf(psm)


# In[]:  filter psm data
# monthly cumulative production for each string
mindate = psm_filtered.index[0]
minmonth = psm_filtered.index.month[0] + psm_filtered.index.year[0] * 12
minday = psm_filtered.index.dayofyear[0] + psm_filtered.index.year[0] * 365
psm_filtered['sequential_month'] = (psm_filtered.index.month + psm_filtered.index.year * 12-minmonth)
psm_filtered['month'] = psm_filtered.index.month
psm_filtered['sequential_day'] = (psm_filtered.index.dayofyear + psm_filtered.index.year * 365-minday)


# Join resampled Teague and Keeton bifacial data to psm

psm_filtered = psm_filtered.join(df_filter['bifacial'].resample('1H').mean()*285)  # normalize Teague to 285W bifacial nameplate power
psm_filtered.rename(columns={'bifacial':'silfab_meas'}, inplace=True)

psm_filtered = psm_filtered.join(data_filtered['bifacial'].resample('1H').mean()*285)  # normalize Keeton to 285W bifacial nameplate power
psm_filtered.rename(columns={'bifacial':'silfab_Keeton_meas'}, inplace=True)

psm_filtered = psm_filtered.join(data_filtered['monofacial'].resample('1H').mean()*300)  # normalize Keeton monofacial to 300W Trina mono nameplate power
psm_filtered.rename(columns={'monofacial':'trina_meas'}, inplace=True)

# In[16]:
# Save combined merged document with psm model and Teague measured.
psm_filtered.to_csv(os.path.join(directory,'Bonanza_Merge.csv'))   
psm_filtered.to_pickle(os.path.join(directory,'Bonanza_Merge.pickle')) 
psm_filtered = pd.read_pickle(os.path.join(directory,'Bonanza_Merge.pickle'))

# In[]   : 

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


# group by daily values - define as separate function. 


#sum cumulative, monthly and daily values
psmmonth = psm_filtered.groupby(psm_filtered.sequential_month).sum()
psmmonth = add_month_values(psmmonth,mindate = psm_filtered.index[0])  #defined above
psmmonth = addPR(psmmonth)

#psmday = psm_filtered.groupby(psm_filtered.sequential_day).sum()
#psmday = add_month_values(psmday)  #defined above
psmday = analyzeDay(psm_filtered)  #defined above
psmday = addPR(psmday)
psmyear = psm_filtered.sum()
psmyear = addPR(psmyear)  # Teague hourly data + modeled




# plot measured PR of Keeton and Teague systems. Use final print size.
# looks like there is a few percent gain for InvE and InvF starting in February.
fig = plt.figure()
fig.set_size_inches(3.5, 2)
ax = fig.add_axes((0.12,0.15,0.75,0.7))
#plt.rc('font', family='sans-serif')
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
plt.rc('axes',labelsize=8)
width = 0.18
plt.bar(psmmonth.index-width*3/2, psmmonth.SilfabPR,width = width, label="Teague", alpha = 0.5)
#plt.bar(datamonth.index-width/2, datamonth.bi1/4/psmmonth.Yr, width = width, label="Keeton_invD", alpha = 0.5)
plt.bar(psmmonth.index+width/2, psmmonth.SilfabKeetonPR, width = width, label="Keeton", alpha = 0.5)
#plt.bar(datamonth.index+width*3/2, datamonth.bi3/4/psmmonth.Yr, width = width, label="Keeton_invF", alpha = 0.5)
plt.plot(psmmonth.index, psmmonth.SilfabModelPR, label="modeled", color = 'black', linewidth = 1)

plt.ylabel('AC Performance Ratio')
plt.ylim([0, 1])
plt.title('Bifacial Teague PR meas:  {:.2f}, model: {:.2f}. Keeton (white fabric): {:.2f}  '.format(psmyear.SilfabPR , psmyear.SilfabModelPR, psmyear.SilfabKeetonPR),fontsize=10)
plt.legend(loc='lower right',fontsize=8)
plt.xticks(psmmonth.index, psmmonth.mmyy, rotation='horizontal')
fig.autofmt_xdate()
plt.show()

# run monofacial comparison based on Teague measured, Keeton monofacial
fig = plt.figure()
fig.set_size_inches(4, 2.5)
ax = fig.add_axes((0.12,0.15,0.75,0.7))
#plt.rc('font', family='sans-serif')
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
plt.rc('axes',labelsize=8)
plt.bar(psmmonth.index, psmmonth.SilfabPR, label="bifacial", color = 'steelblue')
plt.bar(psmmonth.index, psmmonth.TrinaPR, label="monofacial", color = 'silver')
plt.ylabel('AC Performance Ratio')
plt.ylim([0, 1])
plt.title('Cumulative bifacial production:  +{:.1f}%  '.format((psmyear.SilfabPR / psmyear.TrinaPR*100)-100),fontsize=9)
plt.legend(loc='lower right',fontsize=8)
plt.xticks(psmmonth.index, psmmonth.mmyy, rotation='horizontal')
fig.autofmt_xdate()
plt.show()

print('Total measured PR gain: +{:.2f}% .   Bifacial gain: +{:.2f}%'.format((psmyear.SilfabPR / psmyear.TrinaPR*100)-100,
      (psmyear.SilfabPR / psmyear.SilfabModelPR / psmyear.TrinaPR * psmyear.TrinaModelPR*100)-100))

# Add white fabric, starting in February
fig = plt.figure()
fig.set_size_inches(4, 2.5)
ax = fig.add_axes((0.12,0.15,0.75,0.7))
#plt.rc('font', family='sans-serif')
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
plt.rc('axes',labelsize=8)
plt.bar(psmmonth.index[7:], psmmonth.SilfabKeetonPR[7:], label="whiteFabric", color = 'orange')
plt.bar(psmmonth.index, psmmonth.SilfabPR, label="bifacial", color = 'steelblue')
plt.bar(psmmonth.index, psmmonth.TrinaPR, label="monofacial", color = 'silver')
plt.ylabel('AC Performance Ratio')
plt.ylim([0, 1])
plt.title('Cumulative bifacial production:  +{:.1f}%  (+ {:.1f}% w/ fabric)'.format((psmyear.SilfabPR / psmyear.TrinaPR*100)-100,
          (psmyear.SilfabKeetonPR / psmyear.TrinaPR*100)-100),fontsize=9)
plt.legend(loc='lower right',fontsize=8)
plt.xticks(psmmonth.index, psmmonth.mmyy, rotation='horizontal')
fig.autofmt_xdate()
plt.show()

print('Total measured PR gain: +{:.2f}% .   Bifacial gain: +{:.2f}%'.format((psmyear.SilfabPR / psmyear.TrinaPR*100)-100,
      (psmyear.SilfabPR / psmyear.SilfabModelPR / psmyear.TrinaPR * psmyear.TrinaModelPR*100)-100))


# In[17]:
## 
#  Look at modeled bifacial gain for Silfab interior row.  Pass in variables directly from data_filtered (Keeton)
## 

# 
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
writefiletitle = 'Silfab_model_1axis.csv'       
gcr = 0.35;  rtr = 1.0/gcr
C = 0.75  # clearance height: 2m table, 1.5m hub height

# NOTE: this uses simulateTWC below which needs to be read into memory.
os.chdir(r'C:\Users\cdeline\Documents\!Work\Bifacial\Irradiance monitoring\Oregon')
from PVSCpaper import simulateTWC
# Albedo range: 0.15 - 0.2 ?
simulateTWC(TWC, writefiletitle,  beta = 0, sazm = 180, C = C, D = None,
             rowType = 'interior', transFactor = 0.0, cellRows = 6, 
             PVfrontSurface = 'glass', PVbackSurface = 'glass',  albedo = 0.2,  
             tracking = True, backtrack = True, rtr = rtr, Cv = None, offset = 0, max_angle = 45,
             lat = 42.07823, lng = -121.27869, tz = -7, name = 'BONANZA', dataInterval= 60,
             maskStartEnd = False, datestart = '2000-11-2 00:00:00', dateend = '3000-11-24 23:59:59')

from bifacialvf import loadVFresults
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

bifimonth = model.groupby(model.sequential_month).sum()
# add modeled bifi gain and measured PR / model PR difference
bifimonth['BGE_model'] = bifimonth.GTIBackavg / bifimonth.GTIFrontavg * 0.95  # Silfab bifaciality: 0.95?
bifimonth['PR_meas_diff'] = (psmmonth.SilfabPR / psmmonth.SilfabModelPR * psmmonth.TrinaModelPR / psmmonth.TrinaPR)-1
bifimonth = add_month_values(bifimonth,mindate)

bifiyear = model.sum()
bifiyear['BGE_model'] = bifiyear.GTIBackavg / bifiyear.GTIFrontavg * 0.95  # Silfab bifaciality: 0.95?
bifiyear['PR_meas_diff'] = (psmyear.SilfabPR / psmyear.SilfabModelPR * psmyear.TrinaModelPR / psmyear.TrinaPR)-1

# plot the modeled bifacial gain vs. monthly measured PR difference between Teague (psm bifacial) and Keeton (data monofacial)

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
