# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:38:36 2019

@author: sayala
"""
from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.io as sio
import sys, os
#sys.path.insert(0, 'BF_BifacialIrradiances')
from bifacialvf.BF_BifacialIrradiances.PortraitSingleHour import PortraitSingleHour    # For calculateBilInterpol
from bifacialvf.BF_BifacialIrradiances.LandscapeSingleHour import LandscapeSingleHour # For calculateBilInterpol
from pvmismatch import *  # this imports everything we need
import bifacialvf
import os
import pandas as pd

def setupforBilinearInterpolation(portraitorlandscape, sensorsy, BilInterpolParams=None):
    r'''Reads dictionary and assigns values. Also calcualtes the center location
    of the sensors in case the number of sensors is more than the number
    of cells-rows required by the bilinear interpolation routine (6 for landscape)
    
    Example:
    cellCenterBI, interpolA, IVArray, beta_voc_all, m_all, bee_all = setupforBilinearInterpolation(portraitorlandscape='landscape', sensorsy=100, BilInterpolParams=None)
   
    NOTE: portrait mode not working yet, will run everything on landscape!!
    #2DO #TODO #FIX
    '''
    
    try:
        interpolA = BilInterpolParams.interpolA
        IVArray = BilInterpolParams.IVArray
        beta_voc_all = BilInterpolParams.beta_voc_all
        m_all = BilInterpolParams.m_all
        bee_all = BilInterpolParams.bee_all
    except:
        # Try to load bilinear interpolation params dictionary from \BF_BifacialIrradiances\BilinearInterpParams
        print("Warning: BilInterpolParams dictionary is None OR is wrongly defined. Using default values for Bilintear Interpolation routine")
        currentdir =  os.path.abspath(os.path.dirname(__file__))

        mat_contents = sio.loadmat(os.path.join(currentdir,'BF_BifacialIrradiances','BilinearInterpParams','IVArrayYingli.mat'))
        IVArray=mat_contents['IVArray']        
        mat_contents = sio.loadmat(os.path.join(currentdir,'BF_BifacialIrradiances','BilinearInterpParams','newBilinearParamsYingLi.mat'))
        beta_voc_all=mat_contents['beta_voc_all']        
        m_all=mat_contents['m_all']
        bee_all=mat_contents['bee_all']
        interpolA = 0.005  # More accurate interpolation. Do 0.01 as an option.
        
    cellCenterBI=[]
    if portraitorlandscape=='portrait':
        print("Sorry! Bilintear Interpolation is not fully working for portrait mode. It will automatically run landscape!")
        portraitorlandscape = 'landscape'
        
    if portraitorlandscape == 'landscape':
        if sensorsy != 6:
            for i in range (0, 6):
                cellCenterBI.append((i*sensorsy/6.0+(i+1)*sensorsy/6.0)/2)
    
    return cellCenterBI, interpolA, IVArray, beta_voc_all, m_all, bee_all      

def setupforPVMismatch(portraitorlandscape, sensorsy):
    r''' Sets values for calling PVMismatch, for ladscape or portrait modes and 
    
    Example:
    cellCenterPVM, stdpl, cellsx, cellsy = setupforPVMismatch(portraitorlandscape='portrait', sensorsy=100):
    '''
    
    stdpl=np.array([[0,	23,	24,	47,	48,	71,	72,	95],
        [1,	22,	25,	46,	49,	70,	73,	94],
        [2,	21,	26,	45,	50,	69,	74,	93],
        [3,	20,	27,	44,	51,	68,	75,	92],
        [4,	19,	28,	43,	52,	67,	76,	91],
        [5,	18,	29,	42,	53,	66,	77,	90],
        [6,	17,	30,	41,	54,	65,	78,	89],
        [7,	16,	31,	40,	55,	64,	79,	88],
        [8,	15,	32,	39,	56,	63,	80,	87],
        [9,	14,	33,	38,	57,	62,	81,	86],
        [10,	13,	34,	37,	58,	61,	82,	85],
        [11,	12,	35,	36,	59,	60,	83,	84]])

    cellCenterPVM=[]
    
    if portraitorlandscape == 'portrait':                    
        cellsx = 8
        cellsy = 12
    if portraitorlandscape == 'landscape':
        stdpl = stdpl.transpose()
        cellsx = 12
        cellsy = 8
                                
    if sensorsy != cellsy:
        for i in range (0, cellsy):
            cellCenterPVM.append((i*sensorsy/cellsy+(i+1)*sensorsy/cellsy)/2)
    
    return cellCenterPVM, stdpl, cellsx, cellsy


def calculateVFPVMismatch(cellCenterPVM, stdpl, cellsx, cellsy, sensorsy, frontGTIrow, backGTIrow, bififactor=1.0, debug=False):
    r''' calls PVMismatch with all the pre-generated values on view factor.
    
    Inputs:
    bififactor: bifaciality factor of the module. Max 1.0. ALL Rear irradiance values saved include the bifi-factor.

    Example:
    PowerAveraged, PowerDetailed = calculateVFPVMismatch(cellCenterPVM, stdpl, cellsy, cellsx, sensorsy, frontGTIrow, backGTIrow, bififactor)
    
    '''
    
    if np.mean(frontGTIrow) < 1.0:
        PowerAveraged = 0
        PowerDetailed = 0
    else:                
        pvsys = pvsystem.PVsystem(numberStrs=1, numberMods=1)  
        # makes the system  # 1 module, in portrait mode. 
        
        if sensorsy != cellsy:                        
            cellCenterValFront= np.interp(cellCenterPVM, list(range(0,sensorsy)), frontGTIrow)
            cellCenterValBack= np.interp(cellCenterPVM, list(range(0,sensorsy)), backGTIrow)
        else:
            cellCenterValFront = frontGTIrow
            cellCenterValBack = backGTIrow

        cellCenterValues_FrontPlusBack = [(x+y*bififactor)/1000 for x,y in zip(cellCenterValFront,cellCenterValBack)]
        
        # New section from bifacial_radiance optimized
        G=np.array([cellCenterValues_FrontPlusBack]).transpose()
        H = np.ones([1,cellsx]) 
        array_det = np.dot(G,H) 
        array_avg = np.ones([cellsy,cellsx])*(np.mean(cellCenterValues_FrontPlusBack))        
        
        # Actually do calculations
        pvsys.setSuns({0: {0: [array_avg, stdpl]}})
        PowerAveraged=pvsys.Pmp
        
        pvsys.setSuns({0: {0: [array_det, stdpl]}})
        PowerDetailed=pvsys.Pmp        

        '''
        OLD Version
        
        cellCenterValues_FrontPlusBack = cellCenterValFront+cellCenterValBack
        AveFront=cellCenterValFront.mean()                
        AveBack=cellCenterValBack.mean()
        
        # Repeat to create a matrix to pass matrix.
        #Chris: I'm sure you can make this matrix in a prettier way than I can!
        for j in range (0, len(cellCenterValues_FrontPlusBack)):
            sunmatDetailed.append([cellCenterValues_FrontPlusBack[j]/1000]*cellsx)
            
        for j in range (0, len(cellCenterValFront)):
            sunmatAveraged.append([(AveFront+AveBack)/1000]*cellsx)
                        
        # ACtually do calculations
        pvsys.setSuns({0: {0: [sunmatAveraged, stdpl]}})
        PowerAveraged=pvsys.Pmp
        
        pvsys.setSuns({0: {0: [sunmatDetailed, stdpl]}})
        PowerDetailed=pvsys.Pmp
        '''
        
    if debug:          
        return PowerAveraged, PowerDetailed, array_avg, array_det
    else:
        return PowerAveraged, PowerDetailed
    
def calculateVFBilinearInterpolation(portraitorlandscape, sensorsy, cellCenterBI, interpolA, IVArray, beta_voc_all, m_all, bee_all, frontGTIrow, backGTIrow, Tamb, VWind):
    r''' calls BilinearInterpolationRoutine for the view factor results frontGTIrow and backGTIrow 
    
    Example:
    PowerAveraged, PowerDetailed = calculateVFBilinearInterpolation(portraitorlandscape, sensorsy, cellCenterBI, interpolA, IVArray, beta_voc_all, m_all, bee_all, frontGTIrow, backGTIrow, Tamb, VWind):

    '''    
    if portraitorlandscape=='portrait':

        if sensorsy != 6:                        
            cellCenterValFront= np.interp(cellCenterBI, list(range(0,sensorsy)), frontGTIrow)
            cellCenterValBack= np.interp(cellCenterBI, list(range(0,sensorsy)), backGTIrow)
        else:
            cellCenterValFront = frontGTIrow
            cellCenterValBack = backGTIrow
        
        [PmaxIdeal, PmaxUnmatched, PmaxAvg] = PortraitSingleHour(cellCenterValFront, cellCenterValBack, Tamb, VWind, 6, interpolA,IVArray,beta_voc_all,m_all,bee_all)

    if portraitorlandscape=='landscape':
        
        if sensorsy != 6:                        
            cellCenterValFront= np.interp(cellCenterBI, list(range(0,sensorsy)), frontGTIrow)
            cellCenterValBack= np.interp(cellCenterBI, list(range(0,sensorsy)), backGTIrow)
        else:
            cellCenterValFront = frontGTIrow
            cellCenterValBack = backGTIrow
            
        [PmaxIdeal, PmaxUnmatched, PmaxAvg] = LandscapeSingleHour(cellCenterValFront, cellCenterValBack, Tamb, VWind, 6, interpolA,IVArray,beta_voc_all,m_all,bee_all)

        
        return PmaxIdeal, PmaxUnmatched

def analyseVFResultsPVMismatch(filename, portraitorlandscape='portrait', bififactor=1.0, writefilename=None):
    '''
    Opens a finished bifacialVF results file with the metadata and the irradiance
    results for Front and back in format "No_1_RowFrontGTI", detects how many
    points were sampled along the panel ('sensorsy' variable), and given if the
    panel for PVMismatch is on portrait or landscape orientation it writes into
    an outputfile.
    
    If no writefilename is passed it uses the same filename of input adding a 
    '_PVMismatch.csv' ending
    
    Inputs:
    bififactor: bifaciality factor of the module. Max 1.0. ALL Rear irradiance values saved include the bifi-factor.
    portraitorlandscape: 'portrait' or 'landscape', for PVMismatch input
                      which defines the electrical interconnects inside the module. 

    Example:
    analyseVFResultsPVMismatch(filename='Output\test.csv', 
                               portraitorlandscape='portrait',
                               writefiletitle='Output\test.csv') 
                                #This will rewrite the input file!
    '''
    
    filename=r'C:\Users\sayala\Documents\GitHub\bifacialvf\bifacialvf\data\Output\test.csv'
    (data, metadata) = bifacialvf.loadVFresults(filename)

    # Checking to see if PVMismatch has already been run:
    if 'CalculatePVOutput (PVMismatch)' in metadata:
        if metadata['CalculatePVOutput (PVMismatch)'] == True:
            print("Selected File already has a PVMismatch Calculation in it. Exiting calculation!")
            return
        else:
            metadata['CalculatePVOutput (PVMismatch)'] = 'True'
    
    metadata['PortraitorLandscape'] = portraitorlandscape

    frontGTI = [col for col in data if col.endswith('RowFrontGTI')]
    backGTI = [col for col in data if col.endswith('RowBackGTI')]
    sensorsy=len(frontGTI)
    
    frontGTI=data[frontGTI]
    backGTI=data[backGTI]

    cellCenterPVM, stdpl, cellsx, cellsy = bifacialvf.analysis.setupforPVMismatch(portraitorlandscape=portraitorlandscape, sensorsy=sensorsy)

    PowerAveraged_all=[]
    PowerDetailed_all=[]
    
    for i in range (0,len(frontGTI)):
        PowerAveraged, PowerDetailed = bifacialvf.analysis.calculateVFPVMismatch(cellCenterPVM=cellCenterPVM, stdpl=stdpl, cellsx=cellsx, cellsy=cellsy, sensorsy=sensorsy, frontGTIrow=frontGTI.iloc[i], backGTIrow=backGTI.iloc[i], bififactor=bififactor)
        PowerAveraged_all.append(PowerAveraged)
        PowerDetailed_all.append(PowerDetailed)
    
    data['PVMismatch PowerAveraged']=PowerAveraged_all
    data['PVMismatch PowerDetailed']=PowerDetailed_all

    metadata['Orientation'] = portraitorlandscape
    metadata2=pd.Series(metadata).to_frame().T
    
    writefilename=(os.path.splitext(filename)[0])+'_PVMismatch.csv'
    metadata2.to_csv(writefilename,index=False)
    data.to_csv(writefilename, mode='a', index=False, header=True)
    
    print("The DC Power Mismatch loss for the year is of: {:.3f}%".format(100-data['PVMismatch PowerDetailed'].sum()*100/data['PVMismatch PowerAveraged'].sum()))
    

    