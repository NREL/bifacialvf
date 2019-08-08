# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:38:36 2019

@author: sayala
"""
from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.io as sio
import sys, os
sys.path.insert(0, 'BF_BifacialIrradiances')
from .BF_BifacialIrradiances.PortraitSingleHour import PortraitSingleHour    # For calculateBilInterpol
from .BF_BifacialIrradiances.LandscapeSingleHour import LandscapeSingleHour # For calculateBilInterpol
from pvmismatch import *  # this imports everything we need


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
        print("Warning: BilInterpolParams dictionary is None OR is wrongly defined. Using default values for Bilintear Interpolation routine")
        mat_contents = sio.loadmat(os.path.join('BF_BifacialIrradiances','BilinearInterpParams','IVArrayYingli.mat'))
        IVArray=mat_contents['IVArray']        
        mat_contents = sio.loadmat(os.path.join('BF_BifacialIrradiances','BilinearInterpParams','newBilinearParamsYingLi.mat'))
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
        cellsy = 12
        cellsx = 8
    else:
        stdpl = stdpl.transpose()
        cellsy = 8
        cellsx = 12
                                
    if sensorsy != cellsy:
        for i in range (0, cellsy):
            cellCenterPVM.append((i*sensorsy/cellsy+(i+1)*sensorsy/cellsy)/2)
    
    return cellCenterPVM, stdpl, cellsx, cellsy


def calculateVFPVMismatch(cellCenterPVM, stdpl, cellsy, cellsx, sensorsy, frontGTIrow, backGTIrow):
    r''' calls PVMismatch with all the pre-generated values on view factor.
    
    Example:
    PowerAveraged, PowerDetailed = calculateVFPVMismatch(cellCenterPVM, stdpl, cellsy, cellsx, sensorsy, frontGTIrow, backGTIrow)
    
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
            
        sunmatDetailed=[]
        sunmatAveraged=[]
        
        cellCenterValues_FrontPlusBack = cellCenterValFront+cellCenterValBack
        AveFront=cellCenterValFront.mean()                
        AveBack=cellCenterValBack.mean()
                 
        # Repeat to create a matrix to pass matrix.
        for j in range (0, len(cellCenterValues_FrontPlusBack)):
            sunmatDetailed.append([cellCenterValues_FrontPlusBack[j]/1000]*cellsx)
            
        for j in range (0, len(cellCenterValFront)):
            sunmatAveraged.append([(AveFront+AveBack)/1000]*cellsx)
    
            
        # ACtually do calculations
        pvsys.setSuns({0: {0: [sunmatAveraged, stdpl]}})
        PowerAveraged=pvsys.Pmp
        
        pvsys.setSuns({0: {0: [sunmatDetailed, stdpl]}})
        PowerDetailed=pvsys.Pmp

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