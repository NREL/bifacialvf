[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# bifacialvf - Bifacial PV View Factor model with Mismatch routines
python, configuration factor model, electrical model mismatch for bifacial modules.

bifacialvf
Original code by Bill Marion
Python translation by Silvana Ayala
Updates by Chris Deline

Original bilinear interpolation code by Sara MacAlpine
Python translation & Updates by Silvana Ayala

(Forthcoming) PVMismatch add-on
Original code by PVMismatch
For this version of bifacialvf_Mismach to work, PVMismatch must be installed (pip install pvmismatch)

Based on the publication: "A Practical Irradiance Model for Bifacial PV Modules"
B. Marion, S. MacAlpine, C. Deline, A. Asgharzadeh, F. Toor, D. Riley, J. Stein, C. Hansen
2017 IEEE Photovoltaic Specialists Conference, Washington DC, 2017
URL: https://www.nrel.gov/docs/fy17osti/67847.pdf

And:


## Introduction

bifacialvf is a self-contained view factor (or configuration factor)
model which replicates a 5-row PV system of infinite extent perpendicular to the module rows. 
Single-axis tracking is supported, and hourly output files based on TMY inputs 
are saved.  Spatial nonuniformity is reported, with multiple rear-facing irradiances collected
on the back of each module row.

## Pre-requisites
This software is written for Anaconda Python 2.7 which can be downloaded here: https://www.anaconda.com/download/

## Install using pip

1. Clone or download the bifacialvf repository.
2. Navigate to the repository directory where setup.py is located: `cd bifacialvf-master`
3. Install via pip: `pip install .`
4. Alternate installation development mode: `pip install -e .`

## Usage

```
import bifacialvf

bifacialvf.simulate(TMYtoread, writefiletitle, beta, sazm, C, D, 
                     rowType, transFactor, cellRows, 
                     PVfrontSurface, PVbackSurface, albedo, 
                     tracking, backtrack, rtr, max_angle,
                     calculateBilInterpol, interpolA, IVArray, beta_voc_all,
                     m_all, bee_all, calculateDetailedMismatch, portraitorlandscape)
(data, metadata) = bifacialvf.loadVFresults(outputfile)
```
For more usage examples, see the Jupyter notebooks in \docs\

## Prerequisites

*none


## Main Functions
`bifacialvf.simulate(TMYtoread=None, writefiletitle=None, beta=0, sazm=180, C=0.5, D=None,
             rowType='interior', transFactor=0.01, cellRows=6, 
             PVfrontSurface='glass', PVbackSurface='glass', albedo=0.2,  
             tracking=False, backtrack=True, rtr=None, max_angle=45,
             calculateBilInterpol=False, interpolA=0.005, IVArray=None, beta_voc_all=None,
             m_all=None, bee_all=None, calculateDetailedMismatch=False, portraitorlandscape='landscape')`:  

This is the main runfile.  Hourly TMY3 inputs are read, and an outputfile is saved with
a number of irradiance points along the module chord specified by 'cellRows'.


`loadVFresults.loadVFresults(filename = none)` : 
read in saved file from bifacialvf.simulate.  If no filename is passed, a tkinter GUI opens for file selection

## Subroutines

`LandscapeSingleHour.py`: 
BilinearInterpolation calculation for landscape modules

`PortraitSingleHour.py`: 
BilinearInterpolation calculation for Portrait modules 

`sun.py`: 
Solar position and irradiance-related helper files including
hrSolarPos, perezComp, solarPos, and sunIncident

`vf.py`:
View Factor helper files to help with configuration-factor calculation
1-axis tracking and self-shading calculation.
Subroutines include:
getBackSurfaceIrradiances, getFrontSurfaceIrradiances, getGroundShadeFactors,
getSkyConfigurationFactors, trackingBFvaluescalculator, rowSpacing

;
