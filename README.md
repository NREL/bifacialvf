![logo](docs/images_wiki/bifacialVF.png)

# bifacialvf - Bifacial PV View Factor model with Mismatch routines
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Build Status](https://travis-ci.org/NREL/bifacialvf.svg?branch=main)](https://travis-ci.org/NREL/bifacialvf)

python, configuration factor model, electrical model mismatch for bifacial modules.

Original bilinear interpolation code by Sara MacAlpine
Python translation & Updates by Silvana Ayala

bifacialvf
Original code by Bill Marion
Python translation by Silvana Ayala
Updates by Chris Deline

(Forthcoming) PVMismatch add-on
Original code by PVMismatch
For this version of bifacialvf_Mismach to work, PVMismatch must be installed (pip install pvmismatch)

Based on the publication:
Marion, B., MacAlpine, S., Deline, C., Asgharzadeh, A., Toor, F., Riley, D., … Hansen, C. (2017). A Practical Irradiance Model for Bifacial PV Modules: Preprint. In 44th IEEE Photovoltaic Specialists Conference. Washington, DC. https://www.nrel.gov/docs/fy17osti/67847.pdf. NREL/CP-5J00-67847

Bilinear Interpolation based on the publication:
De Soto, W., Klein, S. A., & Beckman, W. A. (2006). Improvement and validation of a model for photovoltaic array performance. Solar Energy, 80(1), 78–88. https://doi.org/10.1016/j.solener.2005.06.010

Marion, B., Rummel, S., & Anderberg, A. (2004). Current--voltage curve translation by bilinear interpolation. Progress in Photovoltaics: Research and Applications, 12(8), 593–607.


## Introduction


bifacialvf is a self-contained view factor (or configuration factor) model which
replicates a 5-row PV system of infinite extent perpendicular to the module
rows. The function returns the irradiance profile along the middle (interior)
row by default, but user interface options include `'first'`, `'interior'`,
`'last'`, and `'single'`. Single-axis tracking is supported, and hourly output
files based on TMY inputs are saved. Spatial nonuniformity is reported, with
multiple rear-facing irradiances collected on the back of each module row.

Bilinear interpolation code add-on to bifacialvf (description below) to pre-generate IV arrays and bifacial coefficients, and to examine the energy production with back side irradiance mismatch for either a portrait or landscape module.   
Included are IV curves and bifacial info for a Yingli (standard) module. 

## Pre-requisites
This software is written for Python 2 or 3. NREL recommends [Anaconda Python](https://www.anaconda.com/download/).

## Install using pip
[bifacialvf](https://pypi.org/project/bifacialvf/) is at the Python Package Index (PyPI). Use pip to install the latest release in your conda environment or virtualenv:

    (myenv)$ pip install bifacialvf

### Install development mode from GitHub
For those interested in contributing to bifacialvf:

1. Clone the bifacialvf repository: `$ git clone https://github.com/NREL/bifacialvf.git bifacialvf-main`
2. Navigate to the repository directory where `setup.py` is located: `$ cd bifacialvf-main`
3. Install via pip in development mode: `$ pip install -e .`

## Usage

    import bifacialvf

    bifacialvf.simulate(inputTMY, outputfile, tilt, azm, clearance, rowspacing)
    (data, metadata) = bifacialvf.loadVFresults(outputfile)
```

For more usage examples, see the Jupyter notebooks in \docs\

## Prerequisites

* [NumPy](https://www.numpy.org/)
* [pvlib python](https://pvlib-python.readthedocs.io/en/stable/)


## Main Functions

    bifacialvf.simulate(
        TMYtoread, writefiletitle,  beta, sazm, C=1, D=0.5,
        rowType = 'interior', transFactor=0.01, cellRows=6,
        PVfrontSurface='glass', PVbackSurface='glass',  albedo=0.62,
        tracking=False, backtrack=False, r2r=1.5, Cv=0.05, offset=0)


This is the main runfile.  Hourly TMY3 inputs are read, and an outputfile is saved with
a number of irradiance points along the module chord specified by `cellRows`.


    loadVFresults.loadVFresults(filename=None)

read in saved file from `bifacialvf.simulate`.  If no filename is passed, a tkinter GUI opens for file selection

## Subroutines

`LandscapeSingleHour.py`: 
BilinearInterpolation calculation for landscape modules

`PortraitSingleHour.py`: 
BilinearInterpolation calculation for Portrait modules 

`sun.py`: 
Solar position and irradiance-related helper files including
`hrSolarPos`, `perezComp`, `solarPos`, and `sunIncident`

`vf.py`:
View Factor helper files to help with configuration-factor calculation
1-axis tracking and self-shading calculation.
Subroutines include:
`getBackSurfaceIrradiances`, `getFrontSurfaceIrradiances`, `getGroundShadeFactors`,
`getSkyConfigurationFactors`, `trackingBFvaluescalculator`, `rowSpacing`
