[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Build Status](https://travis-ci.org/NREL/bifacialvf.svg?branch=master)](https://travis-ci.org/NREL/bifacialvf)

# bifacialvf - Bifacial PV View Factor model
python, configuration factor model

Original code by Bill Marion
Python translation by Silvana Ayala
Updates by Chris Deline

Based on the publication: "A Practical Irradiance Model for Bifacial PV Modules"
B. Marion, S. MacAlpine, C. Deline, A. Asgharzadeh, F. Toor, D. Riley, J. Stein, C. Hansen
2017 IEEE Photovoltaic Specialists Conference, Washington DC, 2017
URL: https://www.nrel.gov/docs/fy17osti/67847.pdf

## Introduction

bifacialvf is a self-contained view factor (or configuration factor) model which
replicates a 5-row PV system of infinite extent perpendicular to the module
rows. The function returns the irradiance profile along the middle (interior)
row by default, but user interface options include `'first'`, `'interior'`,
`'last'`, and `'single'`. Single-axis tracking is supported, and hourly output
files based on TMY inputs are saved. Spatial nonuniformity is reported, with
multiple rear-facing irradiances collected on the back of each module row.

## Pre-requisites
This software is written for Python 2 or 3. NREL recommends [Anaconda Python](https://www.anaconda.com/download/).

## Install using pip
[bifacialvf](https://pypi.org/project/bifacialvf/) is at the Python Package Index (PyPI). Use pip to install the latest release in your conda environment or virtualenv:

    (myenv)$ pip install bifacialvf

### Install development mode from GitHub
For those interested in contributing to bifacialvf:

1. Clone the bifacialvf repository: `$ git clone https://github.com/NREL/bifacialvf.git bifacialvf-master`
2. Navigate to the repository directory where `setup.py` is located: `$ cd bifacialvf-master`
3. Install via pip in development mode: `$ pip install -e .`

## Usage

    import bifacialvf

    bifacialvf.simulate(inputTMY, outputfile, tilt, azm, clearance, rowspacing)
    (data, metadata) = bifacialvf.loadVFresults(outputfile)

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

`sun.py`: 
Solar position and irradiance-related helper files including
`hrSolarPos`, `perezComp`, `solarPos`, and `sunIncident`

`vf.py`:
View Factor helper files to help with configuration-factor calculation
1-axis tracking and self-shading calculation.
Subroutines include:
`getBackSurfaceIrradiances`, `getFrontSurfaceIrradiances`, `getGroundShadeFactors`,
`getSkyConfigurationFactors`, `trackingBFvaluescalculator`, `rowSpacing`
