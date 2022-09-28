# ensure python3 compatible division and printing
from __future__ import division, print_function, absolute_import
from bifacialvf.bifacialvf import simulate, getEPW, readInputTMY  # main program
from bifacialvf.vf import getBackSurfaceIrradiances, getFrontSurfaceIrradiances, getGroundShadeFactors  # main subroutines
from bifacialvf.vf import getSkyConfigurationFactors, trackingBFvaluescalculator, rowSpacing # helper functions
from bifacialvf.sun import hrSolarPos, perezComp, solarPos, sunIncident # solar position and value
from bifacialvf.loadVFresults import loadVFresults # utility for reading result files
from bifacialvf.BF_BifacialIrradiances.LandscapeSingleHour import LandscapeSingleHour # For calculateBilInterpol
from bifacialvf.BF_BifacialIrradiances.PortraitSingleHour import PortraitSingleHour # For calculateBilInterpol

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
