#!/usr/bin/env python
# coding: utf-8

# In[1]:


C=0.2
D=1.5
rowType = 'interior'

[rearSkyConfigFactors2, frontSkyConfigFactors2] = getSkyConfigurationFactors(rowType, tilt, C, D)       ## Sky configuration factors are the same for all times, only based on geometry and row type

# In[3]:

beta = myTMY3['trackingdata_surface_tilt']
C = myTMY3['C']
D = myTMY3['D']

# In[4]:
[rearSkyConfigFactors2, frontSkyConfigFactors2] = getSkyConfigurationFactors(rowType, myTMY3['trackingdata_surface_tilt'], 
                                                                              myTMY3['C'], myTMY3['D'])       ## Sky configuration factors are the same for all times, only based on geometry and row type


# In[ ]:




