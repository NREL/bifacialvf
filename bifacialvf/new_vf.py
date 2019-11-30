

# In[1]:


def getSkyConfigurationFactors(rowType, beta, C, D):
    """
    This method determines the sky configuration factors for points on the
    ground from the leading edge of one row of PV panels to the leading edge of
    the next row of PV panels behind it. This row-to-row dimension is divided
    into 100 ground segments and a sky configuration factor is returned for
    each ground segment. The sky configuration factor represents the fraction
    of the isotropic diffuse sky radiation (unobstructed) that is present on
    the ground when partially obstructed by the rows of PV panels. The
    equations follow that on pages in the notebook dated 8/12/2015. 8/20/2015            

    4/15/2016 Modifed for calculations other than just the interior rows. Row
    type is identified with the string `rowType`, with the possilbe values:

    * first = first row of the array
    * interior = interior row of array
    * last = last row of the array
    * single = a single row array

    Because the sky configuration factors may now be different depending on
    row, they are calculated for the row-to-row dimension to the rear of the
    leading module edge and to the front of the leading edge.

    Parameters
    ----------
    rowType : str
        "first", "interior", "last", or "single"
    beta : float
        Tilt from horizontal of the PV modules/panels (deg)
    C : float
        Ground clearance of PV panel (in PV module/panel slope lengths)            
    D : float
        Horizontal distance between rows of PV panels (in PV module/panel slope
        lengths)

    Returns
    -------
    rearSkyConfigFactors : array of size [100]
        Sky configuration factors to rear of leading PVmodule edge (decimal
        fraction)
    frontSkyConfigFactors : array of size [100]
        Sky configuration factors to rear of leading PVmodule edge (decimal
        fraction)

    Notes
    -----
    The horizontal distance between rows, `D`, is from the back edge of one row
    to the front edge of the next, and it is not the row-to-row spacing.
    """

    rearSkyConfigFactors = []
    frontSkyConfigFactors = []
    
    # Tilt from horizontal of the PV modules/panels, in radians
    beta = beta * DTOR
    # Vertical height of sloped PV panel (in PV panel slope lengths)                    
    h = math.sin(beta)
    # Horizontal distance from front of panel to rear of panel (in PV panel
    # slope lengths)
    x1 = math.cos(beta)
    rtr = D + x1  # Row-to-row distance (in PV panel slope lengths)

    # Forced fix for case of C = 0
    # FIXME: for some reason the Config Factors go from 1 to 2 and not 0 to 1.
    # TODO: investigate why this is happening in the code.
    if C==0:
        C=0.0000000001
        
    if C < 0:
        LOGGER.error(
            "Height is below ground level. Function GetSkyConfigurationFactors"
            " will continue but results might be unreliable")
        
    # Divide the row-to-row spacing into 100 intervals and calculate
    # configuration factors
    delta = rtr / 100.0

    if (rowType == "interior"):
        # Initialize horizontal dimension x to provide midpoint of intervals
        x = -delta / 2.0

        for i in range(0,100):        
            x += delta
            #  <--rtr=x1+D--><--rtr=x1+D--><--rtr=x1+D-->
            # |\            |\            |\            |\ 
            # | \ `         | \           | \          /| \
            # h  \   `      h  \          h  \       /  h  \
            # |   \     `   |   \         |   \    /    |   \
            # |_x1_\____D__`|_x1_\____D___|_x1_\_/_D____|_x1_\_
            # |               `   <------x-----/|
            # C                  `           /
            # |              angA   `      /  angB
            # *------------------------`-/---------------------
            #                          x
            # use ATAN2: 4-quadrant tangent instead of ATAN
            # check 2 rows away
            angA = math.atan2(h + C, (2.0 * rtr + x1 - x))
            angB = math.atan2(C, (2.0 * rtr - x))
            beta1 = max(angA, angB)
            # check 1 rows away
            angA = math.atan2(h + C, (rtr + x1 - x))
            angB = math.atan2(C, (rtr - x))
            beta2 = min(angA, angB)
            # check 0 rows away
            beta3 = max(angA, angB)
            beta4 = math.atan2(h + C, (x1 - x))
            beta5 = math.atan2(C, (-x))
            beta6 = math.atan2(h + C, (-D - x))
            sky1 =0; sky2 =0; sky3 =0
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2))
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4))
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6))
            skyAll = sky1 + sky2 + sky3

            # Save as arrays of values, same for both to the rear and front
            rearSkyConfigFactors.append(skyAll)
            frontSkyConfigFactors.append(skyAll)        
        # End of if "interior"

    elif (rowType == "first"):
    
        # RearSkyConfigFactors don't have a row in front, calculation of sky3
        # changed, beta6 = 180 degrees
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);

            angA = math.atan((h + C) / (rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.pi;

            sky1 = 0.0; sky2 = 0.0; sky3 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky1 + sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);    # Save as arrays of values                   
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        

        #  frontSkyConfigFactors don't have a row in front, calculation of sky3 included as part of revised sky2,
        #  beta 4 set to 180 degrees
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);

            angA = math.atan((h + C) / (rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.pi;

            sky1 = 0.0; sky2 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));

            skyAll = sky1 + sky2;
            frontSkyConfigFactors.append(skyAll);  # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        
        # End of if "first"

    elif (rowType == "last"):
    
        #  RearSkyConfigFactors don't have a row to the rear, combine sky1 into sky 2, set beta 3 = 0.0
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;

            beta3 = 0.0;

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.atan((h + C) / (-D - x));
            if (beta6 < 0.0):
                beta6 += math.pi;

            sky2 = 0.0; sky3 = 0.0;
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);     # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        

        #  FrontSkyConfigFactors have beta1 = 0.0
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);
            beta1 = 0.0;

            angA = math.atan((h + C) / (rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.atan((h + C) / (-D - x));
            if (beta6 < 0.0):
                beta6 += math.pi;

            sky1 = 0.0; sky2 = 0.0; sky3 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky1 + sky2 + sky3;
            frontSkyConfigFactors.append(skyAll);      # Save as arrays of values,
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        
        # End of if "last" row

    elif (rowType == "single"):
    
                    
        #  RearSkyConfigFactors don't have a row to the rear ir front, combine sky1 into sky 2, set beta 3 = 0.0,
        #  for sky3, beta6 = 180.0.
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;

            beta3 = 0.0;

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.pi;

            sky2 = 0.0; sky3 = 0.0;
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);     # Save as arrays of values                    
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        

        #  FrontSkyConfigFactors have only a row to the rear, combine sky3 into sky2, set beta1 = 0, beta4 = 180
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);
            beta1 = 0.0;

            angA = math.atan((h + C) / (rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.pi;

            sky1 = 0.0; sky2 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            skyAll = sky1 + sky2;
            frontSkyConfigFactors.append(skyAll);  # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        
        # End of if "single"
    else:
        print("ERROR: Incorrect row type not passed to function GetSkyConfigurationFactors ");

    return rearSkyConfigFactors, frontSkyConfigFactors;

# In[2]:
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


# In[5]:


def getSkyConfigurationFactors2(rowType, beta, C, D, pitch):
    """
    This method determines the sky configuration factors for points on the
    ground from the leading edge of one row of PV panels to the leading edge of
    the next row of PV panels behind it. This row-to-row dimension is divided
    into 100 ground segments and a sky configuration factor is returned for
    each ground segment. The sky configuration factor represents the fraction
    of the isotropic diffuse sky radiation (unobstructed) that is present on
    the ground when partially obstructed by the rows of PV panels. The
    equations follow that on pages in the notebook dated 8/12/2015. 8/20/2015            

    4/15/2016 Modifed for calculations other than just the interior rows. Row
    type is identified with the string `rowType`, with the possilbe values:

    * first = first row of the array
    * interior = interior row of array
    * last = last row of the array
    * single = a single row array

    Because the sky configuration factors may now be different depending on
    row, they are calculated for the row-to-row dimension to the rear of the
    leading module edge and to the front of the leading edge.

    Parameters
    ----------
    rowType : str
        "first", "interior", "last", or "single"
    beta : float
        Tilt from horizontal of the PV modules/panels (deg)
    C : float
        Ground clearance of PV panel (in PV module/panel slope lengths)            
    D : float
        Horizontal distance between rows of PV panels (in PV module/panel slope
        lengths)

    Returns
    -------
    rearSkyConfigFactors : array of size [100]
        Sky configuration factors to rear of leading PVmodule edge (decimal
        fraction)
    frontSkyConfigFactors : array of size [100]
        Sky configuration factors to rear of leading PVmodule edge (decimal
        fraction)

    Notes
    -----
    The horizontal distance between rows, `D`, is from the back edge of one row
    to the front edge of the next, and it is not the row-to-row spacing.
    """

    rearSkyConfigFactors = []
    frontSkyConfigFactors = []
    
    # Tilt from horizontal of the PV modules/panels, in radians
    beta = beta * DTOR
    # Vertical height of sloped PV panel (in PV panel slope lengths)                    
    h = beta.apply(math.sin)
    # Horizontal distance from front of panel to rear of panel (in PV panel
    # slope lengths)
    x1 = beta.apply(math.cos)

    # Forced fix for case of C = 0
    # FIXME: for some reason the Config Factors go from 1 to 2 and not 0 to 1.
    # TODO: investigate why this is happening in the code. #Sil: Did we fix this?        
    C[C==0] = 0.0000000001
    
    if not C[C<0].empty:    
        LOGGER.error("ERROR: Clearance height is below ground level."
                     " Function GetSkyConfigurationFactors "
                     " will continue but results might be unreliable")
        C[C<0] = 0.0000000001
        
    # Divide the row-to-row spacing into 100 intervals and calculate
    # configuration factors
    delta = pitch / 101.0
    x=np.linspace(delta, pitch-delta, num=100, endpoint=True)
    df_x1 = pd.DataFrame(zip(*[x1 for i in range(100)]))
    
    if (rowType == "interior"):
    # Initialize horizontal dimension x to provide midpoint of intervals
    
        #  <--rtr=x1+D--><--rtr=x1+D--><--rtr=x1+D-->
        # |\            |\            |\            |\ 
        # | \ `         | \           | \          /| \
        # h  \   `      h  \          h  \       /  h  \
        # |   \     `   |   \         |   \    /    |   \
        # |_x1_\____D__`|_x1_\____D___|_x1_\_/_D____|_x1_\_
        # |               `   <------x-----/|
        # C                  `           /
        # |              angA   `      /  angB
        # *------------------------`-/---------------------
        #                          x
        # use ATAN2: 4-quadrant tangent instead of ATAN
        # check 2 rows away

        #angA = math.atan2(h + C, (2.0 * pitch + x1 - x))
        # atan -- y/x .. (y, x)
        clear = h+C
        df_clear = pd.DataFrame(zip(*[clear for i in range(100)]))
        df = df_x1.subtract(x)  #x1-x
        df = df+pitch*2.0
        df_yx = df_clear/df
        angA = df_yx.applymap(math.atan)

        #angB = math.atan2(C, (2.0 * pitch - x))
        df_y = pd.DataFrame(zip(*[C for i in range(100)]))
        df = 2.0*pitch - x
        df_yx = df_y/df
        angB = df_yx.applymap(math.atan)
        
        beta1 = angA.where(angA > angB, angB).fillna(angA) # selects the max of each element-wise between both dataframes.
        
        ## check 1 rows away
        #angA = math.atan2(h + C, (pitch + x1 - x))
        df = df_x1.subtract(x)  #x1-x
        df = df+pitch
        df_yx = df_clear/df
        angA = df_yx.applymap(math.atan)

        #angB = math.atan2(C, (pitch - x))
        df = pitch - x
        df_yx = df_y/df
        angB = df_yx.applymap(math.atan)
        
        # beta2 = min(angA, angB)
        beta2 = angA.where(angA < angB, angB).fillna(angA) # selects the max of each element-wise between both dataframes.

        #beta3 = max(angA, angB)
        beta3 = angA.where(angA > angB, angB).fillna(angA) # selects the max of each element-wise between both dataframes.
        
        # check 0 rows away
        #beta4 = math.atan2(h + C, (x1 - x))
        df_yx = df_clear / df_x1.subtract(x)  #x1-x
        beta4 = df_yx.applymap(math.atan)
        
        #beta5 = math.atan2(C, (-x))
        df_yx = df_y / -x
        beta5 = df_yx.applymap(math.atan) 
        
        #beta6 = math.atan2(h + C, (-D - x))
        df_D = pd.DataFrame(zip(*[D for i in range(100)]))
        df_D = -df_D-x
        df_yx = df_clear / df_D
        beta6 = df_yx.applymap(math.atan)

        #SIL CONTINUE FROM HERE
        sky1 =0; sky2 =0; sky3 =0
        if (beta2 > beta1):
            sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2))
        if (beta4 > beta3):
            sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4))
        if (beta6 > beta5):
            sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6))
        skyAll = sky1 + sky2 + sky3

        # Save as arrays of values, same for both to the rear and front
        rearSkyConfigFactors.append(skyAll)
        frontSkyConfigFactors.append(skyAll)        
    # End of if "interior"

    elif (rowType == "first"):
    
        # RearSkyConfigFactors don't have a row in front, calculation of sky3
        # changed, beta6 = 180 degrees
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);

            angA = math.atan((h + C) / (pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.pi;

            sky1 = 0.0; sky2 = 0.0; sky3 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky1 + sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);    # Save as arrays of values                   
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        

        #  frontSkyConfigFactors don't have a row in front, calculation of sky3 included as part of revised sky2,
        #  beta 4 set to 180 degrees
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);

            angA = math.atan((h + C) / (pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.pi;

            sky1 = 0.0; sky2 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));

            skyAll = sky1 + sky2;
            frontSkyConfigFactors.append(skyAll);  # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        
        # End of if "first"

    elif (rowType == "last"):
    
        #  RearSkyConfigFactors don't have a row to the rear, combine sky1 into sky 2, set beta 3 = 0.0
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;

            beta3 = 0.0;

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.atan((h + C) / (-D - x));
            if (beta6 < 0.0):
                beta6 += math.pi;

            sky2 = 0.0; sky3 = 0.0;
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);     # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        

        #  FrontSkyConfigFactors have beta1 = 0.0
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);
            beta1 = 0.0;

            angA = math.atan((h + C) / (pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.atan((h + C) / (-D - x));
            if (beta6 < 0.0):
                beta6 += math.pi;

            sky1 = 0.0; sky2 = 0.0; sky3 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky1 + sky2 + sky3;
            frontSkyConfigFactors.append(skyAll);      # Save as arrays of values,
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        
        # End of if "last" row

    elif (rowType == "single"):
    
                    
        #  RearSkyConfigFactors don't have a row to the rear ir front, combine sky1 into sky 2, set beta 3 = 0.0,
        #  for sky3, beta6 = 180.0.
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;

            beta3 = 0.0;

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.pi;

            sky2 = 0.0; sky3 = 0.0;
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);     # Save as arrays of values                    
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        

        #  FrontSkyConfigFactors have only a row to the rear, combine sky3 into sky2, set beta1 = 0, beta4 = 180
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);
            beta1 = 0.0;

            angA = math.atan((h + C) / (pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.pi;

            sky1 = 0.0; sky2 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            skyAll = sky1 + sky2;
            frontSkyConfigFactors.append(skyAll);  # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        
        # End of if "single"
    else:
        print("ERROR: Incorrect row type not passed to function GetSkyConfigurationFactors ");

    return rearSkyConfigFactors, frontSkyConfigFactors;