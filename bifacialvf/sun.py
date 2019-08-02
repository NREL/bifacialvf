#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sun module - solar helper files for bifacial-viewfactor
    
@author Bill Marion
@translated to python by sayala 06/09/17

"""

from __future__ import division, print_function # ensure python3 compatible division and printing
import math


def aOIcorrection(n2, inc):
        
            # This method calculates the correction factor for angle of incidence, a factor from 
            # zero to one (no correction) using the air-glass model of Sjerps-Koomen et al in "A Simple
            # Model for PV Module Reflection Losses under Field Conditons in Solar Energy
            # 1997; 57:421-432. 5/27/2014
            #
            # Variables passed to the method:
            #      n2 = index of refraction for glazing (1.526 for glass, 1.4 for Tefzel, 1.3 for AR glass)
            #      inc = incident angle in radians

            cor = -9999.0;   # Only calcutates valid value for 0 <= inc <= 90 degrees
            r0 = math.pow((n2 - 1.0) / (n2 + 1), 2);	# Reflectance at normal incidence, Beckman p217           
            if (inc == 0):
            
                cor = 1.0;		# Relative to normal incidence
            
            elif (inc > 0.0 and inc <= math.pi / 2.0):
            
                refrAng = math.asin(math.sin(inc) / n2);# Refracted angle
                r1 = (math.pow(math.sin(refrAng - inc), 2.0) /
                    math.pow(math.sin(refrAng + inc), 2.0))
                r2 = (math.pow(math.tan(refrAng - inc), 2.0) /
                    math.pow(math.tan(refrAng + inc), 2.0))
                cor = 1.0 - 0.5 * (r1 + r2);
                cor /= 1.0 - r0;		# Relative to normal incidence               
            
            return cor;
           # End of AOIcorrection

def hrSolarPos( year, month, day, hour, lat, lng, tz ):
 
		
            #This method is used to determine the solar position of hourly data by determining 
            #the	appropriate minute of the preceding hour that should be used to determine
            #an average solar position for the hour. Makes calls to method SolarPos.
            #List of Parameters Passed to Method:
            #    year   = year (e.g. 1986)
            #    month  = month of year (e.g. 1=Jan)
            #    day    = day of month
            #    hour   = hour at end of hour, local standard time, (1-24, or 0-23)
            #    at    = latitude in degrees, north positive
            #    ng    = longitude in degrees, east positive
            #    tz     = time zone, west longitudes negative
            #List of Out Parameters:
            #azm = sun azimuth in radians, measured east from north, 0 to 2*pi
            #	zen = sun zenith in radians, 0 to pi
            #	elv = sun elevation in radians, -pi/2 to pi/2
            #	dec = sun declination in radians
            #	sunrise = in local standard time (hrs), not corrected for refraction
            #	sunset = in local standard time (hrs), not corrected for refraction
            #	Eo = eccentricity correction factor
            #	tst = true solar time (hrs)
            #	Returns:
            #	suntime = minutes in hour that sun is up        4/1/03  */
               

   DTOR=math.pi/180;

   # Make call to SolarPos to get sunrise/sunset using hour=12 and minute=0.0
   azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, 12, 0.0, lat, lng, tz )
   
   # Rename sunrise and sunset to "positive" values between 0 and 24
   if( sunrise < 0.0 ):
   	pSunrise = sunrise + 24.0;
   elif( sunrise > 24.0 ):
   	pSunrise = sunrise - 24.0;
   else:
   	pSunrise = sunrise;
   
   if( sunset < 0.0 ):
   	pSunset = sunset + 24.0;
   elif( sunset > 24.0 ):
   	pSunset = sunset - 24.0;
   else:
   	pSunset = sunset;
   	
   # Determine appropriate minute and call SolarPos
   
   if( sunset - sunrise > 23.99 ):
      # Daytime hour, sun never sets
   	suntime = 1.0;       # Fraction of hour with sunup
   	minute = 30.0;       # For determining solar position at midpoint of sunup period
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )
   	
   elif( sunset - sunrise < 0.01 ):
      # Nighttime hour, sun never rises
   	suntime = 0.0;       # Fraction of hour with sunup
   	minute = 30.0;       # For determining solar position at midpoint of sunup period
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )
  	
   elif( hour == int(pSunrise + 1.0) and int(sunset+50) - int(sunrise+50) == 24 ):
      # Apply offset to deal with negative numbers in if statement
   	# Sun sets and rises same hour, summer 
   	suntime = pSunset + 1.0 - pSunrise;
   	minute = 60.0*( 1.0 - 0.5*( float(hour) - pSunrise ) );
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )
   	tmp = zen;
   	if( azm/DTOR < 180.0 ):
   		azm += 360.0*DTOR;
   	tmp2 = azm;
   	minute = 60.0*0.5*( pSunset - float(hour) + 1.0 );
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )
   	tmp += zen;
   	tmp2 += azm;
   	zen = tmp/2.0;    # Zenith angle 
   	azm = tmp2/2.0;   # Azimuth angle 
   	if( azm/DTOR > 360.0 ):
   		azm -= 360.0*DTOR;
   	
   elif( hour == int(pSunrise + 1.0) and int(sunrise) == int(sunset) ):
     # Sun rises and sets same hour, winter
   	# For zenith at mid-height
   	suntime = sunset - sunrise;
   	minute = 60.0*( pSunrise + 0.25*suntime - float(hour) + 1.0 );
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )
   	tmp = zen;        # Save zenith
   	# For azimuth at midpoint
   	minute = 60.0*( pSunrise + 0.5*suntime - float(hour) + 1.0 );
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )
   	zen = tmp;
   	    
   
   elif( hour == int(pSunrise + 1.0) ):
      # Sunrise hour
   	suntime = float(hour) - pSunrise;		# Fraction of hour with sunup
   	minute = 60.0*( 1.0 - 0.5*suntime );	# For determining solar position at midpoint of sunup period
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )    
   	
   elif( hour == int(pSunset + 1.0) ):
      # Sunset hour
   	suntime = pSunset - float(hour) + 1.0;  # Fraction of hour with sunup
   	minute = 60.0*0.5*suntime ;   # For determining solar position at midpoint of sunup period
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )    
   	
   
   elif( hour > int(sunrise+1.0) and hour < int(sunset+1.0) ):
      # Daytime hour
   	suntime = 1.0;       # Fraction of hour with sunup
   	minute = 30.0;       # For determining solar position at midpoint of sunup period
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )    
     
   else:
      # Nighttime hour
   	suntime = 0.0;       # Fraction of hour with sunup
   	minute = 30.0;       # For determining solar position at midpoint of sunup period
   	azm, zen, elv, dec, sunrise, sunset, Eo, tst = solarPos( year, month, day, hour-1, minute, lat, lng, tz )    

   
   return azm, zen, elv, dec, sunrise, sunset, Eo, tst, suntime;   
   # End of hrSolarPos

def iEEERemainder(x,y):
    z = x-y*round(x/y)
    return z;
    
def julian(year, month, day):
    		
	# Returns julian day of year
    i=1; jday=0; k=0;
    nday = [31,28,31,30,31,30,31,31,30,31,30,31];

    if( year % 4 == 0 ):
        k = 1;
    while( i < month ):
        jday += nday[i-1];
        i += 1;
    if( month > 2 ):
        jday += k + day;
    else:
        jday += day;	
    
    return jday;


def perezComp(dn, df, alb, inc, tilt, zen):      
    #Modified version of the Perez model to also return separate values for the 
    #diffuse components - isotropic sky, circumsolar, horizon, and ground reflected;
    #and beam component. 11/18/10
    #
    #Defines the Perez method for calculating values of diffuse + direct
    #solar radiation + ground reflected radiation for a tilted surface
    #and returns the total plane-of-array irradiance(poa).  Method does
    #not check all input for valid entries; consequently, this should be
    #done before calling the method.  (Reference: Perez et al, Solar
    #Energy Vol. 44, No.5, pp.271-289,1990.) Based on original FORTRAN
    #program by Howard Bisner.
    #
    #Modified 6/10/98 so that for zenith angles between 87.5 and 90.0 degrees,
    #the diffuse radiation is treated as isotropic instead of 0.0.
    #			
    #Converted 4/1/03 from C to C#.
    #
    #List of Parameters Passed to Method:
    #dn     = direct normal radiation (W/m2)
    #df     = diffuse horizontal radiation (W/m2)
    #alb    = surface albedo (decimal fraction)
    #inc    = incident angle of direct beam radiation to surface in radians
    #tilt   = surface tilt angle from horizontal in radians
    #zen    = sun zenith angle in radians
    #
    #Variable Returned
    #poa    = plane-of-array irradiance (W/m2), sum of direct beam and sky
    #         and ground-reflected diffuse 
    #iso_dif= diffuse from isotropic sky
    #circ_dif= diffuse from circumsolar region
    #horiz_dif= diffuse from sky region near horizon
    #grd_dif = ground-reflected diffuse
    #beam = direct beam component*/

    # Local variables   
    F11R =  ([-0.0083117, 0.1299457, 0.3296958, 0.5682053,
								 0.8730280, 1.1326077, 1.0601591, 0.6777470]) ;
    F12R =  ([0.5877285, 0.6825954, 0.4868735, 0.1874525,
								 -0.3920403, -1.2367284, -1.5999137, -0.3272588]) ;
    F13R =  ([-0.0620636, -0.1513752, -0.2210958, -0.2951290,
								 -0.3616149, -0.4118494, -0.3589221, -0.2504286]) ;
    F21R =  ([-0.0596012, -0.0189325, 0.0554140, 0.1088631,
								 0.2255647, 0.2877813, 0.2642124, 0.1561313]) ;
    F22R =  ([0.0721249, 0.0659650, -0.0639588, -0.1519229,
								 -0.4620442, -0.8230357, -1.1272340, -1.3765031]) ;
    F23R =  ([-0.0220216, -0.0288748, -0.0260542, -0.0139754,
								 0.0012448, 0.0558651, 0.1310694, 0.2506212]) ;
    EPSBINS =  [1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2] ;
    B2 = 0.000005534; DTOR = 0.01745329; 
    
    iso_dif = 0.0;
    circ_dif = 0.0;
    horiz_dif = 0.0;
    grd_dif = 0.0;
    beam = 0.0;

    if (dn < 0.0):           # Negative values may be measured if cloudy   
        dn = 0.0;
    if (zen < 0.0 or zen > 1.5271631): # Zen not between 0 and 87.5 deg   
    
        if (df < 0.0):
            df = 0.0;
        if (math.cos(inc) > 0.0 and zen < 1.5707963):  # Zen between 87.5 and 90   
                                              # and incident < 90 deg     
            poa = df * (1.0 + math.cos(tilt)) / 2.0 + dn * math.cos(inc);
            iso_dif = df * (1.0 + math.cos(tilt)) / 2.0;
            beam = dn * math.cos(inc);
            return poa, iso_dif, circ_dif, horiz_dif, grd_dif, beam;
        
        else:
            poa = df * (1.0 + math.cos(tilt)) / 2.0;   # Isotropic diffuse only 
            iso_dif = poa;
            return poa, iso_dif, circ_dif, horiz_dif, grd_dif, beam;
        
    
    else:                      # Zen between 0 and 87.5 deg   
    
        CZ = math.cos(zen);
        
        # C# code: ZH = (CZ > 0.0871557) ? CZ : 0.0871557;    # Maximum of 85 deg   
        # (condition) ? [true path] : [false path];
        if (CZ > 0.0871557):
            ZH = CZ
        else:
            ZH = 0.0871557 # Maximum of 85 deg   
        
        D = df;                # Horizontal diffuse radiation   
        if (D <= 0.0):        # Diffuse is zero or less        
        
            if (math.cos(inc) > 0.0):    # Incident < 90 deg   
            
                poa = 0.0 + dn * math.cos(inc);
                beam = poa;
                return poa, iso_dif, circ_dif, horiz_dif, grd_dif, beam;
            
            else:
            
                poa = 0.0;
                return poa, iso_dif, circ_dif, horiz_dif, grd_dif, beam;
            
        
        else:                  # Diffuse is greater than zero   
        
            ZENITH = zen / DTOR;
            AIRMASS = 1.0 / (CZ + 0.15 * math.pow(93.9 - ZENITH, -1.253));
            DELTA = D * AIRMASS / 1367.0;
            T = math.pow(ZENITH, 3.0);
            EPS = (dn + D) / D;
            EPS = (EPS + T * B2) / (1.0 + T * B2);
            i = 0;
            while (i < 7 and EPS > EPSBINS[i]):
                i += 1
                # i++;
            x = F11R[i] + F12R[i] * DELTA + F13R[i] * zen;
            #F1 = (0.0 > x) ? 0.0 : x;
            if (0.0 > x):
                F1 = 0.0
            else:
                F1 = x
            F2 = F21R[i] + F22R[i] * DELTA + F23R[i] * zen;
            COSINC = math.cos(inc);
            if (COSINC < 0.0):
                ZC = 0.0;
            else:
                ZC = COSINC;
            A = D * (1.0 + math.cos(tilt)) / 2.0;
            B = ZC / ZH * D - A;
            C = D * math.sin(tilt);
            # Diffuse on surface by component
            iso_dif = D * (1.0 - F1) * (1.0 + math.cos(tilt)) / 2.0;
            circ_dif = D * F1 * ZC / ZH;
            horiz_dif = F2 * C;
            grd_dif = alb * (dn * CZ + D) * (1.0 - math.cos(tilt)) / 2.0;
            beam = dn * ZC;
            # Total diffuse on surface
            poa = A + F1 * B + F2 * C + alb * (dn * CZ + D) * (1.0 - math.cos(tilt)) / 2.0 + dn * ZC;
            return poa, iso_dif, circ_dif, horiz_dif, grd_dif, beam;
            # End of perezComp
            

def solarPos( year, month, day, hour, minute, lat, lng, tz ): 		
    # This method is based on a paper by Michalsky published in Solar Energy
    # Vol. 40, No. 3, pp. 227-235, 1988. It calculates solar position for the
    # time and location passed to the method based on the Astronomical
    # Almanac's Algorithm for the period 1950-2050. For data averaged over an
    # interval, the appropriate time passed is the midpoint of the interval.
    # (Example: For hourly data averaged from 10 to 11, the time passed to the
    # method should be 10 hours and 30 minutes). The exception is when the time
    # interval includes a sunrise or sunset. For these intervals, the appropriate
    # time should be the midpoint of the portion of the interval when the sun is
    # above the horizon. (Example: For hourly data averaged from 7 to 8 with a
    # sunrise time of 7:30, the time passed to the method should be 7 hours and
    # and 45 minutes).
    #
    # Revised 5/15/98. Replaced algorithm for solar azimuth with one by Iqbal
    # so latitudes below the equator are correctly handled. Also put in checks
    # to allow an elevation of 90 degrees without crashing the program and prevented
    # elevation from exceeding 90 degrees after refraction correction.
    #			
    # Revised 4/1/03. Converted to C# and simplified in a few places. 
    #
    # This method calls the method Julian to get the julian day of year.
    #
    # List of Parameters Passed to Method:
    # year   = year (e.g. 1986)
    # month  = month of year (e.g. 1=Jan)
    # day    = day of month
    # hour   = hour of day, local standard time, (1-24, or 0-23)
    # minute = minutes past the hour, local standard time
    # lat    = latitude in degrees, north positive
    # lng    = longitude in degrees, east positive
    # tz     = time zone, west longitudes negative
    # List of Out Parameters
    # azm = sun azimuth in radians, measured east from north, 0 to 2*pi
    # zen = sun zenith in radians, 0 to pi
    # elv = sun elevation in radians, -pi/2 to pi/2
    # dec = sun declination in radians
    # sunrise = in local standard time (hrs), not corrected for refraction
    # sunset = in local standard time (hrs), not corrected for refraction
    # Eo = eccentricity correction factor
    # tst = true solar time (hrs)                */

			pi=math.pi; DTOR=math.pi/180
			zulu = 0.0; jd = 0.0; time = 0.0; mnlong = 0.0; mnanom = 0.0 
			eclong= 0.0; oblqec = 0.0; num = 0.0; den = 0.0; den = 0.0; ra = 0.0
			gmst = 0.0; lmst = 0.0; ha = 0.0; refrac = 0.0; E = 0.0; ws = 0.0; arg = 0.0
     
			jday = julian(year,month,day);		# Get julian day of year
			zulu = hour + minute/60.0 - tz;		# Convert local time to zulu time
			delta = year - 1949;
			leap = int(delta/4);
			jd = 32916.5 + delta*365 + leap + jday + zulu/24.0;
			time = jd - 51545.0;	# Time in days referenced from noon 1 Jan 2000

			mnlong = 280.46 + 0.9856474*time;
			mnlong = iEEERemainder(mnlong,360.0);	# Finds floating point remainder
			if( mnlong < 0.0 ):
				mnlong += 360.0;    # Mean longitude between 0-360 deg

			mnanom = 357.528 + 0.9856003*time;
			mnanom = iEEERemainder(mnanom,360.0);
			if( mnanom < 0.0 ):
				mnanom += 360.0;
			mnanom = mnanom*DTOR;	# Mean anomaly between 0-2pi radians 

			eclong = mnlong + 1.915*math.sin(mnanom) + 0.020*math.sin(2.0*mnanom);
			eclong = iEEERemainder(eclong,360.0);
			if( eclong < 0.0 ):
				eclong += 360.0;
			eclong = eclong*DTOR;	# Ecliptic longitude between 0-2pi radians

			oblqec = ( 23.439 - 0.0000004*time )*DTOR;   # Obliquity of ecliptic in radians
			num = math.cos(oblqec)*math.sin(eclong);
			den = math.cos(eclong);
			ra  = math.atan(num/den);	# Right ascension in radians
			if( den < 0.0 ):
				ra += pi;
			elif( num < 0.0 ):
				ra += 2.0*pi;

			dec = math.asin( math.sin(oblqec)*math.sin(eclong) );  # Declination in radians

			gmst = 6.697375 + 0.0657098242*time + zulu;
			gmst = iEEERemainder(gmst,24.0);
			if( gmst < 0.0 ):
				gmst += 24.0;			# Greenwich mean sidereal time in hours 

			lmst = gmst + lng/15.0;
			lmst = iEEERemainder(lmst,24.0);
			if( lmst < 0.0 ):
				lmst += 24.0;
			lmst = lmst*15.0*DTOR;		# Local mean sidereal time in radians 

			ha = lmst - ra;
			if( ha < -pi ):
				ha += 2*pi;
			elif( ha > pi ):
				ha -= 2*pi;				# Hour angle in radians between -pi and pi 

			lat = lat*DTOR;				# Change latitude to radians 

			arg = math.sin(dec)*math.sin(lat) + math.cos(dec)*math.cos(lat)*math.cos(ha);  # For elevation in radians
			if( arg > 1.0 ):
				elv = pi/2.0;
			elif( arg < -1.0 ):
				elv = -pi/2.0;
			else:
				elv = math.asin(arg);

			if( math.cos(elv) == 0.0 ):
				azm = pi;		# Assign azimuth = 180 deg if elv = 90 or -90
			else:
						# For solar azimuth in radians per Iqbal
				arg = ((math.sin(elv)*math.sin(lat)-math.sin(dec))/(math.cos(elv)*math.cos(lat))); # for azimuth
				if( arg > 1.0 ):
					azm = 0.0;              # Azimuth(radians)
				elif( arg < -1.0 ):
					azm = pi;
				else:
					azm = math.acos(arg);

				if( ( ha <= 0.0 and ha >= -pi) or ha >= pi ):
					azm = pi - azm;
				else:
					azm = pi + azm;
			

			elv = elv/DTOR;		# Change to degrees for atmospheric correction
			if( elv > -0.56 ):
				refrac = 3.51561*( 0.1594 + 0.0196*elv + 0.00002*elv*elv )/( 1.0 + 0.505*elv + 0.0845*elv*elv );
			else:
				refrac = 0.56;
			if( elv + refrac > 90.0 ):
				elv = 90.0*DTOR;
			else:
				elv = ( elv + refrac )*DTOR ; # Atmospheric corrected elevation(radians)

			E = ( mnlong - ra/DTOR )/15.0;       # Equation of time in hours
			if( E < - 0.33 ):   # Adjust for error occuring if mnlong and ra are in quadrants I and IV
				E += 24.0;
			elif( E > 0.33 ):
				E -= 24.0;

			arg = -math.tan(lat)*math.tan(dec);
			if( arg >= 1.0 ):
				ws = 0.0;						# No sunrise, continuous nights
			elif( arg <= -1.0 ):
				ws = pi;						# No sunset, continuous days
			else:
				ws = math.acos(arg);			# Sunrise hour angle in radians

			# Sunrise and sunset in local standard time
			sunrise = 12.0 - (ws/DTOR)/15.0 - (lng/15.0 - tz) - E;
			sunset  = 12.0 + (ws/DTOR)/15.0 - (lng/15.0 - tz) - E;

			Eo = 1.00014 - 0.01671*math.cos(mnanom) - 0.00014*math.cos(2.0*mnanom);  # Earth-sun distance (AU)   
			Eo = 1.0/(Eo*Eo);					# Eccentricity correction factor
			tst = hour + minute/60.0 + (lng/15.0 - tz) + E;  # True solar time (hr) 
			zen = 0.5*pi - elv;					#  Zenith angle		
				# End of SolarPos method
            
            
			return azm, zen, elv, dec, sunrise, sunset, Eo, tst;
			# End of solarPos


def sunIncident(mode, tilt, sazm, rlim, zen, azm):	
        #This method calculates the incident angle of direct beam radiation to a
        #surface for a given sun position, latitude, and surface orientation. The
        #modes available are fixed tilt, 1-axis tracking, and 2-axis tracking.
        #Azimuth angles are for N=0 or 2pi, E=pi/2, S=pi, and W=3pi/2.  8/13/98
        #
        #Converted 4/1/03 from C to C#.
        #			
        #List of Parameters Passed to Method:
        #mode   = 0 for fixed-tilt, 1 for 1-axis tracking, 2 for 2-axis tracking
        #tilt   = tilt angle of surface from horizontal in degrees (mode 0),
        #        		 or tilt angle of tracker axis from horizontal in degrees (mode 1),
        #        		 MUST BE FROM 0 to 90 degrees.
        #sazm   = surface azimuth in degrees of collector (mode 0), or surface
        #			azimuth of tracker axis (mode 1) with axis azimuth directed from
        #			raised to lowered end of axis if axis tilted.
        #rlim   = plus or minus rotation in degrees permitted by physical constraints
        #			of tracker, range is 0 to 180 degrees.
        #zen    = sun zenith in radians, MUST BE LESS THAN PI/2
        #azm    = sun azimuth in radians, measured east from north
        #
        #List of Out Parameters:
        #inc  = incident angle in radians
        #tiltr = tilt angle of surface from horizontal in radians
        #sazmr = surface azimuth in radians, measured east from north
        #
        #Local variables: rot is the angle that the collector is rotated about the
        #axis when viewed from the raised end of the 1-axis tracker. If rotated
        #counter clockwise the angle is negative. Range is -180 to +180 degrees.
        #When xsazm = azm : rot = 0, tilt = xtilt, and sazm = xsazm = azm  

			pi=3.1415927;DTOR=0.017453293;
			
			inc = 0.0;		
	

			if (mode == 0):               # Fixed-Tilt 
				                  
					tilt = tilt*DTOR;   # Change tilt and surface azimuth to radians   
					sazm = sazm*DTOR;
					arg = math.sin(zen)*math.cos(azm-sazm)*math.sin(tilt) + math.cos(zen)*math.cos(tilt);
					if( arg < -1.0 ):
						inc = pi;
					elif( arg > 1.0  ):
						inc = 0.0;
					else:
						inc = math.acos(arg);
					tiltr=tilt;
					sazmr=sazm;
					return inc, tiltr, sazmr;
                
			if (mode == 1):                # One-Axis Tracking  
			#	case 1:					 
					xtilt = tilt*DTOR;	# Change axis tilt, surface azimuth, and rotation limit to radians   
					xsazm = sazm*DTOR;
					rlim  = rlim*DTOR;
					# Find rotation angle of axis for peak tracking   
					if( abs( math.cos(xtilt) ) < 0.001745 ):    # 89.9 to 90.1 degrees   
										# For vertical axis only   
						if( xsazm <= pi ):
						
							if( azm <= xsazm + pi ):
								rot = azm - xsazm;
							else:
								rot = azm - xsazm - 2.0*pi;
						
						else	:		# For xsazm > pi   
						
							if( azm >= xsazm - pi ):
								rot = azm - xsazm;
							else:
								rot = azm - xsazm + 2.0*pi;
						
					
					else:				# For other than vertical axis   
					
						arg = (math.sin(zen)*math.sin(azm-xsazm)/
							( math.sin(zen)*math.cos(azm-xsazm)*math.sin(xtilt) + math.cos(zen)*math.cos(xtilt) ));
						if( arg < -99999.9 ):
							rot = -pi/2.0;
						elif( arg > 99999.9 ):
							rot = pi/2.0;
						else:
							rot = math.atan(arg);
						# Put rot in II or III quadrant if needed   
						if( xsazm <= pi ):
						
							if( azm > xsazm and azm <= xsazm + pi ):
										# Ensure positive rotation   
								if( rot < 0.0 ):
									rot = pi + rot;   # Put in II quadrant: 90 to 180 deg   
							
							else:
										# Ensure negative rotation    
								if( rot > 0.0 ):
									rot = rot - pi;   # Put in III quadrant: -90 to -180 deg   
							
						
						else	:		# For xsazm > pi   
						
							if( azm < xsazm and azm >= xsazm - pi ):
										# Ensure negative rotation    
								if( rot > 0.0 ):
									rot = rot - pi;   # Put in III quadrant: -90 to -180 deg   
							
							else:
										# Ensure positive rotation   
								if( rot < 0.0 ):
									rot = pi + rot;   # Put in II quadrant: 90 to 180 deg   
							
						
					
					#    printf("rot=%6.1f azm=%6.1f xsazm=%6.1f xtilt=%6.1f zen=%6.1f\n",rot/DTOR,azm/DTOR,xsazm/DTOR,xtilt/DTOR,zen/DTOR);    

					if( rot < -rlim ):	# Do not let rotation exceed physical constraints   
						rot = -rlim;
					elif( rot > rlim ):
						rot = rlim;
					# Find tilt angle for the tracking surface   
					arg = math.cos(xtilt)*math.cos(rot);
					if( arg < -1.0 ):
						tilt = pi;
					elif( arg > 1.0  ):
						tilt = 0.0;
					else:
						tilt = math.acos(arg);
					# Find surface azimuth for the tracking surface   
					if( tilt == 0.0 ):
						sazm = pi;     # Assign any value if tilt is zero   
					else:
					
						arg = math.sin(rot)/math.sin(tilt);
						if( arg < -1.0 ):
							sazm = 1.5*pi + xsazm;
						elif( arg > 1.0  ):
							sazm = 0.5*pi + xsazm;
						elif( rot < -0.5*pi ):
							sazm = xsazm - pi - math.asin(arg);
						elif( rot > 0.5*pi ):
							sazm = xsazm + pi - math.asin(arg);
						else:
							sazm = math.asin(arg) + xsazm;
						if( sazm > 2.0*pi ):       # Keep between 0 and 2pi   
							sazm = sazm - 2.0*pi;
						elif( sazm < 0.0 ):
							sazm = sazm + 2.0*pi;
					
					# printf("zen=%6.1f azm-sazm=%6.1f tilt=%6.1f arg=%7.4f\n",zen/DTOR,(azm-sazm)/DTOR,tilt/DTOR,arg);   
					# Find incident angle   
					arg = math.sin(zen)*math.cos(azm-sazm)*math.sin(tilt) + math.cos(zen)*math.cos(tilt);
					if( arg < -1.0 ):
						inc = pi;
					elif( arg > 1.0  ):
						inc = 0.0;
					else:
						inc = math.acos(arg);
					tiltr=tilt;
					sazmr=sazm;
					return inc, tiltr, sazmr;
             
			if (mode == 2):                # Two-Axis Tracking  
					tilt = zen;
					sazm = azm;
					inc = 0.0;
					tiltr=tilt;
					sazmr=sazm;
					return inc, tiltr, sazmr;
         # End of sunIncident method
