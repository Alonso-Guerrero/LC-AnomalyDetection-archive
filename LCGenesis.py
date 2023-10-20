"""This module contains functions to simulate Light-Curves of transiting objects given a dataset of their parameters. It uses the
Pryngles library, which allows to simulate planets both non-ringed and ringed.
"""

# -*- coding: utf-8 -*-

import pryngles as pr
from pryngles import Consts
import numpy as np
import time as TIME

def LC_nonringed_simulator(dataframe, id_label, label_prad, label_smass, label_srad, label_semi, label_ecc, label_incl, label_peri, time, cadence):
    """
    Simulate the LC of non-ringed transiting planets.
    
    This function runs some Pryngles code to simulate the transiting LC of planets over a given period of time and a given cadence.
    
    Note about limb-darkening: In these simulations, a fixed linear coefficient of 0.65 will be kept for the stars' limb darkening.
    
    Note about time of periapsis: For each simulation, this algorithm will randomly select the moment at which periapsis occurs.
    
    Parameters:
        dataframe (Pandas DataFrame): Dataset of orbital parameters of from which to simulate the LCs. This dataset's values must be in a
    particular set of units, as explained in the other arguments below.
        id_label (str): Name of the dataset column which contains the planets' IDs.
        label_prad (str): Name of the column containing the planet's radius values (unit: Jupiter radius).
        label_smass (str): Name of the column containing the stars' mass values (unit: Solar mass).
        label_srad (str): Name of the column containing the stars' radius values (unit: Solar radius).
        label_semi (str): Name of the column containing the orbital semi-major axis values (unit: AU).
        label_ecc (str): Name of the column containing the orbital eccentricity values.
        label_incl (str): Name of the column containing the orbital inclination (unit: Degrees). Note: In Pryngles, an inclination of 90° means a face-on orbit, but in most datasets 90° means edge-on. Thus, in Pryngles we use (incl - 90) as the inclination value to make up for this incongruence.
        label_peri (str): Name of the column containing the orbital longitude of periastron (unit: Degrees). Note: This value refers to a heliocentric angular position, but Pryngles is a planetocentric model. During these simulations, we will use an observer's angular position of 0°, which means that if the planetocentric longitude of periastron is also 0°, the star's closest position will be directly behind the planet. Thus, during the simulations we will use (peri + 90°) to indicate the planetocentric longitude of periastron.
        time (float): Total time of observation to simulate (unit: Years).
        cadence (float): Simulated photometric cadence (unit: Seconds).
        
    Returns:
        ([list, list, list]): A first list of IDs (str), a second list of orbital periods in seconds (float), and a third list containing each simulated LC (each of them as a 1D Numpy array).
    
    """

    IDS = []
    periods = []
    lightcurves = []

    for j in range(0, len(dataframe[label_prad].values)):

        if j%10==0:
            print('The process is currently starting simulation ',j)

        planet_radius = dataframe.iloc[j][label_prad]
        star_mass = dataframe.iloc[j][label_smass]
        star_radius = dataframe.iloc[j][label_srad]
        semimajor_axis = dataframe.iloc[j][label_semi]
        eccentricity = dataframe.iloc[j][label_ecc]
        inclination = dataframe.iloc[j][label_incl]
        longitude_periastron = dataframe.iloc[j][label_peri]

        sys=pr.System()
        S=sys.add(kind="Star",radius=star_radius*Consts.rsun/sys.ul,limb_coeffs=[0.65])
        P=sys.add(kind="Planet",parent=S,a=semimajor_axis,e=eccentricity,radius=planet_radius*Consts.rjupiter/sys.ul)
        R=sys.add(kind="Ring",parent=P,fi=1.5,fe=1.5,i=30*Consts.deg)
        RP=sys.ensamble_system(lamb=0*Consts.deg,beta=(inclination-90)*Consts.deg)
        RP.Mstar=star_mass*Consts.msun/sys.um
        RP.lambq = (longitude_periastron+90)*Consts.deg

        obs_time = time*Consts.yr
        n_points = int(obs_time/cadence)
        instants = np.linspace(0, obs_time, n_points)
        #seed = int(str(TIME.time())[-6:])
        seed = int(TIME.time()+np.random.randint(100))
        np.random.seed(seed)
        period = RP.T*Consts.yr/(2*np.pi)
        randomizator =np.random.choice(np.linspace(0, period, 100))

        Tps=[]

        for t in instants:
            RP.changeStellarPosition(kepler=True, x=(t+randomizator)/sys.ut)
            #RP.updateOpticalFactors()
            RP.updateTransit()
            Tp=RP.Tip.sum()
            Tps+=[Tp]
  
        Tps = np.array(Tps)

        IDS.append(dataframe.iloc[j][id_label])
        periods.append(period)
        lightcurves.append(1-Tps)

    return [IDS, periods, lightcurves]


def LC_ringed_simulator(dataframe, id_label, label_prad, label_smass, label_srad, label_semi, label_ecc, label_incl, label_peri, label_ext_ring, label_int_ring, label_ring_incl, n_spangles, time, cadence):
    """
    Simulate the LC of non-ringed transiting planets.
    
    This function runs some Pryngles code to simulate the transiting LC of planets over a given period of time and a given cadence.
    
    Note about limb-darkening: In these simulations, a fixed linear coefficient of 0.65 will be kept for the stars' limb darkening.
    
    Note about time of periapsis: For each simulation, this algorithm will randomly select the moment at which periapsis occurs.
    Parameters:
        dataframe (Pandas DataFrame): Dataset of orbital parameters of from which to simulate the LCs. This dataset's values must be in a
    particular set of units, as explained in the other arguments below.
        id_label (str): Name of the dataset column which contains the planets' IDs.
        label_prad (str): Name of the column containing the planet's radius values (unit: Jupiter radius).
        label_smass (str): Name of the column containing the stars' mass values (unit: Solar mass).
        label_srad (str): Name of the column containing the stars' radius values (unit: Solar radius).
        label_semi (str): Name of the column containing the orbital semi-major axis values (unit: AU).
        label_ecc (str): Name of the column containing the orbital eccentricity values.
        label_incl (str): Name of the column containing the orbital inclination (unit: Degrees). Note: In Pryngles, an inclination of 90° eans a face-on orbit, but in most datasets 90° means edge-on. Thus, in Pryngles we use (incl - 90) as the inclination value to make p for this incongruence.
        label_peri (str): Name of the column containing the orbital longitude of periastron (unit: Degrees). Note: This value refers to a heliocentric angular position, but Pryngles is a planetocentric model. During these simulations, we will use an observer's angular position of 0°, which means that if the planetocentric longitude of periastron is also 0°, the star's closest position will be directly behind the planet. Thus, during the simulations we will use (peri + 90°) to indicate the planetocentric longitude of periastron.
        label_ext_ring (str): Name of the column containing the ring's exterior radius values (unit: Planetary radius).
        label_int_ring (str): Name of the column containing the ring's interior radius values (unit: Planetary radius).
        label_ring_incl (str): Name of the column containing the ring's inclination with respect to the orbital plane (unit: Degrees).
        n_spangles (int): Amount of individual "spangles" that conform the ring.        
        time (float): Total time of observation to simulate (unit: Years).
        cadence (float): Simulated photometric cadence (unit: Seconds).

        
    Returns:
        ([list, list, list]): A first list of IDs (str), a second list of orbital periods in seconds (float), and a third list containing each simulated LC (each of them as a 1D Numpy array).
        
    
    """

    IDS = []
    periods = []
    lightcurves = []

    for j in range(0, len(dataframe[label_prad].values)):

        if j%10==0:
            print('The process is currently starting simulation ',j)

        planet_radius = dataframe.iloc[j][label_prad]
        star_mass = dataframe.iloc[j][label_smass]
        star_radius = dataframe.iloc[j][label_srad]
        semimajor_axis = dataframe.iloc[j][label_semi]
        eccentricity = dataframe.iloc[j][label_ecc]
        inclination = dataframe.iloc[j][label_incl]
        longitude_periastron = dataframe.iloc[j][label_peri]
        exterior_ring = dataframe.iloc[j][label_ext_ring]
        interior_ring = dataframe.iloc[j][label_int_ring]
        ring_inclination = dataframe.iloc[j][label_ring_incl]
        
        
        sys=pr.System()
        S=sys.add(kind="Star",radius=star_radius*Consts.rsun/sys.ul,limb_coeffs=[0.65])
        P=sys.add(kind="Planet",parent=S,a=semimajor_axis,e=eccentricity,radius=planet_radius*Consts.rjupiter/sys.ul)
        R=sys.add(kind="Ring",parent=P,fi=interior_ring,fe=exterior_ring,i=ring_inclination*Consts.deg, nspangles=n_spangles)
        RP=sys.ensamble_system(lamb=0*Consts.deg,beta=(inclination-90)*Consts.deg)
        RP.Mstar=star_mass*Consts.msun/sys.um
        RP.lambq = (longitude_periastron+90)*Consts.deg

        obs_time = time*Consts.yr
        n_points = int(obs_time/cadence)
        instants = np.linspace(0, obs_time, n_points)
        #seed = int(str(TIME.time())[-6:])
        seed = int(TIME.time()+np.random.randint(100))
        np.random.seed(seed)
        period = RP.T*Consts.yr/(2*np.pi)
        randomizator =np.random.choice(np.linspace(0, period, 100))
        
        Tps=[]
        Ts=[]

        for t in instants:
            RP.changeStellarPosition(kepler=True, x=(t+randomizator)/sys.ut)
            #RP.updateOpticalFactors()
            RP.updateTransit()
            Tp=RP.Tip.sum()
            T=Tp+RP.Tir.sum()
            Tps+=[Tp]
            Ts+=[T]
  
        Ts = np.array(Ts)

        IDS.append(dataframe.iloc[j][id_label])
        periods.append(period)
        lightcurves.append(1-Ts)

    return [IDS, periods, lightcurves]

