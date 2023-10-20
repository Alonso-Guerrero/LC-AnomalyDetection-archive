"""Synthetic Dataset Genesis functions:

This module contains a few functions to build synthetic datasets of simple ringed or non-ringed planetary systems out of an initial real dataset of real simple non-ringed planetary systems."""

# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import string
import math

def get_longests_periods(original_data, interest_columns, period_label, sample_size, save, filename):
    
    """ 
    Get the systems with longests periods.
    
    This function takes a dataset of real-life simple planetary systems, it selects only the columns we are interested in, then it cleans it, it sorts the dataset from longest to shortest period, and finally it selects the number of elements that we indicate (counting from longest to shortest periods) 
    
    Parameters:
        original_data (Pandas DataFrame): Whole dataset of real-life simple planetary systems.
        interest_columns (list of str): List of names of all the columns that we wish to recover from original_data.
        period_label (str): Name of the columns from original_data which contains the orbital period values.
        sample_size (int): Amount of sorted elements that we wish to recover from original_data.
        save (bool): If True, the resulting dataset will be saved as .csv. If False, it won't.
        filename (str): Name we wish to give to the saved file if save=True. Warning: don't write ".csv" when writing the file's name.
        
    Returns:
        (Pandas DataFrame): Resulting dataset.
    """
    
    resulting_dataset = original_data[interest_columns].dropna().sort_values(period_label, axis=0, ascending=False, ignore_index=True)[:sample_size]
    
    alphabet = list(string.ascii_uppercase)
    amount_first_letters = math.ceil(resulting_dataset.shape[0]/len(alphabet))
    numberless_names = []
  
    for i in alphabet[:amount_first_letters]:
        if i == alphabet[amount_first_letters - 1]:
            for j in alphabet[:resulting_dataset.shape[0]%len(alphabet)]:
                numberless_names.append(i+j)
        else:   
            for j in alphabet:
                numberless_names.append(i+j)
    
    numberless_names = np.array(numberless_names)
    df_names = pd.DataFrame(numberless_names, columns=['ID'])
    resulting_dataset = pd.concat([df_names, resulting_dataset], axis=1)
    
    print('Dataset has been succesfully cleaned, sorted and cropped.')
    if save:
        resulting_dataset.to_csv(filename+'.csv')
    return resulting_dataset



def synthetic_nonringed_generator(dataset, labels, upp, low, eccen_label, exception_labels, draw_size, save, filename):
    
    """
    Generate a synthetic dataset of simple non-ringed planetary systems."
    
    This function builds a synthetic dataset of simple non-ringed planetary systems out of a pre-existing dataset of real simple planetary systems, by performing gaussian random draws for each of its parameters. Each draw is centered at the real value of a parameter, and uses the mean absolute value of their uncertainties as standard deviation.

    After all the process, an  integrity check is ran onto the columns. This check makes sure of two things: one, that the eccentricty values have no forbidden values (since the condition for elliptical orbits is: 0 < e < 1), and two, that there are no negative numbers in the parameters (except for parameters which we want to let contain negative numbers).
    
    **IMPORTANT:** This piece of code REQUIRES that the dataset used as input includes the UPPER AND LOWER UNCERTAINTIES of its values. Said uncertainties must be stored in TWO COLUMNS: One for the UPPER one and another for the LOWER one. Each of those columns must have a MATCHING NAME with the column containing the VALUE of the parameter, PLUS A DISTINCTIVE SUFFIX.

    Example: If the VALUE of the mass of a planet is stored in the column "PlanetMass", then its  UPPER UNCERTAINTY could be stored in a column named "PlanetMass_upper", and its LOWER UNCERTAINTY could be stored in a column named "PlanetMass_lower". These suffixes must be the same for all orbital parameters.

    Parameters:
        dataset (Pandas DataFrame): Dataset of real-life simple planetary systems. MUST INCLUDE UPPER AND LOWER UNCERTAINTIES OF VALUES AS INDIVIDUAL COLUMNS.
        labels (list of str): List of names of columns of dataset containing the parameters' values (NOT its uncertainties. Only the values).
        upp (str): Suffix of the names of the columns storing the upper uncertainties.
        low (str): Suffix of the names of the columns storing the lower uncertainties.
        eccen_label (str): Name of the column storing the eccentricity values in dataset.
        exceptions_labels (list of str): List of names of the columns which we want to let contain negative numbers.
        draw_size (int): Amount of random draws to be performed from the same element of dataset (e.g.: if dataset has 10 elements, and draw_size=3, the code will make 3 random draws for each of the 10 elements. Thus, the resulting synthetic dataset will contain 30 elements).
        save (bool): If True, the resulting dataset will be saved as .csv. If False, it won't.
        filename (str): Name we wish to give to the saved file if save=True. Warning: don't write ".csv" when writing the file's name.
     
     Returns:
        (Pandas DataFrame): Generated synthetic dataset.
         
"""
  
    #First, we use the upper and lower uncertainties of "dataset" to calculate the STANDARD DEVIATIONS for each parameter.
    #Then, those calculations are appended to the original dataframe for an easier manipulation in the following steps.
    extra_columns = []
    for i in labels:
        extra_columns.append(((abs(dataset[i+upp]) + abs(dataset[i+low]))/2).rename(i+'_devstd'))
    working_dataframe = pd.concat([dataset, pd.DataFrame(extra_columns).T], axis=1)

    #Now we have to name each draw with a naming convention. Each draw from the same planet will be named with a string of letters, and each
    #consecutive draw will add '+1' to its name, starting with 0. For example, the third draw of the planet "AB" would be named "AB2".

    #WARNING: THIS NAMING ALGORITHM CAN ONLY INDIVIDUALLY NAME UP TO 26*26=676 PLANETS, BECAUSE WE ARE USING TWO-LETTER NAMES.

    alphabet = list(string.ascii_uppercase)
    amount_first_letters = math.ceil(dataset.shape[0]/len(alphabet))

    numberless_names = []
  
    for i in alphabet[:amount_first_letters]:
        if i == alphabet[amount_first_letters - 1]:
            for j in alphabet[:dataset.shape[0]%len(alphabet)]:
                numberless_names.append(i+j)
        else:   
            for j in alphabet:
                numberless_names.append(i+j)

    #Now, we perform the random draws:
    draws=[]
    for counter in range(0, draw_size):
        drew_dataframe = pd.DataFrame(np.random.normal(loc=working_dataframe[labels].values, scale=working_dataframe[(i + '_devstd' for i in labels)].values), columns=labels)
    
        #We assign each of these iterations their respective names:
        names = pd.DataFrame(np.array([n + str(counter) for n in numberless_names]), columns=['ID'])
        named_drew_dataframe = pd.concat([names, drew_dataframe], ignore_index=False, axis=1)

        #Now an integrity check is performed on the eccentricities. If they got assigned forbidden values, they will be re-rolled.


        print('Draw', counter+1, 'of', draw_size)
        print('')

        while (named_drew_dataframe[eccen_label] <= 0).sum() > 0:
            print('Found ', (named_drew_dataframe[eccen_label] <= 0).sum(), 'eccentricities below 0.')
            print('Re-rolling them...')
            print(' ')
            bad_indexes = np.where(named_drew_dataframe[eccen_label] <0)[0]
            named_drew_dataframe[eccen_label][bad_indexes] = np.random.normal(loc=working_dataframe[eccen_label][bad_indexes].values, scale=working_dataframe[eccen_label+'_devstd'][bad_indexes].values)
            
        while (named_drew_dataframe[eccen_label] >= 1).sum() > 0:
            print('Found ', (named_drew_dataframe[eccen_label] >= 1).sum(), 'eccentricities above 1.')
            print('Re-rolling them...')
            print(' ')
            bad_indexes = np.where(named_drew_dataframe[eccen_label] >= 1)[0]
            named_drew_dataframe[eccen_label][bad_indexes] = np.random.normal(loc=working_dataframe[eccen_label][bad_indexes].values, scale=working_dataframe[eccen_label+'_devstd'][bad_indexes].values)

        print('All forbidden eccentricities have been corrected.')
        print('Procceeding to check the integrity of the remaining columns.')
        print('')

        for L in labels:

            if L not in exception_labels:

                while (named_drew_dataframe[L] < 0).sum() > 0:
                    print('Found ', (named_drew_dataframe[L] < 0).sum(), 'negative values in ', L,'.')
                    print('Re-rolling them...')
                    print('')
                    bad_indexes = np.where(named_drew_dataframe[L] < 0)[0]
                    named_drew_dataframe[L][bad_indexes] = np.random.normal(loc=working_dataframe[L][bad_indexes].values, scale=working_dataframe[L+'_devstd'][bad_indexes].values)
        print('--------------------------------------------------------------------')
        print('')

        draws.append(named_drew_dataframe)
  
    synthetic_dataset = pd.concat(draws, ignore_index=True)

    if save: 
        synthetic_dataset.to_csv(filename+'.csv')
  
    return synthetic_dataset


def synthetic_ringed_generator(dataset, labels, upp, low, eccen_label, plmass_label, plradius_label, plmass_unit_kg, plradius_unit_m, ring_density, max_ring_incl, exception_labels, draw_size, save, filename):

    """
    Generate a synthetic dataset of simple ringed planetary systems.
    
    This function takes a dataset of real-life, then it applies the function synthetic_nonringed_generator() once. After that, it calculates the necessary parameters for placing an artificial ring system of density ring_density around all objects. Once the grid is completely generated, it discards all the elements in which the ring system is not physical (i.e. its Roche Limit ends up being lower than the planet's radius). 
    
    Parameters:
        dataset (Pandas DataFrame): Dataset of real-life simple planetary systems. MUST INCLUDE UPPER AND LOWER UNCERTAINTIES OF VALUES AS INDIVIDUAL COLUMNS.
        labels (list of str): List of names of columns of dataset containing the parameters' values (NOT its uncertainties. Only the values).
        upp (str): Suffix of the names of the columns storing the upper uncertainties.
        low (str): Suffix of the names of the columns storing the lower uncertainties.
        eccen_label (str): Name of the column storing the eccentricity values in dataset.
        plmass_label (str): Name of the column storing the planet mass values in dataset.
        plradius_label (str): Name of the column storing the planet radius values in dataset.
        plmass_unit_kg (float): Corresponding kilograms of measure units in which planet mass is measured in dataset.
        plradius_unit_m (float): Corresponding meters of measure units in which planet radius is measured in dataset.
        ring_density (float): Density of artificial ring for all planets (in kg/m**3).
        max_ring_incl (float): Max ring inclination with respect to orbital plane (deg).    
        exceptions_labels (list of str): List of names of the columns which we want to let contain negative numbers.
        draw_size (int): Amount of random draws to be performed from the same element of dataset before discarding elements by Roche Limit criteria (e.g.: if dataset has 10 elements, and draw_size=3, the code will make 3 random draws for each of the 10 elements. Thus, the resulting synthetic dataset will contain 30 elements before discarding elements by Roche Limit criteria).
        save (bool): If True, the resulting dataset will be saved as .csv. If False, it won't.
        filename (str): Name we wish to give to the saved file if save=True. Warning: don't write ".csv" when writing the file's name.
     
     Returns:
        (Pandas DataFrame): Generated synthetic dataset.
    
    """
    
    
    
    
    df = synthetic_nonringed_generator(dataset, labels, upp, low, eccen_label, exception_labels, draw_size, False, None)
    pl_masses = df[plmass_label].values
    pl_radius = df[plradius_label].values
    pl_volume = (4/3)*np.pi*(pl_radius*plradius_unit_m)**3
    pl_density = (pl_masses*plmass_unit_kg)/pl_volume
    roche_limits = (2*pl_density/ring_density)**(1/3) #in units of planet radius
    exterior_radius = []
    interior_radius = []

    for i in roche_limits:
        exterior_radius.append(np.random.uniform(1, i))

    for i in exterior_radius:
        interior_radius.append(np.random.uniform(1, i))

    exterior_radius = np.array(exterior_radius)
    interior_radius = np.array(interior_radius)
    ring_inclination = np.random.uniform(-max_ring_incl, max_ring_incl, size=exterior_radius.shape[0])
    
    array_of_arrays = np.array([pl_volume, pl_density, roche_limits, exterior_radius, interior_radius, ring_inclination]).T
    list_of_names = ['pl_volume_m**3', 'pl_rho_kg/m**3', 'pl_r_limit_plrad', 'ring_ext_plrad', 'ring_int_plrad', 'ring_incl_deg']
    additional_df = pd.DataFrame(array_of_arrays, columns=list_of_names)
    
    whole_df = pd.concat([df, additional_df], axis=1)
    
    valid_df = whole_df[whole_df['pl_r_limit_plrad']>1].reset_index(drop=True) #Select only those that can indeed present rings around them.
    
    if save: 
        valid_df.to_csv(filename+'.csv')
        
    return valid_df
    
    
    