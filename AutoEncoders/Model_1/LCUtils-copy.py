"""
This module contains a few functions that are useful when analyzing and manipulating light-curves.
"""
import numpy as np
import time as TIME

def get_transit_durations(lcs_array, cadence):
    """
    Get the transit duration from a transiting light-curve.
    
    This function scans multiple light-curves and measures the approximate transit duration in hours.
    
    Parameters:
        lcs_array (2D Numpy Array): Array of light-curves. Each row is a different LC, and each column is a point along the LC.
        cadence (float): Cadence (in seconds) of the light-curves.
        
    Returns:
        (2D Numpy Array): Array in which the first row corresponds to the approximate transit duration of each LC (in hours), the second row indicates whether a transit was found or not found in the corresponding LC (if a transit is not found, the corresponding transit duration will be set to 0 hours), the third row indicates the index at which the transit of each LC begins, and the fourth row indicates the index at which the transit of each LC ends.
    """
    
    cadence_in_hours = cadence/3600
    transit_durations = []
    transit_exists = []
    indices_begining_transits = []
    indices_end_transits = []
    
    for i in range(lcs_array.shape[0]):
        flux = lcs_array[i]
        first_one = False
        ingress = False
        done = False
        counter = 0
        index = -1
        for point in flux:
            index += 1
            if point==1 and first_one == False:
                first_one = True
            if ingress==True:
                if point!=1:
                    counter+=1
                if point==1:
                    index_end_of_transit = index
                    done=True
                    break
            if point!=1 and first_one==True and ingress==False:
                counter+=1
                index_begining_of_transit = index
                ingress=True
        
        if done==False:
            print('No transit found LC in pos: ',i)
            transit_durations.append(0)
            transit_exists.append(0)
            indices_begining_transits.append(0)
            indices_end_transits.append(0)
        
        if done==True:
            transit_durations.append(counter*cadence_in_hours)
            transit_exists.append(1)
            indices_begining_transits.append(index_begining_of_transit)
            indices_end_transits.append(index_end_of_transit)
            
    transit_durations = np.array(transit_durations)
    transit_exists = np.array(transit_exists)
    indices_begining_transits = np.array(indices_begining_transits)
    indices_end_transits = np.array(indices_end_transits)
    
    output = np.vstack((transit_durations, transit_exists, indices_begining_transits, indices_end_transits))
    
    return output
            

    
def LC_cropper(lcs_array, n_total_points, cadence, return_duration):
    """
    Crop given light-curves into a fixed time window.
    
    This function takes an array of light-curves, and returns a cropped copy of each transit (in a fixed time window). While it's guaranteed that each crop will contain the entirety of a transit (if it exists), the position of each transit along its light-curve is random.
    
    Parameters:
        lcs_array (2D Numpy Array): Array of light-curves. Each row is a different LC, and each column is a point along the LC.
        n_total_points (int): Desired amount of individual photometric points in the resulting cropped version of the light-curves.
        cadence (float): Cadence (in seconds) of the light-curves.
        return_duration (bool): If True, the function will also return an array of the light-curves' transit durations.
        
    Returns:
        if return_duration == True:
            ([2D Numpy Array, 1D Numpy Array, 1D Numpy Array]): A first array containing the cropped versions of the light-curves, a second array containing the transit durations (in hours), and a third array indicating whether a transit was found or not found in the corresponding LC (if a transit is not found, the corresponding transit duration will be set to 0 hours).
            
        if return_duration == False:
            ([2D Numpy Array, 1D Numpy Array]): A first array containing the cropped versions of the light-curves, and a second array indicating whether a transit was found or not found in the corresponding LC (if a transit is not found, the corresponding transit duration will be set to 0 hours).
    """
    helper = get_transit_durations(lcs_array, cadence)
    transit_durations = helper[0]
    transit_exists = helper[1]
    transit_indices = helper[2]
    cadence_in_hours = cadence/3600
    n_points_duration = transit_durations/cadence_in_hours
    cropped_lcs = np.zeros((lcs_array.shape[0], n_total_points), dtype=float)
    
    for i in range(len(n_points_duration)):

        #seed = int(TIME.time()+np.random.randint(100))
        #np.random.seed(seed)
        possible_shifts = np.arange(1, (n_total_points-n_points_duration[i]))
        shift = np.random.choice(possible_shifts)
        first_index = int(transit_indices[i]-shift)
        last_index = int(transit_indices[i] - shift + n_total_points)
        if transit_exists[i] == 1:
            if first_index > 1:
                if last_index < lcs_array.shape[1]:
                    cropped_lcs[i, :] = lcs_array[i, first_index:last_index]
                if last_index > lcs_array.shape[1]:
                    cropped_lcs[i, :] = lcs_array[i,-(n_total_points+2):(-2)]
            if first_index < 1:
                cropped_lcs[i, :] = lcs_array[i, 1:1+n_total_points]
        
        if transit_exists[i] == 0:
            cropped_lcs[i, :] = np.linspace(1,1,n_total_points)
        
    if return_duration == True:
        
        return [cropped_lcs, transit_durations, transit_exists]
    
    if return_duration == False:
        
        return [cropped_lcs, transit_exists]
        
        

def get_transit_start(LC_array, cadence):
    """
    Get the hour at which transit starts.
    
    This function takes an array of light-curves and measures the hour at which transit begins.
    
    Parameters:
        LC_array (Numpy 2D Array): Array of light-curves. Each row is a different LC, and each column is a point along the LC.
        cadence (float): Cadence (in seconds) of the light-curves.
    
    Returns:
        (Numpy 1D Array): Array of hours at which each light-curve's transit begins.
    """
    
    cadence_in_hours = cadence/3600
    first_non_one_index = np.argmax(LC_array != 1, axis=1)
    transit_start_hours = first_non_one_index*cadence_in_hours
    return transit_start_hours
    
    
    
    
    
          
            
      
        