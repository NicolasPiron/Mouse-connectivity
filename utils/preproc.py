import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob

def get_grp_mat(pop):
    ''' Load all the matrices in a group and return them as a list of numpy arrays

    Parameters
    ----------
    pop : str
        The name of the group

    Returns
    -------
    mat_list : list
    '''

    path_list = glob.glob(f'data/{pop}/*.txt')
    mat_list = []
    for fname in path_list:
        data = pd.read_csv(fname, sep='\;', header=0, engine='python', index_col=0)
        arr = data.values
        if np.isnan(np.min(arr)):
            mean_val = np.nanmean(arr)
            arr[np.isnan(arr)] = mean_val
        mat_list.append(arr)
    print(f'Out of the {len(path_list)} files in {pop}, {len(mat_list)} were successfully loaded')
    return mat_list
    
def get_nbs_inputs(pop1, pop2):
    ''' Get the input for the NBS function

    Parameters
    ----------
    pop1 : str
        The name of the first group
    pop2 : str
        The name of the second group

    Returns
    -------
    stack : np.ndarray
        A 3D array of the connectivty data from both groups (n_edges x n_edges x n_subjects)
    y_vec : np.ndarray
        A 1D array of the group labels
    '''

    mat_list1 = get_grp_mat(pop1)
    mat_list2 = get_grp_mat(pop2)
    npop1 = len(mat_list1)
    npop2 = len(mat_list2)
    stack = np.stack((mat_list1 + mat_list2), axis=-1)
    y_vec = np.zeros((npop1 + npop2,))  
    y_vec[:npop1] = 1
    y_vec[npop1:] = 2

    return stack, y_vec, npop1, npop2
