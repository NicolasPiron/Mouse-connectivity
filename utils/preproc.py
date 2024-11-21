import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
#from params import groups

def check_tree():
    ''' Create the directory tree for the results (derivative)'''

    if not os.path.exists('derivative/average'):
        os.makedirs('derivative/average')
    if not os.path.exists('derivative/stats'):
        os.makedirs('derivative/stats')
    if not os.path.exists('derivative/stats/null'):
        os.makedirs('derivative/stats/null')
    if not os.path.exists('derivative/stats/pvals'):
        os.makedirs('derivative/stats/pvals')
    if not os.path.exists('derivative/stats/adjacency'):
        os.makedirs('derivative/stats/adjacency')
    if not os.path.exists('derivative/stats/figures'):
        os.makedirs('derivative/stats/figures')
    
def get_single_mat(id):

    path = glob.glob(f'data/*/souris_{id}.csv')[0]
    data = pd.read_csv(path, index_col=0)

    return data.values

def get_grp_mat(pop, females=False):
    ''' Load all the matrices in a group and return them as a list of numpy arrays

    Parameters
    ----------
    pop : str
        The name of the group

    Returns
    -------
    mat_list : list
    '''

    mat_list = []
    if females:
        desc = pd.read_csv('data/all_df.csv')
        ids = desc[(desc['group']==pop) & (desc['sex']=='f')]['id']
        for id in ids:
            path = glob.glob(f'data/*/souris_{id}.csv')[0]
            data = pd.read_csv(path, index_col=0)
            arr = data.values
            if np.isnan(np.min(arr)):
                mean_val = np.nanmean(arr)
                arr[np.isnan(arr)] = mean_val
            mat_list.append(arr)
    else:
        path_list = glob.glob(f'data/{pop}/*.csv')
        mat_list = []
        for fname in path_list:
            data = pd.read_csv(fname, index_col=0)
            arr = data.values
            if np.isnan(np.min(arr)):
                mean_val = np.nanmean(arr)
                arr[np.isnan(arr)] = mean_val
            mat_list.append(arr)
        print(f'Out of the {len(path_list)} files in {pop}, {len(mat_list)} were successfully loaded')

    return mat_list


def get_nbs_inputs(pop1, pop2, females=False):
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

    mat_list1 = get_grp_mat(pop1, females=females)
    mat_list2 = get_grp_mat(pop2, females=females)
    npop1 = len(mat_list1)
    npop2 = len(mat_list2)
    stack = np.stack((mat_list1 + mat_list2), axis=-1)
    y_vec = np.zeros((npop1 + npop2,))  
    y_vec[:npop1] = 1
    y_vec[npop1:] = 2

    return stack, y_vec, npop1, npop2


def all_df():
    ''' Creates a df with average connectivty value, the id,
    the sex and the group for each animal'''

    desc = pd.read_excel('data/code_animaux.xlsx')
    new_rows = []
    no_data = []
    for row in desc.iterrows():
        grp, sex, id = row[1]
        if grp == 'C57BL/6J':
            grp = 'WT'
        elif grp == '3xTg-AD':
            grp = '3xTgAD'
        elif grp == 'Tspo KO':
            grp = 'TSPO_KO'
        elif grp == '3xTg-AD-TSPO':
            grp = '3xTgAD_TSPO_KO'
        print(f'Processing {id}')
        try:
            val = np.nanmean(get_single_mat(id))
            new_row = pd.DataFrame({'id': [id], 'group': [grp], 'sex':[sex], 'average_connectivity': [val]})
            new_rows.append(new_row)
        except:
            print(f'Could not process {id}')
            no_data.append(id)

    df = pd.concat(new_rows)
    df.to_csv('data/all_df.csv', index=False)
    print('Dataframe created')
    print(f'Could not process {len(no_data)} animals: {no_data}')

    return None


def txt_csv(groups=['WT', '3xTgAD', 'TSPO_KO', '3xTgAD_TSPO_KO']):
    ''' Convert all the txt files in a group to csv files.
    Only keeps the id in the filename for easier access.
    '''

    for pop in groups:
        path_list = glob.glob(f'data/{pop}/*.txt')
        for fname in path_list:
            data = pd.read_csv(fname, sep='\;', header=0, engine='python', index_col=0)
            animal_id = fname.split('Souris')[1].split('_')[1]
            new_fname = f'data/{pop}/souris_{animal_id}.csv'
            data.to_csv(new_fname, index=True)
            print(f'{fname} was converted to {new_fname}')

    return None