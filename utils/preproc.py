import pandas as pd
import numpy as np
import os
import glob

def pre_run_check():
    ''' Check if the data is available and preprocessed. If not, preprocess the data.
    Also ensures that the directory tree for the results is created.
    Should be run before any analysis script.'''

    check_tree()
    if not len(glob.glob('data/*/*souris*.csv')) == 0:
        print('Data already preprocessed')
    else:
        try:
            txt_csv()
            print('Successfully transformed .txt to .csv')
        except:
            print('Could not transform .txt to .csv, some analyses may not work')
    if not len(glob.glob('data/*/*souris*zscore.csv')) == 0:
        print('Data already z-scored')
    else:
        try:
            zscore_mat()
            print('Successfully z-scored the matrices')
        except:
            print('Could not z-score the data, some analyses may not work')
    if not os.path.exists('data/all_df.csv'):
        try:
            all_df()
        except:
            print('Could not create the dataframe with the average connectivity values, some analyses may not work')

def check_tree():
    ''' Create the directory tree for the results (derivative)'''

    paths = ['derivative/average/diff',
            'derivative/average/raw',
            'derivative/average/zscored',
            'derivative/average/boxplot',
            'derivative/nbs/null',
            'derivative/nbs/pvals',
            'derivative/nbs/adjacency',
            'derivative/nbs/figures',
            'derivative/ttest/figures/raw_pvals',
            'derivative/ttest/pvals',
            'derivative/ttest/pvals/raw_pvals',
            'derivative/permutations/figures/raw_pvals',
            'derivative/permutations/pvals',
            'derivative/permutations/pvals/raw_pvals',
            'derivative/anova',
            'derivative/individuals/',
    ]
    
    for path in paths:
        if not os.path.exists(path):
            os.makedirs(path)
            print(f'{path} was created')
    
def get_single_mat(id, z=False):

    # not used yet
    if z:
        path = glob.glob(f'data/*/souris_{id}_zscore.csv')[0]
    else:
        path = glob.glob(f'data/*/souris_{id}.csv')[0]
    data = pd.read_csv(path, index_col=0)

    return data.values


def get_av_grp_mat(pop, females=False, z=True):
    ''' '''

    mat_list = get_grp_mat(pop, females=females, z=z)
    stack = np.stack((mat_list), axis=-1)

    return np.mean(stack, axis=-1)

def get_grp_mat(pop, females=False, z=False):
    ''' Load all the matrices in a group and return them as a list of numpy arrays

    Parameters
    ----------
    pop : str
        The name of the group
    females : bool
        If True, only load the matrices of female mice. Default is False.
    z : bool
        If True, load the z-scored matrices. Default is False.

    Returns
    -------
    mat_list : list
    '''

    mat_list = []
    
    desc = pd.read_csv('data/all_df.csv')
    if females:
        mask = (desc['group']==pop) & (desc['sex']=='f')
    else:
        mask = desc['group']==pop
    ids = desc[mask]['id']
    for id in ids:
        if z:
            path = glob.glob(f'data/*/souris_{id}_zscore.csv')[0]
        else:
            path = glob.glob(f'data/*/souris_{id}.csv')[0]
        data = pd.read_csv(path, index_col=0)
        arr = data.values # extract values
        if np.isnan(np.min(arr)): # if NaNs in data, replace with average
            tril_indices = np.tril_indices(arr.shape[0], k=-1) # extract lower triangle values (not counting the diagonal)
            tril = arr[tril_indices]
            mean_val = np.nanmean(tril) # calculate the average of the lower triangle values
            arr[np.isnan(arr)] = mean_val # replace NaNs with the average
            np.fill_diagonal(arr, 1) # set diagonal values to 1 (not to the average)
        mat_list.append(arr)

    print(f'{len(mat_list)} connectivity matrices were successfully loaded')

    return mat_list

def get_ttest_inputs(pop1, pop2, females=False):
    ''' Gets the input for the ttest function. Takes 2 lists of matrices, and for each matrix, 
    extract the lower triangle values. Stack them in a 2D array.
    Returns two 2D arrays.
    
    Parameters
    ----------
    pop1 : str
        The name of the first group
    pop2 : str
        The name of the second group
    females : bool
        If True, only load the matrices of female mice. Default is False.
        
    Returns
    -------
    x1 : np.ndarray
        A 2D array of the matrices of the first group, shape (n_samples, n_edges x n_edges / 2)
    x2 : np.ndarray
        A 2D array of the matrices of the second group, shape (n_samples, n_edges x n_edges / 2)
    '''
    mat_list1 = get_grp_mat(pop1, females=females, z=True)
    mat_list2 = get_grp_mat(pop2, females=females, z=True)

    n_sample1 = len(mat_list1)
    n_sample2 = len(mat_list2)
    n_edges = mat_list1[0].shape[0]

    assert n_edges == mat_list2[0].shape[0], "Matrices have different shapes"

    tril_indices = np.tril_indices(n_edges, k=-1) # indices of the lower triangle (unique edges)
    n_datapoints = len(tril_indices[0]) # number of unique edges (not counting the diagonal)

    x1 = np.zeros((n_sample1, n_datapoints)) # initialize the result arrays
    for i, matrix in enumerate(mat_list1): 
        x1[i] = matrix[tril_indices] # fill the result array with the lower triangle values
    x2 = np.zeros((n_sample2, n_datapoints))
    for i, matrix in enumerate(mat_list2):
        x2[i] = matrix[tril_indices]

    return x1, x2

def back2mat(data, n_edges=26):
    ''' Convert a 1D array of the lower triangle values of a matrix to a full matrix

    Parameters
    ----------
    data : np.ndarray
        A 1D array of the lower triangle values of a matrix
    n_edges : int
        The number of edges of the matrix

    Returns
    -------
    mat : np.ndarray
        A 2D array of the full matrix
    '''

    tril_indices = np.tril_indices(n_edges, k=-1)
    mat = np.zeros((n_edges, n_edges))
    mat[tril_indices] = data
    mat = mat + mat.T
    np.fill_diagonal(mat, 1)

    return mat

def get_nbs_inputs(pop1, pop2, females=False):
    ''' Get the input for the NBS function

    Parameters
    ----------
    pop1 : str
        The name of the first group
    pop2 : str
        The name of the second group
    females : bool
        If True, only load the matrices of female mice. Default is False.

    Returns
    -------
    stack : np.ndarray
        A 3D array of the connectivty data from both groups (n_edges x n_edges x n_subjects)
    y_vec : np.ndarray
        A 1D array of the group labels
    '''

    mat_list1 = get_grp_mat(pop1, females=females, z=True)
    mat_list2 = get_grp_mat(pop2, females=females, z=True)
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

def zscore_mat(groups=['WT', '3xTgAD', 'TSPO_KO', '3xTgAD_TSPO_KO']):
    ''' Z-scores the connectivity matrices of each animal. Saves as souris_id_zscore.csv'''

    for pop in groups:
        path_list = glob.glob(f'data/{pop}/*.csv')
        for fname in path_list:
            if 'zscore' in fname:
                continue
            data = pd.read_csv(fname, index_col=0)
            mat = data.values
            # exctract lower triangle values (because symetric matrix -> redundant values + diagonal 
            # doesn't reflect actual connectivity)
            tril_indices = np.tril_indices(mat.shape[0], k=-1)
            tril = mat[tril_indices]
            z = (tril - np.nanmean(tril)) / np.nanstd(tril)
            # convert back to full matrix
            z = back2mat(z, n_edges=mat.shape[0])
            animal_id = fname.split('souris_')[1].split('.csv')[0]
            new_fname = f'data/{pop}/souris_{animal_id}_zscore.csv'
            pd.DataFrame(z, columns=data.columns, index=data.index).to_csv(new_fname, index=True)    

    return None

def txt_csv(groups=['WT', '3xTgAD', 'TSPO_KO', '3xTgAD_TSPO_KO']):
    ''' Convert all the txt files in a group to csv files.
    Only keeps the id in the filename for easier access.
    '''

    for pop in groups:
        path_list = glob.glob(f'data/{pop}/*.txt')
        for fname in path_list:
            data = pd.read_csv(fname, sep=r'\;', header=0, engine='python', index_col=0)
            animal_id = fname.split('Souris')[1].split('_')[1]
            new_fname = f'data/{pop}/souris_{animal_id}.csv'
            data.to_csv(new_fname, index=True)
            print(f'{fname} was converted to {new_fname}')

    return None