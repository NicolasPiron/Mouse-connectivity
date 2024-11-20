import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils.preproc import *


def plot_averrage():
    ...

####################################################################################################
# for matrix plots
####################################################################################################

def sbj_mat(fname, title, fout):
    ''' Plots the matrix of a single animal'''

    data = pd.read_csv(fname, sep='\;', header=0, engine='python', index_col=0)
    fig = plot_mat(data, title)
    fig.savefig(fout)

def grp_mat(pop):
    ''' Plots the average matrix of a group'''

    mat_list = get_grp_mat(pop)
    stack = np.stack((mat_list), axis=-1)
    av = np.mean(stack, axis=-1)
    fig = plot_mat(av, pop)
    fig.savefig(f'derivative/average/{pop}.png')

def plot_mat(data, title):
    ''' Plots a matrix'''
    
    fig, ax = plt.subplots()
    sns.heatmap(data, ax=ax, cmap='coolwarm', center=0)
    ax.set_title(title)
    
    return fig
