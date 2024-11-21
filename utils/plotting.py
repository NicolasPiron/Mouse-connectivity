import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils.preproc import *


def plot_grp_box(df):
    ''' Plots a boxplot of the average connectivity values of each group'''

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.boxplot(data=df, x='group', y='average_connectivity', color='grey',
                width=0.3, fliersize=0.5, linewidth=1.3, showmeans=True,
                meanprops={"markerfacecolor": "blue", "markeredgecolor": "black"}, ax=ax)
    sns.swarmplot(data=df, x='group', y='average_connectivity', color='black',
                alpha=0.3, size=3, ax=ax)
    ax.set_title('Average connectivity values of each group')
    ax.set_ylabel('Average connectivity value')
    ax.set_xlabel('Group')
    sns.despine()
    plt.tight_layout()

    return fig

def plot_sexdiff_box(df):

    fig, ax = plt.subplots(figsize=(6, 4))  

    sns.boxplot(data=df, x='group', y='average_connectivity', hue='sex', palette='dark:grey',
                width=0.6, fliersize=0.5, linewidth=1.3, showmeans=True,
                meanprops={"markerfacecolor": "blue", "markeredgecolor": "black"}, ax=ax)
    ax.set_title('Average connectivity values of each group')
    ax.set_ylabel('Average connectivity value')
    ax.set_xlabel('Group')
    sns.despine()
    plt.tight_layout()

    return fig

def plot_female_box(df):

    fig, ax = plt.subplots(figsize=(6, 4))  

    sns.boxplot(data=df[df['sex']=='f'], x='group', y='average_connectivity', color='grey',
                width=0.3, fliersize=0.5, linewidth=1.3, showmeans=True,
                meanprops={"markerfacecolor": "blue", "markeredgecolor": "black"}, ax=ax)
    sns.swarmplot(data=df[df['sex']=='f'], x='group', y='average_connectivity', color='black',
                alpha=0.3, size=3, ax=ax)

    ax.set_title('Female average connectivity values of each group')
    ax.set_ylabel('Average connectivity value')
    ax.set_xlabel('Group')
    sns.despine()
    plt.tight_layout()

    return fig
####################################################################################################
# for matrix plots
####################################################################################################

def plot_sgl_mat(fname, title, fout):
    ''' Plots the matrix of a single animal'''

    data = pd.read_csv(fname, sep='\;', header=0, engine='python', index_col=0)
    fig = plot_mat(data, title)
    fig.savefig(fout)

def plot_grp_mat(pop, females=False):
    ''' Plots the average matrix of a group'''

    mat_list = get_grp_mat(pop, females=females)
    stack = np.stack((mat_list), axis=-1)
    av = np.mean(stack, axis=-1)
    if females:
        fig = plot_mat(av,f'{pop} - females')
        fig.savefig(f'derivative/average/females_{pop}.png')
    else:
        fig = plot_mat(av, pop)
        fig.savefig(f'derivative/average/{pop}.png')

def plot_mat(data, title):
    ''' Plots a matrix'''
    
    fig, ax = plt.subplots()
    sns.heatmap(data, ax=ax, cmap='coolwarm', center=0)
    ax.set_title(title)
    
    return fig
