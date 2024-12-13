import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils.preproc import *
from utils.params import acronyms as ac

def plot_diff_group_mat(pop1, pop2, females=False):
    ''' Plots the difference between the average connectivity matrices of two groups'''

    title = f'average difference of {pop1} - {pop2}'
    if females:
        title = f'average difference of {pop1} - {pop2} (females)'

    mat_list1 = get_grp_mat(pop1, females=females, z=True)
    mat_list2 = get_grp_mat(pop2, females=females, z=True)

    def get_av(mat_list):
        stack = np.stack((mat_list), axis=-1)
        av = np.mean(stack, axis=-1)
        return av

    diff = get_av(mat_list1) - get_av(mat_list2)
    fig = plot_mat(diff, title, vmin=None, vmax=None)

    return fig

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
    ''' Plots a boxplot of the average connectivity values of each group, and shows the difference
    between males and females'''

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
    ''' Plots a boxplot of the average connectivity values of each group, with only 
    the females'''

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

def plot_sgl_mat(id, title, fout):
    ''' Plots the matrix of a single animal'''
 
    data = get_single_mat(id, z=True)
    fig = plot_mat(data, title, vmin=None, vmax=None)
    fig.savefig(fout)

def plot_grp_mat(pop, females=False, z=False):
    ''' Plots the average matrix of a group'''

    # define paths and title depending on args. 
    path_map = {
        (False, True): ('raw mean connectivity - {pop} - females', 'derivative/average/raw/females_{pop}.png'),
        (True, True): ('z-scored mean connectivity - {pop} - females', 'derivative/average/zscored/females_{pop}.png'),
        (False, False): ('raw mean connectivity - {pop}', 'derivative/average/raw/{pop}.png'),
        (True, False): ('z-scored mean connectivity - {pop}', 'derivative/average/zscored/{pop}.png'),
    }
    key = (z, females)
    if key in path_map:
        title_template, fout_template = path_map[key]
        title = title_template.format(pop=pop)
        fout = fout_template.format(pop=pop)

    if not z:
        mat_list = get_grp_mat(pop, females=females, z=False)
        vmin, vmax = (-1, 1)
    elif z:
        mat_list = get_grp_mat(pop, females=females, z=True)
        vmin, vmax = (None, None)
    
    stack = np.stack((mat_list), axis=-1)
    vals = np.mean(stack, axis=-1)
    fig = plot_mat(vals,title, vmin=vmin, vmax=vmax)
    fig.savefig(fout, dpi=300)


def plot_mat(data, title, vmin=-1, vmax=1): 
    ''' Plots a matrix'''
    
    fig, ax = plt.subplots(figsize=(7.5, 6))
    sns.heatmap(data, ax=ax, cmap='coolwarm', center=0,
                xticklabels=ac, yticklabels=ac,
                vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    plt.tight_layout()
    
    return fig
