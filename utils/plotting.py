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

    mat_list1 = get_grp_mat(pop1, females=females)
    mat_list2 = get_grp_mat(pop2, females=females)

    def get_av(mat_list):
        stack = np.stack((mat_list), axis=-1)
        av = np.mean(stack, axis=-1)
        return av

    diff = get_av(mat_list1) - get_av(mat_list2)

    fig, ax = plt.subplots()
    sns.heatmap(diff, cmap='coolwarm', center=0, vmin=-1, vmax=1, ax=ax)
    ax.set_title(title)

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

def plot_sgl_mat(fname, title, fout):
    ''' Plots the matrix of a single animal'''

    data = pd.read_csv(fname, sep='\;', header=0, engine='python', index_col=0)
    fig = plot_mat(data, title)
    fig.savefig(fout)

def plot_grp_mat(pop, females=False, metric='mean'):
    ''' Plots the average matrix of a group'''

    # define paths and title depending on args. 
    path_map = {
        ('mean', True): ('mean connectivity - {pop} - females', 'derivative/average/raw/females_{pop}.png'),
        ('sd', True): ('connectivity sd - {pop} - females', 'derivative/average/sd/females_{pop}.png'),
        ('mean', False): ('mean connectivity - {pop}', 'derivative/average/raw/{pop}.png'),
        ('sd', False): ('connectivity sd - {pop}', 'derivative/average/sd/{pop}.png'),
    }
    key = (metric, females)
    if key in path_map:
        title_template, fout_template = path_map[key]
        title = title_template.format(pop=pop)
        fout = fout_template.format(pop=pop)

    mat_list = get_grp_mat(pop, females=females)
    stack = np.stack((mat_list), axis=-1)

    if metric == 'mean':
        vals = np.mean(stack, axis=-1)
        vmin, vmax = (-1, 1)
    elif metric == 'sd':
        vals = np.std(stack, axis=-1)
        vmin, vmax = (0, 0.5)

    fig = plot_mat(vals,title, vmin=vmin, vmax=vmax)
    fig.savefig(fout)


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
