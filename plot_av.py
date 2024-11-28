from utils.preproc import *
from utils.plotting import *
from utils.params import groups, comparisons

# This script compares the average connectivity matrices of all the pairs of 
# groups using three different methods: t-test, permutation test and NBS.
pre_run_check() # check if the data is available and preprocessed

df = pd.read_csv('data/all_df.csv')
fig1 = plot_grp_box(df)
fig2 = plot_sexdiff_box(df)
fig3 = plot_female_box(df)
fig1.savefig('derivative/average/boxplot/average_connectivity.png', dpi=300)    
fig2.savefig('derivative/average/boxplot/average_connectivity_sexdiff.png', dpi=300)
fig3.savefig('derivative/average/boxplot/average_connectivity_female.png', dpi=300)

# plot the average connectivity matrix of each group, and with females only

zscore = [True, False]
females = [True, False]
for z in zscore:
    for female in females:
        for pop in groups:
            plot_grp_mat(pop, females=female, z=z)
            plt.close('all')

for female in females:
    for (pop1, pop2) in comparisons:
        if female == True:
            fname = f'derivative/average/diff/{pop1}_{pop2}.png'
        elif female == False:
            fname = f'derivative/average/diff/females_{pop1}_{pop2}.png'
        fig = plot_diff_group_mat(pop1, pop2, females=female)
        fig.savefig(fname, dpi=300)
        plt.close('all')