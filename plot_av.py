from utils.preproc import *
from utils.plotting import *
from utils.params import groups, comparisons

# This script compares the average connectivity matrices of all the pairs of 
# groups using three different methods: t-test, permutation test and NBS.


# preprocessing
check_tree()
if not len(glob.glob('data/*/*souris*.csv')) == 0:
    print('Data already preprocessed')
else:
    txt_csv()
    print('Successfully transformed .txt to .csv')
if not os.path.exists('data/all_df.csv'):
    all_df()

df = pd.read_csv('data/all_df.csv')
fig1 = plot_grp_box(df)
fig2 = plot_sexdiff_box(df)
fig3 = plot_female_box(df)
fig1.savefig('derivative/average/boxplot/average_connectivity.png', dpi=300)    
fig2.savefig('derivative/average/boxplot/average_connectivity_sexdiff.png', dpi=300)
fig3.savefig('derivative/average/boxplot/average_connectivity_female.png', dpi=300)

# plot the average connectivity matrix of each group, and with females only

metrics = ['mean', 'sd']
females = [True, False]
for metric in metrics:
    for female in females:
        for pop in groups:
            plot_grp_mat(pop, females=female, metric=metric)
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