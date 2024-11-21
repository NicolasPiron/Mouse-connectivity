from utils.preproc import *
from utils.plotting import *
from utils.params import groups

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
fig1.savefig('derivative/average/average_connectivity.png', dpi=300)    
fig2.savefig('derivative/average/average_connectivity_sexdiff.png', dpi=300)
fig3.savefig('derivative/average/average_connectivity_female.png', dpi=300)

# plot the average connectivity matrix of each group, and with females only
pop_vals = []
pop_names = []
for pop in groups:

    plot_grp_mat(pop, females=False)
    plot_grp_mat(pop, females=True)
