from utils.plotting import *
from utils.preproc import pre_run_check, zscore_mat


pre_run_check()
desc = pd.read_csv('data/all_df.csv')
ids = desc['id']
for id in ids:
    plot_sgl_mat(id, f'z-scored matrix mouse {id}', f'derivative/individuals/souris_{id}.png')
    plt.close('all')