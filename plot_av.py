from utils.preproc import *
from utils.plotting import *


# to plot the average connectivity matrix of each group
pop_vals = []
pop_names = []
for pop in os.listdir('data'):
    pop_name = pop.split('/')[0]
    if pop_name == '.DS_Store':
        continue
    grp_mat(pop)

    # mat_list = get_grp_mat(pop_name)
    
    # average_list = [mat.mean() for mat in mat_list]
    # print(average_list)
    # pop_vals.append(average_list)
    # pop_names.append(pop)
