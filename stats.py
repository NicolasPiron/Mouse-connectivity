from utils.conn import nbs_bct_corr_z
from utils.preproc import *
from utils.plotting import *
from params.params import comparisons

################################################################################
# This script compares the average connectivity matrices of all the pairs of
# groups
################################################################################

for pop1, pop2 in comparisons:

    stack, y, npop1, npop2 = get_nbs_inputs(pop1, pop2)
    pval, adj, null = nbs_bct_corr_z(stack, thresh=0.15, y_vec=y, k=1000)

    outdir = f'derivative/stats/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cmp_name = f'{pop1}-vs-{pop2}'

    # group 1 average * adj
    print(f'--- multipliying {pop1} by adj ---')
    mat_list1 = get_grp_mat(pop1)
    stack1 = np.stack((mat_list1), axis=-1)
    av1 = np.mean(stack1, axis=-1)
    av1 = av1 * adj

    # group 2 average * adj
    print(f'--- multipliying {pop2} by adj ---')
    mat_list2 = get_grp_mat(pop2)
    stack2 = np.stack((mat_list2), axis=-1)
    av2 = np.mean(stack2, axis=-1)
    av2 = av2 * adj

    # (average group 1 - average group 2) * adj
    print(f'--- multipliying ({pop1} - {pop2}) by adj ---')
    diff = av1 - av2 

    #fig1 = plot_mat(av1, f'{pop1} - {cmp_name}')
    #fig2 = plot_mat(av2, f'{pop2} - {cmp_name}')
    pval = np.min(pval)
    fig3 = plot_mat(diff, f'{pop1} < {pop2} - pval={pval}')
    #fig1.savefig(os.path.join(outdir, f'{pop1}__{cmp_name}.png'))
    #fig2.savefig(os.path.join(outdir, f'{pop2}__{cmp_name}.png'))
    fig3.savefig(os.path.join(outdir, f'{cmp_name}.png')) # dpi=300

    # print(pval)
    # print(adj)
    # print(null)