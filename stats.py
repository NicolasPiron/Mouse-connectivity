from utils.conn import nbs_bct_corr_z
from utils.preproc import *
from utils.plotting import *
from utils.params import comparisons

################################################################################
# This script compares the average connectivity matrices of all the pairs of
# groups
################################################################################
check_tree()

def run_nbs(comparisons, females=False):

    for pop1, pop2 in comparisons:

        stack, y, _, _ = get_nbs_inputs(pop1, pop2, females=females)
        pval, adj, null = nbs_bct_corr_z(stack, thresh=0.15, y_vec=y, k=1000)

        outdir = f'derivative/stats/'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if females:
            cmp_name = f'fem_{pop1}-vs-{pop2}'
        else:
            cmp_name = f'{pop1}-vs-{pop2}'

        # save the null distribution, p-values and adjacency matrix in .csv files
        np.savetxt(os.path.join(outdir, 'null', f'{cmp_name}_null.csv'), null, delimiter=',')
        np.savetxt(os.path.join(outdir, 'pvals', f'{cmp_name}_pval.csv'), pval, delimiter=',')
        np.savetxt(os.path.join(outdir, 'adjacency', f'{cmp_name}_adj.csv'), adj, delimiter=',')

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
        pval = np.min(pval)
        if females:
            fig3 = plot_mat(diff, f'females - {pop1} < {pop2} - pval={pval}')
        else:
            fig3 = plot_mat(diff, f'{pop1} < {pop2} - pval={pval}')
        fig3.savefig(os.path.join(outdir, 'figures', f'{cmp_name}.png')) # dpi=300


run_nbs(comparisons=comparisons, females=False)
run_nbs(comparisons=comparisons, females=True)

