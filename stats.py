from utils.conn import nbs_bct_corr_z, ttest_with_fdr, permutation_test_with_fdr
from utils.preproc import *
from utils.plotting import *
from utils.params import comparisons

################################################################################
# This script compares the average connectivity matrices of all the pairs of
# groups using three different methods: t-test, permutation test and NBS.
################################################################################
check_tree()

def run_stat_comp(comparisons, test='ttest', females=False):
    ''' Run a statistical comparison between the average connectivity matrices of two groups
    using either a t-test or a permutation test. The results are saved in .csv files and figures.
    '''

    for pop1, pop2 in comparisons:

        if test == 'ttest':
            raw_pvals, fdr_pvals = ttest_with_fdr(pop1, pop2)
            outdir = f'derivative/ttest/'
        elif test == 'permutations':
            raw_pvals, fdr_pvals = permutation_test_with_fdr(pop1, pop2)
            outdir = f'derivative/permutations/'

        if females:
            cmp_name = f'fem_{pop1}-vs-{pop2}'
            title = f'{test} {pop1} - {pop2}, FDR < 0.05 (females)'
        else:
            cmp_name = f'{pop1}-vs-{pop2}'
            title = f'{test} {pop1} - {pop2}, FDR < 0.05'

        np.savetxt(os.path.join(outdir, 'pvals', f'{cmp_name}_pval.csv'), fdr_pvals, delimiter=',')

        # plot (pop1 - pop2) * mask
        mask = fdr_pvals < 0.05
        mask = mask.astype(int)

        av1 = get_av_grp_mat(pop1)
        av2 = get_av_grp_mat(pop2)
        diff = av1 - av2
        diff = diff * mask

        fig = plot_mat(diff, title)
        fig.savefig(os.path.join(outdir, 'figures', f'{cmp_name}.png'), dpi=300) # dpi=300

        # with raw p-values
        np.savetxt(os.path.join(outdir, 'pvals', 'raw_pvals', f'{cmp_name}_raw_pval.csv'), raw_pvals, delimiter=',')
        mask = raw_pvals < 0.05
        mask = mask.astype(int)
        diff = av1 - av2
        diff = diff * mask

        fig = plot_mat(diff, f'{test} {pop1} - {pop2}, p < 0.05')
        fig.savefig(os.path.join(outdir, 'figures', 'raw_pvals', f'{cmp_name}_raw_pval.png'), dpi=300) # dpi=300
    
        
def run_nbs(comparisons, females=False):
    ''' Run a Network Based Statistics comparison between the average connectivity matrices of two groups.'''

    for pop1, pop2 in comparisons:

        stack, y, _, _ = get_nbs_inputs(pop1, pop2, females=females)
        pval, adj, null = nbs_bct_corr_z(stack, thresh=0.15, y_vec=y, k=1000)

        outdir = f'derivative/nbs/'
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
        fig3.savefig(os.path.join(outdir, 'figures', f'{cmp_name}.png'), dpi=300) # dpi=300


run_nbs(comparisons=comparisons, females=False)
run_nbs(comparisons=comparisons, females=True)
run_stat_comp(comparisons=comparisons, test='permutations', females=False)
run_stat_comp(comparisons=comparisons, test='permutations', females=True)
run_stat_comp(comparisons=comparisons, test='ttest', females=False)
run_stat_comp(comparisons=comparisons, test='ttest', females=True)

