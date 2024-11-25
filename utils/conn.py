from __future__ import division
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from utils.preproc import get_ttest_inputs, back2mat

#############################################################################
# Permutation test and t-test with FDR correction
#############################################################################


def ttest_with_fdr(pop1, pop2, alpha=0.05):
    ''' Performs a t-test on each coordinate of two groups of matrices
    and corrects p-values using False Discovery Rate (FDR).
    
    Parameters:
    ----------
    pop1 : str
        Name of the first group.
    pop2 : str
        Name of the second group.
    alpha :float 
        Significance level for FDR correction.
    
    Returns:
    ----------
    raw_pvals : numpy.ndarray
        Raw p-values for each coordinate in shape (n_edges, n_edges).
    adj_pvals : numpy.ndarray
        FDR-adjusted p-values for each coordinate. Same shape as raw_pvals.
    '''

    arr1, arr2 = get_ttest_inputs(pop1, pop2)
    T_stats, raw_pvals = stats.ttest_ind(arr1, arr2, axis=0, equal_var=False, nan_policy='omit')

    # FDR correction using Benjamini-Hochberg
    _, fdr_pvals, _, _ = multipletests(raw_pvals, alpha=alpha, method='fdr_bh')

    raw_pvals = back2mat(raw_pvals) # convert to matrix
    fdr_pvals = back2mat(fdr_pvals)

    return raw_pvals, fdr_pvals

def permutation_test_with_fdr(pop1, pop2, n_permutations=10000, alpha=0.05):
    '''
    Performs a permutation test on each coordinate of two groups of matrices
    and corrects p-values using False Discovery Rate (FDR).
    
    Parameters:
    ----------
    pop1 : str
        Name of the first group.
    pop2 : str
        Name of the second group.
    n_permutations : int 
        Number of permutations for the test.
    alpha :float 
        Significance level for FDR correction.
        
    Returns:
    ----------
    raw_pvals : numpy.ndarray
        Raw p-values for each coordinate in shape (n_edges, n_edges).
    adj_pvals : numpy.ndarray
        FDR-adjusted p-values for each coordinate. Same shape as raw_pvals.
    '''

    arr1, arr2 = get_ttest_inputs(pop1, pop2)
    n_samples1 = arr1.shape[0]
    n_coords = arr1.shape[1]
    
    # Compute observed test statistic (difference in means)
    obs_stat = np.mean(arr1, axis=0) - np.mean(arr2, axis=0)
    
    # Combine data for permutation
    combined_data = np.vstack([arr1, arr2])
    n_combined = combined_data.shape[0]
    
    # Permutation test
    perm_stats = np.zeros((n_permutations, n_coords))
    for i in range(n_permutations):
        # Shuffle labels
        perm_indices = np.random.permutation(n_combined)
        perm_group1 = combined_data[perm_indices[:n_samples1], :]
        perm_group2 = combined_data[perm_indices[n_samples1:], :]
        
        # Compute permuted statistic
        perm_stats[i, :] = np.mean(perm_group1, axis=0) - np.mean(perm_group2, axis=0)
    
    # Calculate p-values
    raw_pvals = np.mean(np.abs(perm_stats) >= np.abs(obs_stat), axis=0)
    
    # FDR correction using Benjamini-Hochberg
    _, fdr_pvals, _, _ = multipletests(raw_pvals, alpha=alpha, method='fdr_bh')

    raw_pvals = back2mat(raw_pvals) # convert to matrix
    fdr_pvals = back2mat(fdr_pvals)

    return raw_pvals, fdr_pvals


#############################################################################
# NBS functions
#############################################################################

######################
# The BCT functions are stolen from https://pypi.org/project/bctpy/
######################

class BCTParamError(RuntimeError):
    pass

def binarize(W, copy=True):
    '''
    Binarizes an input weighted connection matrix.  If copy is not set, this
    function will *modify W in place.*

    Parameters
    ----------
    W : NxN np.ndarray
        weighted connectivity matrix
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : NxN np.ndarray
        binary connectivity matrix
    '''
    if copy:
        W = W.copy()
    W[W != 0] = 1
    return W


def get_components(A, no_depend=False):
    '''
    Returns the components of an undirected graph specified by the binary and
    undirected adjacency matrix adj. Components and their constitutent nodes
    are assigned the same index and stored in the vector, comps. The vector,
    comp_sizes, contains the number of nodes beloning to each component.

    Parameters
    ----------
    A : NxN np.ndarray
        binary undirected adjacency matrix
    no_depend : Any
        Does nothing, included for backwards compatibility

    Returns
    -------
    comps : Nx1 np.ndarray
        vector of component assignments for each node
    comp_sizes : Mx1 np.ndarray
        vector of component sizes

    Notes
    -----
    Note: disconnected nodes will appear as components with a component
    size of 1

    Note: The identity of each component (i.e. its numerical value in the
    result) is not guaranteed to be identical the value returned in BCT,
    matlab code, although the component topology is.

    Many thanks to Nick Cullen for providing this implementation
    '''

    if not np.all(A == A.T):  # ensure matrix is undirected
        raise BCTParamError('get_components can only be computed for undirected'
                            ' matrices.  If your matrix is noisy, correct it with np.around')
    
    A = binarize(A, copy=True)
    n = len(A)
    np.fill_diagonal(A, 1)

    edge_map = [{u,v} for u in range(n) for v in range(n) if A[u,v] == 1]
    union_sets = []
    for item in edge_map:
        temp = []
        for s in union_sets:

            if not s.isdisjoint(item):
                item = s.union(item)
            else:
                temp.append(s)
        temp.append(item)
        union_sets = temp

    comps = np.array([i+1 for v in range(n) for i in 
        range(len(union_sets)) if v in union_sets[i]])
    comp_sizes = np.array([len(s) for s in union_sets])

    return comps, comp_sizes


######################
# The NBS functinon is stolen from https://github.com/GidLev/NBS-correlation
######################

def nbs_bct_corr_z(corr_arr, thresh, y_vec, k=1000, extent=True, verbose=False):

    '''
    Performs the NBS for matrices [corr_arr] and vector [y_vec]  for a Pearson's r-statistic threshold of
    [thresh].

    Parameters
    ----------
    corr_arr : NxNxP np.ndarray
        matrix representing the correlation matrices population with P subjects. must be
        symmetric.

    y_vec : 1xP vector representing the behavioral/physiological values to correlate against

    thresh : float
        minimum Pearson's r-value used as threshold
    k : int
        number of permutations used to estimate the empirical null
        distribution, recommended - 10000
    verbose : bool
        print some extra information each iteration. defaults value = False

    Returns
    -------
    pval : Cx1 np.ndarray
        A vector of corrected p-values for each component of the networks
        identified. If at least one p-value is less than thres, the omnibus
        null hypothesis can be rejected at alpha significance. The null
        hypothesis is that the value of the connectivity from each edge has
        equal mean across the two populations.
    adj : IxIxC np.ndarray
        an adjacency matrix identifying the edges comprising each component.
        edges are assigned indexed values.
    null : Kx1 np.ndarray
        A vector of K sampled from the null distribution of maximal component
        size.

    Notes
    -----
    ALGORITHM DESCRIPTION
    The NBS is a nonparametric statistical test used to isolate the
    components of an N x N undirected connectivity matrix that differ
    significantly between two distinct populations. Each element of the
    connectivity matrix stores a connectivity value and each member of
    the two populations possesses a distinct connectivity matrix. A
    component of a connectivity matrix is defined as a set of
    interconnected edges.

    The NBS is essentially a procedure to control the family-wise error
    rate, in the weak sense, when the null hypothesis is tested
    independently at each of the N(N-1)/2 edges comprising the undirected
    connectivity matrix. The NBS can provide greater statistical power
    than conventional procedures for controlling the family-wise error
    rate, such as the false discovery rate, if the set of edges at which
    the null hypothesis is rejected constitues a large component or
    components.

    The NBS comprises fours steps:
    1. Perform a Pearson r test at each edge indepedently to test the
       hypothesis that the value of connectivity between each edge and an
       external variable, corelates across all nodes.
    2. Threshold the Pearson r-statistic available at each edge to form a set of
       suprathreshold edges.
    3. Identify any components in the adjacency matrix defined by the set
       of suprathreshold edges. These are referred to as observed
       components. Compute the size of each observed component
       identified; that is, the number of edges it comprises.
    4. Repeat K times steps 1-3, each time randomly permuting the extarnal
       variable vector and storing the size of the largest component
       identified for each permutation. This yields an empirical estimate
       of the null distribution of maximal component size. A corrected
       p-value for each observed component is then calculated using this
       null distribution.

    [1] Zalesky A, Fornito A, Bullmore ET (2010) Network-based statistic:
        Identifying differences in brain networks. NeuroImage.
        10.1016/j.neuroimage.2010.06.041

     Adopted from the python implementation of the BCT - https://sites.google.com/site/bctnet/, https://pypi.org/project/bctpy/
     Credit for implementing the vectorized version of the code to Gideon Rosenthal

    '''

    def corr_with_vars(x, y):
        # check correlation X -> M (Sobel's test)
        r, _ = stats.pearsonr(x, y)
        z = 0.5 * np.log((1 + r)/(1 - r))
        return z.item(0)

    ix, jx, nx = corr_arr.shape
    ny, = y_vec.shape

    if not ix == jx:
        raise ValueError('Matrices are not symmetrical')
    else:
        n = ix

    if nx != ny:
        raise ValueError('The [y_vec dimension must match the [corr_arr] third dimension')

    # only consider upper triangular edges
    ixes = np.where(np.triu(np.ones((n, n)), 1))

    # number of edges
    m = np.size(ixes, axis=1)

    # vectorize connectivity matrices for speed
    xmat = np.zeros((m, nx))

    for i in range(nx):
        xmat[:, i] = corr_arr[:, :, i][ixes].squeeze()
    del corr_arr

    # perform pearson corr test at each edge

    z_stat = np.apply_along_axis(corr_with_vars, 1, xmat, y_vec)
    print('z_stat: ', z_stat)

    # threshold
    ind_r, = np.where(z_stat > thresh)

    if len(ind_r) == 0:
        raise ValueError("Unsuitable threshold")

    # suprathreshold adjacency matrix
    adj = np.zeros((n, n))
    adjT = np.zeros((n, n))

    if extent:
        adj[(ixes[0][ind_r], ixes[1][ind_r])] = 1
        adj = adj + adj.T  # make symmetrical
    else:
        adj[(ixes[0][ind_r], ixes[1][ind_r])] = 1
        adj = adj + adj.T  # make symmetrical
        adjT[(ixes[0], ixes[1])] = z_stat
        adjT = adjT + adjT.T  # make symmetrical
        adjT[adjT <= thresh] = 0

    a, sz = get_components(adj)

    # convert size from nodes to number of edges
    # only consider components comprising more than one node (e.g. a/l 1 edge)
    ind_sz, = np.where(sz > 1)
    ind_sz += 1
    nr_components = np.size(ind_sz)
    sz_links = np.zeros((nr_components,))
    for i in range(nr_components):
        nodes, = np.where(ind_sz[i] == a)
        if extent:
            sz_links[i] = np.sum(adj[np.ix_(nodes, nodes)]) / 2
        else:
            sz_links[i] = np.sum(adjT[np.ix_(nodes, nodes)]) / 2

        adj[np.ix_(nodes, nodes)] *= (i + 2)

    # subtract 1 to delete any edges not comprising a component
    adj[np.where(adj)] -= 1

    if np.size(sz_links):
        max_sz = np.max(sz_links)
    else:
        # max_sz=0
        raise ValueError('True matrix is degenerate')
    print('max component size is %i' % max_sz)

    # estimate empirical null distribution of maximum component size by
    # generating k independent permutations
    print('estimating null distribution with %i permutations' % k)

    null = np.zeros((k,))
    hit = 0

    ind_shuff1 = np.array(range(0, y_vec.__len__()))
    ind_shuff2 = np.array(range(0, y_vec.__len__()))

    for u in range(k):
        # randomize
        np.random.shuffle(ind_shuff1)
        np.random.shuffle(ind_shuff2)
        # perform pearson corr test at each edge
        z_stat_perm = np.apply_along_axis(corr_with_vars, 1, xmat, y_vec[ind_shuff1])

        ind_r, = np.where(z_stat_perm > thresh)

        adj_perm = np.zeros((n, n))

        if extent:
            adj_perm[(ixes[0][ind_r], ixes[1][ind_r])] = 1
            adj_perm = adj_perm + adj_perm.T
        else:
            adj_perm[(ixes[0], ixes[1])] = z_stat_perm
            adj_perm = adj_perm + adj_perm.T
            adj_perm[adj_perm <= thresh] = 0

        a, sz = get_components(adj_perm)

        ind_sz, = np.where(sz > 1)
        ind_sz += 1
        nr_components_perm = np.size(ind_sz)
        sz_links_perm = np.zeros((nr_components_perm))
        for i in range(nr_components_perm):
            nodes, = np.where(ind_sz[i] == a)
            sz_links_perm[i] = np.sum(adj_perm[np.ix_(nodes, nodes)]) / 2

        if np.size(sz_links_perm):
            null[u] = np.max(sz_links_perm)
        else:
            null[u] = 0

        # compare to the true dataset
        if null[u] >= max_sz:
            hit += 1
        if verbose:
            print('permutation %i of %i.  Permutation max is %s.  Observed max'
                  ' is %s.  P-val estimate is %.3f') % (
                u, k, null[u], max_sz, hit / (u + 1))
        elif (u % (k / 10) == 0 or u == k - 1):
            print('permutation %i of %i.  p-value so far is %.3f' % (u, k,
                                                                     hit / (u + 1)))
    pvals = np.zeros((nr_components,))
    # calculate p-vals
    for i in range(nr_components):
        pvals[i] = np.size(np.where(null >= sz_links[i])) / k

    return pvals, adj, null