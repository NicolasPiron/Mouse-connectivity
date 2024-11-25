from itertools import combinations

groups = ['WT', '3xTgAD', 'TSPO_KO', '3xTgAD_TSPO_KO']
comparisons = list(combinations(groups, 2))

# ROIs acronyms
acronyms = ['RSplen-L', 'RSplen-R', 'Vis-L', 'Vis-R', 'PPAssoc-L', 'PPAssoc-R', 'Audit-L',
        'Audit-R', 'TAssoc-L', 'TAssoc-R', 'EC-L', 'EC-R', 'Olf-L', 'Olf-R', 'DG-L',
        'DG-R', 'Sub-L', 'Sub-R', 'DHipp-L', 'DHipp-R', 'VHipp-L', 'VHipp-R',
        'SNr-L', 'SNr-R', 'SNc-L', 'SNc-R']
