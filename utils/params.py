from itertools import combinations

groups = ['WT', '3xTgAD', 'TSPO_KO', '3xTgAD_TSPO_KO']
comparisons = list(combinations(groups, 2))
