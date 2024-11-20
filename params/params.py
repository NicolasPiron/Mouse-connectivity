from itertools import combinations

groups = ['WT', '3xTgAD', 'TSPO_KO', '3xTgAD_TSPO_KO']
comparisons = list(combinations(groups, 2))

# comp_dict = {'WT': '3xTgAD',
#             'WT': 'TSPO_KO',
#             'WT': '3xTgAD_TSPO_KO',
#             '3xTgAD': 'TSPO_KO',
#             '3xTgAD': '3xTgAD_TSPO_KO',
#             '3xTgAD_TSPO_KO': 'TSPO_KO',
# }