import pandas as pd

data = pd.read_excel('data/code_animaux.xlsx')

print(data.head())
print(data.columns)
print(data['Sexe'])
print(data['numero'])
print(data['Lignée / Souche (Nom)'])