import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import seaborn as sns

# 加载数据
file_path = "/home/szy/proj/crf/LCK_VS_anal/dockingscore.xlsx"
data = pd.read_excel(file_path)

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol_weight = Descriptors.MolWt(mol)
        logP = Descriptors.MolLogP(mol)
        return mol_weight, logP
    return None, None

# 计算分子量和LogP
data['MolWt'], data['LogP'] = zip(*data['SMILES'].map(calculate_properties))  

# 添加标签以区分活性和装饰化合物
data['Label'] = data['Name'].apply(lambda x: 'Active' if x.startswith('P') else 'Decoy')

# # 箱型图
# fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# sns.boxplot(x=data['Label'], y=data['MolWt'], ax=ax[0], palette="Set2")
# ax[0].set_title('Molecular Weight Boxplot')

# sns.boxplot(x=data['Label'], y=data['LogP'], ax=ax[1], palette="Set2")
# ax[1].set_title('LogP Boxplot')

# plt.tight_layout()
# plt.show()

# 小提琴图
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

sns.violinplot(x=data['Label'], y=data['MolWt'], ax=ax[0], palette="Set2", inner="quartile")
ax[0].set_title('Molecular Weight Violin plot')

sns.violinplot(x=data['Label'], y=data['LogP'], ax=ax[1], palette="Set2", inner="quartile")
ax[1].set_title('LogP Violin plot')

plt.tight_layout()
plt.show()
