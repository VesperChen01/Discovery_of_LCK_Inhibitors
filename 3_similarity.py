from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit.Chem import Draw, rdFMCS
from IPython.display import display 

# 读取excel文件
df = pd.read_excel('/home/szy/proj/crf/LCK_VS_anal/smi/output.xlsx')

# 获取分子smiles
query_smiles = df['Query_smiles'].dropna().values
data_smiles = df['Data_smiles'].dropna().unique()

# 计算所有Data_smiles的分子指纹
data_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 2) for smi in data_smiles]

# 初始化最大相似度和最相似分子列
df['Max_similarity'] = 0.0
df['Max_sim_smiles'] = ''

num_queries = len(query_smiles)
num_data = len(data_smiles)


tanimoto_matrix = np.zeros((num_queries, num_data))

for i, q_smile in enumerate(query_smiles):
    # 检查q_smile是否为字符串
    if not isinstance(q_smile, str):
        print(f"Skipping non-string SMILES at index {i}: {q_smile}")
        continue

    # 创建分子指纹
    query_mol = Chem.MolFromSmiles(q_smile)
    if query_mol is None: 
        print(f"Invalid SMILES at index {i}: {q_smile}")
        continue  # 如果smiles无效，跳过

    query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2)

    # 计算Tanimoto系数
    tcs = DataStructs.BulkTanimotoSimilarity(query_fp, data_fps)
    
    # 找到最大的相似度及其索引
    max_sim_idx = np.argmax(tcs)
    max_similarity = tcs[max_sim_idx]

    df.loc[i, 'Max_similarity'] = max_similarity
    df.loc[i, 'Max_sim_smiles'] = data_smiles[max_sim_idx]
    tanimoto_matrix[i, :] = tcs
# 保存结果
df.to_excel('/home/szy/proj/crf/LCK_VS_anal/smi/LCK-smi.xlsx')

sns.heatmap(tanimoto_matrix, annot=True, fmt=".2f", cmap="YlGnBu")
plt.title('Tanimoto Similarity Heatmap')
plt.xlabel('Data Molecules')
plt.ylabel('Query Molecules')
plt.show()


# 选择一个样本来展示
sample_df = df.sample(n=5)  # 随机选择5个查询分子

# 准备一个图像列表
mols_imgs = []

for idx, row in sample_df.iterrows():
    # 分子和相似度
    query_smi = row['Query_smiles']
    target_smi = row['Max_sim_smiles']
    similarity = row['Max_similarity']

    # 分子对象
    query_mol = Chem.MolFromSmiles(query_smi)
    target_mol = Chem.MolFromSmiles(target_smi)

    # 找到最大公共子结构
    mcs = rdFMCS.FindMCS([query_mol, target_mol])
    common_substructure = Chem.MolFromSmarts(mcs.smartsString)

    # 高亮最大公共子结构
    query_atom_indices = query_mol.GetSubstructMatch(common_substructure)
    target_atom_indices = target_mol.GetSubstructMatch(common_substructure)

    # 生成图像
    img = Draw.MolsToGridImage(
        [query_mol, target_mol],
        molsPerRow=2,
        subImgSize=(300, 300),
        highlightAtomLists=[query_atom_indices, target_atom_indices],
        legends=[f'Query (Sim: {similarity:.2f})', 'Most Similar']
    )
    mols_imgs.append(img)

# 显示图像
for img in mols_imgs:
    display(img)


