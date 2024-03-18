import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from scipy.spatial.distance import cdist

# Function to convert RDKit ExplicitBitVect to numpy array
def bitvector_to_array(bitvector):
    bitlist = list(bitvector.ToBitString())
    return np.array([int(b) for b in bitlist], dtype=np.uint8)

# 读取包含SMILES的Excel文件
excel_file_path = 'LCK_P06239.xlsx'
data = pd.read_excel(excel_file_path)

# 从SMILES创建分子对象
data['Molecule'] = data['SMILES'].apply(lambda x: Chem.MolFromSmiles(x))

# 计算分子指纹
data['Fingerprint'] = data['Molecule'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=512))

# 转换指纹为向量
fingerprints = np.array([bitvector_to_array(fp) for fp in data['Fingerprint']])

# 使用Tanimoto距离计算相似度矩阵
distances = cdist(fingerprints, fingerprints, metric='jaccard')

# 使用层次聚类进行聚类
from sklearn.cluster import AgglomerativeClustering

# 聚类
cluster = AgglomerativeClustering(n_clusters=50, affinity='precomputed', linkage='average') #聚成50类
clusters = cluster.fit_predict(distances)

# 将聚类结果添加到原始数据
data['Cluster'] = clusters

# 保存结果到新的Excel文件
output_file_path = 'LCK_P06239_clustered.xlsx'
data.to_excel(output_file_path)
