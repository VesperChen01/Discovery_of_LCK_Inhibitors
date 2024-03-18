# pip install pandas chembl_webresource_client matplotlib rdkit tqdm numpy

import os 
import requests
import numpy as np
from tqdm.auto import tqdm
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import SaltRemover
import pandas as pd
import matplotlib.pyplot as plt

cwd = os.getcwd()
ids_df = pd.read_csv('ids.csv')

def get_active_compounds_from_chembl(uniprot_id):
    if pd.isnull(uniprot_id):
        print(f"Uniprot ID is nan, skipping.")
        return None
    target = new_client.target
    target_query = target.search(uniprot_id)
    if target_query:
        target_id = target_query[0]['target_chembl_id']
    else:
        print(f"No target found for Uniprot ID: {uniprot_id}")
        return None
    activity = new_client.activity
    activity_query = activity.filter(target_chembl_id=target_id).filter(standard_type__in=["IC50"])

    compounds = []
    for act in tqdm(activity_query):
        compound = {
            'chembl_id': act['molecule_chembl_id'],
            'smiles': act['canonical_smiles'],
        }
        compound[act['standard_type'].lower()] = act['standard_value']
        compound[act['standard_type'].lower() + '_unit'] = act['standard_units'] 
        compounds.append(compound)
    compounds_df = pd.DataFrame(compounds)
    compounds_df = compounds_df.fillna(0)
    return compounds_df

def clean_smiles(smiles_list):
    clean_smiles = []
    remover = SaltRemover.SaltRemover()
    for smiles in tqdm(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            mol = remover.StripMol(mol)
            clean_smiles.append(Chem.MolToSmiles(mol))
    return clean_smiles

def download_pdb_file(pdb_id):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    response = requests.get(url)
    file_path = f'{pdb_id}.pdb'
    with open(file_path, 'w') as f:
        f.write(response.text)
    print(f"The PDB file has been saved to: {os.path.join(os.getcwd(), file_path)}.")

for index, row in ids_df.iterrows():
    uniprot_id = row['uniprot_id']
    pdb_id = row['pdb_id']
    name = row['name']
    print(f"Processing：{uniprot_id} {pdb_id}")
    if pd.isnull(pdb_id):
        print(f"No PDB ID available, using only Uniprot ID: {uniprot_id} to obtain inhibitors.")
        active_compounds = get_active_compounds_from_chembl(uniprot_id)
    else:
        download_pdb_file(pdb_id)
        active_compounds = get_active_compounds_from_chembl(uniprot_id)
    if active_compounds is not None:

        # wash data
        smiles_list = active_compounds['smiles'].tolist()
        clean_smiles_list = clean_smiles(smiles_list)
        clean_compounds = active_compounds[active_compounds['smiles'].isin(clean_smiles_list)][['chembl_id', 'smiles', 'ic50', 'ic50_unit']]
        clean_compounds = clean_compounds.drop_duplicates(subset=['smiles'])
        clean_compounds['ic50'] = pd.to_numeric(clean_compounds['ic50'], errors='coerce')
        clean_compounds = clean_compounds.dropna(subset=['ic50'])
        clean_compounds = clean_compounds.sort_values(by='ic50')

        # calculate pIC50
        clean_compounds['ic50_M'] = clean_compounds['ic50'].astype(float) / 1_000_000_000
        clean_compounds = clean_compounds[clean_compounds['ic50_M'] != 0]
        clean_compounds['pIC50'] = -np.log10(clean_compounds['ic50_M'])
        clean_compounds = clean_compounds.drop(columns=['ic50_M'])

        # Save the DataFrame to a CSV file.
        csv_file_path = os.path.join(cwd, f'{name}.csv')
        clean_compounds.to_csv(csv_file_path, index=False)
        print(f"CSV file has been saved to：{csv_file_path}")

        # create pIC50 distribution chart
        plt.figure(figsize=(10, 5))
        plt.hist(clean_compounds['pIC50'], bins=50, color='blue', alpha=0.7)
        plt.xlabel('pIC50')
        plt.ylabel('number of molecules')
        plt.title(f'The distribution of pIC50 values for {name}')
        plt.savefig(f"{name}_pIC50_distribution.png")
        print(f"The pIC50 distribution chart has been saved to：{name}_pIC50_distribution.png")





 



