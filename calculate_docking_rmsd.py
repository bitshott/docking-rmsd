import os
from rdkit import Chem
import pandas as pd
import glob
from rdkit.Chem import AllChem, rdFMCS, SDWriter, rdMolAlign
from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolops import SanitizeMol
import argparse
import sys
from copy import deepcopy

parser = argparse.ArgumentParser(
    prog='calcDockRmsd'
)

parser.add_argument('--workspace_dir')
parser.add_argument('--ref_ligand')
parser.add_argument('--output_csv')
parser.add_argument('--data_csv')
parser.add_argument('--rec_id')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()
workspace_dir = args.workspace_dir
data_csv = args.data_csv
sdf_files = [file for temp in os.walk(workspace_dir) for file in temp[2] if file.endswith(".sdf")]
sdf_files = sorted(sdf_files, key=lambda x: int(x.split(".")[0]))
print(sdf_files)
ref_csv = [file for temp in os.walk(data_csv) for file in temp[2]]

id_best_rmsd = dict()
ref_mol = Chem.MolFromPDBFile(args.ref_ligand)

csv = f'{workspace_dir.split("/")[-1].split("_")[0]}.csv'
df = pd.read_csv(data_csv + "/" + csv, names=['id'])
print(df)

for sdf in sdf_files:
    print(f'{"-" * 30}\nSTARTED FILE:{sdf}\n{"-" * 30}')
    exp_id = int(sdf.split('.')[0])
    id_best_rmsd[exp_id] = [None, float('inf')]
    sdf_dir = workspace_dir + "/" + sdf
    try:
        suppl = SDMolSupplier(sdf_dir, sanitize=False, removeHs=False)
    except OSError as e:
        print(e)
        continue

    for mol in suppl:
        temp_mol = deepcopy(mol)
        mols = []
        try:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL)
            mols = [ref_mol, mol]
            mcs = rdFMCS.FindMCS(mols)
        except Exception as e:
            print(f'{"-" * 40}\nSANITIZE ERROR, SWITCHING TO NO-SANITIZE MODE {sdf}:{e}\n{"-" * 40}')
            mols = [ref_mol, temp_mol]
            mcs = rdFMCS.FindMCS(mols)

        pattern = Chem.MolFromSmarts(mcs.smartsString)
        ref = mols[0]
        mol_from_list = mols[-1]
        refMatch = ref.GetSubstructMatch(pattern)
        molMatch = mol_from_list.GetSubstructMatch(pattern)
        rms = Chem.rdMolAlign.CalcRMS(mol_from_list, ref_mol, map=[list(zip(molMatch, refMatch))])
        affinity_score = float(mol.GetProp('minimizedAffinity'))

        if rms <= 3.0 and affinity_score <= -8.0:
            if rms < id_best_rmsd.get(exp_id)[1]:
                id_best_rmsd[exp_id] = [mol_from_list, rms]
                print(f'{"-" * 20}\nAdded rmsd\n{"-" * 20}')
                print(Chem.MolToSmiles(id_best_rmsd[exp_id][0]))

    if id_best_rmsd[exp_id][0]:
        id_best_rmsd[exp_id][0].SetProp("chembl_id", df.iloc[exp_id, 0])
        sdf_file = open(f'{workspace_dir}.sdf', "a")
        writer = SDWriter(sdf_file)
        writer.write(mol=id_best_rmsd[exp_id][0])
