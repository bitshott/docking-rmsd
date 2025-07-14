from rdkit import Chem
import os
import seaborn as sns
import shutil
from pathlib import Path
import datetime
import shutil


reference_receptors_cnr2 = ['5ZTY', '6KPF', '6PT0', '8GUS', '8GUR', '8GUT', '8GUQ']
reference_receptors_cnr1 = ['5TGZ', '5U09', '5XR8', '5XRA', '6KQI', '7FEE', '6N4B', '6KPG', '8GHV']
files = ['sdf_filter/' + file for posix in os.walk('sdf_filter') for file in posix[2]]


md_input = Path(f'md_input_{datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")}')
rec_data = Path('data')

for file in files:
    file_receptor_name = file.split('/')[-1].split('_')[0].upper()
    print(f'\n{"-" * 50}\nProcessing receptor: {file_receptor_name}\n{"-" * 50}\n')
    suppl = Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
    for mol in suppl:
        mol_id = mol.GetProp('chembl_id')
        mol_id = mol_id.split('_')[0]
        path_to_mol_id = Path(md_input).joinpath(f'{mol_id}_{file_receptor_name}')

        # if path_to_mol_id.is_dir():
        #     print(f'Duplicate {file_receptor_name}: {mol_id}')
        #     path_to_mol_id = path_to_mol_id.parent / f"{path_to_mol_id.name}_{file_receptor_name}"

        path_to_mol_id.mkdir(parents=True)
        f = path_to_mol_id.joinpath(mol_id + '.sdf')
        print(f)
        with Chem.SDWriter(path_to_mol_id.joinpath(mol_id + '.sdf')) as writer:
            writer.write(mol)
            shutil.copy(rec_data.joinpath(f'{file_receptor_name}.pdb'),
                        path_to_mol_id.joinpath(f'{file_receptor_name}.pdb'))

print(f'{"-" * 50}\nZipping directory: {md_input}\n{"-" * 50}')
shutil.make_archive(md_input, 'zip', md_input)
