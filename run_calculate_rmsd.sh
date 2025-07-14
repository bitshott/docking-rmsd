#!/bin/bash
#RUN FROM RECEPTOR_GROUP DIRECTORY
LIG_DIR="data_lig"
SDF_DIR="sdf_filter"

# for dir in *_package; do
# 	receptor_id="${dir%_package}"
# 	echo "Running: $receptor_id"
# 	echo "${LIG_DIR}/${receptor_id}_ligand.pdb"
# 	python3 calculate_docking_rmsd.py --workspace_dir "${dir}" --rec_id "$receptor_id" --ref_ligand "${LIG_DIR}/${receptor_id}_ligand.pdb" --output_csv "${receptor_id}.csv" --data_csv "data_csv"
# done

mkdir $SDF_DIR
mv *_package.sdf $SDF_DIR

python3 ligand_packing.py