#!/bin/bash
#RUN FROM RECEPTOR_GROUP DIRECTORY
for dir in *_package; do
	receptor_id="${dir%_package}"
	echo "Running: $receptor_id"
	python3 calculate_docking_rmsd.py --workspace_dir "${dir}" --rec_id "$receptor_id" --ref_ligand "${dir}/ligand.pdb" --output_csv "${receptor_id}.csv" --data_csv "data_csv"
done
