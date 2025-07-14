LIG_DIR := data_lig
SDF_DIR := sdf_filter

.PHONY: RMSD_FILTER CREATE_SDF_DIR RMSD_PACKING

ALL: RMSD_FILTER CREATE_SDF_DIR RMSD_PACKING

RMSD_FILTER: 
	@for dir in *_package; do \
		receptor_id="$${dir%_package}"; \
		echo "Running: $${receptor_id}"; \
		echo "${LIG_DIR}/$${receptor_id}_ligand.pdb"; \
		python3 calculate_docking_rmsd.py \
			--workspace_dir "$${dir}" \
			--rec_id "$${receptor_id}" \
			--ref_ligand "${LIG_DIR}/$${receptor_id}_ligand.pdb" \
			--output_csv "$${receptor_id}.csv" \
			--data_csv "data_csv"; \
		done
CREATE_SDF_DIR:
	@mkdir $(SDF_DIR)
	@cp *_package.sdf $(SDF_DIR)

RMSD_PACKING:
	@python3 ligand_packing.py

REMOVE_SDF_DIR:
	@rm -r $(SDF_DIR)

REMOVE_FILTERED_SDF:
	@rm *_package.sdf
