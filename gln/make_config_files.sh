#!/bin/bash

# Define base paths and data types
save_path='/home/panwei/he000176/jinwen/gln/'
declare -a datatype=('LSimp' 'prs' 'real')

# Iterate over each data type
for dt in "${datatype[@]}"; do
  # Iterate from 1 to 9
  for j in {1..9}; do
    # Define and copy the global config file, then modify it
    f="${save_path}configs/base/${dt}/global_config_0.${j}.yaml"
    output_f="${save_path}model_output/${dt}/0.${j}"
    cp "${save_path}configs/global_config.yaml" "$f"
    sed -i "s|output_folder:|output_folder: $output_f|g" "$f"

    # Define and copy the input config file, then modify it
    f="${save_path}configs/base/${dt}/input_config_0.${j}.yaml"
    cov_f="/home/panwei/he000176/jinwen/data/pheno/training_data_for_base_models/pheno_cov_0.${j}.csv"
    cp "${save_path}configs/input_config.yaml" "$f"
    sed -i "s|input_source: /home/panwei/he000176/jinwen/data/pheno/training_data_for_base_models/pheno_cov_0.1.csv|input_source: $cov_f|g" "$f"

    # Define and copy the output config file, then modify it
    f="${save_path}configs/base/${dt}/output_config_0.${j}.yaml"
    pheno_f="/home/panwei/he000176/jinwen/data/pheno/training_data_for_base_models/pheno_${dt}_0.${j}.csv"
    cp "${save_path}configs/output_config.yaml" "$f"
    sed -i "s|output_source:|output_source: $pheno_f|g" "$f"
  done
done

echo "Configuration files have been prepared."
