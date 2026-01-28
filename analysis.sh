./perturbation_assignment/CIRI_assign.sh  --target_directory "/30tb/3tb/data/ratto/AB11_screening_CRISPR/analysis_11/."  --script_directory "/30tb/3tb/data/ratto/CIRI_bioinfo/perturbation_assignment/."   --guide_strategy 1   --matrix_filename "filtered_feature_bc_matrix.h5"   --threshold_a -1   --threshold_i -1

docker run -it  -v /30tb/3tb/data/ratto/CIRI_bioinfo/scripts:/scripts -v /30tb/3tb/data/ratto/AB11_screening_CRISPR/analysis_11/:/data docker.io/hedgelab/ciri_analysis bash

Rscript /scripts/anno_filter.R /data filtered_feature_bc_matrix.h5


./perturbation_assignment/CIRI_assign.sh  --target_directory "/30tb/3tb/data/ratto/AB012/analysis_12/."  --script_directory "/30tb/3tb/data/ratto/CIRI_bioinfo/perturbation_assignment/."   --guide_strategy 2   --matrix_filename "filtered_feature_bc_matrix.h5"   --threshold_a 23.5   --threshold_i -1