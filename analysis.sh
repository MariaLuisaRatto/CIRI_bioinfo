#first step done with baryon (opens docker behind the scenes)
./perturbation_assignment/CIRI_assign.sh  --target_directory "/30tb/3tb/data/ratto/AB11_screening_CRISPR/analysis_11/."  --script_directory "/30tb/3tb/data/ratto/CIRI_bioinfo/perturbation_assignment/."   --guide_strategy 1   --matrix_filename "filtered_feature_bc_matrix.h5"   --threshold_a -1   --threshold_i -1

#run docker before the following steps
docker run -it  -v /30tb/3tb/data/ratto/CIRI_bioinfo/scripts:/scripts -v /30tb/3tb/data/ratto/AB11_screening_CRISPR/analysis_11/:/data docker.io/hedgelab/ciri_analysis bash

Rscript /scripts/anno_filter.R /data/ filtered_feature_bc_matrix.h5

Rscript /scripts/CIRI_load.R /data/ annotated_matrix.csv 0.5e-4 

Rscript /scripts/CIRI_sec.R /data/ 5 "NTCa-NA"

Rscript /scripts/CIRI_sec_subclusters.R /data/ 5 SOX2 muscle 1e-3

Rscript /scripts/CIRI_pseudotime.R /data/ cds_sample_1_ordered.RData pseudotime_muscle_1.csv muscle_1 TRUE "NTCa-NA" 8

Rscript /scripts/CIRI_pseudotime.R /data/ cds_sample_2_ordered.RData pseudotime_muscle_2.csv muscle_2 TRUE "NTCa-NA" 8

Rscript /scripts/guide_genes_expr.R /data/ 


#same for experiment with 2 guides
./perturbation_assignment/CIRI_assign.sh  --target_directory "/30tb/3tb/data/ratto/AB012/analysis_12/."  --script_directory "/30tb/3tb/data/ratto/CIRI_bioinfo/perturbation_assignment/."   --guide_strategy 2   --matrix_filename "filtered_feature_bc_matrix.h5"   --threshold_a 23.5   --threshold_i -1

docker run -it  -v /30tb/3tb/data/ratto/CIRI_bioinfo/scripts:/scripts -v /30tb/3tb/data/ratto/AB012/analysis_12/:/data docker.io/hedgelab/ciri_analysis bash

Rscript /scripts/anno_filter.R /data/ filtered_feature_bc_matrix.h5

Rscript /scripts/CIRI_load.R /data/ annotated_matrix.csv 1e-5

Rscript /scripts/CIRI_sec.R /data/ 3 "NTCa_1A;NTCa_1B-NA"

Rscript /scripts/CIRI_sec_subclusters.R /data/ 3 SOX2 muscle 1e-4

Rscript /scripts/CIRI_pseudotime.R /data/ cds_sample_1_ordered.RData pseudotime_muscle_1.csv muscle_1 TRUE "NTCa_1A;NTCa_1B-NA" 8

Rscript /scripts/CIRI_pseudotime.R /data/ cds_sample_2_ordered.RData pseudotime_muscle_2.csv muscle_2 TRUE "NTCa_1A;NTCa_1B-NA" 8

Rscript /scripts/guide_genes_expr.R /data/ 