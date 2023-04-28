#!/bin/bash
## Code to run SDA analysis on isotypr metrics
## Run for x number of iterations...! 
for VARIABLE in {1..10}
do
mkdir /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/TensorDecomposition/${VARIABLE}
export OPENMP_NUM_THREADS=${NSLOTS:-1}
/apps/well/sda/1.1/./sda --data /gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/isotyper_metrics_filtered_FINAL_METRICS_TENSOR_FORMAT_1500_PRODUCTIVE.txt --N 194 --out /well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/TensorDecomposition/${VARIABLE} --num_openmp_threads ${NSLOTS:-1} --num_blocks ${NSLOTS:-1} --eigen_parallel true --num_comps 100 --remove_zero_comps true --max_iter 3000 --impute_missing true 
done 
