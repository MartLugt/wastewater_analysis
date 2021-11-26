#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-5:00
#SBATCH -p short
#SBATCH --mem=20G
#SBATCH -o logs/run_kallisto_%j.out
#SBATCH -e logs/run_kallisto_%j.err

dataset=$1
num_bootstraps=$2
ref_dir=$3
min_err=$4
ab=$5

HDF5_USE_FILE_LOCKING=FALSE # prevent HDF5 problems (https://github.com/pachterlab/kallisto/issues/197)

outdir=benchmarks/${dataset}/out_${ref_dir}
mkdir -p ${outdir}

for VOC in B.1.1.7_EPI_ISL_1008906; do \
  for err in 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0; do \
    kallisto quant -t 12 -b ${num_bootstraps} -i reference_sets/${ref_dir}/sequences.kallisto_idx -o ${outdir}/${VOC}_ab${ab}_er${err} benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_s0_i${err}_d${err}_1.fastq benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_s0_i${err}_d${err}_2.fastq | tee ${outdir}/${VOC}_ab${ab}_s0_i${err}_d${err}.log # > ${outdir}/${VOC}_ab${ab}.log 2>&1;
    python pipeline/output_abundances.py -m ${min_err} -o ${outdir}/${VOC}_ab${ab}_er${err}/predictions_m${min_err}.tsv --metadata reference_sets/${ref_dir}/metadata.tsv --voc B.1.1.7 ${outdir}/${VOC}_ab${ab}_er${err}/abundance.tsv | tee -a ${outdir}/${VOC}_ab${ab}_s0_i${err}_d${err}.log # >> ${outdir}/${VOC}_ab${ab}.log 2>&1;
  done;
done;
