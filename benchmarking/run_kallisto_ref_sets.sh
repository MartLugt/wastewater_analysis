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
min_ab=$4

HDF5_USE_FILE_LOCKING=FALSE # prevent HDF5 problems (https://github.com/pachterlab/kallisto/issues/197)

outdir=benchmarks/${dataset}/out_${ref_dir}
mkdir -p ${outdir}

for VOC in P.1_EPI_ISL_1194849 B.1.1.7_EPI_ISL_889440 B.1.351_EPI_ISL_1001460 B.1.617.2_EPI_ISL_1924762; do \
  for ab in 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100; do \
    kallisto quant -t 12 -b ${num_bootstraps} -i reference_sets/${ref_dir}/sequences.kallisto_idx -o ${outdir}/${VOC}_ab${ab} benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_1.fastq benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_2.fastq | tee ${outdir}/${VOC}_ab${ab}.log # > ${outdir}/${VOC}_ab${ab}.log 2>&1;
    python pipeline/output_abundances.py -m $min_ab -o ${outdir}/${VOC}_ab${ab}/predictions_m${min_ab}.tsv --metadata reference_sets/${ref_dir}/metadata.tsv --voc B.1.1.7,B.1.351,B.1.617.2,P.1 ${outdir}/${VOC}_ab${ab}/abundance.tsv | tee -a ${outdir}/${VOC}_ab${ab}.log # >> ${outdir}/${VOC}_ab${ab}.log 2>&1;
  done;
done;
