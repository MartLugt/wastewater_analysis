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

outdir=benchmarks/${dataset}/out
mkdir -p ${outdir}

for VOC in P.1_EPI_ISL_1194849 B.1.1.7_EPI_ISL_1008906 B.1.351_EPI_ISL_1001460 B.1.617.2_EPI_ISL_1924762; do \
  for ab in 1 2 5 10 20 50 100; do \
    kallisto quant -t 12 -b ${num_bootstraps} -i reference_sets/${ref_dir}/sequences.kallisto_idx -o ${outdir}/${VOC}_ab${ab} benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_1.fastq benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_2.fastq | tee ${outdir}/${VOC}_ab${ab}.log # > ${outdir}/${VOC}_ab${ab}.log 2>&1;
    python pipeline/output_abundances.py -m $min_ab -o ${outdir}/${VOC}_ab${ab}/predictions_m${min_ab}.tsv --metadata reference_sets/${ref_dir}/metadata.tsv --voc B.1.1.7,B.1.351,B.1.617.2,P.1 ${outdir}/${VOC}_ab${ab}/abundance.tsv | tee -a ${outdir}/${VOC}_ab${ab}.log # >> ${outdir}/${VOC}_ab${ab}.log 2>&1;
  done;
done;
