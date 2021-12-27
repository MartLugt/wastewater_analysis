#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import pandas as pd
from math import floor, log10

from create_error_benchmarks import select_benchmark_genomes

from select_samples import filter_fasta, read_metadata


def main():
    parser = argparse.ArgumentParser(description="Create wastewater benchmarks.")
    parser.add_argument('-m, --metadata', dest='metadata', type=str, required=True, help="metadata tsv file for full sequence database")
    parser.add_argument('-s, --state', dest='state', type=str, default="Connecticut", help="sample location")
    parser.add_argument('-d, --date', dest='date', type=str, default="2021-02-11", help="sample date")
    parser.add_argument('-fr, --fasta_ref', dest='fasta_ref', required=True, type=str, help="fasta file representing full sequence database")
    parser.add_argument('-fv, --fasta_voc', dest='fasta_VOC', required=True, type=str, help="comma-separated list of fasta files for Variants Of Concern (VOC)")
    parser.add_argument('-fc, --fasta_con', dest='fasta_con', required=True, type=str, help="comma-seperated list of fasta files for contaminating viruses.")
    parser.add_argument('-o, --outdir', dest='outdir', required=True, type=str, help="output directory")
    parser.add_argument('--sars2_perc', dest='sars2_perc', required=True, type=str, help="comma-separated list of total SARS-CoV-2 abuncance frequency (%) to be simulated.")
    parser.add_argument('--total_sars2_cov', dest='total_sars2_cov', default=10000, type=int, help="total sequencing depth to be simulated")
    # parser.add_argument('--VOC_perc', dest='voc_perc', default=10, type=int, help="VOC coverage percentage (of total_sars2_cov)")
    parser.add_argument('--spike_only', action='store_true', help="simulate reads for spike region only")
    parser.add_argument('--no_errors', action='store_true', help="disable sequencing error simulation (ART)")
    parser.add_argument('--conts_amount', dest='cont_amount', required=True, type=int, help="divisor for contamination files coverage. (Amount of contaminants).")

    args = parser.parse_args()

    # create output directory
    try:
        os.makedirs(args.outdir)
    except FileExistsError:
        pass

    SARS_CoV_2_frequencies = args.sars2_perc.split(',')

    total_sars2_cov = args.total_sars2_cov
    VOC_files = args.fasta_VOC.split(',')
    VOC_names = [filepath.split('/')[-1] for filepath in VOC_files]
    exclude_list = [name.split('_')[0] for name in VOC_names]

    con_files = args.fasta_con.split(',')

    err_rate = " "
    if (args.no_errors):
        err_rate = " -qs 93 -qs2 93 -ir 0 -ir2 0 -dr 0 -dr2 0 "

    full_df = read_metadata(args.metadata)
    selection_df = select_benchmark_genomes(full_df, args.state, args.date,
                                            exclude_list)
    # filter fasta according to selection and write new fasta
    fasta_selection = args.outdir + "/sequences.fasta"
    filter_fasta(args.fasta_ref, fasta_selection, selection_df)
    print("Selected sequences written to {}".format(fasta_selection))
    # write corresponding metadata to tsv
    metadata_out = args.outdir + "/metadata.tsv"
    selection_df.to_csv(metadata_out, sep='\t', index=False)
    print("Metadata for selected sequences is in {}".format(metadata_out))

    # combine contamination fastas
    fasta_contamination = args.outdir + "/contamination.fasta"
    subprocess.check_call("cat {0} > {1}".format(" ".join(con_files), fasta_contamination), shell=True)


    if args.spike_only:
        # trim sequences to select spike region
        print("\nTrimming genomes around spike region (21063--25884)")
        trimmed_selection = args.outdir + "/sequences.trimmed.fasta"
        subprocess.check_call("reformat.sh in={} out={} fastawrap=0 overwrite=t forcetrimleft=21063 forcetrimright=25884".format(fasta_selection, trimmed_selection), shell=True)
        # also trim VOC sequences
        for filename in VOC_files:
            VOC_name = filename.rstrip('.fasta').split('/')[-1]
            trimmed_file = args.outdir + "/{}.trimmed.fasta".format(VOC_name)
            subprocess.check_call("reformat.sh in={} out={} fastawrap=0 overwrite=t forcetrimleft=21063 forcetrimright=25884".format(filename, trimmed_file), shell=True)
        fasta_selection = trimmed_selection
        print("\nSpike sequences ready\n")
    


    # simulate reads
    for sars2_freq in SARS_CoV_2_frequencies:
        VOC_cov = round(total_sars2_cov * 10 / 100, 2)
        background_cov = round((total_sars2_cov - VOC_cov) / len(selection_df.index), 2)
        contamination_cov = round((total_sars2_cov * (1 / (float(sars2_freq) / 100)) - total_sars2_cov) / args.cont_amount, )# len(con_files), 2)

        # cov = round(total_cov * (float(sars2_freq) / 100), 2)
        # VOC_cov = round(cov * 10 / 100, 2)
        # background_cov = round((cov - VOC_cov) / len(selection_df.index), 2)
        # contamination_cov = round((total_cov - cov) / len(con_files), 2)

        # simulate contamination reads
        print("Simulating contamination reads from {} at {}x coverage ".format(fasta_contamination, contamination_cov))
        subprocess.check_call("art_illumina -ss HS25 -rs 0 -i {0} -l 150 -f {1} -p -o {2}/contamination_ab{1}_ -m 250 -s 10 {3}"
            .format(fasta_contamination, contamination_cov, args.outdir, err_rate), shell=True)

        print("Simulating background reads from {} at {}x coverage ".format(fasta_selection, background_cov))
        subprocess.check_call("art_illumina -ss HS25 -rs 0 -i {0} -l 150 -f {1} -p -o {2}/background_ab{1}_ -m 250 -s 10 {3}"
            .format(fasta_selection, background_cov, args.outdir, err_rate), shell=True)

        for filename in VOC_files:
            VOC_name = filename.rstrip('.fasta').split('/')[-1]
            if args.spike_only:
                voc_fasta = args.outdir + "/{}.trimmed.fasta".format(VOC_name)
            else:
                voc_fasta = filename

            print("Simulating reads from {} at {}x coverage".format(VOC_name, VOC_cov))
            subprocess.check_call("art_illumina -ss HS25 -rs 0 -i {0} -l 150 -f {1} -p -o {2}/{3}_ab{1}_ -m 250 -s 10 {4}"
                .format(voc_fasta, VOC_cov, args.outdir, VOC_name, err_rate), shell=True)

            print("\nMerging fastqs...")
            subprocess.check_call("cat {0}/contamination_ab{4}_1.fq {0}/background_ab{1}_1.fq {0}/{2}_ab{3}_1.fq > {0}/tmp1.fq"
                .format(args.outdir, background_cov, VOC_name, VOC_cov, contamination_cov), shell=True)
            subprocess.check_call("cat {0}/contamination_ab{4}_2.fq {0}/background_ab{1}_2.fq {0}/{2}_ab{3}_2.fq > {0}/tmp2.fq"
                .format(args.outdir, background_cov, VOC_name, VOC_cov, contamination_cov), shell=True)

            print("Shuffling reads...")
            subprocess.check_call("shuffle.sh in={0}/tmp1.fq in2={0}/tmp2.fq out={0}/wwsim_{1}_ab{2}_1.fastq out2={0}/wwsim_{1}_ab{2}_2.fastq overwrite=t fastawrap=0"
                .format(args.outdir, VOC_name, sars2_freq), shell=True)
     
        print("\nBenchmarks with a SARS-CoV-2 frequency of {}% are ready!\n\n".format(sars2_freq))
    # clean up temporary files
    os.remove("{}/tmp1.fq".format(args.outdir))
    os.remove("{}/tmp2.fq".format(args.outdir))
    return


def round_sig(x, sig=2):
    if x == 0: return x
    return round(x, sig-int(floor(log10(abs(x))))-1)


if __name__ == "__main__":
    sys.exit(main())
