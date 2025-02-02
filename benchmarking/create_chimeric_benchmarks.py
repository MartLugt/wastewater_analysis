#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import pandas as pd
from math import floor, log10

from select_samples import filter_fasta, read_metadata

from create_error_benchmarks import select_benchmark_genomes

from rpy2 import robjects
from rpy2.robjects.packages import importr

def main():
    parser = argparse.ArgumentParser(description="Create wastewater benchmarks.")
    parser.add_argument('-m, --metadata', dest='metadata', type=str, required=True, help="metadata tsv file for full sequence database")
    parser.add_argument('-s, --state', dest='state', type=str, default="Connecticut", help="sample location")
    parser.add_argument('-d, --date', dest='date', type=str, default="2021-02-11", help="sample date")
    parser.add_argument('-fr, --fasta_ref', dest='fasta_ref', required=True, type=str, help="fasta file representing full sequence database")
    parser.add_argument('-fv, --fasta_voc', dest='fasta_VOC', required=True, type=str, help="comma-separated list of fasta files for Variants Of Concern (VOC)")
    parser.add_argument('-o, --outdir', dest='outdir', required=True, type=str, help="output directory")
    parser.add_argument('--voc_perc', dest='voc_perc', required=True, type=str, help="comma-separated list of VOC frequencies (%) to be simulated")
    parser.add_argument('--chim_perc', dest='chim_perc', required=True, type=str, help="comma-separated list of chimeric read frequencies (%) to be simulated")
    parser.add_argument('--total_cov', dest='total_cov', default=10000, type=int, help="total sequencing depth to be simulated")
    parser.add_argument('--data_exploration_only', action='store_true', help="exit after sequence selection")
    parser.add_argument('--spike_only', action='store_true', help="simulate reads for spike region only")
    # parser.add_argument('--sub_error_rate', dest='sub_error_rate', default=1.0, type=float, help="substitution error rate for art_illumina")
    # parser.add_argument('--ins_error_rate', dest='ins_error_rate', default=1.0, type=float, help="insertion error rate for art_illumina")
    # parser.add_argument('--del_error_rate', dest='del_error_rate', default=1.0, type=float, help="deletion error rate for art_illumina")

    args = parser.parse_args()

    # create output directory
    try:
        os.makedirs(args.outdir)
    except FileExistsError:
        pass

    VOC_frequencies = args.voc_perc.split(',')
    chim_frequencies = args.chim_perc.split(',')
    total_cov = args.total_cov
    VOC_files = args.fasta_VOC.split(',')
    VOC_names = [filepath.split('/')[-1] for filepath in VOC_files]
    exclude_list = [name.split('_')[0] for name in VOC_names]

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
    # write all sequences to new fasta
    fasta_full = args.outdir + '/full_sequences.fasta'
    subprocess.check_call("cat {0}/sequences.fasta {1} > {2}"
        .format(args.outdir, " ".join(VOC_files), fasta_full), shell=True)
    print("Full fasta for all sequences is in {}/full_sequences.fasta".format(args.outdir))

    # Init R
    SimFFPE = importr('SimFFPE')
    Biostrings = importr('Biostrings')

    robjects.r('''
        ## Get (redundant) phred score profile
        PhredScoreProfilePath <- system.file("extdata", "PhredScoreProfile2.txt", package = "SimFFPE")

        PhredScoreProfile <- as.matrix(read.table(PhredScoreProfilePath, skip = 1))
        colnames(PhredScoreProfile) <- strsplit(readLines(PhredScoreProfilePath)[1], "\t")[[1]]
    ''')

    if args.data_exploration_only:
        sys.exit()

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
    for voc_freq in VOC_frequencies:
        VOC_cov = round(total_cov * float(voc_freq)/100, 2)
        background_cov = round((total_cov - VOC_cov) / len(selection_df.index), 2)
        voc_freq = str(round(float(voc_freq), 3))
        for chim_freq in chim_frequencies:
            err_freq = float(chim_freq)
            # chim_freq = str(round_sig(float(chim_freq) / 100))
            # simulate background sequence read
            print("Simulating background reads from {} at {}x coverage ".format(fasta_selection, background_cov))
            print("Chimeric read rate: {}%".format(chim_freq))
            if background_cov == 0:
                print("Not simulating, coverage is 0")
                subprocess.check_call("touch {2}/background_ab{0}_er{1}_1.fq".format(background_cov, chim_freq, os.path.abspath(args.outdir)), shell=True)
                subprocess.check_call("touch {2}/background_ab{0}_er{1}_2.fq".format(background_cov, chim_freq, os.path.abspath(args.outdir)), shell=True)
            else:
                robjects.r('''
                    sourceSeq = readDNAStringSet("{0}")
                    refp <- "{1}"
                    outFile <- "{4}/background_ab{2}_er{5}"
                    readSimFFPE(sourceSeq, refp, PhredScoreProfile, outFile, coverage={2}, chimericProp={3}, sdInsertLen=10, chimMutRate=0, noiseRate=0, highNoiseRate=0, pairedEnd=TRUE, overWrite=TRUE)
                '''.format(os.path.abspath(fasta_selection), os.path.abspath(fasta_selection), background_cov, str(err_freq/100), os.path.abspath(args.outdir), chim_freq))
            
            # simulate reads for VOC, merge and shuffle
            for filename in VOC_files:
                VOC_name = filename.rstrip('.fasta').split('/')[-1]
                if args.spike_only:
                    voc_fasta = args.outdir + "/{}.trimmed.fasta".format(VOC_name)
                else:
                    voc_fasta = filename
                print("Simulating reads from {} at {}x coverage".format(VOC_name, VOC_cov))
                print("Chimeric read rate: {}%".format(chim_freq))
                robjects.r('''
                    sourceSeq = readDNAStringSet("{0}")
                    refp <- "{1}"
                    outFile <- "{4}/{5}_ab{2}_er{6}"
                    readSimFFPE(sourceSeq, refp, PhredScoreProfile, outFile, coverage={2}, chimericProp={3}, sdInsertLen=10, chimMutRate=0, noiseRate=0, highNoiseRate=0, pairedEnd=TRUE, overWrite=TRUE)
                '''.format(os.path.abspath(voc_fasta), os.path.abspath(fasta_selection), VOC_cov, str(err_freq/100), os.path.abspath(args.outdir), VOC_name, chim_freq))


                print("\nMerging fastqs...")
                subprocess.check_call("cat {0}/background_ab{3}_er{2}_1.fq {0}/{1}_ab{4}_er{2}_1.fq > {0}/tmp1.fq"
                    .format(args.outdir, VOC_name, chim_freq, background_cov, VOC_cov), shell=True)
                subprocess.check_call("cat {0}/background_ab{3}_er{2}_2.fq {0}/{1}_ab{4}_er{2}_2.fq > {0}/tmp2.fq"
                    .format(args.outdir, VOC_name, chim_freq, background_cov, VOC_cov), shell=True)

                print("Shuffling reads...")
                subprocess.check_call("shuffle.sh in={0}/tmp1.fq in2={0}/tmp2.fq out={0}/wwsim_{1}_ab{3}_er{2}_1.fastq out2={0}/wwsim_{1}_ab{3}_er{2}_2.fastq overwrite=t fastawrap=0 ignorebadquality"
                    .format(args.outdir, VOC_name, chim_freq, voc_freq), shell=True)
                
            print("\nBenchmarks with a chimeric read frequency of {}% are ready!\n\n".format(chim_freq))
    # clean up temporary files
    os.remove("{}/tmp1.fq".format(args.outdir))
    os.remove("{}/tmp2.fq".format(args.outdir))
    return


def round_sig(x, sig=2):
    if x == 0: return x
    return round(x, sig-int(floor(log10(abs(x))))-1)


if __name__ == "__main__":
    sys.exit(main())
