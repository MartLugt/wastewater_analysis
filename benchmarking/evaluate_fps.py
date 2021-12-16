#!/usr/bin/env python3

import argparse
import os
import sys

import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt

class Error:
    def __init__(self, voc_freq, err_freq, sum):
        self.voc_freq = voc_freq
        self.err_freq = err_freq
        self.sum = sum

def main():
    parser = argparse.ArgumentParser(description="Plot false positive percentages between datasets.")
    parser.add_argument('predictions', type=str, nargs='+', help="prediction files")
    parser.add_argument('--voc', dest='voc', type=str, required=True, help="comma-separated list of VoCs")
    parser.add_argument('-o,--outdir', dest='outdir', required=True)
    parser.add_argument('--plot_abundance_value', dest='plot_ab_val', default=1, type=float, help="abundance value for plots in which the abundance is not plotted. (Value is set to closest datapoint). Default: 1")
    args = parser.parse_args()

    err_dict = dict()
    dataset_set = set()
    ab_set = set()
    err_set = set()
    voc_list = args.voc.split(',')

    os.makedirs(args.outdir, exist_ok=True)

    for filename in args.predictions:
        dataset  = filename.split('/')[-4]
        dir_name = filename.split('/')[-2]
        voc_name = dir_name.split('_')[0]
        err_freq = float(dir_name.split('_')[-1].lstrip('er'))
        voc_freq = float(dir_name.split('_')[-2].lstrip('ab'))
        if voc_name not in voc_list:
            continue

        ab_set.add(voc_freq)
        err_set.add(err_freq)
        dataset_set.add(dataset)
        print(dataset)

        with open(filename, 'r') as f:
            # number of false pos and sum of false pos
            num_voc = 0
            sum_ab = 0
            for line in f:
                if line[0] == '#':
                    continue
                [variant, tpm, ab, corrected_ab] = line.rstrip('\n').split('\t')
                # TODO: RELATIVE ab to actual voc ab?
                if variant == voc_name:
                    continue
                
                ab = float(ab)

                if ab == 0:
                    continue

                num_voc += 1
                sum_ab += float(ab)

            if dataset not in err_dict:
                err_dict[dataset] = []
            err_dict[dataset].append(Error(voc_freq, err_freq, sum_ab))

    colors = {dataset : cm.tab10((i))
         for i, dataset in enumerate(dataset_set)}

    plot_single_frq_val = min(ab_set, key=lambda x:abs(x-args.plot_ab_val))

    plt.figure()
    for dataset in dataset_set:

        sum_of_sums = dict()
        for e in [x for x in err_dict[dataset] if x.voc_freq == plot_single_frq_val]:
            if e.err_freq not in sum_of_sums:
                sum_of_sums[e.err_freq] = 0
            sum_of_sums[e.err_freq] += e.sum

        plot_tups = [(k,v) for k, v in sum_of_sums.items()]

        plot_tups.sort(key= lambda x: x[0])

        errs, sums = zip(*plot_tups)

        # # Sort lists on error frequency
        # err_dict[dataset].sort(key = lambda x : x.err_freq)
        # err_freqs = [x.err_freq for x in err_dict[dataset] if x.voc_freq == plot_single_frq_val]
        # sums = [x.sum for x in err_dict[dataset] if x.voc_freq == plot_single_frq_val]

        print(len(errs), len(sums))
        plt.plot(errs, sums, label=dataset, color=colors[dataset])
    
    plt.legend()
    plt.grid(which="both", alpha=0.2)
    plt.xlabel("Induced error frequency (%)")
    plt.ylabel("Sum of false positive percentages")
    
    plt.xlim(0.5, 10)
    plt.xscale('log')
    plt.tight_layout()

    plt.show()     
    

if __name__ == "__main__":
    sys.exit(main())
