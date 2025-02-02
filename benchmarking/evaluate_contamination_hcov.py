#!/usr/bin/env python3

import sys
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import math

def main():
    parser = argparse.ArgumentParser(description="Evaluate predicted frequencies.")
    parser.add_argument('predictions', type=str, nargs='+', help="prediction files")
    parser.add_argument('--voc', dest='voc', type=str, required=True, help="comma-separated list of strains of interest")
    parser.add_argument('-o,--outdir', dest='outdir', required=True)
    parser.add_argument('--suffix', dest='suffix', default="", help="add suffix to output figure names")
    parser.add_argument('-v,--verbose', dest='verbose', action='store_true')
    parser.add_argument('-m', dest='min_ab', default=0, type=float, help="minimal abundance (any samples with true abundance below this threshold are skipped; any predictions below this threshold are considered absent)")
    parser.add_argument('--no_plots', action='store_true')
    parser.add_argument('--output_format', dest='output_format', default='png', help="comma-separated list of desired output formats")
    parser.add_argument('--font_size', dest='font_size', default=12, type=int, help="set font size for the plots")
    parser.add_argument('--conts_in_meta', action='store_true', help="Enable if contaminants are present in metadata and kallisto index")
    args = parser.parse_args()

    false_pos_count = 0
    false_neg_count = 0
    true_pos_count = 0
    true_neg_count = 0
    err_list = []
    variant_set = set()
    voc_list = args.voc.split(',')
    output_formats = args.output_format.split(',')

    # read predictions
    for filename in args.predictions:
        dir_name = filename.split('/')[-2]
        dataset = filename.split('/')[-4]
        voc_name = dir_name.split('_')[0]
        hcov = dataset.split('_')[-2]
        voc_freq = float(dir_name.split('_')[-1].lstrip('ab'))
        if voc_name not in voc_list:
            continue
        elif voc_freq < args.min_ab:
            continue
        variant_set.add(voc_name)
        with open(filename, 'r') as f:
            variant_found = False
            err_tups = []
            positives = []
            for line in f:
                if line[0] == "#":
                    continue
                [variant, tpm, ab, corrected_ab] = line.rstrip('\n').split('\t')
                if variant not in voc_list:
                    continue
                ab = float(ab)
                if (args.conts_in_meta):
                    ab = ab * (100 / voc_freq)
                abs_err = abs(ab - 10)
                if ab < args.min_ab:
                    continue
                positives.append(variant)
                if variant == voc_name:
                    variant_found = True
                    err_tups.append((voc_name, voc_freq, abs_err, ab, hcov))
                else:
                        false_pos_count += 1
                        if args.verbose:
                            print("False positive: {} predicted at {}% in {}".format(
                                    variant, ab, filename))
            if variant_found:
                true_pos_count += 1
                if len(err_tups) == 1:
                    err_list.append(err_tups[0])
                else:
                    voc_name = err_tups[0][0]
                    voc_freq = err_tups[0][1]
                    ab = sum([x[3] for x in err_tups])
                    abs_err = abs(ab - voc_freq)
                    err_list.append((voc_name, voc_freq, abs_err, ab, hcov))
            else:
                false_neg_count += 1
                if args.verbose:
                    print("VOC not found in {}".format(filename))
                # add zero estimate to error list?
                # err_list.append((voc_name, voc_freq, voc_freq, 0))
            for variant in voc_list:
                if variant not in positives and variant != voc_name:
                    # true negative
                    true_neg_count += 1
            true_neg_count += len([x for x in voc_list if
                                    x not in positives and x != voc_name ])

    _, f, e, _, h = zip(*err_list)
    unique_err_vals = list(set(e))
    unique_freq_vals = list(set(f))
    unique_hvoc_vals = list(set(h))

    unique_hvoc_vals.remove("Other")
    unique_hvoc_vals.append("Other")

    # compute averages
    av_err_list = []
    for hcov in unique_hvoc_vals:
        for freq in unique_freq_vals:
            f = list(filter(lambda x: x[4] == hcov and x[1] == freq, err_list))
            _, _, err, ab, _ = zip(*f)
            av_err = sum(err) / len(f)
            av_ab  = sum(ab)  / len(f)
            av_err_list.append((hcov, freq, av_err, av_ab))

    for hcov in unique_hvoc_vals:
        f = list(filter(lambda x: x[0] == hcov, av_err_list))
        _, _, err, _ = zip(*f)
        pct_err = map(lambda x: x/10*100, err)
        av = sum(pct_err) / len(f)
        print("Average error for {}: {}%".format(hcov, av))

    # compute stats
    average_rel_err = sum([x[2]/x[1]*100 for x in err_list]) / len(err_list)
    average_rel_err_tp = (sum([x[2]/x[1]*100 for x in err_list if x[3] > 0])
                            / len(err_list))
    # print("average relative error: {}%".format(average_rel_err))
    print("average relative error of true positives: {}%".format(
                                                            average_rel_err_tp))
    print("total # true positives: {}".format(true_pos_count))
    print("total # true negatives: {}".format(true_neg_count))
    print("total # false positives: {}".format(false_pos_count))
    print("total # false negatives: {}".format(false_neg_count))

    fpr = save_dev(false_pos_count, (false_pos_count + true_neg_count))
    fnr = save_dev(false_neg_count, (false_neg_count + true_pos_count))
    recall = save_dev(true_pos_count, (true_pos_count + false_neg_count))
    precision = save_dev(true_pos_count, (true_pos_count + false_pos_count))
    print("FPR = {}".format(fpr))
    print("FNR = {}".format(fnr))
    print("Precision = {}".format(precision))
    print("Recall = {}\n".format(recall)) # sensitivity

    if args.no_plots:
        sys.exit()

    # sort error tuples by voc frequency
    av_err_list.sort(key = lambda x : x[1])
    # err_list.sort(key = lambda x : x[1])
    # variant_list = sorted(list(variant_set))

    # fix color per voc
    # colormap = cm.get_cmap('Accent', len(variant_list))
    # colors = {voc : colormap((i)/len(variant_list))
    #             for i, voc in enumerate(variant_list)}

    colors = {hcov : cm.tab10((i))
            for i, hcov in enumerate(unique_hvoc_vals)}

    plt.rcParams.update({'font.size': args.font_size}) # increase font size
    plt.figure()
    for hcov in unique_hvoc_vals:
        freq_values = [x[1] for x in av_err_list if x[0] == hcov] # and x[2] < 10]
        err_values = [x[2]/10*100 for x in av_err_list if x[0] == hcov] # and x[2] < 10]
        plt.plot(freq_values, err_values, label=hcov, color=colors[hcov])
        if (freq_values[0] > min(unique_freq_vals)):
            plt.plot(freq_values[0], err_values[0], marker="s", color=colors[hcov], markersize=6)
    plt.legend()
    plt.grid(which="both", alpha=0.2)
    plt.ylim(-5, 105)
    plt.xlabel("Total SARS-CoV-2 frequency (%)")
    plt.ylabel("Relative prediction error (%)")
    # plt.gcf().set_size_inches(4, 3)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_error_plot{}.{}".format(args.outdir,
                                                     args.suffix,
                                                     format))

    # also plot on log scale
    plt.xscale('log')
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_error_plot_logscale{}.{}".format(args.outdir,
                                                              args.suffix,
                                                              format))

    # plot true vs estimated frequencies on a scatterplot
    plt.figure()
    for hcov in unique_hvoc_vals:
        freq_values = [x[1] for x in av_err_list if x[0] == hcov]
        est_values = [x[3] for x in av_err_list if x[0] == hcov]
        plt.scatter(freq_values, est_values, label=hcov, alpha=0.7,
                    color=colors[hcov], s=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.7, 150)
    plt.ylim(0.7, 150)
    # plt.plot([0, 100], [0, 100], 'k-', lw=0.75)
    plt.hlines(10, 0, 150, 'k', lw=0.75)
    plt.legend(prop={'size': args.font_size}) #ncol=len(variants_list),
    plt.grid(which="both", alpha=0.2)
    plt.xlabel("Total SARS-CoV-2 frequency (%)")
    plt.ylabel("Estimated VoC frequency (%)")
    # # Hide the right and top spines
    # ax = plt.gca()
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_scatter_loglog{}.{}".format(args.outdir,
                                                         args.suffix,
                                                         format))

    return

def save_dev(a: int, b: int):
    if a == 0:
        return 0
    else:
        return a / b

if __name__ == "__main__":
    sys.exit(main())
