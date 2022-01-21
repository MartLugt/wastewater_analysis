#!/usr/bin/env python3

import sys
import os
import argparse
import matplotlib
import matplotlib.pyplot as plt
import json
import matplotlib.ticker
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from math import log10

superscript = str.maketrans("-0123456789.", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹·")

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
    parser.add_argument('--correct', action='store_true', help="Correct for pseudoaligned reads")
    parser.add_argument('--full_cont_count', dest='cont_count', default=1, type=float, help="Contamination count in full dataset if evaluating a sub-set")

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
        run_info = filename.split('/')
        run_info[-1] = "run_info.json"
        run_info = "/".join(run_info)

        dir_name = filename.split('/')[-2]
        voc_name = dir_name.split('_')[0]
        sars_freq = 100 - float(dir_name.split('_')[-1].lstrip('ab'))
        if voc_name not in voc_list:
            continue
        if sars_freq == 0:
            continue
        # elif voc_freq < args.min_ab:
        #     continue
        variant_set.add(voc_name)
        with open(filename, 'r') as f, open(run_info) as r:
            run_data = json.load(r)
            n = run_data["n_processed"]
            n_p = run_data["n_pseudoaligned"]
            n_sars2 = n * (sars_freq / 100)

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
                    ab = ab * (100 / sars_freq)
                if args.correct:
                    ab = ab * (n_p / n_sars2)   
                abs_err = abs(ab - 10)
                if ab < args.min_ab:
                    continue

                # voc_freq = (-1.0 * sars_freq / args.cont_count)
                voc_freq = sars_freq

                positives.append(variant)
                if variant == voc_name:
                    variant_found = True
                    err_tups.append((voc_name, voc_freq, abs_err, ab))
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
                    err_list.append((voc_name, voc_freq, abs_err, ab))
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

    _, f, e, _ = zip(*err_list)
    unique_err_vals = list(set(e))
    unique_freq_vals = list(set(f))

    for voc in voc_list:
        f = list(filter(lambda x: x[0] == voc , err_list))
        _, _, err, _ = zip(*f)
        pct_err = map(lambda x: x/10*100, err)
        av = sum(pct_err) / len(f)
        print("Average error for {}: {}%".format(voc, av))

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
    err_list.sort(key = lambda x : x[1])
    variant_list = sorted(list(variant_set))

    # fix color per voc
    # colormap = cm.get_cmap('Accent', len(variant_list))
    # colors = {voc : colormap((i)/len(variant_list))
    #             for i, voc in enumerate(variant_list)}

    colors = {voc : cm.tab10((i))
            for i, voc in enumerate(variant_list)}

    _, f, e, _ = zip(*err_list)
    unique_err_vals = list(set(e))
    unique_freq_vals = list(set(f))

    plt.rcParams.update({'font.size': args.font_size}) # increase font size
    plt.figure()
    for voc in variant_list:
        freq_values = [x[1] for x in err_list if x[0] == voc]# and x[2] < 10]
        err_values = [x[2]/10*100 for x in err_list if x[0] == voc]# and x[2] < 10]
        plt.plot(freq_values, err_values, label=voc, color=colors[voc])
        # if (freq_values[0] > min(unique_freq_vals)):
        #     plt.plot(freq_values[0], err_values[0], marker="s", color=colors[voc], markersize=6)
    ax = plt.gca()
    plt.legend()
    plt.grid(which="both", alpha=0.2)
    plt.ylim(-5, 105)
    plt.xlabel("Relative contamination (%)")
    plt.ylabel("Relative prediction error (%)")
    # plt.gcf().set_size_inches(4, 3)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_error_plot{}.{}".format(args.outdir,
                                                     args.suffix,
                                                     format))

    # also plot on log scale
    f = lambda a: 10**a
    g = lambda b: np.log10(b)
    # plt.gca().set_xscale('function', functions=(g, f))
    # plt.xscale('symlog', linthresh=0.005)
    # plt.grid(True)
    # plt.xlim(1.2 * min(unique_freq_vals), 0.8 * max(unique_freq_vals))

    plt.xscale('log')
    # plt.gca().invert_xaxis()

    # print(plt.xticks()[0])
    # plt.gca().set_xticklabels(plt.xticks()[0][::-1])

    # ax.xaxis.set_minor_locator(matplotlib.ticker.SymmetricalLogLocator(subs=np.arange(1, 10), linthresh=0.005, base=10.0)) 
    # ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(symlog_tick_formatter))
    # ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(1,10) / -10))
    # plt.gca().xaxis.grid(True, which='minor')  # minor grid on too
    # plt.gca().invert_xaxis()
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_error_plot_logscale{}.{}".format(args.outdir,
                                                              args.suffix,
                                                              format))

    # plt.show()    

    # plot true vs estimated frequencies on a scatterplot
    plt.figure()
    for voc in variant_list:
        freq_values = [x[1] for x in err_list if x[0] == voc]
        est_values = [x[3] for x in err_list if x[0] == voc]
        plt.scatter(freq_values, est_values, label=voc, alpha=0.7,
                    color=colors[voc], s=20)
    
    # plt.xscale('symlog', linthresh=0.005)
    # plt.gca().xaxis.set_minor_locator(matplotlib.ticker.SymmetricalLogLocator(subs=np.arange(1, 10), linthresh=0.005, base=10.0)) 
    # plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(symlog_tick_formatter))

    plt.grid(True)

    # plt.xscale('function', functions=(f, g))
    plt.xscale('log')
    plt.yscale('log')
    # plt.xlim(1.5 * min(unique_freq_vals), 0.7 * max(unique_freq_vals))
    plt.xlim(0.7, 150)
    plt.ylim(0.7, 150)
    # plt.plot([0, 100], [0, 100], 'k-', lw=0.75)
    plt.hlines(10, -150, 150, 'k', lw=0.75)
    plt.legend(prop={'size': args.font_size}) #ncol=len(variants_list),
    plt.grid(which="both", alpha=0.2)
    plt.xlabel("Relative contamination (%)")
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

    # plt.show()
    return

def save_dev(a: int, b: int):
    if a == 0:
        return 0
    else:
        return a / b

def symlog_tick_formatter(val, pos=None):
    return "10" + str(int(log10(1 * val))).translate(superscript)

if __name__ == "__main__":
    sys.exit(main())
