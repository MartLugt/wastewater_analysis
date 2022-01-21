#!/usr/bin/env python3

import sys
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--freqs', dest='freqs', type=str, required=True)
    parser.add_argument('--voc_freq', dest='voc_freq', default=10, type=int)
    parser.add_argument('--sars_fold', dest='sars_fold', default=1000)
    parser.add_argument('--contaminants', dest='contaminants', type=str, required=True)
    parser.add_argument('--conts_amount', dest='conts_amount', default=0)
    parser.add_argument('--font_size', dest='font_size', default=12)
    parser.add_argument('--outdir', dest='outdir', required=True)
    args = parser.parse_args()

    contaminant_names = args.contaminants.split(',')
    freqs = [float(i) for i in args.freqs.split(',')]

    conts_amount = float(args.conts_amount) if args.conts_amount != 0 else len(contaminant_names)

    voc_ct = args.sars_fold * args.voc_freq / 100
    background_ct = args.sars_fold - voc_ct

    voc = [voc_ct for _ in freqs]
    background = [background_ct for _ in freqs]

    contaminants = []
    for c in contaminant_names:
        contaminants.append(
            [(args.sars_fold * (1 / (f / 100)) - args.sars_fold) / conts_amount for f in freqs]
        )

    labels=['10⁰', '10¹', '10²']

    plt.rcParams.update({'font.size': args.font_size}) # increase font size
    plt.figure()
    plt.bar(labels, voc, 0.5, label='VoC')
    plt.bar(labels, background, 0.5, bottom=voc, label='Background')
    bottom = np.add(voc, background)

    for i, c in enumerate(contaminant_names):

        plt.bar(labels, contaminants[i], 0.5, bottom=bottom, label=c)
        bottom = np.add(bottom, contaminants[i])

    plt.legend(prop={'size': 12})
    plt.ylim(0, (voc_ct + background_ct + ((args.sars_fold * (1 / (min(freqs) / 100)) - args.sars_fold) / conts_amount)) * 1.05)
    plt.ylabel('Absolute read count')
    plt.xlabel('Total SARS-CoV-2 frequency (%)')

    plt.tight_layout()

    plt.savefig("{}/count_bar_plot.png".format(args.outdir))



if __name__ == "__main__":
    sys.exit(main())