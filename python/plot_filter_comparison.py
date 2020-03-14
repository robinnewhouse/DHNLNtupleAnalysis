import sys
if sys.version_info < (3, 0):
    sys.stdout.write("Sorry, requires Python 3.x, not Python 2.x\n")
    sys.exit(1)

import argparse, os, math, ROOT, glob, uproot, time
import numpy as np
import matplotlib.pyplot as plt
import random, string
from IPython.display import Image

from ROOT import gROOT
# from pylab import *

gROOT.SetStyle("ATLAS");

tmp_path = "/tmp/rnewhous/"


def make_cutflow_dict(hcutflow):
    cutflow_dict = {}
    for i in range(hcutflow.GetXaxis().GetNbins()):
        label = str(hcutflow.GetXaxis().GetLabels()[i])
        if "filter" in label: label = "filter mismatch"
        value = hcutflow[i + 1]
        cutflow_dict[label] = value
    return cutflow_dict


def make_cutflow(ch_name='4filter_FP', color=ROOT.kAzure, histogram_path="", labels=[]):
    # Open file and get cutflow
    file = ROOT.TFile(histogram_path)
    hcutflow = file.Get('CutFlow_' + ch_name)

    # Create canvas
    # hcutflow.GetXaxis().SetBinLabel(3, ch_name)
    canvas_name = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
    canvas = ROOT.TCanvas(canvas_name, "cutflow", 1000, 800)
    canvas.cd(1)

    # Prepare histogram settings for drawing
    ymax_cutflow = hcutflow.GetMaximum()
    hcutflow.GetYaxis().SetRangeUser(0, ymax_cutflow * 1.05)
    hcutflow.SetFillColor(color)
    hcutflow.SetLineWidth(0)
    hcutflow.GetXaxis().SetTickLength(0)
    hcutflow.SetMarkerSize(1.5)
    hcutflow.Draw("HIST TEXT0 SAME")

    # Add text information
    draw_text(*labels)

    # Display and save image
    # MyC02.Draw(canvas_name)
    image_file = tmp_path + 'Cutflow_' + ch_name + '.png'
    canvas.SaveAs(image_file)
    #     display(Image(image_file))

    return make_cutflow_dict(hcutflow)


def print_cutflow_to_file(cutflow_dict_tp, cutflow_dict_fp, name=""):
    # Print cutflows to text file
    with open(tmp_path + 'Cutflow.txt', 'w') as f:
        # Print title
        print("cut".rjust(20), "True Positive".rjust(15), "False Positive".rjust(15), "Scale Factor".rjust(15), sep=" -", file=f)
        print("                 ------------------------------------------------------", file=f)
        # Print table
        for key in cutflow_dict_tp.keys():
            val_tp = cutflow_dict_tp[key]
            val_fp = cutflow_dict_fp[key]
            scale_factor = val_tp / (val_fp + val_tp)
            if key in ["all", "trigger"]: scale_factor = 1  # don't show scale factor before filter cut
            print(
            key.rjust(20), str(val_tp).rjust(15), str(val_fp).rjust(15), "{0:.4f}".format(scale_factor).rjust(15), sep="  ", file=f)


def plot_stacked_cutflow(cutflow_dict_tp, cutflow_dict_fp, labels=[], name=""):
    # Prepare inputs
    fig, ax1 = plt.subplots(figsize=(13, 8))
    ind = np.arange(len(cutflow_dict_tp))
    vals_tp = np.array(list(cutflow_dict_tp.values()))
    vals_fp = np.array(list(cutflow_dict_fp.values()))
    vals_fp[0:2] = [0, 0]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][0:2]  # get colors for labels later

    # Plot stacked cutflow
    p_tp = ax1.bar(ind, vals_tp, width=1, label="True Positive", color=colors[0])  # plot bar graph
    p_fp = ax1.bar(ind[2:], vals_fp[2:], width=1, label="False Positive", bottom=vals_tp[2:],
                   color=colors[1])  # plot bar graph
    plt.ylim((0, vals_tp[0] * 1.2))
    plt.xticks(ind, cutflow_dict_tp.keys(), rotation=60, ha='right', fontsize=15)
    plt.ylabel('Number of passed events', fontsize=20)
    ax1.legend(loc=(.65, .4), fontsize=20)

    # Print text labels on top of cutflow
    cmap = plt.cm.get_cmap('Oranges')
    for v, i in enumerate(ind):
        plt.text(i, vals_tp[i] + vals_fp[i] + 2800, "{:.0f}".format(vals_tp[i]),
                 color=colors[0], fontweight='bold', fontsize=12, ha='center')
        if i < 2: continue
        plt.text(i, vals_tp[i] + vals_fp[i] + 800, "+{:.0f}".format(vals_fp[i]),
                 color=colors[1], fontweight='bold', fontsize=12, ha='center')

        # Plot scale factor
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    vals_sf = np.divide(vals_tp, (vals_tp + vals_fp))
    vals_sf[1] = vals_sf[2]  # trick for plotting left side properly
    p_sf = ax2.step(ind[1:] + .5, vals_sf[1:], label="Scale Factor", linewidth=3, color="red")  # plot bar graph
    plt.text(ind[-1], vals_sf[-1] - .04, "{:.4f}".format(vals_sf[-1]), color="red", fontweight='bold', fontsize=15,
             ha='center')
    plt.ylim((0, 1))
    plt.ylabel('Scale Factor TP/(TP+FP) ', fontsize=20, labelpad=15)
    ax2.legend(loc=(.65, .6), fontsize=20)

    # Print information on plot
    plt.text(.65, .75, "ATLAS", fontweight='bold', fontstyle='italic', fontsize=28, transform=ax1.transAxes)
    plt.text(.80, .75, "Internal", fontsize=28, transform=ax1.transAxes)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ypos = 0.75
    for label in labels:
        plt.text(.35, ypos, label, fontsize=20, transform=ax2.transAxes)
        ypos = ypos - 0.075

    plt.title(name, fontsize=20)
    image_file = tmp_path + 'Cutflow_' + 'stacked' + '.png'
    plt.savefig(image_file)
    plt.show()


def make_cutflows(hist_path, labels, name):
    cutflow_dict_tp = make_cutflow(ch_name='4filter_TP', color=ROOT.kAzure - 4, histogram_path=hist_path,
                                   labels=[*labels, "True Positives"])
    cutflow_dict_fp = make_cutflow(ch_name='4filter_FP', color=ROOT.kRed - 4, histogram_path=hist_path,
                                   labels=[*labels, "False Positives"])

    print_cutflow_to_file(cutflow_dict_tp, cutflow_dict_fp, name)

    plot_stacked_cutflow(cutflow_dict_tp, cutflow_dict_fp, labels, name)


make_cutflows(hist_path="/eos/home-r/rnewhous/public/HNL/filter_test/WmuHNL50_10G_lt10dd_el/histograms.root",
              labels=['Mass: 10 GeV', 'Lifetime: 10 mm', 'DV type: e-mu'], name="WmuHNL50_10G_lt10dd_el")