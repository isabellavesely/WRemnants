import argparse
from utilities.io_tools import combinetf_input
import pandas as pd
import matplotlib.pyplot as plt
import os
import mplhep as hep
import re

"""
For plotting impact uncertainties over various luminosities, comparing high pile-up and combined fits.
Additionally, plots the percent decrease in uncertainty from simple high pile-up.
"""

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--ungroup", action='store_true', help="Use ungrouped nuisances")
    parser.add_argument("-s", "--sort", action='store_true', help="Sort nuisances by impact")
    # parser.add_argument("-ih", "--inputHighPU", default = 'fitresults_123456789_highPU.root', type=str, help="HighPU fitresults file, output ROOT file from combinetf")
    parser.add_argument("-ic", "--inputfolder", type=str, help="Folder with combined fitresults, output ROOT file from combinetf")
    parser.add_argument("-o", "--outfile", default = 'combine_lumi_projection.png', type=str, help="name of output file")
    parser.add_argument("--addPDF", default=0, type=int, nargs="*", help="Additionally outputs in PDF if == 1")
    parser.add_argument("-lsc", "--lumiscales", default = [1.0], type = float, nargs='+', help="Additional luminosity scales to examine (default 1.0)")
    return parser.parse_args()

def get_file_names(folder_path):
    return os.listdir(folder_path)

def updateImpactsDict(args, fitresult, df, lumiscale, poi='Wmass'):
    """
    Update dictionary for the labels and values, for one lumiscale
    Returns nuisance labels if plotAll set (for columns)
    """
    # h_impacts, h_labels,_ = combinetf_input.read_impacts_poi(args.inputHighPU, not args.ungroup, sort=args.sort, poi=poi, normalize = False)
    impacts,labels,_ = combinetf_input.read_impacts_poi(fitresult, not args.ungroup, sort=args.sort, poi=poi, normalize = False)
    lumiVal = lumiscale * 0.200181002 / 100

    plot_labels = ['HighPU Total', 'Combine Fit Total', 'Luminosity']
    plot_vals = [0.09, impacts[list(labels).index('Total')], lumiVal]

    new_row = {k: v*100 for k, v in zip(plot_labels, plot_vals)}
    df.loc[len(df)] = new_row

def updatePercentDict(args, fitresult, df_percents, lumiscale, poi='Wmass'):
    impacts,labels,_ = combinetf_input.read_impacts_poi(fitresult, not args.ungroup, sort=args.sort, poi=poi, normalize = False)
    lumiVal = lumiscale * 0.200181002 / 100
    combine_total = impacts[list(labels).index('Total')]
    hPU_total = 0.09

    plot_labels = ['Percent Decreases in Uncertainty', 'Luminosity']
    plot_vals = [(hPU_total - combine_total)/(hPU_total), lumiVal]
    # print((combine_total - hPU_total)/(hPU_total))

    new_row = {k: v*100 for k, v in zip(plot_labels, plot_vals)}
    df_percents.loc[len(df_percents)] = new_row


def sortFileNames(inputFiles): # natural sorting for fit result files (if lumi inputted in any order)
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(inputFiles, key=alphanum_key)

if __name__ == '__main__':
    args = parseArgs()
    lumiscales = sorted(args.lumiscales)
    inFolder = args.inputfolder
    inputFiles = sortFileNames(get_file_names(inFolder))

    df = pd.DataFrame(columns=['HighPU Total', 'Combine Fit Total', 'Luminosity'])
    df_percents = pd.DataFrame(columns=['Percent Decreases in Uncertainty', 'Luminosity'])

    for i in range(len(lumiscales)):
        inputFile = inFolder + '/' + inputFiles[i]
        fitresult = combinetf_input.get_fitresult(inputFile)
        for poi in combinetf_input.get_poi_names(fitresult):
            updateImpactsDict(args, fitresult, df, lumiscales[i], poi)
            updatePercentDict(args, fitresult, df_percents, lumiscales[i], poi)

    # print(df_percents)

    plt.figure(figsize=(8, 8))
    hep.cms.label(fontsize=20, data=False, label="Projection", com=13.6)
    colors = ['royalblue', 'crimson']
    for color, column in zip(colors, df.columns):
        if column != 'Luminosity':  # Exclude the Luminosity column from plotting
            plt.plot(df['Luminosity'], df[column], label=column, marker='o', color = color)
    plt.xlabel("Integrated luminosity (fb$^{-1})$", fontsize=14)
    plt.ylabel("Uncertainty in $m_{W}$ (MeV)", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0, 2.25)
    plt.ylim(0, 15)
    plt.grid(True)
    plt.legend(fontsize=14)
    plt.savefig(args.outfile)

    plt.figure(figsize=(8, 8))
    hep.cms.label(fontsize=20, data=False, label = "Projection", com=13.6)
    plt.plot(df_percents['Luminosity'], df_percents['Percent Decreases in Uncertainty'], label='Percent Decreases in Uncertainty', marker='o', color = 'crimson')
    plt.xlabel("Integrated luminosity (fb$^{-1})$", fontsize=14)
    plt.ylabel("Percent Decrease in Uncertainty (%)", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0, 2.25)
    plt.ylim(0, 20)
    plt.grid(True)
    plt.legend(fontsize=14)
    plt.savefig('combine_lumi_projection_percent.png')

    if args.addPDF:
        pdfFileName = (args.outfile).split('.')[0] + ".pdf"
        plt.savefig(pdfFileName, format = 'pdf')
