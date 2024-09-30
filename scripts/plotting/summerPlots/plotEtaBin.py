import argparse
from utilities.io_tools import combinetf_input #intput_tools, output_tools?
import pandas as pd
import matplotlib.pyplot as plt
import os
import mplhep as hep
import re

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--ungroup", action='store_true', help="Use ungrouped nuisances")
    parser.add_argument("-n", "--nuisance", type=str, help="Only print value for specific nuisance")
    parser.add_argument("-s", "--sort", action='store_true', help="Sort nuisances by impact")
    parser.add_argument("-i", "--inputfolder", type=str, help="Folder with fitresults, output ROOT file from combinetf")
    parser.add_argument("-o", "--outfile", default = 'uncertainty_from_eta.png', type=str, help="name of output file")
    parser.add_argument("-o2", "--addPercentPlot", default = 'percent_eta_bin_projection.png', type=str, help="name of output file, if uncertainty plot is requested")
    parser.add_argument("--addPDF", default=0, type=int, nargs="*", help="Also output in PDF if == 1")
    parser.add_argument("-eta", "--etaRebins", default = [1.0], type = float, nargs='+', help="eta rebinning values")
    return parser.parse_args()

def get_file_names(folder_path):
    return os.listdir(folder_path)

def updateImpactsDict(args, fitresult, df, etaRebin, poi='Wmass'):
    """
    Update dictionary for the labels and values, for one lumiscale
    """
    impacts,labels,_ = combinetf_input.read_impacts_poi(fitresult, not args.ungroup, sort=args.sort, poi=poi, normalize = False)
    
    # plot_labels = ['Total', 'Theory', 'Experiment', 'PDF', 'Data Stat.', 'Fakes', 'Number of Eta Bins']
    plot_labels = ['Total', 'PDF', 'Data Stat.', 'Number of Eta Bins']
    largestRebin = args.etaRebins[-1] / 100
    plot_vals = [impacts[list(labels).index('Total')],\
                # impacts[list(labels).index('theory')], impacts[list(labels).index('experiment')], \
                impacts[list(labels).index('pdfCT18Z')], \
                impacts[list(labels).index('stat')], \
                # impacts[list(labels).index('Fake')], 
                largestRebin / etaRebin]
                # otherwise, pdfNNPDF31, pdfNNPDF40, pdfCT18Z
                
    new_row = {k: v*100 for k, v in zip(plot_labels, plot_vals)}
    df.loc[len(df)] = new_row
    
def natural_key(filename):
    match = re.search(r'(\d+)\.root$', filename)
    if match:
        return int(match.group(1))
    return filename

def updatePercentDict(args, fitresult, df_percents, original_bin_unc, etaRebin, poi='Wmass'):
    impacts,labels,_ = combinetf_input.read_impacts_poi(fitresult, not args.ungroup, sort=args.sort, poi=poi, normalize = False)
    new_unc = impacts[list(labels).index('Total')]
    print(f"new_unc: {new_unc}")
    print(f"original: {original_bin_unc}")
    largestRebin = args.etaRebins[-1]
    
    plot_vals = [100*(original_bin_unc-new_unc)/original_bin_unc, largestRebin / etaRebin]
    df_percents.loc[len(df_percents)] = plot_vals

if __name__ == '__main__':
    args = parseArgs()
    etaRebins = sorted(args.etaRebins)
    inFolder = args.inputfolder
    inputFiles = sorted((get_file_names(inFolder)), key = natural_key) 
    
    df = pd.DataFrame(columns=['Total', 'PDF', 'Data Stat.', 'Number of Eta Bins'])
    df_percents = pd.DataFrame(columns=['Percent Decreases in Uncertainty', 'Number of Eta Bins'])
    
    for i in range(len(etaRebins)):
        inputFile = inFolder + '/' + inputFiles[i]
        fitresult = combinetf_input.get_fitresult(inputFile)
        for poi in combinetf_input.get_poi_names(fitresult):
            print(poi)
            updateImpactsDict(args, fitresult, df, etaRebins[i], poi)
            print(df) 
    orginal_bin_unc = df['Total'].iloc[len(etaRebins) - 1] / 100
    
    if args.addPercentPlot:
        for i in range(len(etaRebins)):
            inputFile = inFolder + '/' + inputFiles[i]
            fitresult = combinetf_input.get_fitresult(inputFile)
            for poi in combinetf_input.get_poi_names(fitresult):
                updatePercentDict(args, fitresult, df_percents, orginal_bin_unc, etaRebins[i], poi)
                
    plt.figure(figsize=(8, 8))
    hep.cms.label(fontsize=20, data=False, label="Projection", com=13.6)
    # colors = ['steelblue', 'darkorange', 'olivedrab', 'mediumpurple', 'crimson', 'hotpink']
    colors = ['steelblue', 'mediumpurple', 'crimson', 'blue']
    for i, column in enumerate(df.columns):
        print(column)
        if column != 'Number of Eta Bins':  # Exclude the Luminosity column from plotting
            plt.plot(df['Number of Eta Bins'], df[column], label=column, marker='o', color = colors[i])
    plt.xlim(0, 50)
    plt.ylim(0, 55)
    plt.title("Scaled pt-eta-charge, onlyPDF lowPU CT18Z eta bin study. \n Lumi: 1.0 fb^-1", y = 1.05, fontsize = 16)
    plt.xlabel("Number of Eta Bins", fontsize = 16)
    plt.ylabel("Uncertainty in $m_{W}$ (MeV)", fontsize = 16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.legend(fontsize=16)
    plt.savefig(args.outfile)
    if args.addPDF:
        pdfFileName = (args.outfile).split('.')[0] + ".pdf"
        plt.savefig(pdfFileName, format = 'pdf')
    
    if args.addPercentPlot:
        plt.figure(figsize=(8, 8))
        hep.cms.label(fontsize=20, data=False, label = "Projection", com=13.6)
        plt.title("Scaled pt-eta-charge, onlyPDF lowPU eta bin study, Percent change.\n Lumi: 1.0 fb^-1", y = 1.05, fontsize = 16)
        plt.plot(df_percents['Number of Eta Bins'], df_percents['Percent Decreases in Uncertainty'], label='Percent Decreases in Uncertainty', marker='o', color = 'crimson')
        plt.xlabel("Number of Eta Bins", fontsize=14)
        plt.ylabel("Percent Decrease in Uncertainty (%)", fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlim(0, 50)
        plt.ylim(-5, 30)
        plt.grid(True)
        plt.legend(fontsize=16)
        plt.savefig(args.addPercentPlot)
