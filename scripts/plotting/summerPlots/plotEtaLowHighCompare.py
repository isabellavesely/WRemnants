import argparse
from utilities.io_tools import combinetf_input #intput_tools, output_tools?
import pandas as pd
import matplotlib.pyplot as plt
import os
import mplhep as hep
import re
import numpy as np

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--ungroup", action='store_true', help="Use ungrouped nuisances")
    parser.add_argument("-n", "--nuisance", type=str, help="Only print value for specific nuisance")
    parser.add_argument("-s", "--sort", action='store_true', help="Sort nuisances by impact")
    parser.add_argument("-i1", "--inputfolder1", type=str, help="Folder with fitresults, output ROOT file from combinetf")
    parser.add_argument("-i2", "--inputfolder2", type=str, help="Folder with second fitresults, output ROOT file from combinetf")
    
    parser.add_argument("-o", "--outfile", default = 'compare_eta_uncertainty.png', type=str, help="name of output file")
    parser.add_argument("-o2", "--outfile_PercentPlot", default = 'compare_percent_eta_uncertainty.png', type=str, help="name of output file")
    parser.add_argument("-o3", "--outfile_CumulativePercentPlot", default = 'compare_cumulative_percent_eta_uncertainty.png', type=str, help="name of output file")
    parser.add_argument("-eta1", "--etaRebins1", default = [1.0], type = float, nargs='+', help="eta rebinning values")
    parser.add_argument("-eta2", "--etaRebins2", default = [1.0], type = float, nargs='+', help="eta rebinning values")

    parser.add_argument("--addPercentPlot", default = False)
    parser.add_argument("--cumulativePercents", default = False)
    parser.add_argument("--addPDF", default=0, type=int, nargs="*", help="Also output in PDF if == 1")
    return parser.parse_args()


# eta_bin_fits6_CT18Z_lumiPDF
# eta_bin_fits6_CT18Z_lowPU_pdf

def get_file_names(folder_path):
    return os.listdir(folder_path)

def updateImpactsDict(args, fitresult, df, etaRebin, poi='Wmass', lowPU = True):
    """
    Update dictionary for the labels and values, for one lumiscale
    """
    impacts,labels,_ = combinetf_input.read_impacts_poi(fitresult, not args.ungroup, sort=args.sort, poi=poi, normalize = False)
    if lowPU: 
        rebinToBin = 0.24
        plot_labels = ['LowPU (mt-eta-charge) Total', 'LowPU PDF', 'Number of Eta Bins']
    else: 
        rebinToBin = 0.24 ##FIXXXXXX!!
        # plot_labels = ['HighPU (pt-eta-charge) Total', 'HighPU PDF', 'Number of Eta Bins']
        plot_labels = ['HighPU (pt-eta-charge) Total', 'highPU PDF', 'Number of Eta Bins']
    print(f"Possible labels: {labels} and length {len(impacts)}")
    plot_vals = [impacts[list(labels).index('Total')],\
                # impacts[list(labels).index('stat')], \
                # impacts[list(labels).index('binByBinStat')], \
                impacts[list(labels).index('pdfCT18Z')], \
                rebinToBin/etaRebin] #pdfCT18Z, pdfNNPDF31
                
    new_row = {k: v*100 for k, v in zip(plot_labels, plot_vals)}
    df.loc[len(df)] = new_row
    
def natural_key(filename):
    match = re.search(r'(\d+)\.root$', filename)
    if match:
        return int(match.group(1))
    return filename

def updatePercentDict(args, fitresult, df_percents, original_bin_unc, etaRebin, poi='Wmass', lowPU = True):
    impacts,labels,_ = combinetf_input.read_impacts_poi(fitresult, not args.ungroup, sort=args.sort, poi=poi, normalize = False)
    new_unc = impacts[list(labels).index('Total')]
    print("updating Percent Dict now")
    
    if lowPU:
        largestRebin = args.etaRebins1[-1]
    else:
        largestRebin = args.etaRebins2[-1]
    
    plot_vals = [100*(original_bin_unc-new_unc)/original_bin_unc, largestRebin/etaRebin]
    df_percents.loc[len(df_percents)] = plot_vals

if __name__ == '__main__':
    args = parseArgs()
    etaRebins1 = sorted(args.etaRebins1)
    inFolder1 = args.inputfolder1
    inputFiles1 = sorted((get_file_names(inFolder1)), key = natural_key)
    
    etaRebins2 = sorted(args.etaRebins2)
    inFolder2 = args.inputfolder2 
    inputFiles2 = sorted((get_file_names(inFolder2)), key = natural_key) 
    
    df1 = pd.DataFrame(columns=['LowPU (mt-eta-charge) Total', 'LowPU PDF','Number of Eta Bins'])
    df2 = pd.DataFrame(columns=['HighPU (pt-eta-charge) Total', 'highPU PDF', 'Number of Eta Bins'])
    
    dfpercents1 = pd.DataFrame(columns=['LowPU (mt-eta-charge) Percent Decreases in Unc', 'Number of Eta Bins'])
    dfpercents2 = pd.DataFrame(columns=['HighPU (pt-eta-charge) Percent Decreases in Unc', 'Number of Eta Bins2'])
    
    # for first file
    for i in range(len(etaRebins1)):
        inputFile = inFolder1 + '/' + inputFiles1[i]
        fitresult = combinetf_input.get_fitresult(inputFile)
        for poi in combinetf_input.get_poi_names(fitresult):
            updateImpactsDict(args, fitresult, df1, etaRebins1[i], poi)
    df1 = df1.iloc[::-1].reset_index(drop=True)
    print(df1)
    orginal_bin_unc_1 = df1['LowPU (mt-eta-charge) Total'].iloc[0] / 100
    
    for i in range(len(etaRebins2)):
        inputFile = inFolder2 + '/' + inputFiles2[i]
        fitresult = combinetf_input.get_fitresult(inputFile)
        for poi in combinetf_input.get_poi_names(fitresult):
            updateImpactsDict(args, fitresult, df2, etaRebins2[i], poi, lowPU = False)
    df2 = df2.iloc[::-1].reset_index(drop=True)
    print(df2)
    orginal_bin_unc_2 = df2['HighPU (pt-eta-charge) Total'].iloc[0] / 100
    
    df_combined = pd.merge(df1, df2, on='Number of Eta Bins', how='outer')
    print(df_combined)
    
    for i in range(len(etaRebins1)):
        inputFile = inFolder1 + '/' + inputFiles1[i]
        fitresult = combinetf_input.get_fitresult(inputFile)
        for poi in combinetf_input.get_poi_names(fitresult):
            updatePercentDict(args, fitresult, dfpercents1, orginal_bin_unc_1, etaRebins1[i], poi)
    print(dfpercents1)
    
    if args.addPercentPlot or args.cumulativePercents:
        for i in range(len(etaRebins2)):
            inputFile = inFolder2 + '/' + inputFiles2[i]
            fitresult = combinetf_input.get_fitresult(inputFile)
            for poi in combinetf_input.get_poi_names(fitresult):
                updatePercentDict(args, fitresult, dfpercents2, orginal_bin_unc_2, etaRebins2[i], poi, lowPU = False)
        print(dfpercents2)

        df_percent_combined = pd.concat([dfpercents1, dfpercents2], axis=1)
        df_percent_combined.sort_values(by='Number of Eta Bins2', inplace=True)
        print(df_percent_combined)
    
    # plt.figure(figsize=(8, 8))
    # hep.cms.label(fontsize=20, data=False, label="Projection", com=13.6)
    
    # colors = ['steelblue', 'steelblue', 'black', 'crimson', 'crimson']
    
    # for i, column in enumerate(df_combined.columns):
    #     if column != 'Number of Eta Bins':  # Exclude the Eta bin column from plotting
    #         if "PDF" not in column:
    #             plt.plot(df_combined['Number of Eta Bins'], df_combined[column], label=column, marker='o', color = colors[i])
    #         else:
    #             plt.plot(df_combined['Number of Eta Bins'], df_combined[column], label=column, marker='o', linestyle = "--", color = colors[i])
    # plt.xlim(0, 50)
    # plt.ylim(0, 45)
    # plt.title("pt-eta-charge, scaled lumi, onlyPDF: highPU vs lowPU. \n Lumi: 1.0 fb^-1, PDF: CT18Z", y = 1.05, fontsize = 16)
    # plt.xlabel("Number of Eta Bins", fontsize = 16)
    # plt.ylabel("Uncertainty in $m_{W}$ (MeV)", fontsize = 16)
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    # plt.grid(True)
    # plt.legend(fontsize=16)
    # plt.savefig(args.outfile)
    # if args.addPDF:
    #     pdfFileName = (args.outfile).split('.')[0] + ".pdf"
    #     plt.savefig(pdfFileName, format = 'pdf')
    
    # if args.addPercentPlot:
    #     plt.figure(figsize=(8, 8))
    #     hep.cms.label(fontsize=20, data=False, label = "Projection", com=13.6)
    #     plt.title("pt-eta-charge, scaled lumi, onlyPDF: Percent change from total.\n Lumi: 1.0 fb^-1, PDF: CT18Z", y = 1.05, fontsize = 16)
        
    #     print(f"all columns: {df_percent_combined.columns}")
        
    #     plt.plot(df_percent_combined['Number of Eta Bins'], df_percent_combined['LowPU (mt-eta-charge) Percent Decreases in Unc'], label='LowPU Percent Decreases in Unc', marker='o', color = "steelblue")
    #     plt.plot(df_percent_combined['Number of Eta Bins2'], df_percent_combined['HighPU (pt-eta-charge) Percent Decreases in Unc'], label='HighPU Percent Decreases in Unc', marker='o', color = "crimson")

    #     plt.xlabel("Number of Eta Bins", fontsize=14)
    #     plt.ylabel("Percent Decrease in Uncertainty (%)", fontsize=14)
    #     plt.xticks(fontsize=14)
    #     plt.yticks(fontsize=14)
    #     plt.xlim(0, 50)
    #     plt.ylim(-5, 25)
    #     plt.grid(True)
    #     plt.legend(fontsize=16)
    #     plt.savefig(args.outfile_PercentPlot)
