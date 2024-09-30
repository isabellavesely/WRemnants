import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
from matplotlib import cm
import numpy as np

if __name__ == "__main__":
    
    x_pu_vals = [1, 2, 5, 10]
    y_pu_unc = [9.03, 9.42, 10.55, 12.26]
    
    fig = plt.figure()
    ax1 = fig.add_subplot()
    hep.cms.label(ax=ax1, lumi="1", label="Preliminary", com=13.6)
    
    ax1.plot(x_pu_vals, y_pu_unc)
    ax1.scatter(x_pu_vals, y_pu_unc, marker='o')
    
    ax1.set_title("2023 Wminus sample: Data Stat. Unc vs PU", y = 1.07, fontsize = 16)
    ax1.set_xlabel("Pile-up "+ r'$\langle\mu\rangle$', fontsize='14')
    ax1.set_ylabel('Data Statistics Uncertainty (MeV)', fontsize='14')
    ax1.grid()

    ax1.set_xlim([0, 12])
    ax1.set_ylim([6, 14])
    
    fname = "/home/submit/ivesely/public_html/studies_lowpu/plots/pileUp_statUnc_2023.pdf"
    fig.savefig(fname, bbox_inches='tight')
    fig.savefig(fname.replace(".pdf", ".png"), bbox_inches='tight')
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    hep.cms.label(ax=ax2, lumi="1", label="Preliminary", com=13.6)
    y_percent_unc = []
    for unc in y_pu_unc:
        y_percent_unc.append((unc - 9.03)/(9.03)*100)
        
    ax2.grid(True, which='major', linewidth=1)
    ax2.grid(True, which='minor', linestyle='--', linewidth=0.5)
    ax2.minorticks_on()
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(10))  # Major ticks every 10 MeV
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(1))
        
    ax2.plot(x_pu_vals, y_percent_unc, zorder = 5, color = "crimson")
    ax2.scatter(x_pu_vals, y_percent_unc, marker='o', zorder = 5, color = "crimson")
    
    ax2.set_title("2023 Wminus sample: Percent Increases in Data Stat. Unc from PU1", y = 1.07, fontsize = 16)
    ax2.set_xlabel("Pile-up "+ r'$\langle\mu\rangle$', fontsize='14')
    ax2.set_ylabel('Data Statistics Uncertainty (MeV)', fontsize='14')
    
    # ax2.grid()
    ax2.set_xlim([0, 12])
    ax2.set_ylim([-5, 45])

    fname = "/home/submit/ivesely/public_html/studies_lowpu/plots/pileUp_statUnc_percents_2023.pdf"
    fig2.savefig(fname, bbox_inches='tight')
    fig2.savefig(fname.replace(".pdf", ".png"), bbox_inches='tight')
    
    # fig3 = plt.figure()
    # ax3 = fig3.add_subplot()
    # hep.cms.label(ax=ax2, lumi="1", label="Preliminary", com=13.6)
    # y_percent_unc = []
    # for i, unc in enumerate(y_pu_unc):
    #     if i != 0:
    #         y_percent_unc.append((unc - y_pu_unc[i-1])/(y_pu_unc[i-1])*100)
    #     else:
    #         y_percent_unc.append(0)
        
    # ax3.plot(x_pu_vals, y_percent_unc)
    # ax3.scatter(x_pu_vals, y_percent_unc, marker='o')
    
    # ax3.set_title("2023 Wminus sample: Percent Increases in Data Stat. Unc vs PU", y = 1.07, fontsize = 16)
    # ax3.set_xlabel("Pile-up "+ r'$\langle\mu\rangle$', fontsize='14')
    # ax3.set_ylabel('Data Statistics Uncertainty (MeV)', fontsize='14')
    # ax3.grid()
    # ax3.set_xlim([0, 12])
    # ax3.set_ylim([-5, 46])

    # fname = "/home/submit/ivesely/public_html/studies_lowpu/plots/pileUp_statUnc_percents_2023_2.pdf"
    # fig3.savefig(fname, bbox_inches='tight')
    # fig3.savefig(fname.replace(".pdf", ".png"), bbox_inches='tight')
    
    
    
