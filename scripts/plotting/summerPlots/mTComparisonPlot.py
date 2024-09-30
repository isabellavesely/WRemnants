from wremnants import plot_tools,theory_tools
from utilities import boostHistHelpers as hh
from utilities.io_tools import input_tools, output_tools
import hist
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib import cm
import numpy as np
import boost_histogram as bh

s = hist.tag.Slicer()

if __name__ == "__main__":

    procs = ["Wminusmunu"]
    histName = "transverseMass"

    # pu_files = ["Wminusmunu/mw_lowPU_mu_nnpdf31_PUAVE1.hdf5", "Wminusmunu/mw_lowPU_mu_nnpdf31_PUAVE2.hdf5", "Wminusmunu_poisson/mw_lowPU_mu_nnpdf31_PUAVE5_target3.hdf5", "Wminusmunu/mw_lowPU_mu_nnpdf31_2017_Wminus.hdf5", "Wminusmunu/mw_lowPU_mu_nnpdf31_PUAVE5.hdf5", "Wminusmunu/mw_lowPU_mu_nnpdf31_PUAVE10.hdf5"]
    # pu_labels = [r"$\langle\mu\rangle$ = 1", r"$\langle\mu\rangle$ = 2", r"$\langle\mu\rangle$ = 3, poisson weighted. Unc: 9.74 MeV", r"$\langle\mu\rangle$ = 3, 2017. Unc: 4.32 MeV", r"$\langle\mu\rangle$ = 5", r"$\langle\mu\rangle$ = 10"]
    # pu_files = ["Wminusmunu_poisson/mw_lowPU_mu_nnpdf31_PUAVE5_target3.hdf5", "Wminusmunu/mw_lowPU_mu_nnpdf31_2017_Wminus.hdf5", "mw_lowPU_mu_nnpdf31_2017.hdf5"]
    # pu_labels = [r"$\langle\mu\rangle$ = 3, poisson from 5", r"$\langle\mu\rangle$ = 3, 2017 Wminus", r"$\langle\mu\rangle$ = 3, 2017"]
    pu_files = ["Wminusmunu/mw_lowPU_mu_nnpdf31_2017_Wminus.hdf5", "Wminusmunu_poisson/mw_lowPU_mu_nnpdf31_PUAVE5_target3.hdf5"]
    pu_labels = [r"$\langle\mu\rangle$ = 3 (2017). Stat unc: 4.32 MeV", r"$\langle\mu\rangle$ = 3 (2023 poisson-scaled). Stat unc: 9.66 MeV"]
    pu_lumis = [0.2, 1.0] 
    
    target_lumi = 1.0

    # folder_str = "outputFolder/new_lowPU_samples/Wminusmunu/"
    folder_str = "outputFolder/new_lowPU_samples/"
    pu_files = [folder_str + file for file in pu_files]
    
    hists = []
    for pu, lumi in zip(pu_files, pu_lumis):
        hist = input_tools.read_all_and_scale(pu, procs, [histName], False)[0]
        hist = hist[{'passIso' : True}].project("mt")
        hist = hist[s[40j:120j]]
        
        # scaling_factor = target_lumi / lumi
        # hist = hist * scaling_factor
        
        print(hist)
        hists.append(hist)

    fig = plt.figure()

    ax1 = fig.add_subplot() 
    hep.cms.label(ax=ax1, lumi="1", label="Preliminary", com=13.6)

    ax1.set_xlabel(r"m$_{T}^{\ell\nu}$ (GeV)")
    ax1.set_ylabel("Events")
    ax1.set_xlim([40, 120])
    ax1.set_ylim([0, 500])
    #ax1.autoscale(axis='y')

    for i,hist in enumerate(hists):
        hep.histplot(
            hist,
            histtype="step",
            label=pu_labels[i],
            ax=ax1,
            zorder=1,
        )

    plt.legend(loc="upper left")

    fname = "/home/submit/ivesely/public_html/studies_lowpu/plots/mt_comparison_overlay_mu3_minus_new.pdf"
    fig.savefig(fname, bbox_inches='tight')
    fig.savefig(fname.replace(".pdf", ".png"), bbox_inches='tight')
