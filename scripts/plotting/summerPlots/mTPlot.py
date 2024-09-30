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

    procs = ["Wplusmunu", "Wminusmunu"]
    histName = "transverseMass"
    data = input_tools.read_all_and_scale("outputFolder/new_lowPU_samples/mw_lowPU_mu_nnpdf31_2017.hdf5", procs, [histName], lumi=True)[0]
    data = data[{'passIso' : True}].project("mt")
    data = data[s[40j:120j]]
    # data = data
    print(data)


    fig = plt.figure()

    ax1 = fig.add_subplot()
    hep.cms.label(ax=ax1, lumi="0.2", label="Preliminary", com=13.6)

    ax1.set_xlabel(r"m$_{T}^{\ell\nu}$ (GeV)")
    ax1.set_ylabel("Events")
    ax1.set_xlim([40, 120])
    ax1.set_ylim([0, 60000])
    # ax1.autoscale(axis='y')


    hep.histplot(
        data,
        histtype="step",
        label=r"$\langle\mu\rangle$ = 1",
        ax=ax1,
        zorder=1,
    )

    plt.legend(loc="upper left")

    fname = "/home/submit/ivesely/public_html/studies_lowpu/plots/mt_comparison_2017.pdf"
    fig.savefig(fname, bbox_inches='tight')
    fig.savefig(fname.replace(".pdf", ".png"), bbox_inches='tight')
