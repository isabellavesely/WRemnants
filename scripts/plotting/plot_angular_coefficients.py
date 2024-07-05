from utilities import logging, common
from utilities.io_tools import input_tools, output_tools
from utilities.styles import styles
from wremnants import theory_tools, plot_tools

import numpy as np
import matplotlib as mpl
import hist
import math
import mplhep as hep
import h5py

def a0(theta, phi):
    return 1. + np.cos(theta)**2

def a1(theta, phi):
    return 0.5 * (1. - 3. * np.cos(theta)**2)

def a2(theta, phi):
    return np.sin(2*theta) * np.cos(phi)
    
def a3(theta, phi):
    return 0.5 * np.sin(theta)**2 * np.cos(2*phi)

def a4(theta, phi):
    return np.sin(theta) * np.cos(phi)

def a5(theta, phi):
    return np.cos(theta)

def a6(theta, phi):
    return np.sin(theta)**2 * np.sin(2*phi)

def a7(theta, phi):
    return np.sin(2*theta) * np.sin(phi)

def a8(theta, phi):
    return np.sin(theta) * np.sin(phi)

def load_moments_coefficients(filename, process="Z"):
    with h5py.File(filename, "r") as ff:
        out = input_tools.load_results_h5py(ff)

    moments = out[process]

    corrh = theory_tools.moments_to_angular_coeffs(moments)

    if 'muRfact' in corrh.axes.name:
        corrh = corrh[{'muRfact' : 1.j,}]
    if 'muFfact' in corrh.axes.name:
        corrh = corrh[{'muFfact' : 1.j,}]
    
    axes_names = ['massVgen','absYVgen','ptVgen','chargeVgen', 'helicity']
    if not list(corrh.axes.name) == axes_names:
        raise ValueError (f"Axes [{corrh.axes.name}] are not the ones this functions expects ({axes_names})")
    
    if np.count_nonzero(corrh[{"helicity" : -1.j}] == 0):
        logger.warning("Zeros in sigma UL for the angular coefficients will give undefined behaviour!")
    # histogram has to be without errors to load the tensor directly
    corrh_noerrs = hist.Hist(*corrh.axes, storage=hist.storage.Double())
    corrh_noerrs.values(flow=True)[...] = corrh.values(flow=True)

    return corrh_noerrs

if __name__ == '__main__':
    parser = common.plot_parser()
    parser.add_argument("infile", help="Moments file `w_z_moments.hdf` with coefficients produced in w_z_gen_dists.py histmaker")

    args = parser.parse_args()
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    if args.infile is not None:
        raise NotImplementedError("Using moments file from histmaker output is not yet supported")
        # hcoeff = load_moments_coefficients(args.infile)

    outdir = output_tools.make_plot_dir(args.outpath, args.outfolder, eoscp=args.eoscp)

    colors = mpl.colormaps["tab10"]
    linestyles = ["solid","dotted","dashed","dashdot"]

    npoints = 100

    theta_edges = np.linspace(-np.pi, 0, npoints+1)
    cosTheta_edges = np.cos(theta_edges)
    theta = theta_edges[:-1] + np.diff(theta_edges)/2.
    cosTheta = cosTheta_edges[:-1] + np.diff(cosTheta_edges)/2.

    phi = np.linspace(-np.pi, np.pi, npoints)
    x, y = np.meshgrid(theta, phi)
    
    histos = []
    for ai in [a0, a1, a2, a3, a4, a5, a6, a7, a8]:
        histo = hist.Hist(
            hist.axis.Variable(cosTheta_edges, name = "cosTheta", underflow=False, overflow=False),
            hist.axis.Regular(npoints, -math.pi, math.pi, circular = True, name = "phi"), 
            storage=hist.storage.Double()
        )

        histo.values()[...] = ai(x,y).T

        histos.append(histo)

    scales = [0., 20./3., 5., 20., 4., 4., 5., 5., 4.]
    offsets = [1., 2./3., 0., 0., 0., 0., 0., 0., 0.]

    for moments in (True, False):

        # 2D plots
        xlabel = styles.xlabels.get("costhetastarll", "cosTheta")
        ylabel = styles.xlabels.get("phistarll", "phi")
        for i, histo in enumerate(histos):
            fig = plot_tools.makeHistPlot2D(histo, cms_label=args.cmsDecor, xlabel=xlabel, ylabel=ylabel, scaleleg=args.scaleleg, has_data=False)

            if i == 0:
                idx = "UL"
            else:
                idx = str(i-1)

            outfile = "angular_"
            if moments:
                outfile += "moment"
            else:
                outfile += "coefficient"
            outfile += f"_{idx}_cosTheta_phi"
            if args.postfix:
                outfile += f"_{args.postfix}"
            plot_tools.save_pdf_and_png(outdir, outfile)
            plot_tools.write_index_and_log(outdir, outfile, args=args)

        # 1D plots
        for axis_name, ais in [("cosTheta", [0,1,5]), ("phi", [3,4,6,8])]:
            h1ds = [histo.project(axis_name)/np.product([histo.axes[n].size for n in histo.axes.name if n != axis_name]) for histo in histos]

            fig, ax1 = plot_tools.figure(h1ds[0], xlabel=styles.xlabels.get(f"{axis_name.lower()}starll", axis_name), ylabel="Frequency",
                grid=True, automatic_scale=False, width_scale=1.2, logy=False)    
            
            j=0
            for i, h1d in enumerate(h1ds):
                if i not in ais:
                    continue

                val_x = h1d.axes[0].centers
                val_y = h1d.values()
                if i == 0:
                    idx = "\mathrm{UL}"
                else:
                    idx = str(i-1)
                if moments:
                    val_y = val_y * scales[i] + offsets[i]
                    label=f"$\mathrm{{M}}_{idx}$"
                else:
                    label=f"$\mathrm{{A}}_{idx}$"

                ax1.plot(val_x, val_y, color=colors(i), linestyle=linestyles[j], label=label)
                j += 1

            plot_tools.addLegend(ax1, ncols=2, text_size=12, loc="upper left")
            plot_tools.fix_axes(ax1, logy=False)

            scale = max(1, np.divide(*ax1.get_figure().get_size_inches())*0.3)
            hep.cms.label(ax=ax1, lumi=None, fontsize=20*args.scaleleg*scale, 
                label=args.cmsDecor, data=False)

            outfile = "angular_"
            if moments:
                outfile += "moments"
            else:
                outfile += "coefficients"
            outfile += f"_{axis_name}"
            if args.postfix:
                outfile += f"_{args.postfix}"
            plot_tools.save_pdf_and_png(outdir, outfile)
            plot_tools.write_index_and_log(outdir, outfile, args=args)

    if output_tools.is_eosuser_path(args.outpath) and args.eoscp:
        output_tools.copy_to_eos(outdir, args.outpath, args.outfolder)

