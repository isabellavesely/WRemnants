#!/usr/bin/env python3

# example
# python scripts/analysisTools/w_mass_13TeV/validateVetoDilepton.py /scratch/mciprian/CombineStudies/testZmumuVeto/fromDavide/mz_dilepton_scetlib_dyturboCorr_maxFiles_m1_vetoGlobal.hdf5 scripts/analysisTools/plots/fromMyWremnants/testZmumuVeto/validateVetoSF_fromDavide/global/ -n nominal_muonsonly --plotNonTrig

# histogram template:
#      nominal_muonsonly
# Axes = ('trigMuons_eta0', 'trigMuons_pt0', 'trigMuons_charge0', 'nonTrigMuons_eta0', 'nonTrigMuons_pt0', 'nonTrigMuons_charge0')
# Regular(24, -2.4, 2.4, name='trigMuons_eta0')
# Regular(34, 26, 60, name='trigMuons_pt0')
# Regular(2, -2, 2, underflow=False, overflow=False, name='trigMuons_charge0')
# Regular(24, -2.4, 2.4, name='nonTrigMuons_eta0')
# Regular(50, 15, 65, name='nonTrigMuons_pt0')
# Regular(2, -2, 2, underflow=False, overflow=False, name='nonTrigMuons_charge0')

from wremnants.datasets.datagroups import Datagroups
from wremnants import histselections as sel
#from wremnants import plot_tools,theory_tools,syst_tools
from utilities import boostHistHelpers as hh
from utilities import common, logging
from utilities.io_tools import input_tools, output_tools

import narf
import wremnants
from wremnants import theory_tools,syst_tools
import hist

import numpy as np

import pickle
import lz4.frame

import argparse
import os
import shutil
import re

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

from scripts.analysisTools.plotUtils.utility import *
from scripts.analysisTools.w_mass_13TeV.plotPrefitTemplatesWRemnants import plotPrefitHistograms

if __name__ == "__main__":
    parser = common_plot_parser()
    parser.add_argument("inputfile", type=str, nargs=1, help="Input file with histograms (pkl.lz4 or hdf5 file)")
    parser.add_argument("outdir",   type=str, nargs=1, help="Output folder")
    parser.add_argument("-n", "--baseName", type=str, help="Histogram name in the file (it depends on what study you run)", default="nominal_muonsonly")
    parser.add_argument('-p','--processes', default=None, nargs='*', type=str,
                        help='Choose what processes to plot, otherwise all are done')
    #parser.add_argument("--mtRange", type=float, nargs=2, default=[0,40], choices=[-1.0, 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 120.], help="Apply mT cut, if upper edge is negative integrate the overflow")
    parser.add_argument("-c", "--charges", type=int, default=[-1, 1], nargs='+', choices=[-1, 1], help="Charge selection for chosen muon")
    parser.add_argument("--rr", "--ratio-range", dest="ratioRange", default=(0.95,1.05), type=float, nargs=2, help="Range for ratio plot")
    parser.add_argument('--plotNonTrig', action='store_true', help='Plot non triggering muon (with the veto selection), otherwise plot triggering muon')
    parser.add_argument('--scaleProc', default=None, nargs='*', type=str,
                        help='Apply scaling factor to process by name, with syntax proc=scale=charge (=charge can be omitted, if given it must be plus or minus). Can specify multiple times')
    args = parser.parse_args()

    logger = logging.setup_logger(os.path.basename(__file__), args.verbose)

    allCharges = { -1 : "minus", 1 : "plus" }
    charges = allCharges if len(args.charges) == 2 else {args.charges[0]: allCharges[args.charges[0]]}
    muonTag = "nonTrigMuon" if args. plotNonTrig else "trigMuons"
    muonTitle = "Non triggering" if args. plotNonTrig else "Triggering"
    ###############################################################################################
    fname = args.inputfile[0]
    
    ROOT.TH1.SetDefaultSumw2()

    outdir_original = f"{args.outdir[0]}/"
    outdir = createPlotDirAndCopyPhp(outdir_original, eoscp=args.eoscp)
    
    canvas = ROOT.TCanvas("canvas", "", 800, 700)
    cwide = ROOT.TCanvas("cwide","",2400,600)                      
    adjustSettings_CMS_lumi()
    canvas1D = ROOT.TCanvas("canvas1D", "", 800, 900)

    groups = Datagroups(fname, mode="z_dilepton")
    datasets = groups.getNames()
    if args.processes is not None and len(args.processes):
        datasets = list(filter(lambda x: x in args.processes, datasets))
    else:
        datasets = list(filter(lambda x: x != "QCD", datasets))

    logger.debug(f"Using these processes: {datasets}")
    inputHistName = args.baseName
    groups.setNominalName(inputHistName)
    groups.loadHistsForDatagroups(inputHistName, syst="", procsToRead=datasets, applySelection=False)
    histInfo = groups.getDatagroups() # keys are same as returned by groups.getNames()
    s = hist.tag.Slicer()

    scaleDict = {"plus": {},
                 "minus" : {}}
    if args.scaleProc:
        for sp in args.scaleProc:
            tokens = sp.split("=")
            proc,scale = tokens[0],float(tokens[1])
            if len(tokens) == 3:
                ch = tokens[2]
                scaleDict[ch][proc] = scale
            else:
                scaleDict["plus"][proc] = scale
                scaleDict["minus"][proc] = scale
    
    for charge in charges.keys():
        chargeTag = charges[charge]
        chargeBin = 0 if charge == -1 else 1
        outdirTag = f"{muonTag}_{chargeTag}/"
        outdirCharge = f"{outdir}/{outdirTag}/"
        createPlotDirAndCopyPhp(outdirCharge, eoscp=args.eoscp)

        scaleDictByCharge = scaleDict[charges[charge]]

        hdata2D = None
        hmc2D = []
        for d in datasets:
            logger.info(f"Running on process {d}")
            hin = histInfo[d].hists[inputHistName]
            logger.debug(hin.axes)
            ###
            # select charge and integrate other muon
            if args.plotNonTrig:
                h = hin[{"nonTrigMuons_charge0": s[chargeBin],
                         "trigMuons_pt0" : s[::hist.sum],
                         "trigMuons_eta0" : s[::hist.sum],
                         "trigMuons_charge0" : s[::hist.sum]}]
            else:
                h = hin[{"trigMuons_charge0": s[chargeBin],
                         "nonTrigMuons_pt0" : s[::hist.sum],
                         "nonTrigMuons_eta0" : s[::hist.sum],
                         "nonTrigMuons_charge0" : s[::hist.sum]}]

            ###
            logger.debug(h.axes)
            if d =="Data":
                hdata2D = narf.hist_to_root(h)
                hdata2D.SetName(f"{d}_{chargeTag}")
                hdata2D.SetTitle(f"{d} {chargeTag}")
            else:
                hmc2D.append(narf.hist_to_root(h))
                hmc2D[-1].SetName(f"{d}_{chargeTag}")
                hmc2D[-1].SetTitle(f"{d}_{chargeTag}")
                if d in scaleDictByCharge:
                    hmc2D[-1].Scale(scaleDictByCharge[d])
                    logger.info(f"Scaling process {d} for charge {chargeTag} by {scaleDictByCharge[d]}")
        # end of process loop
        plotPrefitHistograms(hdata2D, hmc2D, outdirCharge, xAxisName=f"{muonTitle} muon #eta", yAxisName=f"{muonTitle} muon p_{{T}} (GeV)",
                             chargeLabel=chargeTag, canvas=canvas, canvasWide=cwide, canvas1D=canvas1D,
                             ratioRange=args.ratioRange, lumi=16.8)

    copyOutputToEos(outdir, outdir_original, eoscp=args.eoscp)
