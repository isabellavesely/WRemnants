#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools,combine_helpers,combine_theory_helper
from wremnants import histselections as sel
from wremnants.datasets.datagroups import Datagroups
from utilities import common, logging, input_tools
import itertools
import argparse
import os
import pathlib
import hist
import copy
import math
import time

data_dir = common.data_dir

def make_parser(parser=None):
    if not parser:
        parser = common.common_parser_combine()
    parser.add_argument("--fitvar", nargs="+", help="Variable to fit", default=["pt", "eta"])
    parser.add_argument("-n", "--baseName", type=str, help="Histogram name in the file (e.g., 'nominal')", default="nominal")
    parser.add_argument("--noEfficiencyUnc", action='store_true', help="Skip efficiency uncertainty (useful for tests, because it's slow). Equivalent to --excludeNuisances '.*effSystTnP|.*effStatTnP' ")
    parser.add_argument("--ewUnc", action='store_true', help="Include EW uncertainty")
    parser.add_argument("--widthUnc", action='store_true', help="Include uncertainty on W and Z width")
    parser.add_argument("--pseudoData", type=str, help="Hist to use as pseudodata")
    parser.add_argument("--pseudoDataIdx", type=str, default="0", help="Variation index to use as pseudodata")
    parser.add_argument("--pseudoDataFile", type=str, help="Input file for pseudodata (if it should be read from a different file)", default=None)
    parser.add_argument("--pseudoDataProcsRegexp", type=str, default=".*", help="Regular expression for processes taken from pseudodata file (all other processes are automatically got from the nominal file). Data is excluded automatically as usual")
    parser.add_argument("-x",  "--excludeNuisances", type=str, default="", help="Regular expression to exclude some systematics from the datacard")
    parser.add_argument("-k",  "--keepNuisances", type=str, default="", help="Regular expression to keep some systematics, overriding --excludeNuisances. Can be used to keep only some systs while excluding all the others with '.*'")
    parser.add_argument("--scaleMuonCorr", type=float, default=1.0, help="Scale up/down dummy muon scale uncertainty by this factor")
    parser.add_argument("--noHist", action='store_true', help="Skip the making of 2D histograms (root file is left untouched if existing)")
    parser.add_argument("--effStatLumiScale", type=float, default=None, help="Rescale equivalent luminosity for efficiency stat uncertainty by this value (e.g. 10 means ten times more data from tag and probe)")
    parser.add_argument("--binnedScaleFactors", action='store_true', help="Use binned scale factors (different helpers and nuisances)")
    parser.add_argument("--isoEfficiencySmoothing", action='store_true', help="If isolation SF was derived from smooth efficiencies instead of direct smoothing")
    parser.add_argument("--axlim", type=float, default=[], nargs='*', help="Restrict axis to this range (assumes pairs of values by axis, with trailing axes optional)")
    parser.add_argument("--unfolding", action='store_true', help="Prepare datacard for unfolding")
    parser.add_argument("--genAxes", type=str, default=None, nargs="+", help="Specify which gen axis should be used in unfolding, if 'None', use all (inferred from metadata).")
    parser.add_argument("--fitXsec", action='store_true', help="Fit signal inclusive cross section")
    parser.add_argument("--correlatedNonClosureNuisances", action='store_true', help="get systematics from histograms for the Z non-closure nuisances without decorrelation in eta and pt")
    parser.add_argument("--sepImpactForNC", action="store_true", help="use a dedicated impact gropu for non closure nuisances, instead of putting them in muonScale")
    parser.add_argument("--genModel", action="store_true", help="Produce datacard with the xnorm as model (binned according to axes defined in --fitvar)")
    parser.add_argument("--simultaneousABCD", action="store_true", help="Produce datacard for simultaneous fit of ABCD regions")
    return parser


def setup(args,xnorm=False):   

    # NOTE: args.filterProcGroups and args.excludeProcGroups should in principle not be used together
    #       (because filtering is equivalent to exclude something), however the exclusion is also meant to skip
    #       processes which are defined in the original process dictionary but are not supposed to be (always) run on
    if args.addQCDMC or "QCD" in args.filterProcGroups:
        logger.warning("Adding QCD MC to list of processes for the fit setup")
    else:
        if "QCD" not in args.excludeProcGroups:
            logger.warning("Automatic removal of QCD MC from list of processes. Use --filterProcGroups 'QCD' or --addQCDMC to keep it")
            args.excludeProcGroups.append("QCD")
    filterGroup = args.filterProcGroups if args.filterProcGroups else None
    excludeGroup = args.excludeProcGroups if args.excludeProcGroups else None
    if args.simultaneousABCD and (excludeGroup is None or "Fake" not in excludeGroup):
        excludeGroup.append("Fake")
    logger.debug(f"Filtering these groups of processes: {args.filterProcGroups}")
    logger.debug(f"Excluding these groups of processes: {args.excludeProcGroups}")

    datagroups = Datagroups(args.inputFile, excludeGroups=excludeGroup, filterGroups=filterGroup, applySelection= not xnorm and not args.simultaneousABCD)

    if args.axlim or args.rebin:
        if len(args.axlim) % 2 or len(args.axlim)/2 > len(args.fitvar) or len(args.rebin) > len(args.fitvar):
            raise ValueError("Inconsistent rebin or axlim arguments. axlim must be at most two entries per axis, and rebin at most one")

        sel = {}
        for var,low,high,rebin in itertools.zip_longest(args.fitvar, args.axlim[::2], args.axlim[1::2], args.rebin):
            s = hist.tag.Slicer()
            if low is not None and high is not None:
                logger.info(f"Restricting the axis '{var}' to range [{low}, {high}]")
                sel[var] = s[complex(0, low):complex(0, high):hist.rebin(rebin) if rebin else None]
            elif rebin:
                sel[var] = s[hist.rebin(rebin)]
            if rebin:
                logger.info(f"Rebinning the axis '{var}' by [{rebin}]")

        logger.info(f"Will apply the global selection {sel}")
        datagroups.setGlobalAction(lambda h: h[sel])

    wmass = datagroups.wmass
    wlike = datagroups.wlike
    lowPU = datagroups.lowPU
    dilepton = datagroups.dilepton

    constrainMass = dilepton and not "mll" in args.fitvar

    if wmass:
        base_group = "Wenu" if datagroups.flavor == "e" else "Wmunu"
    else:
        base_group = "Zee" if datagroups.flavor == "ee" else "Zmumu"

    if args.unfolding and args.fitXsec:
        raise ValueError("Options --unfolding and --fitXsec are incompatible. Please choose one or the other")
    elif args.fitXsec:
        datagroups.unconstrainedProcesses.append(base_group)
    elif args.unfolding:
        constrainMass = True
        datagroups.setGenAxes(args.genAxes)
        
        if wmass:
            # gen level bins, split by charge
            datagroups.defineSignalBinsUnfolding(base_group, "W_qGen0", member_filter=lambda x: x.name.startswith("Wminus"))
            datagroups.defineSignalBinsUnfolding(base_group, "W_qGen1", member_filter=lambda x: x.name.startswith("Wplus"))
            # out of acceptance contribution
            datagroups.groups[base_group].deleteMembers([m for m in datagroups.groups[base_group].members if not m.name.startswith("Bkg")])
        else:
            datagroups.defineSignalBinsUnfolding(base_group, "Z", member_filter=lambda x: x.name.startswith(base_group))
            # out of acceptance contribution
            datagroups.groups[base_group].deleteMembers([m for m in datagroups.groups[base_group].members if not m.name.startswith("Bkg")])

    if args.noHist and args.noStatUncFakes:
        raise ValueError("Option --noHist would override --noStatUncFakes. Please select only one of them")

    # Start to create the CardTool object, customizing everything
    cardTool = CardTool.CardTool(xnorm=xnorm)
    cardTool.setDatagroups(datagroups)
    logger.debug(f"Making datacards with these processes: {cardTool.getProcesses()}")
    if args.absolutePathInCard:
        cardTool.setAbsolutePathShapeInCard()
    cardTool.setProjectionAxes(args.fitvar)
    if wmass and args.simultaneousABCD:
        cardTool.setChannels(["inclusive"])
        cardTool.setWriteByCharge(False)
        fitvars = ["passIso", "passMT", *args.fitvar]
        cardTool.setProjectionAxes(fitvars)
        cardTool.unroll=True
    if args.sumChannels or xnorm or dilepton:
        cardTool.setChannels(["inclusive"])
        cardTool.setWriteByCharge(False)
    if xnorm:
        datagroups.select_xnorm_groups() # only keep processes where xnorm is defined
        datagroups.deleteGroup("Fake")
        if args.unfolding:
            cardTool.setProjectionAxes(["count"])
        else:
            if wmass:
                # add gen charge as additional axis
                datagroups.groups[base_group].add_member_axis("qGen", datagroups.results, 
                    member_filters={-1: lambda x: x.name.startswith("Wminus"), 1: lambda x: x.name.startswith("Wplus")}, 
                    hist_filter=lambda x: x.startswith("xnorm"))
            cardTool.unroll = True
            # remove projection axes from gen axes, otherwise they will be integrated before
            datagroups.setGenAxes([a for a in datagroups.gen_axes if a not in cardTool.project])
    else:
        cardTool.setHistName(args.baseName)
        cardTool.setNominalName(args.baseName)
    if args.unfolding:
        cardTool.addPOISumGroups(additional_axes=["qGen"] if base_group[0] == "W" else None)
    if args.noHist:
        cardTool.skipHistograms()
    cardTool.setFakeName(args.qcdProcessName)
    cardTool.setSpacing(52)
    if args.noStatUncFakes:
        cardTool.setProcsNoStatUnc(procs=args.qcdProcessName, resetList=False)
    cardTool.setCustomSystForCard(args.excludeNuisances, args.keepNuisances)
    if args.pseudoData:
        cardTool.setPseudodata(args.pseudoData, args.pseudoDataIdx, args.pseudoDataProcsRegexp)
        if args.pseudoDataFile:
            cardTool.setPseudodataDatagroups(make_datagroup(args.pseudoDataFile,
                                                                  excludeGroups=excludeGroup,
                                                                  filterGroups=filterGroup,
                                                                  applySelection= not xnorm)
            )
    cardTool.setLumiScale(args.lumiScale)

    logger.info(f"cardTool.allMCProcesses(): {cardTool.allMCProcesses()}")
        
    passSystToFakes = wmass and not args.skipSignalSystOnFakes and args.qcdProcessName not in excludeGroup and (filterGroup == None or args.qcdProcessName in filterGroup) and not xnorm

    cardTool.addProcessGroup("single_v_samples", lambda x: x[0] in ["W", "Z"] and ("mu" in x or "tau" in x))
    if wmass:
        cardTool.addProcessGroup("single_v_nonsig_samples", lambda x: x[0] == "Z"and ("mu" in x or "tau" in x))

    cardTool.addProcessGroup("single_vmu_samples", lambda x: x[0] in ["W", "Z"] and "mu" in x)
    cardTool.addProcessGroup("signal_samples", lambda x: ((x[0] == "W" and wmass) or (x[0] == "Z" and not wmass)) and "mu" in x)
    cardTool.addProcessGroup("signal_samples_inctau", lambda x: ((x[0] == "W" and wmass) or (x[0] == "Z" and not wmass)) and ("mu" in x or "tau" in x))
    cardTool.addProcessGroup("MCnoQCD", lambda x: x not in ["QCD", "Data"])

    logger.info(f"All MC processes {cardTool.procGroups['MCnoQCD']}")
    logger.info(f"Single V samples: {cardTool.procGroups['single_v_samples']}")
    if wmass:
        logger.info(f"Single V no signal samples: {cardTool.procGroups['single_v_nonsig_samples']}")
    logger.info(f"Signal samples: {cardTool.procGroups['signal_samples']}")

    constrainedZ = constrainMass and not wmass
    label = 'W' if wmass else 'Z'
    massSkip = [(f"^massShift[W|Z]{i}MeV.*",) for i in range(0, 110 if constrainedZ else 100, 10)]
    if wmass:
        cardTool.addSystematic(f"massWeightZ",
                                processes=['single_v_nonsig_samples'],
                                group=f"massShiftZ",
                                skipEntries=massSkip[:]+[("^massShiftZ100MeV.*",)],
                                mirror=False,
                                noConstraint=False,
                                systAxes=["massShift"],
                                passToFakes=passSystToFakes,
        )

    if not (constrainMass or wmass):
        massSkip.append(("^massShift.*2p1MeV.*",))

    cardTool.addSystematic(f"massWeight{label}", 
                            processes=["signal_samples_inctau"],
                            group=f"massShift{label}",
                            noiGroup=not constrainMass,
                            skipEntries=massSkip,
                            mirror=False,
                            #TODO: Name this
                            noConstraint=not constrainMass,
                            systAxes=["massShift"],
                            passToFakes=passSystToFakes,
    )

    if args.doStatOnly:
        # print a card with only mass weights
        logger.info("Using option --doStatOnly: the card was created with only mass nuisance parameter")
        return cardTool
    
    if args.widthUnc:
        widthSkipZ = [("widthZ2p49333GeV",), ("widthZ2p49493GeV",), ("widthZ2p4952GeV",)] 
        widthSkipW = [("widthW2p09053GeV",), ("widthW2p09173GeV",), ("widthW2p085GeV",)]
        if wmass and not xnorm:
            cardTool.addSystematic(f"widthWeightZ",
                                    processes=single_v_nonsig_samples,
                                    group=f"widthZ",
                                    skipEntries=widthSkipZ[:],
                                    mirror=False,
                                    systAxes=["width"],
                                    passToFakes=passSystToFakes,
            )
        cardTool.addSystematic(f"widthWeight{label}",
                                processes=signal_samples_inctau,
                                skipEntries=widthSkipZ[:] if label=="Z" else widthSkipW[:],
                                group=f"width{label}",
                                mirror=False,
                                #TODO: Name this
                                systAxes=["width"],
                                passToFakes=passSystToFakes,
        )

    if not xnorm:
        if wmass:
            cardTool.addSystematic("luminosity",
                                   processes=['MCnoQCD'],
                                   outNames=["lumiDown", "lumiUp"],
                                   group="luminosity",
                                   systAxes=["downUpVar"],
                                   labelsByAxis=["downUpVar"],
                                   passToFakes=passSystToFakes)

        else:
            # TOCHECK: no fakes here, most likely
            cardTool.addLnNSystematic("luminosity", processes=['MCnoQCD'], size=1.012, group="luminosity")
    else:
        pass

    if args.ewUnc:
        cardTool.addSystematic(f"horacenloewCorr", 
            processes=['single_v_samples'],
            mirror=True,
            group="theory_ew",
            systAxes=["systIdx"],
            labelsByAxis=["horacenloewCorr"],
            skipEntries=[(0, -1), (2, -1)],
            passToFakes=passSystToFakes,
        )

    to_fakes = passSystToFakes and not args.noQCDscaleFakes and not xnorm
    combine_helpers.add_pdf_uncertainty(cardTool, single_v_samples, passSystToFakes, from_corr=args.pdfUncFromCorr, scale=args.scalePdf)
    combine_helpers.add_modeling_uncertainty(cardTool, args.minnloScaleUnc, signal_samples_inctau, 
        single_v_nonsig_samples if not xnorm else [], to_fakes, args.resumUnc, wmass, scaleTNP=args.scaleTNP)

    if xnorm:
        return cardTool

    # Below: experimental uncertainties

    if wmass:
        #cardTool.addLnNSystematic("CMS_Fakes", processes=[args.qcdProcessName], size=1.05, group="MultijetBkg")
        cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06)
        cardTool.addLnNSystematic("CMS_VV", processes=["Diboson"], size=1.16)
        cardTool.addSystematic("luminosity",
                                processes=allMCprocesses_noQCDMC,
                                outNames=["lumiDown", "lumiUp"],
                                group="luminosity",
                                systAxes=["downUpVar"],
                                labelsByAxis=["downUpVar"],
                                passToFakes=passSystToFakes)
    else:
        cardTool.addLnNSystematic("CMS_background", processes=["Other"], size=1.15)
        cardTool.addLnNSystematic("luminosity", processes=allMCprocesses_noQCDMC, size=1.017 if lowPU else 1.012, group="luminosity")

    if not args.noEfficiencyUnc:

        ## this is only needed when using 2D SF from 3D with ut-integration, let's comment for now
        # if wmass:
        #     cardTool.addSystematic("sf2d", 
        #         processes=allMCprocesses_noQCDMC,
        #         outNames=["sf2dDown","sf2dUp"],
        #         group="SF3Dvs2D",
        #         scale = 1.0,
        #         mirror = True,
        #         mirrorDownVarEqualToNomi=False, # keep False, True is pathological
        #         noConstraint=False,
        #         systAxes=[],
        #         #labelsByAxis=["downUpVar"],
        #         passToFakes=passSystToFakes,
        #     )

        chargeDependentSteps = common.muonEfficiency_chargeDependentSteps
        effTypesNoIso = ["reco", "tracking", "idip", "trigger"]
        effStatTypes = [x for x in effTypesNoIso]
        if args.binnedScaleFactors or not args.isoEfficiencySmoothing:
            effStatTypes.extend(["iso"])
        else:
            effStatTypes.extend(["iso_effData", "iso_effMC"])
        allEffTnP = [f"effStatTnP_sf_{eff}" for eff in effStatTypes] + ["effSystTnP"]
        for name in allEffTnP:
            if "Syst" in name:
                axes = ["reco-tracking-idip-trigger-iso", "n_syst_variations"]
                axlabels = ["WPSYST", "_etaDecorr"]
                nameReplace = [("WPSYST0", "reco"), ("WPSYST1", "tracking"), ("WPSYST2", "idip"), ("WPSYST3", "trigger"), ("WPSYST4", "iso"), ("effSystTnP", "effSyst"), ("etaDecorr0", "fullyCorr") ]
                scale = 1.0
                mirror = True
                mirrorDownVarEqualToNomi=False
                groupName = "muon_eff_syst"
                splitGroupDict = {f"{groupName}_{x}" : f".*effSyst.*{x}" for x in list(effTypesNoIso + ["iso"])}
                splitGroupDict[groupName] = ".*effSyst.*" # add also the group with everything
                # decorrDictEff = {                        
                #     "x" : {
                #         "label" : "eta",
                #         "edges": [round(-2.4+i*0.1,1) for i in range(49)]
                #     }
                # }

            else:
                nameReplace = [] if any(x in name for x in chargeDependentSteps) else [("q0", "qall")] # for iso change the tag id with another sensible label
                mirror = True
                mirrorDownVarEqualToNomi=False
                if args.binnedScaleFactors:
                    axes = ["SF eta", "nPtBins", "SF charge"]
                else:
                    axes = ["SF eta", "nPtEigenBins", "SF charge"]
                axlabels = ["eta", "pt", "q"]
                nameReplace = nameReplace + [("effStatTnP_sf_", "effStat_")]           
                scale = 1.0
                groupName = "muon_eff_stat"
                splitGroupDict = {f"{groupName}_{x}" : f".*effStat.*{x}" for x in effStatTypes}
                splitGroupDict[groupName] = ".*effStat.*" # add also the group with everything
            if args.effStatLumiScale and "Syst" not in name:
                scale /= math.sqrt(args.effStatLumiScale)

            cardTool.addSystematic(
                name, 
                mirror=mirror,
                mirrorDownVarEqualToNomi=mirrorDownVarEqualToNomi,
                group=groupName,
                systAxes=axes,
                labelsByAxis=axlabels,
                baseName=name+"_",
                processes=['MCnoQCD'],
                passToFakes=passSystToFakes,
                systNameReplace=nameReplace,
                scale=scale,
                splitGroup=splitGroupDict,
                decorrelateByBin = {}
            )
            # if "Syst" in name and decorrDictEff != {}:
            #     # add fully correlated version again
            #     cardTool.addSystematic(
            #         name,
            #         rename=f"{name}_EtaDecorr",
            #         mirror=mirror,
            #         mirrorDownVarEqualToNomi=mirrorDownVarEqualToNomi,
            #         group=groupName,
            #         systAxes=axes,
            #         labelsByAxis=axlabels,
            #         baseName=name+"_",
            #         processes=allMCprocesses_noQCDMC,
            #         passToFakes=passSystToFakes,
            #         systNameReplace=nameReplace,
            #         scale=scale,
            #         splitGroup=splitGroupDict,
            #         decorrelateByBin = decorrDictEff
            #     )

    to_fakes = passSystToFakes and not args.noQCDscaleFakes and not xnorm

    theory_helper = combine_theory_helper.TheoryHelper(cardTool)
    theory_helper.configure(resumUnc=args.resumUnc, 
        propagate_to_fakes=to_fakes,
        np_model=args.npUnc,
        tnp_magnitude=args.tnpMagnitude,
        mirror_tnp=True,
        pdf_from_corr=args.pdfUncFromCorr,
        scale_pdf_unc=args.scalePdf,
    )
    theory_helper.add_all_theory_unc()

    msv_config_dict = {
        "smearingWeights":{
            "hist_name": "muonScaleSyst_responseWeights",
            "syst_axes": ["unc", "downUpVar"],
            "syst_axes_labels": ["unc", "downUpVar"]
        },
        "massWeights":{
            "hist_name": "muonScaleSyst",
            "syst_axes": ["downUpVar", "scaleEtaSlice"],
            "syst_axes_labels": ["downUpVar", "ieta"]
        },
        "manualShift":{
            "hist_name": "muonScaleSyst_manualShift",
            "syst_axes": ["downUpVar"],
            "syst_axes_labels": ["downUpVar"]
        }
    }

    if not lowPU:
        # FIXME: remove this once msv from smearing weights is implemented for the Z
        msv_config = msv_config_dict[args.muonScaleVariation] if wmass else msv_config_dict["massWeights"]

        cardTool.addSystematic(msv_config['hist_name'], 
            processes=['single_v_samples' if wmass else 'single_vmu_samples'],
            group="muonScale",
            baseName="CMS_scale_m_",
            systAxes=msv_config['syst_axes'],
            labelsByAxis=msv_config['syst_axes_labels'],
            passToFakes=passSystToFakes,
            scale = args.scaleMuonCorr
        )
        cardTool.addSystematic("muonL1PrefireSyst", 
            processes=['MCnoQCD'],
            group="muonPrefire",
            baseName="CMS_prefire_syst_m",
            systAxes=["downUpVar"],
            labelsByAxis=["downUpVar"],
            passToFakes=passSystToFakes,
        )
        cardTool.addSystematic("muonL1PrefireStat", 
            processes=['MCnoQCD'],
            group="muonPrefire",
            baseName="CMS_prefire_stat_m_",
            systAxes=["downUpVar", "etaPhiRegion"],
            labelsByAxis=["downUpVar", "etaPhiReg"],
            passToFakes=passSystToFakes,
        )
        cardTool.addSystematic("ecalL1Prefire", 
            processes=['MCnoQCD'],
            group="ecalPrefire",
            baseName="CMS_prefire_ecal",
            systAxes=["downUpVar"],
            labelsByAxis=["downUpVar"],
            passToFakes=passSystToFakes,
        )
        if wmass or wlike:
            combine_helpers.add_recoil_uncertainty(cardTool, ['signal_samples'], passSystToFakes=passSystToFakes, flavor="mu")

    if (wmass or wlike) and not input_tools.args_from_metadata(cardTool, "noRecoil"):
        combine_helpers.add_recoil_uncertainty(cardTool, signal_samples, passSystToFakes=passSystToFakes, flavor=datagroups.flavor if datagroups.flavor else "mu")

    if wmass:
        if not lowPU:
            non_closure_scheme = input_tools.args_from_metadata(cardTool, "nonClosureScheme")
            if non_closure_scheme == "A-M-separated":
                cardTool.addSystematic("Z_non_closure_parametrized_A", 
                    processes=['single_v_samples'],
                    group="nonClosure" if args.sepImpactForNC else "muonScale",
                    baseName="Z_nonClosure_parametrized_A_",
                    systAxes=["unc", "downUpVar"] if not (args.correlatedNonClosureNuisances) else ["downUpVar"],
                    labelsByAxis=["unc", "downUpVar"] if not (args.correlatedNonClosureNuisances) else ["downUpVar"],
                    passToFakes=passSystToFakes
                )
            if non_closure_scheme in ["A-M-separated", "binned-plus-M"]:
                cardTool.addSystematic("Z_non_closure_parametrized_M", 
                    processes=['single_v_samples'],
                    group="nonClosure" if args.sepImpactForNC else "muonScale",
                    baseName="Z_nonClosure_parametrized_M_",
                    systAxes=["unc", "downUpVar"] if not (args.correlatedNonClosureNuisances) else ["downUpVar"],
                    labelsByAxis=["unc", "downUpVar"] if not (args.correlatedNonClosureNuisances) else ["downUpVar"],
                    passToFakes=passSystToFakes
                )            
            if non_closure_scheme == "A-M-combined":
                cardTool.addSystematic("Z_non_closure_parametrized", 
                    processes=['single_v_samples'],
                    group="nonClosure" if args.sepImpactForNC else "muonScale",
                    baseName="Z_nonClosure_parametrized_",
                    systAxes=["unc", "downUpVar"] if not (args.correlatedNonClosureNuisances) else ["downUpVar"],
                    labelsByAxis=["unc", "downUpVar"] if not (args.correlatedNonClosureNuisances) else ["downUpVar"],
                    passToFakes=passSystToFakes
                )
            if non_closure_scheme in ["binned", "binned-plus-M"]:
                cardTool.addSystematic("Z_non_closure_binned", 
                    processes=['single_v_samples'],
                    group="nonClosure" if args.sepImpactForNC else "muonScale",
                    baseName="Z_nonClosure_binned_",
                    systAxes=["unc_ieta", "unc_ipt", "downUpVar"] if not (args.correlatedNonClosureNuisances) else ["downUpVar"],
                    labelsByAxis=["unc_ieta", "unc_ipt", "downUpVar"] if not (args.correlatedNonClosureNuisances) else ["downUpVar"],
                    passToFakes=passSystToFakes
                )
    
    # Previously we had a QCD uncertainty for the mt dependence on the fakes, see: https://github.com/WMass/WRemnants/blob/f757c2c8137a720403b64d4c83b5463a2b27e80f/scripts/combine/setupCombineWMass.py#L359

    return cardTool

def main(args,xnorm=False):
    cardTool = setup(args, xnorm)
    cardTool.setOutput(args.outfolder, fitvars=args.fitvar, doStatOnly=args.doStatOnly, postfix=args.postfix)
    cardTool.writeOutput(args=args, forceNonzero=not args.unfolding, check_systs=not args.unfolding, simultaneousABCD=args.simultaneousABCD)
    return

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()

    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
    
    time0 = time.time()

    if args.genModel:
        main(args, xnorm=True)
    else:
        main(args)
        if args.unfolding:
            main(args, xnorm=True)

    logger.info(f"Running time: {time.time()-time0}")
