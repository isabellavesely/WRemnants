from utilities import boostHistHelpers as hh, common, logging
from utilities.io_tools import input_tools

import numpy as np
import hist
import re

logger = logging.child_logger(__name__)


def add_mass_diff_variations(cardTool, **kwargs):
    mass_diff_args = kwargs.copy()
    if args.fitMassDiff == "charge":
        cardTool.addSystematic(**mass_diff_args,
            # # on gen level based on the sample, only possible for mW
            # preOpMap={m.name: (lambda h, swap=swap_bins: swap(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown"))
            #     for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members if "minus" in m.name},
            # on reco level based on reco charge
            preOpMap={m.name: (lambda h:
                hh.swap_histogram_bins(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", "charge", 0)
                ) for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members},
        )
    elif args.fitMassDiff == "cosThetaStarll":
        cardTool.addSystematic(**mass_diff_args,
            preOpMap={m.name: (lambda h:
                hh.swap_histogram_bins(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", "cosThetaStarll", hist.tag.Slicer()[0:complex(0,0):])
                ) for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members},
        )
    elif args.fitMassDiff == "eta-sign":
        cardTool.addSystematic(**mass_diff_args,
            preOpMap={m.name: (lambda h:
                hh.swap_histogram_bins(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", "eta", hist.tag.Slicer()[0:complex(0,0):])
                ) for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members},
        )
    elif args.fitMassDiff == "eta-range":
        cardTool.addSystematic(**mass_diff_args,
            preOpMap={m.name: (lambda h:
                hh.swap_histogram_bins(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", "eta", hist.tag.Slicer()[complex(0,-0.9):complex(0,0.9):])
                ) for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members},
        )
    elif args.fitMassDiff.startswith("etaRegion"):
        # 3 bins, use 3 unconstrained parameters: mass; mass0 - mass2; mass0 + mass2 - mass1
        mass_diff_args["rename"] = f"massDiff1{suffix}{label}"
        mass_diff_args["systNameReplace"] = [("Shift",f"Diff1{suffix}")]
        cardTool.addSystematic(**mass_diff_args,
            preOpMap={m.name: (lambda h: hh.swap_histogram_bins(
                hh.swap_histogram_bins(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", args.fitMassDiff, 2), # invert for mass2
                "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", args.fitMassDiff, 1, axis1_replace=f"massShift{label}0MeV") # set mass1 to nominal
                ) for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members},
        )
        mass_diff_args["rename"] = f"massDiff2{suffix}{label}"
        mass_diff_args["systNameReplace"] = [("Shift",f"Diff2{suffix}")]
        cardTool.addSystematic(**mass_diff_args,
            preOpMap={m.name: (lambda h:
                hh.swap_histogram_bins(h, "massShift", f"massShift{label}50MeVUp", f"massShift{label}50MeVDown", args.fitMassDiff, 1)
                ) for g in cardTool.procGroups[signal_samples_forMass[0]] for m in cardTool.datagroups.groups[g].members},
        )


def add_recoil_uncertainty(card_tool, samples, passSystToFakes=False, pu_type="highPU", flavor="", group_compact=True):
    met = input_tools.args_from_metadata(card_tool, "met")
    if flavor == "":
        flavor = input_tools.args_from_metadata(card_tool, "flavor")
    if pu_type == "highPU" and (met in ["RawPFMET", "DeepMETReso", "DeepMETPVRobust", "DeepMETPVRobustNoPUPPI"]):

        '''
        card_tool.addSystematic("recoil_syst",
            processes=samples,
            mirror = True,
            group = "recoil" if group_compact else "recoil_syst",
            splitGroup={"experiment": f".*"},
            systAxes = ["recoil_unc"],
            passToFakes=passSystToFakes,
        )
        '''

        card_tool.addSystematic("recoil_stat",
            processes=samples,
            mirror = True,
            group = "recoil" if group_compact else "recoil_stat",
            splitGroup={"experiment": f".*"},
            systAxes = ["recoil_unc"],
            passToFakes=passSystToFakes,
        )


def add_explicit_MCstat(cardTool, recovar, samples='signal_samples', wmass=False, source=None):
    """
    add explicit bin by bin stat uncertainties 
    Parameters:
    source (tuple of str): take variations from histogram with name given by f"{source[0]}_{source[1]}" (E.g. used to correlate between masked channels). 
        If None, use variations from nominal histogram 
    """
    datagroups = cardTool.datagroups

    recovar_syst = [f"_{n}" for n in recovar]
    info=dict(
        baseName="MCstat_"+"_".join(cardTool.procGroups[samples])+"_",
        rename="statMC",
        group=f"statMC",
        passToFakes=False,
        processes=[samples],
        mirror=True,
        labelsByAxis=[f"_{p}" if p != recovar[0] else p for p in recovar],
    )
    cardTool.setProcsNoStatUnc(cardTool.procGroups[samples])    
    if source is not None:
        # signal region selection
        if wmass:
            action_sel = lambda h, x: histselections.SignalSelectorABCD(h[x]).get_hist(h[x])  
        else:
            action_sel = lambda h, x: h[x]

        integration_var = {a:hist.sum for a in datagroups.gen_axes_names} # integrate out gen axes for bin by bin uncertainties
        cardTool.addSystematic(**info,
            nominalName=source[0],
            name=source[1],
            systAxes=recovar,
            actionRequiresNomi=True,
            action=lambda hv, hn:
                hh.addHists(
                    hn[{"count":hist.sum, "acceptance":hist.sum}].project(*datagroups.gen_axes_names),
                    action_sel(hv, {"acceptance":True}).project(*recovar, *datagroups.gen_axes_names),
                    scale2=(
                        np.sqrt(action_sel(hv, {"acceptance":hist.sum, **integration_var}).variances(flow=True))
                        / action_sel(hv, {"acceptance":hist.sum, **integration_var}).values(flow=True)
                    )[...,*[np.newaxis] * len(datagroups.gen_axes_names)]
                )
        )
    else:
        if args.fitresult:
            info["group"] = "binByBinStat"
        cardTool.addSystematic(**info,
            name=cardTool.nominalName,
            systAxes=recovar_syst,
            action=lambda h: 
                hh.addHists(h.project(*recovar),
                    hh.expand_hist_by_duplicate_axes(h.project(*recovar), recovar, recovar_syst),
                    scale2=np.sqrt(h.project(*recovar).variances(flow=True))/h.project(*recovar).values(flow=True))
        )


def add_electroweak_uncertainty(card_tool, ewUncs, flavor="mu", samples="single_v_samples", passSystToFakes=True, wlike=False):
    info = dict(
        systAxes=["systIdx"],
        mirror=True,
        splitGroup={"theory_ew" : f".*", "theory" : f".*"},
        passToFakes=passSystToFakes,
    )
    # different uncertainty for W and Z samples
    all_samples = card_tool.procGroups[samples]
    z_samples = [p for p in all_samples if p[0]=="Z"]
    w_samples = [p for p in all_samples if p[0]=="W"]
    mode = card_tool.datagroups.mode
    
    for ewUnc in ewUncs:
        if "renesanceEW" in ewUnc:
            pass
            if w_samples:
                # add renesance (virtual EW) uncertainty on W samples
                card_tool.addSystematic(f"{ewUnc}Corr",
                    processes=w_samples,
                    preOp = lambda h : h[{"var": ["nlo_ew_virtual"]}],
                    labelsByAxis=[f"renesanceEWCorr"],
                    scale=1.,
                    systAxes=["var"],
                    group = "theory_ew_virtW_corr",
                    splitGroup={"theory_ew" : f".*", "theory" : f".*"},
                    passToFakes=passSystToFakes,
                    mirror = True,
                )
        elif ewUnc == "powhegFOEW":
            if z_samples:
                card_tool.addSystematic(f"{ewUnc}Corr",
                    preOp = lambda h : h[{"weak": ["weak_ps", "weak_aem"]}],
                    processes=z_samples,
                    labelsByAxis=[f"{ewUnc}Corr"],
                    scale=1.,
                    systAxes=["weak"],
                    mirror=True,
                    group="theory_ew_virtZ_scheme",
                    splitGroup={"theory_ew" : f".*", "theory" : f".*"},
                    passToFakes=passSystToFakes,
                    rename = "ewScheme",
                )
                card_tool.addSystematic(f"{ewUnc}Corr",
                    preOp = lambda h : h[{"weak": ["weak_default"]}],
                    processes=z_samples,
                    labelsByAxis=[f"{ewUnc}Corr"],
                    scale=1.,
                    systAxes=["weak"],
                    mirror=True,
                    group="theory_ew_virtZ_corr",
                    splitGroup={"theory_ew" : f".*", "theory" : f".*"},
                    passToFakes=passSystToFakes,
                    rename = "ew",
                )
        else:
            if "FSR" in ewUnc:
                if flavor == "e":
                    logger.warning("ISR/FSR EW uncertainties are not implemented for electrons, proceed w/o")
                    continue
                scale=1
            if "ISR" in ewUnc:
                scale=2
            else:
                scale=1

            if "winhac" in ewUnc:
                if not w_samples:
                    logger.warning("Winhac is not implemented for any other process than W, proceed w/o winhac EW uncertainty")
                    continue
                elif all_samples != w_samples:
                    logger.warning("Winhac is only implemented for W samples, proceed w/o winhac EW uncertainty for other samples")
                samples = w_samples
            else:
                samples = all_samples

            s = hist.tag.Slicer()
            if ewUnc.startswith("virtual_ew"):
                preOp = lambda h : h[{"systIdx" : s[0:1]}]
            else:
                preOp = lambda h : h[{"systIdx" : s[1:2]}]

            card_tool.addSystematic(f"{ewUnc}Corr", **info,
                processes=samples,
                labelsByAxis=[f"{ewUnc}Corr"],
                scale=scale,
                preOp = preOp,
                group = f"theory_ew_{ewUnc}",
            )  

def projectABCD(cardTool, h, return_variances=False, dtype="float64"):
    # in case the desired axes are different at low MT and high MT we need to project each seperately, and then concatenate

    if any(ax not in h.axes.name for ax in cardTool.getFakerateAxes()):
        logger.warning(f"Not all desired fakerate axes found in histogram. Fakerate axes are {cardTool.getFakerateAxes()}, and histogram axes are {h.axes.name}")

    fakerate_axes = [n for n in h.axes.name if n in cardTool.getFakerateAxes()]

    lowMT_axes = [n for n in h.axes.name if n in fakerate_axes]
    highMT_failIso_axes = [n for n in h.axes.name if n in [*fakerate_axes, *cardTool.fit_axes]]
    highMT_passIso_axes = [n for n in h.axes.name if n in cardTool.fit_axes]

    hist_lowMT = h[{cardTool.nameMT : cardTool.failMT}].project(*[*lowMT_axes, common.passIsoName])
    hist_highMT_failIso = h[{cardTool.nameMT : cardTool.passMT, **common.failIso}].project(*[*highMT_failIso_axes])
    hist_highMT_passIso = h[{cardTool.nameMT : cardTool.passMT, **common.passIso}].project(*[*highMT_passIso_axes])

    flat_lowMT = hist_lowMT.values(flow=False).flatten().astype(dtype)
    flat_highMT_failIso = hist_highMT_failIso.values(flow=False).flatten().astype(dtype)
    flat_highMT_passIso = hist_highMT_passIso.values(flow=False).flatten().astype(dtype)

    flat = np.append(flat_lowMT, flat_highMT_failIso)
    flat = np.append(flat, flat_highMT_passIso)

    if not return_variances:
        return flat

    flat_variances_lowMT = hist_lowMT.variances(flow=False).flatten().astype(dtype)
    flat_variances_highMT_failIso = hist_highMT_failIso.variances(flow=False).flatten().astype(dtype)
    flat_variances_highMT_passIso = hist_highMT_passIso.variances(flow=False).flatten().astype(dtype)

    flat_variances = np.append(flat_variances_lowMT, flat_variances_highMT_failIso)
    flat_variances = np.append(flat_variances, flat_variances_highMT_passIso)

    return flat, flat_variances


