
import itertools
import hist
import numpy as np

from narf import ioutils

from utilities import logging
from utilities.io_tools import combinetf_input

logger = logging.child_logger(__name__)

def fitresult_pois_to_hist(infile, poi_types = ["mu", "pmaskedexp", "pmaskedexpnorm", "sumpois", "sumpoisnorm", ],
    merge_channels=True, grouped=True, uncertainties=None, 
):
    # convert POIs in fitresult into histograms
    # uncertainties, use None to get all, use [] to get none
    # grouped=True for grouped uncertainties 
    # Different channels can have different year, flavor final state, particle final state, sqrt(s), 
    #   if merge_channels=True the lumi is added up for final states with different flavors or eras with same sqrt(s)
    channel_energy={
        "2017G": "5TeV",
        "2017H": "13TeV",
        "2016preVFP": "13TeV", 
        "2016postVFP":"13TeV", 
        "2017": "13TeV", 
        "2018": "13TeV",
    }
    channel_flavor ={
        "e": "l", 
        "mu": "l",
        "ee": "ll", 
        "mumu": "ll",
    }

    fitresult = combinetf_input.get_fitresult(infile)
    meta = ioutils.pickle_load_h5py(fitresult["meta"])
    meta_info = meta["meta_info"]
    result = {}
    for poi_type in poi_types:
        logger.debug(f"Now at POI type {poi_type}")

        scale = 1 
        if poi_type in ["nois"]:
            scale = 1./(imeta["args"]["scaleNormXsecHistYields"]*imeta["args"]["priorNormXsec"])

        df = combinetf_input.read_impacts_pois(fitresult, poi_type=poi_type, group=grouped, uncertainties=uncertainties)

        result[poi_type] = {}

        if merge_channels:
            channel_info = {}
            for chan, info in meta["channel_info"].items():
                if chan.endswith("masked"):
                    continue
                channel = f"chan_{channel_energy[info['era']]}"
                if info['flavor'] in channel_flavor:
                    channel += f"_{channel_flavor[info['flavor']]}"

                logger.debug(f"Merge channel {chan} into {channel}")

                lumi = info["lumi"]    
                gen_axes = info["gen_axes"]    

                if channel not in channel_info:
                    channel_info[channel] = {
                        "gen_axes": gen_axes,
                        "lumi": lumi,
                    }
                else:
                    if gen_axes != channel_info[channel]["gen_axes"]:
                        raise RuntimeError(f"The gen axes are different among channels {channel_info}, so they can't be merged")
                    channel_info[channel]["lumi"] += lumi
        else:
            channel_info = meta["channel_info"]

        for channel, info in channel_info.items():
            logger.debug(f"Now at channel {channel}")

            channel_scale = scale
            if poi_type in ["pmaskedexp", "sumpois"]:
                channel_scale = info["lumi"]

            result[poi_type][channel] = {}
            for proc, gen_axes_proc in info["gen_axes"].items():
                logger.debug(f"Now at proc {proc}")

                if poi_type.startswith("sum"):
                    # make all possible lower dimensional gen axes combinations; wmass only combinations including qGen
                    gen_axes_permutations = [list(k) for n in range(1, len(gen_axes_proc)) for k in itertools.combinations(gen_axes_proc, n)]
                else:
                    gen_axes_permutations = [gen_axes_proc[:],]

                result[poi_type][channel][proc] = {}
                for axes in gen_axes_permutations:
                    shape = [a.extent for a in axes]
                    axes_names = [a.name for a in axes]

                    data = combinetf_input.select_pois(df, axes_names, base_processes=proc, flow=True)

                    values = np.reshape(data["value"].values/channel_scale, shape)
                    variances = np.reshape( (data["err_total"].values/channel_scale)**2, shape)

                    h_ = hist.Hist(*axes, storage=hist.storage.Weight())
                    h_.view(flow=True)[...] = np.stack([values, variances], axis=-1)

                    hist_name = "hist_" + "_".join(axes_names)
                    logger.debug(f"Save histogram {hist_name}")
                    result[poi_type][channel][proc][hist_name] = h_

                    if "err_stat" in data.keys():
                        # save stat only hist
                        variances = np.reshape( (data["err_stat"].values/channel_scale)**2, shape)
                        h_stat = hist.Hist(*axes, storage=hist.storage.Weight())
                        h_stat.view(flow=True)[...] = np.stack([values, variances], axis=-1)
                        result[poi_type][channel][proc][f"{hist_name}_stat"] = h_stat

                    # save other systematic uncertainties as separately varied histograms
                    labels = [u.replace("err_","") for u in filter(lambda x: x.startswith("err_") and x not in ["err_total", "err_stat"], data.keys())]
                    if labels:
                        # string category always has an overflow bin, we set it to 0
                        systs = np.stack([values, *[values + np.reshape(data[f"err_{u}"].values/channel_scale, shape) for u in labels], np.zeros_like(values)], axis=-1)
                        h_syst = hist.Hist(*axes, hist.axis.StrCategory(["nominal", *labels], name="syst"), storage=hist.storage.Double())
                        h_syst.values(flow=True)[...] = systs
                        result[poi_type][channel][proc][f"{hist_name}_syst"] = h_syst

    return result, meta