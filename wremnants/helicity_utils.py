import ROOT
import pathlib
import hist
import narf.clingutils
import uproot
import pathlib
import hist
import pickle
import lz4.frame
from wremnants.correctionsTensor_helper import makeCorrectionsTensor
from wremnants.theory_tools import helicity_xsec_to_angular_coeffs
from utilities import common, logging
from utilities import boostHistHelpers as hh
from utilities.io_tools import input_tools
import numpy as np
import h5py
import hdf5plugin
import narf

logger = logging.child_logger(__name__)

narf.clingutils.Declare('#include "syst_helicity_utils.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

#UL, A0...A4
axis_helicity = hist.axis.Integer(-1, 8, name="helicity", overflow=False, underflow=False)
axis_helicity_multidim = hist.axis.Integer(-1, 8, name="helicitySig", overflow=False, underflow=False)

#creates the helicity weight tensor
def makehelicityWeightHelper(is_w_like = False, filename=None):
    if filename is None:
        filename = f"{common.data_dir}/angularCoefficients//w_z_moments_theoryAgnosticBinning.hdf5"
    with h5py.File(filename, "r") as ff:
        out = input_tools.load_results_h5py(ff)

    hist_helicity_xsec_scales = out["Z"] if is_w_like else out["W"]

    corrh = helicity_xsec_to_angular_coeffs(hist_helicity_xsec_scales)

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

    return makeCorrectionsTensor(corrh_noerrs, ROOT.wrem.WeightByHelicityHelper, tensor_rank=1)

def make_helper_helicity(axes, nhelicity=6):
    """
    Converts axes into tensor axes expanded by the helicity tensor
    """
    ndim = len(axes)
    shape = [a.size for a in axes]
    try:
        helper = ROOT.wrem.tensor_helper_helicity[ndim, nhelicity, *shape]()
    except Exception as e:
        logger.warning(f"An error occurred while trying to create a helicity tensor helper: {e}")
    tensor_axes = [axis_helicity_multidim, *axes]
    return helper, tensor_axes
