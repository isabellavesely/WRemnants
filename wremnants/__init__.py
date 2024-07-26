import ROOT
import narf
import pathlib

ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")

narf.clingutils.Declare('#include "muonCorr.h"')
narf.clingutils.Declare('#include "histoScaling.h"')
narf.clingutils.Declare('#include "histHelpers.h"')
narf.clingutils.Declare('#include "utils.h"')
narf.clingutils.Declare('#include "csVariables.h"')
narf.clingutils.Declare('#include "EtaPtCorrelatedEfficiency.h"')
narf.clingutils.Declare('#include "theoryTools.h"')
narf.clingutils.Declare('#include "syst_helicity_utils_polvar.h"')
