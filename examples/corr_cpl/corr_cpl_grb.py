from cosmogrb.instruments.gbm import GBMGRB
from cosmogrb.grb import SourceParameter

from .corr_cpl_source import CorrCPLSourceFunction


class GBMGRB_CORR_CPL(GBMGRB):
    peak_flux = SourceParameter()
    alpha = SourceParameter()
    ep_start = SourceParameter()
    ep_tau = SourceParameter()
    Nrest = SourceParameter()
    gamma = SourceParameter()

    def __init__(self, **kwargs):

        # pass up
        super(GBMGRB_CORR_CPL, self).__init__(
            source_function_class=CorrCPLSourceFunction, **kwargs,
        )
