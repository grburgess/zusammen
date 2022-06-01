import numpy as np

from .corr_cpl_grb import GBMGRB_CORR_CPL
from cosmogrb.universe.universe import GRBWrapper, ParameterServer, Universe


class GBM_CORR_CPL_Universe(Universe):
    """
    Documentation for GBM_CORR_CPL_Universe
    """

    def __init__(self, population, save_path="."):

        super(GBM_CORR_CPL_Universe, self).__init__(
            population, save_path=save_path
        )

    def _grb_wrapper(self, parameter_server, serial=False):
        return GBM_CORR_CPL_Wrapper(parameter_server, serial=serial)

    def _process_populations(self):

        # get the Ra and Dec
        super(GBM_CORR_CPL_Universe, self)._process_populations()

        self._local_parameters["ep_start"] = self._population.ep_obs
        self._local_parameters["alpha"] = self._population.alpha
        self._local_parameters["peak_flux"] = self._population.fluxes_latent
        self._local_parameters["ep_tau"] = self._population.ep_tau
        self._local_parameters["Nrest"] = self._population.nrest
        self._local_parameters["gamma"] = self._population.gamma

    def _parameter_server_type(self, **kwargs):

        return GBM_CORR_CPL_ParameterServer(**kwargs)


class GBM_CORR_CPL_Wrapper(GRBWrapper):
    """
    Documentation for GBM_CPL_Wrapper
    """

    def __init__(self, parameter_server, serial=False):
        super(GBM_CORR_CPL_Wrapper, self).__init__(
            parameter_server=parameter_server, serial=serial
        )

    def _grb_type(self, **kwargs):
        return GBMGRB_CORR_CPL(**kwargs)


class GBM_CORR_CPL_ParameterServer(ParameterServer):
    """
    Documentation for GBM_CPL_ParameterServer
    """

    def __init__(
        self,
        name,
        ra,
        dec,
        z,
        duration,
        T0,
        peak_flux,
        alpha,
        ep_start,
        ep_tau,
        Nrest,
        gamma,
    ):
        """FIXME! briefly describe function

        :param name:
        :param ra:
        :param dec:
        :param z:
        :param duration:
        :param T0:
        :param peak_flux:
        :param alpha:
        :param ep:
        :param ep_tau:
        :returns:
        :rtype:

        """

        super(GBM_CORR_CPL_ParameterServer, self).__init__(
            name=name,
            ra=ra,
            dec=dec,
            z=z,
            duration=duration,
            T0=T0,
            peak_flux=peak_flux,
            alpha=alpha,
            ep_start=ep_start,
            ep_tau=ep_tau,
            Nrest=Nrest,
            gamma=gamma,
        )
