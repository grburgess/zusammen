import numpy as np

import popsynth


class TDecaySampler(popsynth.AuxiliarySampler):

    sigma = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1)

    def __init__(self):
        """
        samples the decay of the of the pulse
        """

        super(TDecaySampler, self).__init__(name="tdecay", observed=False)

    def true_sampler(self, size):

        t90 = 10 ** self._secondary_samplers["log_t90"].true_values
        trise = self._secondary_samplers["trise"].true_values

        self._true_values = (
            1.0 / 50.0 * (10 * t90 + trise + np.sqrt(trise) * np.sqrt(20 * t90 + trise))
        )


class DurationSampler(popsynth.AuxiliarySampler):

    sigma = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1)

    def __init__(self):
        """
        samples how long the pulse lasts
        """

        super(DurationSampler, self).__init__(name="duration", observed=False)

    def true_sampler(self, size):

        t90 = 10 ** self._secondary_samplers["log_t90"].true_values

        self._true_values = 1.5 * t90


class EpeakObsSampler(popsynth.AuxiliarySampler):
    def __init__(self):
        """
        Samples Epeak in the observed frame
        """

        super(EpeakObsSampler, self).__init__(
            name="log_ep_obs", observed=False, uses_distance=True
        )

    def true_sampler(self, size):

        secondary = self._secondary_samplers["log_ep"]

        self._true_values = secondary.true_values - np.log10(1 + self._distance)


class LumSampler(popsynth.DerivedLumAuxSampler):
    """
    Sample luminosity from Epeak
    """

    s_scat = popsynth.auxiliary_sampler.AuxiliaryParameter(default=0.3)

    def __init__(self):

        super(LumSampler, self).__init__(name="obs_lum")

    def compute_luminosity(self):

        return self._true_values

    def true_sampler(self, size):

        log_ep = self._secondary_samplers["log_ep"].true_values
        log_nrest = self._secondary_samplers["log_nrest"].true_values
        gamma = self._secondary_samplers["gamma"].true_values

        ep = np.power(10, log_ep)  # keV

        lum = np.power(10, log_nrest) * np.power(ep / 100, gamma)  # erg s^-1

        tmp = np.random.normal(0, self.s_scat * lum)

        self._true_values = lum + tmp


class DerivedEpeakSampler(popsynth.AuxiliarySampler):
    """
    Samples Epeak for a given L - probably not the way to go
    """

    Nrest = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1e52)
    gamma = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1.5, vmin=0)
    s_scat = popsynth.auxiliary_sampler.AuxiliaryParameter(default=0.1)

    s_det = popsynth.auxiliary_sampler.AuxiliaryParameter(default=0.1)

    def __init__(self):

        super(DerivedEpeakSampler, self).__init__(
            "derived_Epeak", observed=True, uses_luminosity=True, uses_distance=True
        )

    def true_sampler(self, size):

        index = (np.log10(self._luminosity) - np.log10(self.Nrest)) / self.gamma
        Ep = np.power(10, index) * 100  # keV

        s = np.random.normal(0, self.s_scat * Ep, size)

        self._true_values = Ep + s

    def observation_sampler(self, size):

        Ep_obs = self._true_values / (1 + self._distance)

        s = np.random.normal(0, self.s_det * self._true_values, size)

        self._obs_values = Ep_obs + s
