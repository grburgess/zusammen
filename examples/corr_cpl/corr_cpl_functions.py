import numpy as np
import numba as nb
from interpolation import interp

from cosmogrb.sampler.cpl_functions import cpl
from cosmogrb.utils.numba_array import VectorFloat64

from popsynth.distributions.cosmological_distribution import luminosity_distance


@nb.njit(fastmath=True, cache=False)
def ep_decay(t, t_start, ep_start, ep_tau):

    if t >= t_start:
        return ep_start * np.power((1 + t) / (1 + t_start), ep_tau)

    else:
        return ep_start


@nb.njit(fastmath=True, cache=False)
def golenetskii_corr(ep, Nrest, gamma, z):

    dl = luminosity_distance(z)  # cm

    logepr = np.log10((1 + z) * ep / 100)

    logF = np.log10(Nrest) - np.log10(4 * np.pi * dl * dl) + gamma * logepr

    return np.power(10, logF)


@nb.njit(fastmath=True, cache=False)
def corr_cpl_evolution(
    energy,
    time,
    peak_flux,
    ep_start,
    ep_tau,
    emin,
    emax,
    alpha,
    redshift,
    Nrest,
    gamma,
):

    N = time.shape[0]
    M = energy.shape[0]

    a = 10.0
    b = 1e4

    out = np.empty((N, M))

    for n in range(N):

        ep = ep_decay(time[n], t_start=0.0, ep_start=ep_start, ep_tau=ep_tau)

        F = golenetskii_corr(ep, Nrest=Nrest, gamma=gamma, z=redshift)

        for m in range(M):
            out[n, m] = cpl(energy[m], alpha=alpha, xp=ep, F=F, a=a, b=b)

    return out


@nb.njit(fastmath=True, cache=False)
def folded_cpl_evolution(
    energy,
    time,
    peak_flux,
    ep_start,
    ep_tau,
    emin,
    emax,
    alpha,
    redshift,
    Nrest,
    gamma,
    response,
):

    return interp(response[0], response[1], energy) * corr_cpl_evolution(
        energy,
        time,
        peak_flux,
        ep_start,
        ep_tau,
        emin,
        emax,
        alpha,
        redshift,
        Nrest,
        gamma,
    )


@nb.njit(fastmath=True, cache=False)
def energy_integrated_evolution(
    emin,
    emax,
    time,
    peak_flux,
    ep_start,
    ep_tau,
    alpha,
    redshift,
    Nrest,
    gamma,
    effective_area,
):

    energy_grid = np.power(10, np.linspace(np.log10(emin), np.log10(emax), 75))

    energy_slice = folded_cpl_evolution(
        energy_grid,
        time,
        peak_flux,
        ep_start,
        ep_tau,
        emin,
        emax,
        alpha,
        redshift,
        Nrest,
        gamma,
        effective_area,
    )

    return np.trapz(energy_slice[0, :], energy_grid)


@nb.njit(fastmath=True, cache=False)
def time_integrated_evolution(
    tmin,
    tmax,
    energy,
    peak_flux,
    ep_start,
    ep_tau,
    alpha,
    redshift,
    Nrest,
    gamma,
    emin,
    emax,
    effective_area,
):

    n_times = 50

    time_grid = np.linspace(tmin, tmax, n_times)

    time_slice = folded_cpl_evolution(
        energy,
        time_grid,
        peak_flux,
        ep_start,
        ep_tau,
        emin,
        emax,
        alpha,
        redshift,
        Nrest,
        gamma,
        effective_area,
    )

    return np.trapz(time_slice[:, 0], time_grid)


@nb.njit(fastmath=True, cache=False)
def sample_events(
    emin,
    emax,
    tstart,
    tstop,
    peak_flux,
    ep_start,
    ep_tau,
    alpha,
    redshift,
    Nrest,
    gamma,
    effective_area,
    fmax,
):

    time = tstart

    arrival_times = VectorFloat64(0)
    arrival_times.append(time)

    vtime = np.empty(1)

    while True:

        time = time - (1.0 / fmax) * np.log(np.random.rand())
        if time > tstop:
            break

        test = np.random.rand()

        vtime[0] = time

        p_test = (
            energy_integrated_evolution(
                emin,
                emax,
                vtime,
                peak_flux,
                ep_start,
                ep_tau,
                alpha,
                redshift,
                Nrest,
                gamma,
                effective_area,
            )
            / fmax
        )

        if test <= p_test:
            arrival_times.append(time)

    return arrival_times.arr


@nb.njit(fastmath=True, cache=False)
def sample_energy(
    times,
    peak_flux,
    ep_start,
    ep_tau,
    alpha,
    emin,
    emax,
    redshift,
    Nrest,
    gamma,
    effective_area,
):

    N = times.shape[0]

    egrid = np.power(10, np.linspace(np.log10(emin), np.log10(emax), 500))

    out = np.zeros(N)

    tmps = folded_cpl_evolution(
        egrid,
        times,
        peak_flux,
        ep_start,
        ep_tau,
        emin,
        emax,
        alpha,
        redshift,
        Nrest,
        gamma,
        effective_area,
    )

    x = np.empty(1)
    vtime = np.empty(1)
    for i in range(N):

        # the maximum is either at the lower bound or the max effective area

        tmp = tmps[i, :]

        idx = np.argmax(tmp)

        # bump up C just in case

        C = tmp[idx] * 5

        # so this scheme for dealing with the effective area
        # is likely very fragile. The idea is that power law
        # envelope needs to be greater than the function every where

        if alpha == -1.0:
            alpha = -1 + 1e-20

        while True:

            # sample from a power law
            u = np.random.uniform(0, 1)
            x[0] = np.power(
                (np.power(emax, alpha + 1) - np.power(emin, alpha + 1)) * u
                + np.power(emin, alpha + 1),
                1.0 / (alpha + 1.0),
            )

            y = np.random.uniform(0, 1) * C * np.power(x[0] / egrid[idx], alpha)

            # here the vtime is just to trick this into being an array

            vtime[0] = times[i]

            if (
                y
                <= (
                    folded_cpl_evolution(
                        x,
                        vtime,
                        peak_flux,
                        ep_start,
                        ep_tau,
                        emin,
                        emax,
                        alpha,
                        redshift,
                        Nrest,
                        gamma,
                        effective_area,
                    )
                )[0, 0]
            ):

                out[i] = x[0]
                break

    return out
