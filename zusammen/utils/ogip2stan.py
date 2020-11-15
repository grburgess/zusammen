import collections
import os
from glob import glob

import h5py
import numpy as np
import yaml
from astropy.cosmology import WMAP9 as cosmo

from threeML import OGIPLike


def sanitize_filename(filename, abspath=False):

    sanitized = os.path.expandvars(os.path.expanduser(filename))

    if abspath:

        return os.path.abspath(sanitized)

    else:

        return sanitized


class GRBDatum(object):
    def __init__(
        self,
        name,
        observation,
        background,
        background_error,
        response,
        mc_energies,
        ebounds,
        exposure,
        mask,
        significance,
    ):

        self._n_chans = len(observation)

        assert self._n_chans == len(background)
        assert self._n_chans == response.shape[0]
        assert self._n_chans == len(background_error)
        assert self._n_chans == len(mask)

        assert sum(mask) > 0

        assert exposure > 0

        self._name = name
        self._observation = np.array(observation).astype(int)
        self._background = np.array(background)
        self._background_error = background_error
        self._response = response

        self._response = response
        self._ebounds = ebounds
        self._mc_energies = mc_energies

        self._exposure = exposure
        self._mask = mask
        self._significance = significance

    @property
    def name(self):
        return self._name

    @property
    def response(self):
        return self._response

    @property
    def response_transpose(self):
        return self._response.T

    @property
    def observation(self):
        return self._observation

    @property
    def background(self):
        return self._background

    @property
    def background_error(self):
        return self._background_error

    @property
    def mask(self):
        return self._mask

    @property
    def mask_stan(self):
        return np.where(self._mask)[0] + 1

    @property
    def n_channels_used(self):
        return sum(self._mask)  # CHECK THIS MOTHER FUCKER

    @property
    def idx_background_zero(self):
        return np.where(self._background[self._mask] == 0)[0]

    @property
    def idx_background_nonzero(self):
        return np.where(self._background[self._mask] > 0)[0]

    @property
    def n_bkg_zero(self):
        return sum(self._background[self._mask] == 0)

    @property
    def n_bkg_nonzero(self):
        return sum(self._background[self._mask] > 0)

    @property
    def n_chans(self):
        return self._n_chans

    @property
    def n_echans(self):
        return self._response.shape[1]

    @property
    def exposure(self):
        return self._exposure

    @property
    def significance(self):
        return self._significance

    @property
    def ebounds(self):
        return self._mc_energies

    @property
    def ebounds_lo(self):
        return self._mc_energies[:-1]

    @property
    def ebounds_hi(self):
        return self._mc_energies[1:]

    @property
    def cbounds(self):
        return self._ebounds

    @property
    def cbounds_lo(self):
        return self._ebounds[:-1]

    @property
    def cbounds_hi(self):
        return self._ebounds[1:]

    @classmethod
    def from_ogip(cls, name, obs_file, bkg_file, rsp, selection, spectrum_number=1):
        """
        Create the base data from FITS files.

        :param cls:
        :param name:
        :param obs_file:
        :param bkg_file:
        :param rsp:
        :param selection:
        :param spectrum_number:
        :returns:
        :rtype:

        """

        # use 3ml to read the
        # FITS files
        ogip = OGIPLike(
            name,
            observation=obs_file,
            background=bkg_file,
            response=rsp,
            spectrum_number=spectrum_number,
            verbose=False,
        )

        # set the mask

        ogip.set_active_measurements(selection)

        return cls(
            name,
            ogip.observed_counts,
            ogip.background_counts,
            ogip.background_count_errors,
            ogip.response.matrix,
            ogip.response.monte_carlo_energies,
            ogip.response.ebounds,
            ogip.exposure,
            ogip.mask,
            ogip.significance,
        )

    def to_hdf5_file_or_group(self, name):
        """
        write the data to and HDF5 file or group

        :param name:
        :returns:
        :rtype:

        """

        if isinstance(name, h5py.File) or isinstance(name, h5py.Group):

            is_file = False
            f = name

        else:

            if_file = True
            f = h5py.File(name, "w")

        f.attrs["name"] = self._name
        f.attrs["exposure"] = self._exposure
        f.attrs["significance"] = self._significance

        f.create_dataset(
            "observation", data=self._observation, compression="lzf", shuffle=True
        )
        f.create_dataset(
            "background", data=self._background, compression="lzf", shuffle=True
        )
        f.create_dataset(
            "background_error",
            data=self._background_error,
            compression="lzf",
            shuffle=True,
        )
        f.create_dataset(
            "response", data=self._response, compression="lzf", shuffle=True
        )
        f.create_dataset(
            "mc_energies", data=self._mc_energies, compression="lzf", shuffle=True
        )
        f.create_dataset("ebounds", data=self._ebounds, compression="lzf", shuffle=True)
        f.create_dataset("mask", data=self._mask, compression="lzf", shuffle=True)

        if is_file:

            f.close()

    @classmethod
    def from_hdf5_file_or_group(cls, name):
        """
        read in the data from an HDF5 file

        :param cls:
        :param name:
        :returns:
        :rtype:

        """

        # check if we are from a bigger file or opening one
        if isinstance(name, h5py.File) or isinstance(name, h5py.Group):

            # keep track if we will need to close
            is_file = False
            f = name

        else:

            if_file = True
            f = h5py.File(name, "r")

        # extract all the shit

        name = f.attrs["name"]
        exposure = f.attrs["exposure"]
        significance = f.attrs["significance"]

        observation = f["observation"][()]
        background = f["background"][()]
        background_error = f["background_error"][()]
        response = f["response"][()]
        mc_energies = f["mc_energies"][()]
        ebounds = f["ebounds"][()]
        mask = f["mask"][()]

        if is_file:
            f.close()

        return cls(
            name,
            observation,
            background,
            background_error,
            response,
            mc_energies,
            ebounds,
            exposure,
            mask,
            significance,
        )


class GRBInterval(object):
    def __init__(self, grb_name, *data):
        """
        A time interval consisting of all the detectors

        :param name:
        :param observations:
        :param backgrounds:
        :param responses:
        :returns:
        :rtype:

        """

        self._data = collections.OrderedDict()
        self._n_dets = len(data)

        n_chans = []
        n_echans = []

        for datum in data:

            assert isinstance(datum, GRBDatum)

            self._data[datum.name] = datum
            n_echans.append(datum.n_echans)
            n_chans.append(datum.n_chans)

        self._name = grb_name

        self._max_n_chans = max(n_chans)
        self._max_n_echans = max(n_echans)

    @property
    def name(self):
        return self._name

    @property
    def data(self):

        return self._data

    @property
    def n_detectors(self):

        return self._n_dets

    @property
    def max_n_echans(self):
        return self._max_n_echans

    @property
    def max_n_chans(self):
        return self._max_n_chans

    @classmethod
    def from_dict(cls, d, grb_name, spectrum_number=1):
        """
        create from a dictionary that is initially
        from the yaml file and will trigger the reading
        of the OGIP files

        :param cls:
        :param d:
        :param spectrum_number:
        :returns:
        :rtype:

        """

        dets = d["detectors"]

        directory = sanitize_filename(d["dir"])

        data = []

        # we always want the sorted
        sorted_dets = list(dets.keys())
        sorted_dets.sort()

        # this depends on the file
        # structure. This is not great
        # but not bad

        for det in sorted_dets:

            # match pha

            phas = glob(os.path.join(directory, f"*{det}*.pha"))

            obs_file = [f for f in phas if "bak" not in f][0]

            # match bak
            bak_file = glob(os.path.join(directory, f"*{det}*bak.pha"))[0]

            # match rsp
            rsp_file = glob(os.path.join(directory, f"*{det}*.rsp"))[0]

            datum = GRBDatum.from_ogip(
                det, obs_file, bak_file, rsp_file, dets[det], spectrum_number
            )

            data.append(datum)
        return cls(grb_name, *data)

    def to_hdf5_file_or_group(self, name):
        """
        write the interval to HDF5 which triggers the recursive
        writers

        :param name:
        :returns:
        :rtype:

        """

        if isinstance(name, h5py.File) or isinstance(name, h5py.Group):

            is_file = False
            f = name

        else:

            if_file = True
            f = h5py.File(name, "w")

        for k, v in self._data.items():

            grp = f.create_group(k)
            grp.attrs["grb_name"] = self._name

            v.to_hdf5_file_or_group(name=grp)

        if is_file:

            f.close()

    @classmethod
    def from_hdf5_file_or_group(cls, name):
        """FIXME! briefly describe function

        :param cls:
        :param name:
        :returns:
        :rtype:

        """

        if isinstance(name, h5py.File) or isinstance(name, h5py.Group):

            is_file = False
            f = name

        else:

            if_file = True
            f = h5py.File(name, "r")

        data = []

        grb_name = f.attrs["grb_name"]

        for det in f.keys():

            datum = GRBDatum.from_hdf5_file_or_group(f[det])

            data.append(datum)

        if is_file:

            f.close()

        return cls(grb_name, *data)


class GRBData(object):
    def __init__(self, grb_name, *intervals, z=None):
        """
        Contains all the intervals from a single
        GRB

        :param grb_name:
        :returns:
        :rtype:

        """

        self._intervals = []

        n_chans = []
        n_echans = []
        self._n_detectors = []

        for interval in intervals:
            assert isinstance(interval, GRBInterval)
            assert grb_name == interval.name

            n_chans.append(interval.max_n_chans)
            n_echans.append(interval.max_n_echans)
            self._n_detectors.append(interval.n_detectors)
            self._intervals.append(interval)

        self._n_intervals = len(intervals)
        self._name = grb_name

        self._max_n_chans = max(n_chans)
        self._max_n_echans = max(n_echans)

        self._z = z

    @property
    def intervals(self):
        return self._intervals

    @property
    def n_intervals(self):
        return self._n_intervals

    @property
    def n_detectors(self):
        return self._n_detectors

    @property
    def max_n_detectors(self):
        return max(self._n_detectors)

    @property
    def name(self):
        return self._name

    @property
    def max_n_echans(self):
        return self._max_n_echans

    @property
    def max_n_chans(self):
        return self._max_n_chans

    @property
    def z(self):
        return self._z

    @property
    def luminosity_distance(self):
        return cosmo.luminosity_distance(self._z).to("cm").value

    @classmethod
    def from_dict(cls, grb_name, d):
        """

        construct from dictionary via
        recursively loading from FITS
        files

        :param cls:
        :param grb_name:
        :param d:
        :returns:
        :rtype:

        """

        intervals = []

        n_intervals = d["n_intervals"]
        z = d["z"]
        for i in range(n_intervals):

            interval = GRBInterval.from_dict(d, grb_name, spectrum_number=i + 1)

            intervals.append(interval)

        return cls(grb_name, *intervals, z=z)

    @classmethod
    def from_hdf5_file_or_group(cls, name):
        """

        construct from an HDF5 file or group

        :param cls:
        :param name:
        :returns:
        :rtype:

        """

        if isinstance(name, h5py.File) or isinstance(name, h5py.Group):

            is_file = False
            f = name

        else:

            if_file = True
            f = h5py.File(name, "r")

        intervals = []

        grb_name = f.attrs["grb_name"]

        n_intervals = f.attrs["n_intervals"]

        z = f.attrs["z"]

        for i in range(n_intervals):

            interval = GRBInterval.from_hdf5_file_or_group(f[f"interval_{i}"])

            intervals.append(interval)

        if is_file:

            f.close()

        return cls(grb_name, *intervals, z=z)

    def to_hdf5_file_or_group(self, name):
        """
        write this GRB to a file or group

        :param name:
        :returns:
        :rtype:

        """

        # check if this is the endpoint, i.e.,
        # this is a file or if we are writing to a group

        if isinstance(name, h5py.File) or isinstance(name, h5py.Group):

            is_file = False
            f = name

        else:

            if_file = True
            f = h5py.File(name, "w")

        f.attrs["z"] = self._z
        f.attrs["n_intervals"] = self._n_intervals

        for i in range(self._n_intervals):

            grp = f.create_group(f"interval_{i}")
            grp.attrs["grb_name"] = self._name

            # this will recurse down the intervals

            self._intervals[i].to_hdf5_file_or_group(grp)

        if is_file:

            f.close()


class DataSet(object):
    def __init__(self, *grbs):
        """
        A data set of GRBs with sub intervals. Serves as an
        interface between FITS/HDF5/Stan to keep things organized

        :returns:
        :rtype:

        """

        self._grbs = collections.OrderedDict()
        self._n_intervals = 0

        n_echans = []
        n_chans = []
        n_dets = []

        self._grb_id = {}

        for i, grb in enumerate(grbs):

            assert isinstance(grb, GRBData)

            self._grbs[grb.name] = grb
            self._n_intervals += grb.n_intervals
            n_echans.append(grb.max_n_echans)
            n_chans.append(grb.max_n_chans)
            n_dets.append(grb.max_n_detectors)
            # tag for stan
            self._grb_id[grb.name] = i + 1

        self._n_grbs = len(grbs)

        self._max_n_chans = max(n_chans)
        self._max_n_echans = max(n_echans)
        self._max_n_detectors = max(n_dets)

    @classmethod
    def from_dict(cls, d):
        """
        construct from a dictionary with the layout
        specified in the from_yaml function

        :param cls:
        :param d:
        :returns:
        :rtype:

        """

        # create and empty list

        grbs = []

        for grb_name, d2 in d.items():

            # creat a GRB data with the dict loader
            # which will recurse and build things
            # from fits files

            grb = GRBData.from_dict(grb_name, d2)

            grbs.append(grb)

        return cls(*grbs)

    @classmethod
    def from_yaml(cls, file_name):
        """

        Construct from a yaml file that specifies
        where the PHA/FITS files are. This is to
        initially build the data sets and then dump
        them to HDF5. The yaml file should look like

        grb_name:
           n_intervals: 4
           dir: ~/home/location
           z=1
           detectors:
             n1:
               - 10-900
             b0:
               - 250-30000


        :param cls:
        :param file_name:
        :returns:
        :rtype:

        """

        # open the file

        with open(file_name, "r") as f:

            # read the shit in from the yaml file

            d = yaml.load(f, Loader=yaml.SafeLoader)

            # call the dict construction method

        return cls.from_dict(d)

    @classmethod
    def from_hdf5_file(cls, file_name):
        """
        Read the data set from an HDF5 file

        :param cls:
        :param file_name:
        :returns:
        :rtype:

        """

        # create an empty list
        grbs = []

        # open the file

        with h5py.File(file_name, "r") as f:

            # see the names of the GRBs

            grb_names = f.keys()

            # for each grb name there will
            # be a group that should have all
            # the recursive info

            for name in grb_names:

                # create a GRB data which will recurse down the
                # file and gather the interval info

                grb = GRBData.from_hdf5_file_or_group(f[name])

                grbs.append(grb)

        return cls(*grbs)

    def to_hdf5_file(self, file_name):
        """
        Write the entire GRB data set to an HDF5 file

        :param file_name:
        :returns:
        :rtype:

        """

        # open the file

        with h5py.File(file_name, "w") as f:

            # for each grb

            for grb_name, grb in self._grbs.items():

                # create a group names have the
                # the GRB

                grp = f.create_group(grb_name)

                # save the grb name to the grp
                # for fidelity

                grp.attrs["grb_name"] = grb_name

                # call the GRB's own HDF5 creator
                # which will recursively write all
                # the information to the file

                grb.to_hdf5_file_or_group(grp)

    def to_stan_dict(self):

        stan_data = collections.OrderedDict()

        stan_data["N_intervals"] = int(self._n_intervals)
        stan_data["N_grbs"] = self._n_grbs

        stan_data["max_n_echan"] = self._max_n_echans
        stan_data["max_n_chan"] = self._max_n_chans

        observed_counts = np.zeros(
            (self._n_intervals, self._max_n_detectors, self._max_n_chans)
        )
        background_counts = np.zeros(
            (self._n_intervals, self._max_n_detectors, self._max_n_chans)
        )
        background_errors = np.zeros(
            (self._n_intervals, self._max_n_detectors, self._max_n_chans)
        )

        idx_background_zero = np.zeros(
            (self._n_intervals, self._max_n_detectors, self._max_n_chans)
        )
        idx_background_nonzero = np.zeros(
            (self._n_intervals, self._max_n_detectors, self._max_n_chans)
        )
        n_bkg_zero = np.zeros((self._n_intervals, self._max_n_detectors))
        n_bkg_nonzero = np.zeros((self._n_intervals, self._max_n_detectors))

        responses = np.zeros(
            (
                self._n_intervals,
                self._max_n_detectors,
                
                self._max_n_chans,
                self._max_n_echans
            )
        )

        exposures = np.zeros((self._n_intervals, self._max_n_detectors))
        n_echan = np.zeros((self._n_intervals, self._max_n_detectors))
        n_chan = np.zeros((self._n_intervals, self._max_n_detectors))

        masks = np.zeros((self._n_intervals, self._max_n_detectors, self._max_n_chans))
        n_channels_used = np.zeros((self._n_intervals, self._max_n_detectors))
        grb_id = np.zeros(self._n_intervals)
        ebounds_lo = np.zeros(
            (self._n_intervals, self._max_n_detectors, self._max_n_echans)
        )
        ebounds_hi = np.zeros(
            (self._n_intervals, self._max_n_detectors, self._max_n_echans)
        )

        cbounds_lo = np.zeros(
            (self._n_intervals, self._max_n_detectors, self._max_n_chans)
        )
        cbounds_hi = np.zeros(
            (self._n_intervals, self._max_n_detectors, self._max_n_chans)
        )

        grb_id = []

        n_dets = []

        z = []
        dl = []

        total_number_of_channels_used = 0

        i = 0

        for name, grb in self._grbs.items():

            for interval in grb.intervals:

                z.append(grb.z)

                dl.append(grb.luminosity_distance)

                n_dets.append(interval.n_detectors)

                grb_id.append(self._grb_id[grb.name])

                for j, (det, datum) in enumerate(interval.data.items()):

                    observed_counts[i, j, : datum.n_chans] = datum.observation
                    background_counts[i, j, : datum.n_chans] = datum.background

                    idx_background_zero[i, j, : datum.n_bkg_zero] = (
                        datum.idx_background_zero + 1
                    )
                    idx_background_nonzero[i, j, : datum.n_bkg_nonzero] = (
                        datum.idx_background_nonzero + 1
                    )

                    n_bkg_zero[i, j] = datum.n_bkg_zero
                    n_bkg_nonzero[i, j] = datum.n_bkg_nonzero

                    background_errors[i, j, : datum.n_chans] = datum.background_error

                    # responses[
                    #     i, j, : datum.n_echans, : datum.n_chans
                    # ] = datum.response_transpose

                    responses[
                        i, j, :datum.n_chans, :  datum.n_echans
                    ] = datum.response


                    
                    this_mask = datum.mask_stan

                    masks[i, j, : len(this_mask)] = this_mask

                    # this could be a bug!
                    n_channels_used[i, j] = datum.n_channels_used

                    n_chan[i, j] = datum.n_chans
                    n_echan[i, j] = datum.n_echans

                    ebounds_lo[i, j, : datum.n_echans] = datum.ebounds_lo
                    ebounds_hi[i, j, : datum.n_echans] = datum.ebounds_hi

                    cbounds_lo[i, j, : datum.n_chans] = datum.cbounds_lo
                    cbounds_hi[i, j, : datum.n_chans] = datum.cbounds_hi

                    exposures[i, j] = datum.exposure

                    total_number_of_channels_used += datum.n_channels_used

                # iterate the interval
                i += 1

        stan_data["object_idx"] = np.array(grb_id).astype(int)
        stan_data["grb_id"] = np.array(grb_id).astype(int)

        stan_data["N_dets"] = n_dets

        stan_data["observed_counts"] = observed_counts
        stan_data["background_counts"] = background_counts
        stan_data["background_errors"] = background_errors
        stan_data["idx_background_nonzero"] = idx_background_nonzero.astype(int)
        stan_data["idx_background_zero"] = idx_background_zero.astype(int)
        stan_data["N_bkg_zero"] = n_bkg_zero.astype(int)
        stan_data["N_bkg_nonzero"] = n_bkg_nonzero.astype(int)
        stan_data["N_chan"] = n_chan.astype(int)
        stan_data["N_echan"] = n_echan.astype(int)

        stan_data["ebounds_lo"] = ebounds_lo
        stan_data["ebounds_hi"] = ebounds_hi
        stan_data["cbounds_lo"] = cbounds_lo
        stan_data["cbounds_hi"] = cbounds_hi
        stan_data["exposure"] = exposures
        stan_data["response"] = responses
        stan_data["mask"] = masks.astype(int)
        stan_data["N_channels_used"] = n_channels_used.astype(int)

        stan_data["dl"] = dl
        stan_data["z"] = z

        return stan_data
