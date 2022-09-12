import collections

from pathlib import Path

import numpy as np
import yaml

from cosmogrb.io.gbm_fits import grbsave_to_gbm_fits
from cosmogrb.universe.survey import Survey
from cosmogrb.utils.file_utils import if_directory_not_existing_then_make
from threeML import TimeSeriesBuilder


class GRBProcessor(object):
    def __init__(self, gbm_grb, n_nai_to_use: int = 3, use_bb: bool = False):
        """
        :param gbm_grb:
        :param n_nai_to_use:
        :returns:
        :rtype:

        """

        self._grb_save = gbm_grb
        assert n_nai_to_use > 0, "yo use some detectors"

        self._n_nai_to_use: int = int(n_nai_to_use)

        self._use_bb: bool = use_bb

        self._config_dict = collections.OrderedDict()

        self._config_dict["z"] = float(self._grb_save.z)

        if_directory_not_existing_then_make(self._grb_save.name)

        # gets the light curves we want
        self._setup_order_by_distance()

        self._create_fits_files()

        self._threeml_process()

    def _setup_order_by_distance(self):

        # now we will go through the lightcurves
        # collect

        angular_distances = []
        bgo_anglular_distance = 1000
        lc_names = []

        for name, det in self._grb_save.items():

            lc = det["lightcurve"]

            if name.startswith("n"):

                lc_names.append(str(name))

                angular_distances.append(lc.extra_info["angle"])

            else:
                if lc.extra_info["angle"] < bgo_anglular_distance:
                    bgo_anglular_distance = lc.extra_info["angle"]
                    bgo_det = str(name)

        angular_distances = np.array(angular_distances)
        lc_names = np.array(lc_names)

        # now attach the lc_names sorted by
        # the detector distance to the GRB

        idx = angular_distances.argsort()

        self._lc_names = list(lc_names[idx][: self._n_nai_to_use])
        self._lc_names.append(bgo_det)
        self._lc_names = [str(x) for x in self._lc_names]

    def _create_fits_files(self):

        self._fits_files = grbsave_to_gbm_fits(
            self._grb_save,
            destination=self._grb_save.name,
            detectors=self._lc_names,
        )

    def _threeml_process(self):

        self._config_dict["dir"] = str(Path(self._grb_save.name).absolute())

        det_dic = {}

        for i, name in enumerate(self._lc_names):

            if name.startswith("n"):

                selection = "10-900"

            else:

                selection = "250-30000"

            det_dic[name] = selection

            ts = TimeSeriesBuilder.from_gbm_tte(
                name=name,
                tte_file=self._fits_files[name]["tte"],
                rsp_file=self._fits_files[name]["rsp"],
                poly_order=0,
                verbose=False,
            )

            ts.set_background_interval(
                "-20--5",
                f"{self._grb_save.duration + 5}- {self._grb_save.duration + 20}",
            )

            # for now do nothing else

            if self._use_bb:

                if i < 1:

                    ts.create_time_bins(
                        -25,
                        self._grb_save.duration + 1,
                        method="bayesblocks",
                        p0=0.1,
                    )

                    bins_to_use = ts

                else:

                    ts.read_bins(bins_to_use)

                intervals = ts.bins.containing_interval(0,
                                                        self._grb_save.duration,
                                                        inner=False)

                n_intervals = len(
                    intervals
                )
                first_interval_num = 0

                # check for the first bin in intervals if it is mostly
                # before the GRB. If more than 50% of the time interval
                # is before GRB we want to jump this interval for the fits
                if ((0 - intervals[0].start) / (intervals[0].stop - intervals[0].start)) > 0.5:
                    first_interval_num = 1
                    n_intervals -= 1

                # check for the last bin in intervals if it is mostly
                # after the GRB. If more than 50% of the time interval
                # is after the GRB we want to skip this interval for the fits
                if ((intervals[-1].stop-self._grb_save.duration) / (intervals[-1].stop - intervals[-1].start)) > 0.5:
                    n_intervals -= 1

                if n_intervals > 1:

                    ts.write_pha_from_binner(
                        file_name=Path(self._grb_save.name) / name,
                        start=0.0,
                        stop=self._grb_save.duration,
                        #inner=True,
                        force_rsp_write=True,
                        overwrite=True,
                    )

                self._config_dict["n_intervals"] = n_intervals
                self._config_dict["first_interval_num"] = first_interval_num

                if n_intervals > 1:

                    fig = ts.view_lightcurve(use_binner=True)

                    fig.savefig(
                        f"{Path(self._grb_save.name) / name}_lc.png",
                        bbox_inches="tight",
                    )

            else:

                self._config_dict["n_intervals"] = 1

                ts.set_active_time_interval(f"0-{self._grb_save.duration}")

                plugin = ts.to_spectrumlike()

                plugin.write_pha(
                    filename=Path(self._grb_save.name) / name,
                    force_rsp_write=True,
                    overwrite=True,
                )

            self._config_dict["detectors"] = det_dic

            for k, v in self._fits_files[name].items():
                Path(v).unlink()

    @property
    def yaml_params(self):

        return self._config_dict


class AnalysisBuilder(object):
    def __init__(self, survey_file, use_all=False, use_bb=False):

        if isinstance(survey_file, str):

            self._survey = Survey.from_file(survey_file)

        else:

            assert isinstance(survey_file, Survey)

            self._survey = survey_file

        self._config_dict = collections.OrderedDict()
        for k, v in self._survey.items():

            print(k)

            process = GRBProcessor(v.grb, use_bb=use_bb)

            self._config_dict[k] = process.yaml_params

    def write_yaml(self, file_name: str) -> None:
        """TODO describe function

        :param file_name:
        :type file_name: str
        :returns:

        """

        with open(file_name, "w") as f:
            yaml.dump(self.yaml_params, stream=f, default_flow_style=False)

    @property
    def yaml_params(self):
        return self._config_dict
