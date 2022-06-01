from threeML.io.plotting.data_residual_plot import ResidualPlot
import warnings

warnings.simplefilter("ignore")
# from threeML_utils.colors import Colors

NO_REBIN = 1e-99


# def scale_colour(self, colour, scalefactor):    # pragma: no cover
#     if isinstance(colour, np.ndarray):
#         r, g, b = colour[:3] * 255.0
#     else:
#         hexx = colour.strip('#')
#         if scalefactor < 0 or len(hexx) != 6:
#             return hexx
#         r, g, b = int(hexx[:2], 16), int(hexx[2:4], 16), int(hexx[4:], 16)
#     r = self._clamp(int(r * scalefactor))
#     g = self._clamp(int(g * scalefactor))
#     b = self._clamp(int(b * scalefactor))
#     return "#%02x%02x%02x" % (r, g, b)


def display_posterior_model_counts(
    plugin,
    model,
    samples,
    data_color="k",
    model_color="r",
    thin=100,
    min_rate=1,
    shade=True,
    q_level=68,
    gradient=0.6,
    axes=None,
    **kwargs
):

    show_residuals = False

    if axes != None:
        residual_plot = ResidualPlot(
            show_residuals=show_residuals, model_subplot=axes
        )
    else:
        residual_plot = ResidualPlot(show_residuals=show_residuals)
        axes = residual_plot.data_axis

    plugin.set_model(model)

    show_legend = False

    for params in samples:

        model.set_free_parameters(params)

        # first with no data

        plugin.display_model(
            data_color=data_color,
            model_color=model_color,
            min_rate=min_rate,
            step=False,
            show_residuals=False,
            show_data=False,
            show_legend=show_legend,
            ratio_residuals=False,
            model_label=None,
            model_subplot=axes,
            model_kwargs=dict(alpha=0.1),
            **kwargs
            #                model_subplot=axes,
            #                data_kwargs=data_kwargs,
        )

    plugin.display_model(
        data_color=data_color,
        model_color=model_color,
        min_rate=min_rate,
        step=True,
        show_residuals=False,
        show_data=True,
        show_legend=show_legend,
        ratio_residuals=False,
        #        model_label=model_label,
        model_subplot=axes,
        model_kwargs=dict(alpha=0.0),  # no model
        **kwargs
        #    data_kwargs=data_kwargs
    )

    return residual_plot.figure
