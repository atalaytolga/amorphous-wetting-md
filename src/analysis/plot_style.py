import matplotlib as mpl
from cycler import cycler


def get_plot_style(
    mode="single_column",
    nrows=1,
    ncols=1,
    for_notebook=False,
):
    if mode == "single_column":
        fig_width = 3.37
        base_font = 9
        label_font = 10
        title_font = 10
        tick_font = 8
        legend_font = 8
    elif mode == "double_column":
        fig_width = 6.69
        base_font = 9
        label_font = 10
        title_font = 10
        tick_font = 8
        legend_font = 8
    else:
        raise ValueError("mode must be 'single_column' or 'double_column'")

    if for_notebook:
        fig_width = 10.0

    panel_aspect = 0.72
    fig_height = fig_width * panel_aspect * (nrows / ncols)

    prop_cycle = cycler(
        color=[
            "#4C72B0",
            "#DD8452",
            "#55A868",
            "#C44E52",
            "#8172B3",
            "#937860",
            "#DA8BC3",
        ]
    ) + cycler(
        linestyle=["-", "--", "-.", ":", "-", "--", "-."],
    )

    return {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "font.size": base_font,
        "axes.labelsize": label_font,
        "axes.titlesize": title_font,
        "xtick.labelsize": tick_font,
        "ytick.labelsize": tick_font,
        "legend.fontsize": legend_font,
        "figure.figsize": (fig_width, fig_height),
        "figure.dpi": 150,
        "figure.constrained_layout.use": True,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "axes.prop_cycle": prop_cycle,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.minor.size": 2,
        "ytick.minor.size": 2,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "lines.linewidth": 1.5,
        "lines.markersize": 4.5,
        "lines.markeredgewidth": 0.7,
        "legend.frameon": False,
        "legend.handlelength": 2.0,
        "svg.fonttype": "none",
    }


def set_plot_style(**kwargs):
    mpl.rcParams.update(get_plot_style(**kwargs))
