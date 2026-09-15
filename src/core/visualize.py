from ultraplot import Figure
from matplotlib import rcParams
from pathlib import Path
from src.io.model import read_vmodel
from numpy import array, asarray, arange, linspace
from pandas import read_csv, concat
from obspy.geodetics.base import degrees2kilometers as d2k

rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 9,
})


def plot_events(config):

    fig_outfile = Path("figures") / "events.png"

    catalog_ini_file = Path("inputs") / "catalog.csv"
    catalog_sel_file = Path("outputs/stage_01") / "select.csv"

    cat_ini = read_csv(catalog_ini_file)
    cat_sel = read_csv(catalog_sel_file)

    min_depth = config.FIGURE_SETTINGS.min_sei_depth
    max_depth = config.FIGURE_SETTINGS.max_sei_depth

    axes_shape = array(
        [
            [1, 1],
            [2, 2],
        ]
    )

    figsize = (5.5 / 2.54, 6.0 / 2.54)
    dpi = 300

    fig = Figure(
        figsize=figsize,
        dpi=dpi,
        share=True,
    )

    axes = fig.add_subplots(axes_shape)

    # ---------------------------------------------------------
    # Propertices
    # ---------------------------------------------------------
    for ax in axes:

        ax.format(
            xlabel=r"Longitude ($\degree$)",
            xlabelweight="bold",
            ylabel=r"Latitude ($\degree$)",
            ylabelweight="bold",
            xformatter="%.2f",
            yformatter="%.2f",
            xlocator=("maxn", 4),
            ylocator=("maxn", 4),
            fontsize=8,
            tickdir="in",
            ticklen=0.5,
            ticklenratio=0.5,
            tickwidth=0.2,
            tickwidthratio=0.5,
            tickminor=True,
        )

        ax.grid(
            which="major",
            alpha=0.2,
            linewidth=0.3,
            linestyle=":",
        )

    # ---------------------------------------------------------
    # Seismicity
    # ---------------------------------------------------------

    ax = axes[0]

    ax.scatter(
        cat_ini.lon.values,
        cat_ini.lat.values,
        s=cat_ini.mag.values * 0.01,
        c=cat_ini.dep.values,
        lw=0.1,
        ec="k",
        label="Initial",
        vmin=min_depth,
        vmax=max_depth,
    )

    ax.format(ultitle=f"N={cat_ini.shape[0]}")

    ax = axes[1]

    points = ax.scatter(
        cat_sel.lon.values,
        cat_sel.lat.values,
        s=cat_sel.mag.values * 0.01,
        c=cat_sel.dep.values,
        lw=0.1,
        ec="k",
        label="Initial",
        vmin=min_depth,
        vmax=max_depth,
    )

    ax.format(ultitle=f"N={cat_sel.shape[0]}")

    # ---------------------------------------------------------
    # Common geographic limits
    # ---------------------------------------------------------

    xmin = min(cat_ini.lon.min(), cat_sel.lon.min())
    xmax = max(cat_ini.lon.max(), cat_sel.lon.max())

    ymin = min(cat_ini.lat.min(), cat_sel.lat.min())
    ymax = max(cat_ini.lat.max(), cat_sel.lat.max())

    for ax in axes:
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

    # ---------------------------------------------------------
    # Colorbar
    # ---------------------------------------------------------

    cbar = fig.colorbar(
        points,
        loc="right",
        label="Depth (km)",
        pad=0.1,
        width=0.1,
    )
    cbar.ax.invert_yaxis()

    fig.save(
        fig_outfile,
        bbox_inches="tight",
    )


def plot_models(config, stage_n, run_n):

    if stage_n == 1:

        root = Path("outputs/stage_01")
        run_dir = root / f"run_{run_n:02d}"

        fig_outfile = Path("figures") / f"models_st01_run_{run_n:02d}.png"

    elif stage_n == 3:

        root = Path("outputs/stage_03")
        run_dir = root / f"run_{run_n:02d}"

        fig_outfile = Path("figures") / f"models_st03_run_{run_n:02d}.png"

    imodel = config.INPUTS.vmodel
    imodel = read_vmodel(imodel)

    imodel_files = sorted(Path(run_dir).rglob("model.mod"))
    omodels_files = sorted(Path(run_dir).rglob("velocity_summary.csv"))
    amodel_file = run_dir / "summary" / "averaged.csv"
    smodel_file = run_dir / "summary" / "simplified.csv"

    imodels = [read_vmodel(m) for m in imodel_files]
    omodels = []
    for omodel_file in omodels_files:
        omodel = read_csv(omodel_file)
        state_type = (
            "final_state"
            if "final_state" in omodel.state_type.values
            else "iteration"
        )
        omodel = omodel[
            (omodel.model == 1)
            & (omodel.iteration == omodel.iteration.max())
            & (omodel.state_type == state_type)
        ]
        omodels.append(omodel[["velocity", "depth"]])

    amodel = read_csv(amodel_file)
    amodel = amodel[
        (amodel.model == 1) & (amodel.iteration == amodel.iteration.max())
    ]
    smodel = read_csv(smodel_file)
    smodel = smodel[
        (smodel.model == 1) & (smodel.iteration == smodel.iteration.max())
    ]

    figsize = (4.5 / 2.54, 6.0 / 2.54)
    dpi = 300
    show_initial = True
    xlim = (
        config.FIGURE_SETTINGS.min_vel_vp,
        config.FIGURE_SETTINGS.max_vel_vp,
    )
    ylim = (
        config.FIGURE_SETTINGS.max_vel_depth,
        config.FIGURE_SETTINGS.min_vel_depth,
    )

    # ---------------------------------------------------------
    # Select velocity type
    # ---------------------------------------------------------

    layer_attribute = "layers_vp"
    initial_layers = imodel.layers_vp

    # ---------------------------------------------------------
    # Create figure
    # ---------------------------------------------------------

    fig = Figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot()

    # ---------------------------------------------------------
    # Plot initial/reference model
    # ---------------------------------------------------------

    if show_initial:

        velocity = asarray(
            [layer.vel for layer in initial_layers],
            dtype=float,
        )

        depth = asarray(
            [layer.depth for layer in initial_layers],
            dtype=float,
        )

        ax.plot(
            velocity,
            depth,
            drawstyle="steps-pre",
            linewidth=0.6,
            label="Initial model",
            color="violet5",
            zorder=70,
        )

    # ---------------------------------------------------------
    # Plot random models
    # ---------------------------------------------------------

    for i, model in enumerate(imodels):

        layers = getattr(model, layer_attribute)

        velocity = asarray(
            [layer.vel for layer in layers],
            dtype=float,
        )

        depth = asarray(
            [layer.depth for layer in layers],
            dtype=float,
        )

        # Stepwise layered velocity model
        ax.plot(
            velocity,
            depth,
            drawstyle="steps-pre",
            linewidth=0.4,
            label="Random models" if i == 0 else None,
            color="gray5",
            zorder=50,
        )

    # ---------------------------------------------------------
    # Plot inverted models
    # ---------------------------------------------------------

    for i, model in enumerate(omodels):

        velocity = model.velocity.values
        depth = model.depth.values

        # Stepwise layered velocity model
        ax.plot(
            velocity,
            depth,
            drawstyle="steps-pre",
            linewidth=0.6,
            label="Inverted models" if i == 0 else None,
            alpha=0.50,
            color="teal5",
            zorder=60,
        )

    # ---------------------------------------------------------
    # Plot averaged model
    # ---------------------------------------------------------

    velocity = amodel.velocity.values
    depth = amodel.depth.values

    # Stepwise layered velocity model
    ax.plot(
        velocity,
        depth,
        drawstyle="steps-pre",
        linewidth=0.6,
        label="Averaged model",
        color="orange5",
        zorder=80,
    )

    # ---------------------------------------------------------
    # Plot simplified model
    # ---------------------------------------------------------

    velocity = smodel.velocity.values
    depth = smodel.depth.values

    # Stepwise layered velocity model
    ax.plot(
        velocity,
        depth,
        drawstyle="steps-pre",
        linewidth=1.0,
        label="Simplified model",
        color="darkred",
        zorder=90,
    )

    # ---------------------------------------------------------
    # Axes
    # ---------------------------------------------------------

    # Depth increases downward
    ax.invert_yaxis()

    # ---------------------------------------------------------
    # Ticks
    # ---------------------------------------------------------

    ax.tick_params(
        which="major",
        direction="in",
        length=2,
        width=0.4,
        labelsize=5,
        top=True,
        right=True,
    )

    ax.tick_params(
        which="minor",
        direction="in",
        length=1,
        width=0.3,
        top=True,
        right=True,
    )

    # ---------------------------------------------------------
    # Minor ticks
    # ---------------------------------------------------------

    ax.minorticks_on()

    # ---------------------------------------------------------
    # Formatter
    # ---------------------------------------------------------

    ax.format(
        xlabel="Vp (km/s)",
        xlabelweight="bold",
        xformatter="%.1f",
        xlim=xlim,
        ylabel="Depth (km)",
        ylabelweight="bold",
        ylim=ylim,
    )

    # ---------------------------------------------------------
    # Grid
    # ---------------------------------------------------------

    ax.grid(which="major", alpha=0.2, linewidth=0.3, linestyle=":")

    # ---------------------------------------------------------
    # Legend
    # ---------------------------------------------------------

    ax.legend(fontsize=4, frameon=False, loc="lower left", ncol=1)

    fig.save(fig_outfile)


def plot_statistics(config, stage_n, run_n):
    if stage_n == 1:

        root = Path("outputs/stage_01")
        run_dir = root / f"run_{run_n:02d}"

        fig_outfile = Path("figures") / f"statistics_st01_run_{run_n:02d}.png"

    elif stage_n == 3:

        root = Path("outputs/stage_03")
        run_dir = root / f"run_{run_n:02d}"

        fig_outfile = Path("figures") / f"statistics_st03_run_{run_n:02d}.png"

    attempts_files = sorted(Path(run_dir).rglob("attempts.csv"))

    attempts = [read_csv(file) for file in attempts_files]
    attempts = concat(attempts, ignore_index=True)

    # ---------------------------------------------------------
    # Statistics
    # ---------------------------------------------------------

    datavar = (
        attempts.groupby("iteration")["data_variance"]
        .agg(
            mean="mean",
            std="std",
            q05=lambda x: x.quantile(0.05),
            q10=lambda x: x.quantile(0.10),
            q25=lambda x: x.quantile(0.25),
            q75=lambda x: x.quantile(0.75),
            q90=lambda x: x.quantile(0.90),
            q95=lambda x: x.quantile(0.95),
        )
        .sort_index()
    )

    rms = (
        attempts.groupby("iteration")["rms_residual"]
        .agg(
            mean="mean",
            std="std",
            q05=lambda x: x.quantile(0.05),
            q10=lambda x: x.quantile(0.10),
            q25=lambda x: x.quantile(0.25),
            q75=lambda x: x.quantile(0.75),
            q90=lambda x: x.quantile(0.90),
            q95=lambda x: x.quantile(0.95),
        )
        .sort_index()
    )

    # ---------------------------------------------------------
    # Figure
    # ---------------------------------------------------------

    figsize = (9.0 / 2.54, 6.0 / 2.54)
    dpi = 300

    fig = Figure(figsize=figsize, dpi=dpi, share=False)

    # Two genuinely independent, adjacent axes
    ax_left, ax_right = fig.add_subplots(
        ncols=2,
    )

    # ---------------------------------------------------------
    # Common x
    # ---------------------------------------------------------

    x = datavar.index

    # =========================================================
    # LEFT AXIS — Data Variance
    # =========================================================

    # 5–95%
    ax_left.fill_between(
        x,
        datavar["q05"],
        datavar["q95"],
        alpha=0.10,
        color="blue",
    )

    # 25–75%
    ax_left.fill_between(
        x,
        datavar["q25"],
        datavar["q75"],
        alpha=0.25,
        color="blue",
    )

    # Mean
    ax_left.plot(
        x,
        datavar["mean"],
        linewidth=1,
        label="Mean",
        color="blue",
    )

    # Mean + std
    ax_left.plot(
        x,
        datavar["mean"] + datavar["std"],
        linestyle="--",
        linewidth=0.5,
        label="Mean ± std",
        color="blue",
    )

    # Mean - std
    ax_left.plot(
        x,
        datavar["mean"] - datavar["std"],
        linestyle="--",
        linewidth=0.5,
        color="blue",
    )

    ax_left.set_ylabel(
        r"Data Variance ($s^2$)",
        fontweight="bold",
    )

    # =========================================================
    # RIGHT AXIS — RMS Residual
    # =========================================================

    # 5–95%
    ax_right.fill_between(
        x,
        rms["q05"],
        rms["q95"],
        alpha=0.10,
        color="red",
    )

    # 25–75%
    ax_right.fill_between(
        x,
        rms["q25"],
        rms["q75"],
        alpha=0.25,
        color="red",
    )

    # Mean
    ax_right.plot(
        x,
        rms["mean"],
        linewidth=1,
        label="Mean",
        color="red",
    )

    # Mean + std
    ax_right.plot(
        x,
        rms["mean"] + rms["std"],
        linestyle="--",
        linewidth=0.5,
        label="Mean ± std",
        color="red",
    )

    # Mean - std
    ax_right.plot(
        x,
        rms["mean"] - rms["std"],
        linestyle="--",
        linewidth=0.5,
        color="red",
    )

    ax_right.set_ylabel(
        r"RMS Residual ($s$)",
        fontweight="bold",
    )

    # ---------------------------------------------------------
    # X labels
    # ---------------------------------------------------------

    ax_left.set_xlabel(
        "Iteration",
        fontweight="bold",
    )

    ax_right.set_xlabel(
        "Iteration",
        fontweight="bold",
    )

    # ---------------------------------------------------------
    # Integer x ticks
    # ---------------------------------------------------------

    xticks = range(
        int(x.min()),
        int(x.max()) + 1,
    )

    ax_left.format(
        xlocator=xticks,
        fontsize=6,
    )

    ax_right.format(
        xlocator=xticks,
        fontsize=6,
    )

    # ---------------------------------------------------------
    # Ticks
    # ---------------------------------------------------------

    for ax in (ax_left, ax_right):

        ax.tick_params(
            which="major",
            direction="in",
            length=1,
            width=0.5,
            labelsize=5,
            top=True,
            right=True,
        )

        ax.tick_params(
            which="minor",
            direction="in",
            length=0.5,
            width=0.15,
            top=True,
            right=True,
        )

        ax.minorticks_on()

    # ---------------------------------------------------------
    # Grid
    # ---------------------------------------------------------

    ax_left.grid(
        which="major",
        alpha=0.2,
        linewidth=0.3,
        linestyle=":",
    )

    ax_right.grid(
        which="major",
        alpha=0.2,
        linewidth=0.3,
        linestyle=":",
    )

    # ---------------------------------------------------------
    # Legends
    # ---------------------------------------------------------

    ax_left.legend(
        fontsize=4,
        frameon=False,
        loc="upper right",
        ncol=1,
    )

    ax_right.legend(
        fontsize=4,
        frameon=False,
        loc="upper right",
        ncol=1,
    )

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------

    fig.save(fig_outfile)


def plot_dislocations(config, stage_n, run_n):

    fig_outfile = Path("figures") / f"dislocation_{run_n:02d}.png"

    catalog_ini_file = Path("outputs/stage_02") / "select.csv"
    catalog_prt_file = (
        Path("outputs/stage_04") / f"run_{run_n:02d}" / "perturb.csv"
    )
    catalog_fin_file = (
        Path("outputs/stage_04") / f"run_{run_n:02d}" / "relocated.csv"
    )

    cat_ini = read_csv(catalog_ini_file)
    cat_prt = read_csv(catalog_prt_file)
    cat_fin = read_csv(catalog_fin_file)

    delta = [
        config.STAGE_04.dx_km*1.5,
        config.STAGE_04.dy_km*1.5,
        config.STAGE_04.dz_km*1.5,
    ]

    # ------------------------------------------------------------------
    # Dislocations
    # ------------------------------------------------------------------
    components = ["lat", "lon", "dep"]

    dis_prt = cat_prt[components] - cat_ini[components]
    dis_fin = cat_fin[components] - cat_ini[components]

    axes_shape = array(
        [
            [1, 2, 3],
            [4, 5, 6],
        ]
    )

    figsize = (5.5 / 2.54, 6.0 / 2.54)
    dpi = 300

    fig = Figure(
        figsize=figsize,
        dpi=dpi,
        share=False,
    )

    axs = fig.add_subplots(axes_shape)
    axs.grid(
        which="major",
        alpha=0.2,
        linewidth=0.3,
        linestyle=":",
    )
    # ---------------------------------------------------------
    # Dislocations
    # ---------------------------------------------------------

    for j, component in enumerate(components):

        # ==============================================================
        # Row 1: dislocations for individual events
        # ==============================================================
        ax = axs[0, j]
        ax.format(
            ylim=(-delta[j], delta[j]),
            xlabel="Event",
            ylabel=f"$\\Delta${component} (km)",
        )

        x = arange(len(cat_ini))

        ax.scatter(
            x,
            d2k(dis_prt[component].values) if j < 2 else dis_prt[component].values,
            s=cat_ini.mag.values * 0.01,
            c="#E07A3F",
            lw=0.1,
            ec="k",
            # label="Initial → Perturbed",
        )
        ax.scatter(
            x,
            d2k(dis_fin[component].values) if j < 2 else dis_fin[component].values,
            s=cat_ini.mag.values * 0.01,
            c="#4C78A8",
            lw=0.1,
            ec="k",
            # label="Initial → Final",
        )

        ax.axhline(0, lw=0.8, ls="--")

        # ==============================================================
        # Row 2: histograms
        # ==============================================================
        ax = axs[1, j]

        ax.format(
            xlim=(-delta[j], delta[j]),
            ylabel="Number of events",
            xlabel=f"$\\Delta${component} (km)",
        )

        data = (
            d2k(dis_prt[component].dropna())
            if j < 2
            else dis_prt[component].dropna()
        ).values
        mean_prt = data.mean()
        std_prt = data.std()
        ax.hist(
            data,
            linspace(-delta[j], delta[j], 11),
            lw=0.3,
            histtype="bar",
            filled=True,
            alpha=0.7,
            edgecolor="k",
            color="#E07A3F",
            label="Initial → Perturbed" if j ==0 else None,
        )
        data = (
            d2k(dis_fin[component].dropna())
            if j < 2
            else dis_fin[component].dropna()
        ).values
        mean_fin = data.mean()
        std_fin = data.std()
        ax.hist(
            data,
            linspace(-delta[j], delta[j], 11),
            lw=0.3,
            histtype="bar",
            filled=True,
            alpha=0.7,
            edgecolor="k",
            color="#4C78A8",
            label="Initial → Final" if j ==0 else None,
        )

        ax.text(
            0.97,
            0.97,
            (
                f"Perturbed: μ = {mean_prt:.2f}\n"
                f"             σ = {std_prt:.2f}\n"
                f"Final:       μ = {mean_fin:.2f}\n"
                f"             σ = {std_fin:.2f}"
            ),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
        )
        ax.axvline(0, lw=0.8, ls="--")

    fig.legend(loc="t")

    fig.save(
        fig_outfile,
        bbox_inches="tight",
    )
