from pathlib import Path
from src.model.Models import LAYER, VMODEL
from pandas import read_csv, DataFrame, concat
from numpy import array_equal, allclose, arange
from shutil import copy


def read_vmodel(vmodel_file: str | Path) -> VMODEL:
    vmodel = VMODEL.model_construct()
    with open(vmodel_file) as f:
        for l, line in enumerate(f):
            if l == 0:
                vmodel.model_name = line.strip()
            elif l == 1:
                vmodel.n_layers, header = [
                    int(v) if i == 0 else v for i, v in enumerate(line.split())
                ]
            elif 1 < l < vmodel.n_layers + 2:
                vel, depth, vdamp = [float(i) for i in line.split()[:3]]
                layer = LAYER.model_validate(
                    {"vel": vel, "depth": depth, "vdamp": vdamp}
                )
                vmodel.layers_vp.append(layer)
            elif vmodel.n_layers + 2 < l < 2 * (vmodel.n_layers + 2):
                vel, depth, vdamp = [float(i) for i in line.split()[:3]]
                layer = LAYER.model_validate(
                    {"vel": vel, "depth": depth, "vdamp": vdamp}
                )
                vmodel.layers_vs.append(layer)
    return vmodel


def simplify_model(
    df,
    averaged_modelfile,
    simplified_modelfile,
    velocity_threshold=0.2,
):
    """
    Calculate the average velocity model from the last inversion
    iteration and simplify it by merging adjacent layers.

    The first layer is always preserved and is never merged.

    Vs is optional. If no Vs model (model=2) is present, the
    simplification is performed using Vp only.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing all inverted velocity models.

        Required columns:

            iteration
            model
            layer
            depth
            velocity

        model=1 -> Vp
        model=2 -> Vs

    averaged_modelfile : str or path-like
        Output file for the averaged model.

    simplified_modelfile : str or path-like
        Output file for the simplified model.

    velocity_threshold : float, default=0.2
        Maximum absolute velocity difference between two adjacent
        layers for merging.

        If Vs exists, a merge occurs only when both conditions
        are satisfied:

            abs(dVp) <= velocity_threshold
            abs(dVs) <= velocity_threshold

        If Vs does not exist, only:

            abs(dVp) <= velocity_threshold

        is required.

    Returns
    -------
    average_df : pandas.DataFrame
        Averaged model from the last inversion iteration.

    final_df : pandas.DataFrame
        Final averaged and simplified velocity model.
    """

    # =============================================================
    # 1. Validation
    # =============================================================

    required_columns = {
        "iteration",
        "model",
        "layer",
        "depth",
        "velocity",
    }

    missing = required_columns.difference(df.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if df.empty:
        raise ValueError("Input DataFrame is empty.")

    if velocity_threshold < 0:
        raise ValueError("velocity_threshold must be >= 0.")

    data = df.copy()

    # =============================================================
    # 2. Select last inversion iteration
    # =============================================================

    state_type = (
        "final_state"
        if "final_state" in data.state_type.values
        else "iteration"
    )

    last_iteration = data[data.state_type == state_type].iteration.max()

    data = data[
        (data["iteration"] == last_iteration) & (data.state_type == state_type)
    ].copy()

    # =============================================================
    # 3. Check Vp and Vs availability
    # =============================================================

    models = set(data["model"].unique())

    # Vp is mandatory
    if 1 not in models:
        raise ValueError("No Vp model found (model=1).")

    # Vs is optional
    has_vs = 2 in models

    # =============================================================
    # 4. Average Vp
    # =============================================================

    vp = (
        data[data["model"] == 1]
        .groupby("layer", as_index=False)
        .agg(
            depth=("depth", "mean"),
            velocity=("velocity", "mean"),
        )
        .sort_values("layer")
        .reset_index(drop=True)
    )

    # =============================================================
    # 5. Average Vs -- only if available
    # =============================================================

    if has_vs:

        vs = (
            data[data["model"] == 2]
            .groupby("layer", as_index=False)
            .agg(
                depth=("depth", "mean"),
                velocity=("velocity", "mean"),
            )
            .sort_values("layer")
            .reset_index(drop=True)
        )

    else:

        vs = None

    # =============================================================
    # 6. Check Vp / Vs structure
    # =============================================================

    if has_vs:

        if not array_equal(
            vp["layer"].to_numpy(),
            vs["layer"].to_numpy(),
        ):
            raise ValueError("Vp and Vs have different layer structures.")

        if not allclose(
            vp["depth"].to_numpy(),
            vs["depth"].to_numpy(),
        ):
            raise ValueError("Vp and Vs have different averaged layer depths.")

    # =============================================================
    # 7. Construct averaged DataFrame
    # =============================================================

    average_vp = DataFrame(
        {
            "iteration": last_iteration,
            "model": 1,
            "layer": vp["layer"],
            "depth": vp["depth"],
            "velocity": vp["velocity"],
        }
    )

    if has_vs:

        average_vs = DataFrame(
            {
                "iteration": last_iteration,
                "model": 2,
                "layer": vs["layer"],
                "depth": vs["depth"],
                "velocity": vs["velocity"],
            }
        )

        average_df = concat(
            [
                average_vp,
                average_vs,
            ],
            ignore_index=True,
        )

    else:

        average_df = average_vp.copy()

    # =============================================================
    # 8. Create common layer representation
    # =============================================================

    layers = DataFrame(
        {
            "layer": vp["layer"].to_numpy(),
            "depth": vp["depth"].to_numpy(),
            "vp": vp["velocity"].to_numpy(dtype=float),
        }
    )

    if has_vs:

        layers["vs"] = vs["velocity"].to_numpy(dtype=float)

    # =============================================================
    # 9. Simplify the model
    #
    # IMPORTANT:
    #
    # Index 0 = first layer.
    #
    # It is permanently fixed and is NEVER merged.
    #
    # Therefore simplification starts from index 1.
    # =============================================================

    changed = True

    while changed:

        changed = False

        # Start from layer 2.
        #
        # Index 0 is the fixed first layer.

        i = 1

        while i < len(layers) - 1:

            current = layers.iloc[i]
            deeper = layers.iloc[i + 1]

            # -----------------------------------------------------
            # Vp velocity difference
            # -----------------------------------------------------

            dvp = abs(deeper["vp"] - current["vp"])

            # -----------------------------------------------------
            # Merge criterion
            # -----------------------------------------------------

            if has_vs:

                dvs = abs(deeper["vs"] - current["vs"])

                merge = dvp <= velocity_threshold and dvs <= velocity_threshold

            else:

                merge = dvp <= velocity_threshold

            # -------------------------------------------------
            # Merge current layer into deeper layer
            # -------------------------------------------------

            if merge:

                # -------------------------------------------------
                # Average Vp
                # -------------------------------------------------

                layers.at[
                    layers.index[i + 1],
                    "vp",
                ] = (current["vp"] + deeper["vp"]) / 2.0

                # -------------------------------------------------
                # Average Vs if available
                # -------------------------------------------------

                if has_vs:

                    layers.at[
                        layers.index[i + 1],
                        "vs",
                    ] = (current["vs"] + deeper["vs"]) / 2.0

                # -------------------------------------------------
                # Keep the TOP depth of the first layer
                # -------------------------------------------------

                layers.at[
                    layers.index[i + 1],
                    "depth",
                ] = current["depth"]

                # -------------------------------------------------
                # Remove current layer
                # -------------------------------------------------

                layers = layers.drop(index=layers.index[i]).reset_index(
                    drop=True
                )

                changed = True

                # Do NOT increment i
                #
                # The newly merged layer must be tested against
                # the next deeper layer.
                # -------------------------------------------------

            else:

                i += 1

    # =============================================================
    # 10. Renumber layers
    # =============================================================

    layers["layer"] = arange(
        1,
        len(layers) + 1,
    )

    # =============================================================
    # 11. Convert back to long DataFrame
    # =============================================================

    final_vp = DataFrame(
        {
            "iteration": last_iteration,
            "model": 1,
            "layer": layers["layer"],
            "depth": layers["depth"],
            "velocity": layers["vp"],
        }
    )

    if has_vs:

        final_vs = DataFrame(
            {
                "iteration": last_iteration,
                "model": 2,
                "layer": layers["layer"],
                "depth": layers["depth"],
                "velocity": layers["vs"],
            }
        )

        final_df = concat(
            [
                final_vp,
                final_vs,
            ],
            ignore_index=True,
        )

    else:

        final_df = final_vp.copy()

    # =============================================================
    # 12. Save results
    # =============================================================

    average_df.to_csv(
        averaged_modelfile,
        index=False,
    )

    final_df.to_csv(
        simplified_modelfile,
        index=False,
    )

    return average_df, final_df


def write_vmodel(vmodel: VMODEL, vmodelout: str | Path):
    with open(vmodelout, "w") as fp:
        fp.write(f"{vmodel.model_name[:40]:40}\n")
        fp.write(
            f"{vmodel.n_layers:3d}        vel,depth,vdamp,phase(f5.2,5x,f7.2,2x,f7.3,3x,a1)\n"
        )
        for i, layer in enumerate(vmodel.layers_vp):
            line = f"{layer.vel:5.2f}     {layer.depth:5.2f}    {layer.vdamp:5.2f}"
            if i == 0:
                fp.write(line + "           P-VELOCITY MODEL\n")
            else:
                fp.write(line + "\n" if i + 2 <= vmodel.n_layers else line)
        if bool(vmodel.layers_vs):
            fp.write(f"\n{vmodel.n_layers:3d}\n")
            for i, layer in enumerate(vmodel.layers_vs):
                line = f"{layer.vel:5.2f}     {layer.depth:5.2f}    {layer.vdamp:5.2f}"
                if i == 0:
                    fp.write(line + "           S-VELOCITY MODEL\n")
                else:
                    fp.write(
                        line + "\n" if i + 2 <= 2 * vmodel.n_layers else line
                    )


def simplified2vmodel(config, simplified_modelfile, final_modelfile):
    first_layer_dmp = config.SYNTHETIC_MODELS.first_layer_dmp
    last_layer_dmp = config.SYNTHETIC_MODELS.last_layer_dmp
    vmodel = VMODEL.model_construct()
    model_df = read_csv(simplified_modelfile)
    for r, row in model_df.iterrows():
        if row.model == 1:
            vmodel.model_name = "Simplified"
            vmodel.n_layers = (
                int(model_df.shape[0] / 2)
                if 2 in model_df.model.values
                else model_df.shape[0]
            )
            layer = LAYER.model_validate(
                {
                    "vel": row.velocity,
                    "depth": row.depth,
                    "vdamp": (
                        first_layer_dmp
                        if r == 0
                        else (
                            last_layer_dmp
                            if r == (vmodel.n_layers - 1)
                            else 1.0
                        )
                    ),
                }
            )
            vmodel.layers_vp.append(layer)
        elif row.model == 2:
            vmodel.n_layers /= 2
            layer = LAYER.model_validate(
                {
                    "vel": row.velocity,
                    "depth": row.depth,
                    "vdamp": (
                        first_layer_dmp
                        if r == vmodel.n_layers
                        else (
                            last_layer_dmp
                            if r == (vmodel.n_layers - 1) * 2
                            else 1.0
                        )
                    ),
                }
            )
            vmodel.layers_vs.append(layer)
    write_vmodel(vmodel, final_modelfile)


def prepare_model(config, stage_n, run_n=None):
    if stage_n == 1:

        # outputs
        rootout = Path("outputs/stage_01")
        run_dir = rootout / f"run_{run_n:02d}"

        # Read / Implement / Write
        if run_n == 1:
            vmodel_input = Path("inputs") / "vmodel.dat"
        else:
            vmodel_input = (
                rootout / f"run_{run_n-1:02d}" / "summary" / "model.mod"
            )
        n_models = config.SYNTHETIC_MODELS.n_models
        velocity_std = config.SYNTHETIC_MODELS.velocity_std
        depth_std = config.SYNTHETIC_MODELS.depth_std
        seed = config.SYNTHETIC_MODELS.seed
        vmodel = read_vmodel(vmodel_input)
        vmodels = vmodel.generate_random_models(
            n_models, velocity_std, depth_std, seed
        )
        for n_model, vmodel in enumerate(vmodels, 1):
            model_dir = run_dir / f"model_{n_model:02d}"
            vmodel_out = model_dir / "model.mod"
            write_vmodel(vmodel, vmodel_out)

    elif stage_n == 2:

        # inputs
        run_n = config.STAGE_01.n_runs
        model_input = (
            Path("outputs/stage_01")
            / f"run_{run_n:02d}"
            / "summary"
            / "model.mod"
        )

        # outputs
        run_dir = Path("outputs/stage_02")

        # Read / Implement / Write
        copy(model_input, run_dir)

    if stage_n == 3:

        # outputs
        rootout = Path("outputs/stage_03")
        run_dir = rootout / f"run_{run_n:02d}"

        # Read / Implement / Write
        if run_n == 1:
            vmodel_input = Path("outputs/stage_02") / "model.mod"
        else:
            vmodel_input = (
                rootout / f"run_{run_n-1:02d}" / "summary" / "model.mod"
            )
        n_models = config.STAGE_03.n_models
        velocity_std = config.SYNTHETIC_MODELS.velocity_std
        depth_std = config.SYNTHETIC_MODELS.depth_std
        seed = config.SYNTHETIC_MODELS.seed
        vmodel = read_vmodel(vmodel_input)
        vmodels = vmodel.generate_random_models(
            n_models, velocity_std, depth_std, seed
        )
        for n_model, vmodel in enumerate(vmodels, 1):
            model_dir = run_dir / f"model_{n_model:02d}"
            vmodel_out = model_dir / "model.mod"
            write_vmodel(vmodel, vmodel_out)

    if stage_n == 4:

        stg_03_run_n = config.STAGE_03.n_runs

        # outputs
        rootout = Path("outputs/stage_04")
        run_dir = rootout / f"run_{run_n:02d}"

        # Read / Implement / Write
        vmodel_input = (
            Path(f"outputs/stage_03/run_{stg_03_run_n:02d}/summary")
            / "model.mod"
        )
        vmodel = read_vmodel(vmodel_input)
        vmodel_out = run_dir / "model.mod"
        write_vmodel(vmodel, vmodel_out)

    elif stage_n == 5:

        # inputs
        run_n = config.STAGE_03.n_runs
        model_input = (
            Path("outputs/stage_03")
            / f"run_{run_n:02d}"
            / "summary"
            / "model.mod"
        )

        # outputs
        run_dir = Path("outputs/stage_05")

        # Read / Implement / Write
        copy(model_input, run_dir)


def finalized_inverted_models(config, stage_n, run_n):

    if stage_n == 1:
        root = Path("outputs/stage_01")
    elif stage_n == 3:
        root = Path("outputs/stage_03")

    run_dir = root / f"run_{run_n:02d}"

    # inputs
    omodel_files = sorted(run_dir.rglob("velocity_summary.csv"))

    # outputs
    averaged_modelfile = run_dir / "summary" / "averaged.csv"
    simplified_modelfile = run_dir / "summary" / "simplified.csv"
    final_modelfile = run_dir / "summary" / "model.mod"

    # processing
    omodels = [read_csv(i) for i in omodel_files]
    omodels = concat(omodels)

    merge_thr = config.SYNTHETIC_MODELS.merging_threshold_km_s
    simplify_model(
        omodels, averaged_modelfile, simplified_modelfile, merge_thr
    )
    simplified2vmodel(config, simplified_modelfile, final_modelfile)
