from src.model.Models import VELEST_CONFIG
from pathlib import Path
from pandas import read_csv
import yaml


def _bool(x: bool) -> int:
    """Convert Python bool to VELEST integer."""
    return int(bool(x))


def write_cmn(config: VELEST_CONFIG, outfile: Path):

    lines = []

    # ------------------------------------------------------------------
    # Line 1
    # ------------------------------------------------------------------

    lines.append(config.TITLE)

    # ------------------------------------------------------------------
    # Line 2
    # ------------------------------------------------------------------

    c = config.COORDINATE_SYSTEM

    lines.append(
        f"{c.origin_latitude:8.4f} "
        f"{c.origin_longitude:8.4f} "
        f"{c.system:d} "
        f"{c.elevation_shift:.2f} "
        f"{c.trial_mode:d} "
        f"{c.trial_depth:.2f} "
        f"{c.input_format:d}"
    )

    # ------------------------------------------------------------------
    # Line 3
    # ------------------------------------------------------------------

    d = config.DATASET

    lines.append(
        f"{d.number_of_earthquakes:d} "
        f"{d.number_of_shots:d} "
        f"{d.rotation:.2f}"
    )

    # ------------------------------------------------------------------
    # Line 4
    # ------------------------------------------------------------------

    inv = config.INVERSION

    lines.append(
        f"{_bool(inv.single_event)} " f"{_bool(inv.compute_resolution)}"
    )

    # ------------------------------------------------------------------
    # Line 5
    # ------------------------------------------------------------------

    m = config.RELOCATION

    lines.append(
        f"{m.maximum_distance:.2f} "
        f"{_bool(m.use_topography)} "
        f"{m.minimum_depth:.2f} "
        f"{m.velocity_adjustment:.2f} "
        f"{m.depth_adjustment:.2f} "
        f"{_bool(m.allow_low_velocity_layers)}"
    )

    # ------------------------------------------------------------------
    # Line 6
    # ------------------------------------------------------------------

    v = config.VELOCITY

    lines.append(
        f"{v.phase_mode:d} "
        f"{v.s_weight_factor:.2f} "
        f"{v.vp_vs_ratio:.3f} "
        f"{v.number_of_models:d}"
    )

    # ------------------------------------------------------------------
    # Line 7
    # ------------------------------------------------------------------

    d = config.DAMPING

    lines.append(
        f"{d.origin_time:.4f} "
        f"{d.hypocenter_xy:.4f} "
        f"{d.hypocenter_depth:.4f} "
        f"{d.velocity:.4f} "
        f"{d.station_correction:.4f}"
    )

    # ------------------------------------------------------------------
    # Line 8
    # ------------------------------------------------------------------

    s = config.STATION_CORRECTION

    lines.append(
        f"{_bool(s.invert_station_corrections)} "
        f"{_bool(s.use_shot_corrections)} "
        f"{_bool(s.fix_shot_corrections)} "
        f"{_bool(s.use_station_elevation)} "
        f"{_bool(s.use_station_corrections)}"
    )

    # ------------------------------------------------------------------
    # Line 9
    # ------------------------------------------------------------------

    o = config.OUTPUT_SWITCH

    lines.append(
        f"{_bool(o.turbo_mode)} "
        f"{_bool(o.write_cnv)} "
        f"{_bool(o.write_station_file)} "
        f"{_bool(o.write_summary_file)}"
    )

    # ------------------------------------------------------------------
    # Line 10
    # ------------------------------------------------------------------

    d = config.DIAGNOSTIC

    lines.append(
        f"{_bool(d.raypaths)} "
        f"{_bool(d.derivatives)} "
        f"{_bool(d.averaging_kernel)} "
        f"{_bool(d.dirichlet_spread)} "
        f"{_bool(d.reflection_points)} "
        f"{_bool(d.refraction_points)} "
        f"{_bool(d.residuals)}"
    )

    # ------------------------------------------------------------------
    # Line 11
    # ------------------------------------------------------------------

    i = config.ITERATION

    lines.append(
        f"{i.minimum_rms_change:.5f} "
        f"{i.maximum_iterations:d} "
        f"{i.invert_ratio_every:d}"
    )

    # ------------------------------------------------------------------
    # Line 12
    # ------------------------------------------------------------------

    q = config.QUALITY_CLASS

    lines.append(
        f"{q.mode:d} "
        f"{q.class0:d} "
        f"{q.class1:d} "
        f"{q.class2:d} "
        f"{q.class3:d} "
        f"{q.class4:d} "
        f"{q.class5:d}"
    )

    # ------------------------------------------------------------------
    # Line 13
    # ------------------------------------------------------------------

    q = config.QUALITY_FACTOR

    lines.append(
        f"{q.mode:d} "
        f"{q.class0:.3f} "
        f"{q.class1:.3f} "
        f"{q.class2:.3f} "
        f"{q.class3:.3f} "
        f"{q.class4:.3f}"
    )

    # ------------------------------------------------------------------
    # Lines 14-22
    # ------------------------------------------------------------------

    inp = config.INPUT_FILE

    lines.extend(
        [
            str(inp.velocity_model),
            str(inp.stations),
            str(inp.seismo),
            str(inp.region_names),
            str(inp.region_coordinates),
            str(inp.topography_primary),
            str(inp.topography_secondary),
            str(inp.earthquakes),
            str(inp.shots),
        ]
    )

    # ------------------------------------------------------------------
    # Lines 23-34
    # ------------------------------------------------------------------

    out = config.OUTPUT_FILE

    lines.extend(
        [
            str(out.main_output),
            str(out.single_event_locations),
            str(out.cnv_output),
            str(out.station_output),
            str(out.summary_output),
            str(out.ray_output),
            str(out.derivative_output),
            str(out.averaging_kernel_output),
            str(out.dirichlet_output),
            str(out.reflection_output),
            str(out.refraction_output),
            str(out.residual_output),
        ]
    )

    with open(outfile, "w") as fp:
        fp.write("\n".join(lines))
        fp.write("\n")


def prepare_cmn(config, stage_n, run_n=None):

    if stage_n == 1:
        
        root = Path("outputs/stage_01")

        # inputs
        events_input = root / "select.csv"
        cmn_input = Path("configs") / "stage_01" / "velest.yaml"

        # outputs
        run_dir = root / f"run_{run_n:02d}"

        # process
        events_df = read_csv(events_input)
        for n_model in range(1, config.SYNTHETIC_MODELS.n_models + 1):
            model_dir = run_dir / f"model_{n_model:02d}"
            cmnout = model_dir / "velest.cmn"
            with open(cmn_input) as f:
                data = yaml.safe_load(f)
            VelestConfig = VELEST_CONFIG.model_validate(data)
            VelestConfig.DATASET.number_of_earthquakes = len(events_df)
            VelestConfig.INPUT_FILE.velocity_model = model_dir / "model.mod"
            VelestConfig.INPUT_FILE.stations = root / "stations.sta"
            VelestConfig.INPUT_FILE.earthquakes = root / "select.cnv"
            VelestConfig.OUTPUT_FILE.main_output = model_dir / "velout.inv"
            VelestConfig.OUTPUT_FILE.cnv_output = model_dir / "velout.cnv"
            VelestConfig.OUTPUT_FILE.station_output = model_dir / "velout.sta"
            VelestConfig.OUTPUT_FILE.summary_output = model_dir / "velout.smp"
            VelestConfig.OUTPUT_FILE.ray_output = model_dir / "velout.ray"
            VelestConfig.OUTPUT_FILE.residual_output = model_dir / "velout.res"
            write_cmn(VelestConfig, cmnout)

    if stage_n == 2:
        
        rootout = Path("outputs/stage_02")
        
        # inputs
        events_input = rootout / "original.csv"
        cmn_input = Path("configs") / "stage_02" / "velest.yaml"

        # outputs
        cmnout = rootout / "velest.cmn"
        events_df = read_csv(events_input)
        with open(cmn_input) as f:
            data = yaml.safe_load(f)
        VelestConfig = VELEST_CONFIG.model_validate(data)
        VelestConfig.DATASET.number_of_earthquakes = len(events_df)
        VelestConfig.INPUT_FILE.velocity_model = rootout / "model.mod"
        VelestConfig.INPUT_FILE.stations = rootout / "stations.sta"
        VelestConfig.INPUT_FILE.earthquakes = rootout / "original.cnv"
        VelestConfig.OUTPUT_FILE.main_output = rootout / "velout.inv"
        VelestConfig.OUTPUT_FILE.single_event_locations = rootout / "single_events.out"
        VelestConfig.OUTPUT_FILE.cnv_output = rootout / "velout.cnv"
        VelestConfig.OUTPUT_FILE.station_output = rootout / "velout.sta"
        VelestConfig.OUTPUT_FILE.summary_output = rootout / "velout.sum"
        write_cmn(VelestConfig, cmnout)

    if stage_n == 3:
        
        root = Path("outputs/stage_03")

        # inputs
        
        cmn_input = Path("configs") / "stage_03" / "velest.yaml"

        # outputs

        run_dir = root / f"run_{run_n:02d}"

        for n_model in range(1, config.STAGE_03.n_models + 1):
            model_dir = run_dir / f"model_{n_model:02d}"
            cmnout = model_dir / "velest.cmn"
            events_input =  Path("outputs/stage_02") / "select.csv"
            events_df = read_csv(events_input)
            with open(cmn_input) as f:
                data = yaml.safe_load(f)
            VelestConfig = VELEST_CONFIG.model_validate(data)
            VelestConfig.DATASET.number_of_earthquakes = len(events_df)
            VelestConfig.INPUT_FILE.velocity_model = model_dir / "model.mod"
            VelestConfig.INPUT_FILE.stations = model_dir / "stations.sta"
            VelestConfig.INPUT_FILE.earthquakes = model_dir / "select.cnv"
            VelestConfig.OUTPUT_FILE.main_output = model_dir / "velout.inv"
            VelestConfig.OUTPUT_FILE.cnv_output = model_dir / "velout.cnv"
            VelestConfig.OUTPUT_FILE.station_output = model_dir / "velout.sta"
            VelestConfig.OUTPUT_FILE.summary_output = model_dir / "velout.smp"
            VelestConfig.OUTPUT_FILE.ray_output = model_dir / "velout.ray"
            VelestConfig.OUTPUT_FILE.residual_output = model_dir / "velout.res"
            write_cmn(VelestConfig, cmnout)
            
    if stage_n == 4:
        
        root = Path("outputs/stage_04")

        # inputs
        
        cmn_input = Path("configs") / "stage_04" / "velest.yaml"

        # outputs

        run_dir = root / f"run_{run_n:02d}"

        cmnout = run_dir / "velest.cmn"
        events_input =  run_dir / "perturb.csv"
        events_df = read_csv(events_input)
        with open(cmn_input) as f:
            data = yaml.safe_load(f)
        VelestConfig = VELEST_CONFIG.model_validate(data)
        VelestConfig.DATASET.number_of_earthquakes = len(events_df)
        VelestConfig.INPUT_FILE.velocity_model = run_dir / "model.mod"
        VelestConfig.INPUT_FILE.stations = run_dir / "stations.sta"
        VelestConfig.INPUT_FILE.earthquakes = run_dir / "perturb.cnv"
        VelestConfig.OUTPUT_FILE.main_output = run_dir / "velout.inv"
        VelestConfig.OUTPUT_FILE.cnv_output = run_dir / "velout.cnv"
        VelestConfig.OUTPUT_FILE.station_output = run_dir / "velout.sta"
        VelestConfig.OUTPUT_FILE.summary_output = run_dir / "velout.smp"
        VelestConfig.OUTPUT_FILE.ray_output = run_dir / "velout.ray"
        VelestConfig.OUTPUT_FILE.residual_output = run_dir / "velout.res"
        write_cmn(VelestConfig, cmnout)            
        
    if stage_n == 5:

        rootout = Path("outputs/stage_05")
        
        # inputs
        events_input = rootout / "original.csv"
        cmn_input = Path("configs") / "stage_05" / "velest.yaml"

        # outputs
        cmnout = rootout / "velest.cmn"
        events_df = read_csv(events_input)
        with open(cmn_input) as f:
            data = yaml.safe_load(f)
        VelestConfig = VELEST_CONFIG.model_validate(data)
        VelestConfig.DATASET.number_of_earthquakes = len(events_df)
        VelestConfig.INPUT_FILE.velocity_model = rootout / "model.mod"
        VelestConfig.INPUT_FILE.stations = rootout / "stations.sta"
        VelestConfig.INPUT_FILE.earthquakes = rootout / "original.cnv"
        VelestConfig.OUTPUT_FILE.main_output = rootout / "velout.inv"
        VelestConfig.OUTPUT_FILE.single_event_locations = rootout / "single_events.out"
        VelestConfig.OUTPUT_FILE.cnv_output = rootout / "velout.cnv"
        VelestConfig.OUTPUT_FILE.station_output = rootout / "velout.sta"
        VelestConfig.OUTPUT_FILE.summary_output = rootout / "velout.sum"
        write_cmn(VelestConfig, cmnout)