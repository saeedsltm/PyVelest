from __future__ import annotations

from typing import Optional
from pydantic import BaseModel, Field, model_validator
from pathlib import Path
from numpy import random, asarray, all, diff, any, empty, cumsum, sum, isclose

# ============================================================================== USER CONFIG


class INPUTS_(BaseModel):
    """Quality thresholds used to select catalog events."""

    catalog: str = Field(default=Path("inputs/select.out"))
    stations: str = Field(default=Path("inputs/stations.csv"))
    vmodel: str = Field(default=Path("inputs/vmodel.dat"))


class STUDY_AREA_(BaseModel):
    """Geographic bounds for events retained in a study."""

    min_lat: float = Field(ge=-90.0, le=90.0)
    max_lat: float = Field(ge=-90.0, le=90.0)
    min_lon: float = Field(ge=-180.0, le=180.0)
    max_lon: float = Field(ge=-180.0, le=180.0)
    lat_hemisphere: str = "N"
    lon_hemisphere: str = "E"

    @model_validator(mode="after")
    def validate_bounds(self) -> "STUDY_AREA_":
        if self.min_lat > self.max_lat:
            raise ValueError("min_lat must be less than or equal to max_lat")
        if self.min_lon > self.max_lon:
            raise ValueError("min_lon must be less than or equal to max_lon")
        return self

    def contains(self, lat: float, lon: float) -> bool:
        """Return True when a coordinate falls inside the region."""

        return (
            self.min_lat <= lat <= self.max_lat
            and self.min_lon <= lon <= self.max_lon
        )


class EVENT_SELECTION_(BaseModel):
    """Quality thresholds used to select catalog events."""

    max_events: Optional[int] = Field(default=100, ge=1)
    cell_count_lat: int = Field(default=4, ge=1)
    cell_count_lon: int = Field(default=4, ge=1)
    cell_count_depth: int = Field(default=4, ge=1)
    max_horizontal_error_km: Optional[float] = Field(default=None, ge=0.0)
    max_depth_error_km: Optional[float] = Field(default=None, ge=0.0)
    max_rms_sec: Optional[float] = Field(default=None, ge=0.0)
    max_gap_deg: Optional[float] = Field(default=None, ge=0.0)
    min_used_stations: Optional[int] = Field(default=None, ge=0)


class REFERENCE_STATION_(BaseModel):
    """Weights used to rank candidate reference stations."""

    phase_count_weight: float = Field(default=0.45, ge=0.0)
    residual_weight: float = Field(default=0.35, ge=0.0)
    distance_weight: float = Field(default=0.20, ge=0.0)


class SYNTHETIC_MODELS_(BaseModel):
    """Settings for generating a simple synthetic 1D velocity model."""

    n_models: int = Field(default=12, ge=1)
    velocity_std: float = Field(default=0.2, ge=0.0)
    depth_std: float = Field(default=2.0, ge=0.0)
    seed: int = Field(default=1, ge=0.0)
    layer_count: int = Field(default=12, ge=2)
    merging_threshold_km_s: float = Field(default=0.2, ge=0.0)
    first_layer_dmp: int = Field(default=50, ge=0)
    last_layer_dmp: int = Field(default=10, ge=0)

class STAGE_01_(BaseModel):
    """Settings for stage_01"""

    n_runs: int = Field(default=5, ge=0)


class STAGE_02_(BaseModel):
    """Settings for stage_02"""

    n_runs: int = Field(default=5, ge=0)


class STAGE_03_(BaseModel):
    """Settings for stage_03"""

    n_runs: int = Field(default=5, ge=0)
    n_models: int = Field(default=10, ge=0)


class STAGE_04_(BaseModel):
    """Settings for stage_04"""

    n_runs: int = Field(default=5, ge=0)
    dx_km: float = Field(default=2.0, ge=0)
    dy_km: float = Field(default=2.0, ge=0)
    dz_km: float = Field(default=2.0, ge=0)
    seed: int = Field(default=1, ge=0.0)


class REPORT_SETTINGS_(BaseModel):
    """Settings for generating a simple report."""

    report_on: bool = Field(default=False)


class FIGURE_SETTINGS_(BaseModel):
    """Settings for generating figures."""

    min_sei_lat: float = Field(ge=-90.0, le=90.0)
    max_sei_lat: float = Field(ge=-90.0, le=90.0)
    min_sei_lon: float = Field(ge=-180.0, le=180.0)
    max_sei_lon: float = Field(ge=-180.0, le=180.0)
    min_sei_depth: float = Field(default=0.0, ge=0.0)
    max_sei_depth: float = Field(default=20.0, ge=0.0)
    min_vel_vp: float = Field(default=4.0, ge=2.0)
    max_vel_vp: float = Field(default=8.5, le=9.0)
    min_vel_depth: float = Field(default=0.0, ge=0.0)
    max_vel_depth: float = Field(default=45.0, le=99.9)


class GLOBAL_SETTINGS_(BaseModel):
    """Settings for global"""

    clean_workspace: bool = Field(default=True)


class USER_CONFIG(BaseModel):
    """Top-level PyVELEST preprocessing configuration."""

    INPUTS: INPUTS_
    STUDY_AREA: STUDY_AREA_
    EVENT_SELECTION: EVENT_SELECTION_
    REFERENCE_STATION: REFERENCE_STATION_
    SYNTHETIC_MODELS: SYNTHETIC_MODELS_
    STAGE_01: STAGE_01_
    STAGE_02: STAGE_02_
    STAGE_03: STAGE_03_
    STAGE_04: STAGE_04_
    REPORT_SETTINGS: REPORT_SETTINGS_
    FIGURE_SETTINGS: FIGURE_SETTINGS_
    GLOBAL_SETTINGS: GLOBAL_SETTINGS_


# ============================================================================== VELOCITY MODEL
class LAYER(BaseModel):
    """a VELEST velocity model layer including vel,depth,vdamp"""

    vel: float = Field(ge=1.0, le=8.5)
    depth: float = Field(ge=-5.0, le=60.0)
    vdamp: float = Field(default=1.0, ge=0.0)


class VMODEL(BaseModel):
    """VELEST velocity model"""

    model_name: Optional[str] = Field(default="VELEST MODEL")
    n_layers: int = Field(ge=2)
    layers_vp: list[LAYER] = Field(default_factory=list)
    layers_vs: Optional[list[LAYER]] = None

    def generate_random_models(
        self,
        n_models,
        velocity_std,
        depth_std,
        seed=None,
        keep_boundaries=True,
        min_velocity_increment=0.01,
        max_attempts=10000,
    ):
        """
        Generate random 1-D velocity models around the current model.

        The current model is treated as the mean/reference model.

        Velocity perturbations are applied to velocity increments
        (velocity gradients between adjacent layers), rather than to
        absolute velocities. This guarantees that the generated model
        remains monotonically increasing with depth.

        Both layer depths and velocities are perturbed.

        Vs is optional. If the reference model does not contain Vs layers,
        only Vp is randomized and the generated models will also contain
        no Vs model.
        """

        # ---------------------------------------------------------
        # Validation
        # ---------------------------------------------------------

        if n_models <= 0:
            raise ValueError("n_models must be greater than zero.")

        if velocity_std < 0:
            raise ValueError("velocity_std must be >= 0.")

        if depth_std < 0:
            raise ValueError("depth_std must be >= 0.")

        if min_velocity_increment < 0:
            raise ValueError("min_velocity_increment must be >= 0.")

        # ---------------------------------------------------------
        # Random number generator
        # ---------------------------------------------------------

        rng = random.default_rng(seed)

        # ---------------------------------------------------------
        # Extract reference model
        # ---------------------------------------------------------

        vp0 = asarray(
            [layer.vel for layer in self.layers_vp],
            dtype=float,
        )

        depth0 = asarray(
            [layer.depth for layer in self.layers_vp],
            dtype=float,
        )

        # Vs is optional
        has_vs = bool(self.layers_vs)

        if has_vs:
            vs0 = asarray(
                [layer.vel for layer in self.layers_vs],
                dtype=float,
            )
        else:
            vs0 = None

        n_layers = len(vp0)

        # ---------------------------------------------------------
        # Basic consistency checks
        # ---------------------------------------------------------

        if len(depth0) != n_layers:
            raise ValueError(
                "Velocity and depth arrays must have the same length."
            )

        if n_layers < 2:
            raise ValueError("At least two layers are required.")

        if not all(diff(depth0) > 0):
            raise ValueError(
                "Reference model depths must be strictly increasing."
            )

        if not all(diff(vp0) >= 0):
            raise ValueError(
                "Reference Vp model must be monotonically increasing."
            )

        # ---------------------------------------------------------
        # Vs consistency checks -- only if Vs exists
        # ---------------------------------------------------------

        if has_vs:

            if len(vs0) != n_layers:
                raise ValueError(
                    "Vp and Vs models must have the same number of layers."
                )

            if not all(diff(vs0) >= 0):
                raise ValueError(
                    "Reference Vs model must be monotonically increasing."
                )

        # ---------------------------------------------------------
        # Reference velocity increments
        # ---------------------------------------------------------

        dvp0 = diff(vp0)

        if has_vs:
            dvs0 = diff(vs0)
        else:
            dvs0 = None

        # ---------------------------------------------------------
        # Generate models
        # ---------------------------------------------------------

        models = []

        attempts = 0

        while len(models) < n_models:

            attempts += 1

            if attempts > n_models * max_attempts:
                raise RuntimeError(
                    "Unable to generate the requested number of "
                    "valid random velocity models. "
                    "Try reducing velocity_std/depth_std or "
                    "reducing min_velocity_increment."
                )

            # =====================================================
            # 1. Perturb layer depths
            # =====================================================

            depth = depth0.copy()

            if depth_std > 0:

                if keep_boundaries:

                    perturbation = rng.normal(
                        loc=0.0,
                        scale=depth_std,
                        size=n_layers - 2,
                    )

                    depth[1:-1] += perturbation

                else:

                    perturbation = rng.normal(
                        loc=0.0,
                        scale=depth_std,
                        size=n_layers,
                    )

                    depth += perturbation

            # -----------------------------------------------------
            # Make sure depths are strictly increasing
            # -----------------------------------------------------

            if not all(diff(depth) > 0):
                continue

            # -----------------------------------------------------
            # Keep boundaries exactly fixed
            # -----------------------------------------------------

            if keep_boundaries:

                depth[0] = depth0[0]
                depth[-1] = depth0[-1]

            # =====================================================
            # 2. Perturb Vp velocity increments
            # =====================================================

            if velocity_std > 0:

                dvp = dvp0 + rng.normal(
                    loc=0.0,
                    scale=velocity_std,
                    size=n_layers - 1,
                )

            else:

                dvp = dvp0.copy()

            # -----------------------------------------------------
            # All Vp velocity increments must be positive
            # -----------------------------------------------------

            if any(dvp < min_velocity_increment):
                continue

            # =====================================================
            # 3. Perturb Vs velocity increments
            # =====================================================

            if has_vs:

                if velocity_std > 0:

                    dvs = dvs0 + rng.normal(
                        loc=0.0,
                        scale=velocity_std,
                        size=n_layers - 1,
                    )

                else:

                    dvs = dvs0.copy()

                # -------------------------------------------------
                # All Vs velocity increments must be positive
                # -------------------------------------------------

                if any(dvs < min_velocity_increment):
                    continue

            # =====================================================
            # 4. Reconstruct Vp
            # =====================================================

            vp = empty(n_layers)

            vp[0] = vp0[0]

            vp[1:] = vp[0] + cumsum(dvp)

            # =====================================================
            # 5. Reconstruct Vs -- only if Vs exists
            # =====================================================

            if has_vs:

                vs = empty(n_layers)

                vs[0] = vs0[0]

                vs[1:] = vs[0] + cumsum(dvs)

            # =====================================================
            # 6. Keep boundary velocities fixed
            # =====================================================

            if keep_boundaries:

                # -------------------------------------------------
                # Normalize Vp increments
                # -------------------------------------------------

                total_dvp0 = sum(dvp0)
                total_dvp = sum(dvp)

                if total_dvp <= 0:
                    continue

                dvp = dvp * (total_dvp0 / total_dvp)

                if any(dvp < min_velocity_increment):
                    continue

                vp = empty(n_layers)

                vp[0] = vp0[0]
                vp[1:] = vp[0] + cumsum(dvp)

                # -------------------------------------------------
                # Normalize Vs increments only if Vs exists
                # -------------------------------------------------

                if has_vs:

                    total_dvs0 = sum(dvs0)
                    total_dvs = sum(dvs)

                    if total_dvs <= 0:
                        continue

                    dvs = dvs * (total_dvs0 / total_dvs)

                    if any(dvs < min_velocity_increment):
                        continue

                    vs = empty(n_layers)

                    vs[0] = vs0[0]
                    vs[1:] = vs[0] + cumsum(dvs)

                # -------------------------------------------------
                # Explicitly enforce boundaries
                # -------------------------------------------------

                vp[-1] = vp0[-1]

                if has_vs:
                    vs[-1] = vs0[-1]

            # =====================================================
            # 7. Final validity checks
            # =====================================================

            if not all(diff(depth) > 0):
                continue

            if not all(diff(vp) >= min_velocity_increment):
                continue

            if has_vs:
                if not all(diff(vs) >= min_velocity_increment):
                    continue

            # -----------------------------------------------------
            # Check boundaries
            # -----------------------------------------------------

            if keep_boundaries:

                if not isclose(depth[0], depth0[0]):
                    continue

                if not isclose(depth[-1], depth0[-1]):
                    continue

                if not isclose(vp[0], vp0[0]):
                    continue

                if not isclose(vp[-1], vp0[-1]):
                    continue

                if has_vs:

                    if not isclose(vs[0], vs0[0]):
                        continue

                    if not isclose(vs[-1], vs0[-1]):
                        continue

            # =====================================================
            # 8. Construct LAYER objects
            # =====================================================

            layers_vp = []

            for i in range(n_layers):

                layers_vp.append(
                    LAYER(
                        vel=float(vp[i]),
                        depth=float(depth[i]),
                        vdamp=self.layers_vp[i].vdamp,
                    )
                )

            # -----------------------------------------------------
            # Construct Vs layers only if Vs exists
            # -----------------------------------------------------

            if has_vs:

                layers_vs = []

                for i in range(n_layers):

                    layers_vs.append(
                        LAYER(
                            vel=float(vs[i]),
                            depth=float(depth[i]),
                            vdamp=self.layers_vs[i].vdamp,
                        )
                    )

            else:

                layers_vs = None

            # =====================================================
            # 9. Construct VMODEL
            # =====================================================

            model = VMODEL(
                model_name=f"{self.model_name}_R{len(models) + 1:04d}",
                layers_vp=layers_vp,
                layers_vs=layers_vs,
                n_layers=n_layers,
            )

            models.append(model)

        return models


# ============================================================================== STATION MODEL
class STATION(BaseModel):
    """a VELEST station info including code,lat,lon,elv,icc,ptcorr,stcorr"""

    code: str = Field(max_length=4)
    lat: float = Field(ge=-90.0, le=90.0)
    lat_hem: str = Field(max_length=1)
    lon: float = Field(ge=-180.0, le=180.0)
    lon_hem: str = Field(max_length=1)
    elv: float = Field(ge=0.0)
    mn: int = Field(ge=1)
    icc: int = Field(ge=0)
    ptcorr: float = Field(ge=-9.9, le=9.9)
    stcorr: float = Field(ge=-9.9, le=9.9)


class STATIONS(BaseModel):
    """VELEST stations"""

    n_stations: int = Field(ge=3)
    stations: list[STATION] = Field(default_factory=list)


# ============================================================================== VELEST CONFIG

# =============================================================================
# Coordinate system (line 2)
# =============================================================================


class COORDINATE_SYSTEM_(BaseModel):
    origin_latitude: float = Field(ge=-90.0, le=90.0)  # olat
    origin_longitude: float = Field(ge=-180.0, le=180.0)  # olon
    system: int = Field(ge=0, le=2)  # icoordsystem
    elevation_shift: float = Field(ge=0.0)  # zshift
    trial_mode: int = Field(ge=0)  # itrial
    trial_depth: float = Field(ge=0.0, le=99.9)  # ztrial
    input_format: int = Field(ge=0, le=2)  # ised


# =============================================================================
# Dataset (line 3)
# =============================================================================


class DATASET_(BaseModel):
    number_of_earthquakes: Optional[int]  # neqs
    number_of_shots: Optional[int]  # nshot
    rotation: Optional[float] = Field(default=0.0, ge=0.0)  # rotate


# =============================================================================
# Inversion mode (line 4)
# =============================================================================


class INVERSION_(BaseModel):
    single_event: bool = False  # isingle
    compute_resolution: bool = False  # iresolcalc


# =============================================================================
# Relocation constraints (line 5)
# =============================================================================


class RELOCATION_(BaseModel):
    maximum_distance: float = Field(ge=0, le=1000)  # dmax
    use_topography: bool = False  # itopo
    minimum_depth: float = Field(ge=-5.0, le=50.0)  # zmininput
    velocity_adjustment: float = Field(ge=0.0, le=2.0)  # veladj
    depth_adjustment: float = Field(ge=0.0, le=5.0)  # zadj
    allow_low_velocity_layers: bool = False  # lowveloclay


# =============================================================================
# Velocity model (line 6)
# =============================================================================


class VELOCITY_(BaseModel):
    phase_mode: int = Field(ge=1, le=3)  # nsp
    s_weight_factor: float = 1.0  # swtfac
    vp_vs_ratio: float = 1.73  # vpvs
    number_of_models: int = Field(ge=1, le=2)  # nmod


# =============================================================================
# Damping (line 7)
# =============================================================================


class DAMPING_(BaseModel):
    origin_time: float = Field(ge=0.01, le=100)  # othet
    hypocenter_xy: float = Field(ge=0.01, le=100)  # xythet
    hypocenter_depth: float = Field(ge=0.01, le=100)  # zthet
    velocity: float = Field(ge=0.01, le=100)  # vthet
    station_correction: float = Field(ge=0.01, le=100)  # stathet


# =============================================================================
# Station corrections (line 8)
# =============================================================================


class STATION_CORRECTION_(BaseModel):
    invert_station_corrections: bool = True  # nsinv
    use_shot_corrections: bool = False  # nshcor
    fix_shot_corrections: bool = False  # nshfix
    use_station_elevation: bool = True  # iuseelev
    use_station_corrections: bool = True  # iusestacorr


# =============================================================================
# Output switches (line 9)
# =============================================================================


class OUTPUT_SWITCH_(BaseModel):
    turbo_mode: bool = False  # iturbo
    write_cnv: bool = True  # icnvout
    write_station_file: bool = True  # istaout
    write_summary_file: bool = False  # ismpout


# =============================================================================
# Diagnostic outputs (line 10)
# =============================================================================


class DIAGNOSTIC_(BaseModel):
    raypaths: bool = False  # irayout
    derivatives: bool = False  # idrvout
    averaging_kernel: bool = False  # ialeout
    dirichlet_spread: bool = False  # idspout
    reflection_points: bool = False  # irflout
    refraction_points: bool = False  # irfrout
    residuals: bool = False  # iresout


# =============================================================================
# Iteration control (line 11)
# =============================================================================


class ITERATION_(BaseModel):
    minimum_rms_change: float = Field(ge=0.001, le=10.0)  # delmin
    maximum_iterations: int = Field(ge=1, le=100)  # ittmax
    invert_ratio_every: int = Field(ge=0, le=1000)  # inverterratio


# =============================================================================
# Quality classes (line 12)
# =============================================================================


class QUALITY_CLASS_(BaseModel):
    mode: int = 0
    class0: int = 0
    class1: int = 1
    class2: int = 2
    class3: int = 3
    class4: int = 4
    class5: int = 4


# =============================================================================
# Quality factors (line 13)
# =============================================================================


class QUALITY_FACTOR_(BaseModel):
    mode: int = 0
    class0: float = 1.00
    class1: float = 0.50
    class2: float = 0.25
    class3: float = 0.125
    class4: float = 0.00


# =============================================================================
# Input files (lines 14-22)
# =============================================================================


class INPUT_FILE_(BaseModel):
    velocity_model: Optional[Path] | str
    stations: Optional[Path] | str
    seismo: Optional[Path] | str
    region_names: Optional[Path] | str
    region_coordinates: Optional[Path] | str
    topography_primary: Optional[Path] | str
    topography_secondary: Optional[Path] | str
    earthquakes: Optional[Path] | str
    shots: Optional[Path] | str


# =============================================================================
# Output files (lines 23-34)
# =============================================================================


class OUTPUT_FILE_(BaseModel):
    main_output: Optional[Path] | str
    single_event_locations: Optional[Path] | str
    cnv_output: Optional[Path] | str
    station_output: Optional[Path] | str
    summary_output: Optional[Path] | str
    ray_output: Optional[Path] | str
    derivative_output: Optional[Path] | str
    averaging_kernel_output: Optional[Path] | str
    dirichlet_output: Optional[Path] | str
    reflection_output: Optional[Path] | str
    refraction_output: Optional[Path] | str
    residual_output: Optional[Path] | str


# =============================================================================
# Root configuration
# =============================================================================


class VELEST_CONFIG(BaseModel):
    TITLE: str = "VELEST inversion"
    COORDINATE_SYSTEM: COORDINATE_SYSTEM_ = Field(
        default_factory=COORDINATE_SYSTEM_
    )
    DATASET: DATASET_ = Field(default_factory=DATASET_)
    INVERSION: INVERSION_ = Field(default_factory=INVERSION_)
    RELOCATION: RELOCATION_ = Field(default_factory=RELOCATION_)
    VELOCITY: VELOCITY_ = Field(default_factory=VELOCITY_)
    DAMPING: DAMPING_ = Field(default_factory=DAMPING_)
    STATION_CORRECTION: STATION_CORRECTION_ = Field(
        default_factory=STATION_CORRECTION_
    )
    OUTPUT_SWITCH: OUTPUT_SWITCH_ = Field(default_factory=OUTPUT_SWITCH_)
    DIAGNOSTIC: DIAGNOSTIC_ = Field(default_factory=DIAGNOSTIC_)
    ITERATION: ITERATION_ = Field(default_factory=ITERATION_)
    QUALITY_CLASS: QUALITY_CLASS_ = Field(default_factory=QUALITY_CLASS_)
    QUALITY_FACTOR: QUALITY_FACTOR_ = Field(default_factory=QUALITY_FACTOR_)
    INPUT_FILE: INPUT_FILE_ = Field(default_factory=INPUT_FILE_)
    OUTPUT_FILE: OUTPUT_FILE_ = Field(default_factory=OUTPUT_FILE_)
