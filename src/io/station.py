from __future__ import annotations


from pandas import read_csv, DataFrame, concat, to_numeric, notna
from pathlib import Path
from src.model.Models import USER_CONFIG
from src.model.Models import STATION, STATIONS

from numpy import nan, sqrt, mean, arange
from shutil import copy


def read_station(config: USER_CONFIG, infile: Path | str):
    min_lat = config.STUDY_AREA.min_lat
    max_lat = config.STUDY_AREA.max_lat
    min_lon = config.STUDY_AREA.min_lon
    max_lon = config.STUDY_AREA.max_lon

    station_df = read_csv(infile)

    mask = station_df.lat.between(min_lat, max_lat) & station_df.lon.between(
        min_lon, max_lon
    )

    filtered_df = station_df[mask]

    return filtered_df


def read_vstation(config: USER_CONFIG, infile: Path | str):
    vstations = STATIONS.model_construct()
    with open(infile) as fo:
        for n, line in enumerate(fo):
            if n > 0 and line.strip():
                data = {
                    "code": line[:4].strip(),
                    "lat": float(line[4:11]),
                    "lat_hem": line[11],
                    "lon": float(line[12:21]),
                    "lon_hem": line[21],
                    "elv": float(line[23:27]),
                    "mn": int(line[28]),
                    "icc": int(line[30:33]),
                    "ptcorr": float(line[34:39]),
                    "stcorr": float(line[41:46]),
                }
                station = STATION.model_validate(data)
                vstations.stations.append(station)
    return vstations


def calculate_average_corr(stations_df, outfile):
    stations_av_df = stations_df.groupby("code", as_index=False).agg(
        {
            "lat": "first",
            "lat_hem": "first",
            "lon": "first",
            "lon_hem": "first",
            "elv": "first",
            "mn": "first",
            "icc": "first",
            "ptcorr": "mean",
            "stcorr": "mean",
        }
    )
    stations_av_df = stations_av_df.sort_values(
        "icc",
        ascending=False,
        ignore_index=True,
    )

    fmt = "(a4,f7.4,a1,1x,f8.4,a1,1x," "i4,1x,i1,1x,i3,1x,f5.2,2x,f5.2,3x,i1)"

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    with outfile.open("w") as f:

        # First line is the FORTRAN format string.
        f.write(fmt + "\n")

        for _, row in stations_av_df.iterrows():

            code = row.code
            lat = float(row["lat"])
            lon = float(row["lon"])
            elv = int(round(float(row["elv"])))

            # VELEST expects positive coordinates + hemisphere character.
            lat_abs = abs(lat)
            lon_abs = abs(lon)

            lat_hemisphere = row.lat_hem
            lon_hemisphere = row.lon_hem

            # One velocity model.
            model = row.mn

            # Station delays.
            ptcor = float(row.ptcorr)
            stcor = float(row.stcorr)

            icc = int(row["icc"])

            line = (
                f"{code:<4}"
                f"{lat_abs:7.4f}"
                f"{lat_hemisphere}"
                f" "
                f"{lon_abs:8.4f}"
                f"{lon_hemisphere}"
                f" "
                f"{elv:4d}"
                f" "
                f"{model:1d}"
                f" "
                f"{icc:3d}"
                f" "
                f"{ptcor:5.2f}"
                f"  "
                f"{stcor:5.2f}"
                f"   "
                f"{model:1d}"
            )

            f.write(line + "\n")
        f.write("\n")


def calculate_reference_station_scores(
    df: DataFrame,
    observation_weight: float = 0.50,
    residual_weight: float = 0.35,
    distance_weight: float = 0.15,
) -> DataFrame:
    """
    Calculate a reference-station score for every station.

    The input dataframe is in long format, with one row per phase
    observation.

    Expected columns
    ----------------
    eid : event identifier
    sta : station code
    dis : epicentral distance (km)
    azm : azimuth (degrees)
    pha : phase code, e.g. "PG" or "SG"
    wet : observation weight
    ptt : observed/calculated travel time (s)
    res : travel-time residual (s)

    Event-summary rows have an empty ``pha`` and are ignored when
    calculating station statistics.

    Parameters
    ----------
    df
        Long-format phase observation dataframe.

    observation_weight
        Weight assigned to the number of observations / event coverage.

    residual_weight
        Weight assigned to the travel-time residual.

    distance_weight
        Weight assigned to the average epicentral distance.

    Returns
    -------
    pandas.DataFrame
        Columns:

        code
        nP
        nS
        nObs
        coverage
        median_abs_residual
        mean_abs_residual
        rms_residual
        mean_distance
        observation_score
        residual_score
        distance_score
        score
        icc

        ``icc`` is the rank of the station according to the final
        reference-station score. The best station has the highest
        score but the lowest ICC.
    """

    # ------------------------------------------------------------
    # Validate required columns
    # ------------------------------------------------------------

    required_columns = {
        "eid",
        "sta",
        "dis",
        "pha",
        "res",
    }

    missing = required_columns.difference(df.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    # ------------------------------------------------------------
    # Keep only phase-observation rows
    # ------------------------------------------------------------
    #
    # Event-summary rows have pha == NaN.
    #
    # We also require a station code because an observation without
    # a station cannot contribute to station statistics.
    # ------------------------------------------------------------

    observations = df.loc[df["pha"].notna() & df["sta"].notna()].copy()

    if observations.empty:
        return DataFrame(
            columns=[
                "code",
                "nP",
                "nS",
                "nObs",
                "coverage",
                "median_abs_residual",
                "mean_abs_residual",
                "rms_residual",
                "mean_distance",
                "observation_score",
                "residual_score",
                "distance_score",
                "score",
                "icc",
            ]
        )

    # ------------------------------------------------------------
    # Normalize data types
    # ------------------------------------------------------------

    observations["pha"] = observations["pha"].astype(str).str.upper()

    observations["res"] = to_numeric(
        observations["res"],
        errors="coerce",
    )

    observations["dis"] = to_numeric(
        observations["dis"],
        errors="coerce",
    )

    # ------------------------------------------------------------
    # Discover stations
    # ------------------------------------------------------------

    stations = sorted(observations["sta"].dropna().astype(str).unique())

    # Number of events
    n_events = observations["eid"].nunique()

    if n_events == 0:
        raise ValueError("No valid event IDs found in dataframe.")

    rows = []

    # ------------------------------------------------------------
    # Calculate station statistics
    # ------------------------------------------------------------

    for sta in stations:

        station = observations.loc[observations["sta"].astype(str) == sta]

        # --------------------------------------------------------
        # P and S observations
        # --------------------------------------------------------

        P = station.loc[
            station["pha"] == "PG",
            "res",
        ].dropna()

        S = station.loc[
            station["pha"] == "SG",
            "res",
        ].dropna()

        # All valid residuals
        residuals = concat(
            [
                P.abs(),
                S.abs(),
            ],
            ignore_index=True,
        )

        nP = len(P)
        nS = len(S)
        nObs = nP + nS

        # --------------------------------------------------------
        # Event coverage
        # --------------------------------------------------------
        #
        # Count DISTINCT events, not rows.
        #
        # This is important because an event may contain multiple
        # observations from the same station/phase.
        # --------------------------------------------------------

        coverage_events = station["eid"].nunique()

        coverage = coverage_events / n_events

        # --------------------------------------------------------
        # Distance
        # --------------------------------------------------------

        D = station["dis"].dropna()

        # --------------------------------------------------------
        # Statistics
        # --------------------------------------------------------

        if nObs:
            median_abs_residual = residuals.median()
            mean_abs_residual = residuals.mean()
            rms_residual = sqrt(mean(residuals.to_numpy() ** 2))
        else:
            median_abs_residual = nan
            mean_abs_residual = nan
            rms_residual = nan

        mean_distance = D.mean() if len(D) else nan

        rows.append(
            {
                "code": sta,
                "nP": nP,
                "nS": nS,
                "nObs": nObs,
                "coverage": coverage,
                "median_abs_residual": median_abs_residual,
                "mean_abs_residual": mean_abs_residual,
                "rms_residual": rms_residual,
                "mean_distance": mean_distance,
            }
        )

    stats = DataFrame(rows)

    # ------------------------------------------------------------
    # Normalize the three criteria
    # ------------------------------------------------------------

    # More observations / higher coverage = better
    max_coverage = stats["coverage"].max()

    if max_coverage > 0:
        stats["observation_score"] = stats["coverage"] / max_coverage
    else:
        stats["observation_score"] = 0.0

    # Smaller residual = better
    max_residual = stats["median_abs_residual"].max()

    if notna(max_residual) and max_residual > 0:
        stats["residual_score"] = (
            1.0 - stats["median_abs_residual"] / max_residual
        )
    else:
        stats["residual_score"] = 0.0

    # Smaller distance = better
    #
    # Use inverse normalization:
    #
    #     score = 1 - D / Dmax
    #
    # so the closest station receives the highest score.
    max_distance = stats["mean_distance"].max()

    if notna(max_distance) and max_distance > 0:
        stats["distance_score"] = 1.0 - stats["mean_distance"] / max_distance
    else:
        stats["distance_score"] = 0.0

    # ------------------------------------------------------------
    # Handle stations with missing statistics
    # ------------------------------------------------------------

    stats = stats.fillna(0.0)

    # Numerical protection
    stats["observation_score"] = stats["observation_score"].clip(0.0, 1.0)

    stats["residual_score"] = stats["residual_score"].clip(0.0, 1.0)

    stats["distance_score"] = stats["distance_score"].clip(0.0, 1.0)

    # ------------------------------------------------------------
    # Final weighted score
    # ------------------------------------------------------------

    stats["score"] = (
        observation_weight * stats["observation_score"]
        + residual_weight * stats["residual_score"]
        + distance_weight * stats["distance_score"]
    )

    # ------------------------------------------------------------
    # Sort by final score
    # ------------------------------------------------------------

    stats = stats.sort_values(
        "score",
        ascending=False,
        ignore_index=True,
    )

    # ------------------------------------------------------------
    # ICC
    # ------------------------------------------------------------
    #
    # Best station:
    #     score highest
    #     icc = N
    #
    # Worst station:
    #     score lowest
    #     icc = 1
    #
    # Therefore ICC is simply the position in the sorted table.
    # ------------------------------------------------------------

    stats["icc"] = arange(
        len(stats),
        0,
        -1,
        dtype=int,
    )

    return stats

def finalized_inverted_stations(config, stage_n, run_n):

    if stage_n == 1:
        root = Path("outputs/stage_01")
    elif stage_n == 3:
        root = Path("outputs/stage_03")        
    
    run_dir = root / f"run_{run_n:02d}"
    osta_files = sorted(run_dir.rglob("velout.sta"))

    # outputs
    averaged_stafile = run_dir / "summary" / "stations.sta"

    # processing
    ostations = [read_vstation(config, i) for i in osta_files]
    stations = []
    for vstations in ostations:
        stations.append(DataFrame([i.dict() for i in vstations.stations]))
    ostations = concat(stations)
    calculate_average_corr(ostations, averaged_stafile)


def wrire_station(
    config: USER_CONFIG,
    infile: Path | str,
    outfile: Path | str,
    station_out_csv: Path | str,
):
    cat_inpfile = Path("inputs/catalog.csv")
    stations_sel = read_station(config, infile)
    select_cat_df = read_csv(cat_inpfile)
    sorted_station_df = calculate_reference_station_scores(select_cat_df)
    stations_df = sorted_station_df.merge(stations_sel, on="code", how="inner")
    p_delay = 0
    s_delay = 0
    fmt = "(a4,f7.4,a1,1x,f8.4,a1,1x," "i4,1x,i1,1x,i3,1x,f5.2,2x,f5.2,3x,i1)"

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    stations_sel.to_csv(station_out_csv, index=False)

    with outfile.open("w") as f:

        # First line is the FORTRAN format string.
        f.write(fmt + "\n")

        for _, row in stations_df.iterrows():

            station = str(row["code"]).upper()

            if len(station) > 4:
                raise ValueError(
                    f"VELEST station code must be <= 4 characters: "
                    f"{station}"
                )

            lat = float(row["lat"])
            lon = float(row["lon"])
            elv = int(round(float(row["elv"])))

            # VELEST expects positive coordinates + hemisphere character.
            lat_abs = abs(lat)
            lon_abs = abs(lon)

            lat_hemisphere = config.STUDY_AREA.lat_hemisphere
            lon_hemisphere = config.STUDY_AREA.lon_hemisphere

            # One velocity model.
            model = 1

            # Station delays.
            ptcor = float(p_delay)
            stcor = float(s_delay)

            icc = int(row["icc"])

            line = (
                f"{station:<4}"
                f"{lat_abs:7.4f}"
                f"{lat_hemisphere}"
                f" "
                f"{lon_abs:8.4f}"
                f"{lon_hemisphere}"
                f" "
                f"{elv:4d}"
                f" "
                f"{model:1d}"
                f" "
                f"{icc:3d}"
                f" "
                f"{ptcor:5.2f}"
                f"  "
                f"{stcor:5.2f}"
                f"   "
                f"{model:1d}"
            )

            f.write(line + "\n")
        f.write("\n")


def prepare_station(config, stage_n, run_n=None):

    if stage_n == 1:

        # inputs
        station_input = Path("inputs") / "stations.csv"

        # outputs
        root_dir = Path("outputs/stage_01")
        # run_dir = root_dir / f"run_{run_n:02d}"
        station_csv_outfile = root_dir / "stations_sel.csv"

        # Read / Implement / Write
        # for n_model in range(1, config.SYNTHETIC_MODELS.n_models + 1):
        #     model_dir = run_dir / f"model_{n_model:02d}"
        station_out = root_dir / "stations.sta"
        wrire_station(config, station_input, station_out, station_csv_outfile)

    elif stage_n == 2:

        # inputs
        run_n = config.STAGE_01.n_runs
        station_input = (
            Path("outputs/stage_01")
            / f"run_{run_n:02d}"
            / "summary"
            / "stations.sta"
        )

        station_sel = Path("outputs/stage_01") / "stations_sel.csv"

        # outputs
        run_dir = Path("outputs/stage_02")

        # Read / Implement / Write
        copy(station_input, run_dir)
        copy(station_sel, run_dir)

    if stage_n == 3:

        # inputs
        station_input = Path("outputs/stage_02") / "stations.sta"

        # outputs
        root_dir = Path("outputs/stage_03")
        run_dir = root_dir / f"run_{run_n:02d}"

        # Read / Implement / Write
        for n_model in range(1, config.STAGE_03.n_models + 1):
            model_dir = run_dir / f"model_{n_model:02d}"
            copy(station_input, model_dir)

    if stage_n == 4:
        
        stg_03_run_n = config.STAGE_03.n_runs

        # inputs
        station_input = Path(f"outputs/stage_03/run_{stg_03_run_n:02d}/summary") / "stations.sta"

        # outputs
        root_dir = Path("outputs/stage_04")
        run_dir = root_dir / f"run_{run_n:02d}"

        # Read / Implement / Write
        copy(station_input, run_dir)

    elif stage_n == 5:

        # inputs
        run_n = config.STAGE_03.n_runs
        station_input = (
            Path("outputs/stage_03")
            / f"run_{run_n:02d}"
            / "summary"
            / "stations.sta"
        )
        station_sel = Path("outputs/stage_01") / "stations_sel.csv"

        # outputs
        run_dir = Path("outputs/stage_05")

        # Read / Implement / Write
        copy(station_input, run_dir)
        copy(station_sel, run_dir)
