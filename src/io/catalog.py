from __future__ import annotations

import re
from pathlib import Path

from obspy import UTCDateTime as utc
from obspy.core.event import (
    Catalog,
    Event,
    Origin,
    Magnitude,
    Pick,
    Arrival,
    WaveformStreamID,
)

from obspy import read_events
from pandas import DataFrame, read_csv
from src.io.select3D import select_events_3d
from datetime import timezone
from obspy.geodetics.base import gps2dist_azimuth as gps
from obspy.geodetics.base import kilometer2degrees as k2d
from shutil import copy


from datetime import datetime as dt
from math import sqrt
from numpy import nan, random

import warnings

warnings.filterwarnings("ignore")


def get_picks_info(event, eid, station_inpfile=Path("inputs/stations.csv")):
    """Extract VELEST-ready pick rows from one ObsPy event.

    Distances are returned in kilometers and travel times in seconds. Picks
    whose phase or station cannot be resolved are ignored because they cannot
    be represented reliably in the catalog summary.
    """
    origin = event.preferred_origin()
    if origin is None:
        if not event.origins:
            raise ValueError(f"Event {eid} has no origin")
        origin = event.origins[0]

    stations_df = read_csv(station_inpfile)
    station_info = {}
    for row in stations_df.itertuples(index=False):
        distance_m, azimuth, _ = gps(
            row.lat, row.lon, origin.latitude, origin.longitude
        )
        station_info[row.code] = (distance_m * 1e-3, azimuth)

    phases = []
    for arrival in origin.arrivals:
        phase = (arrival.phase or "").upper()
        if not phase or phase[0] not in "PS" or arrival.pick_id is None:
            continue

        pick = arrival.pick_id.get_referred_object()
        if pick is None or pick.waveform_id is None:
            continue
        sta = pick.waveform_id.station_code
        if sta not in station_info:
            continue

        dist, azimuth = station_info[sta]
        if "extra" not in pick.__dict__:
            weight = 0
        else:
            weight_info = pick.extra.get("nordic_pick_weight", {"value": 0})
            weight = weight_info.get("value", 0)
        phases.append(
            {
                "eid": eid,
                "sta": sta,
                "dis": dist,
                "azm": azimuth,
                "pha": phase,
                "wet": weight,
                "ptt": pick.time - origin.time,
                "res": arrival.time_residual,
            }
        )
    return phases


def summarize_catalog(config):

    # inputs
    cat_inpfile = config.INPUTS.catalog

    # outputs
    cat_outfile = Path("inputs/catalog.csv")

    # processing
    cat = read_events(cat_inpfile)
    data = []
    for eid, event in enumerate(cat, 1):
        eid = eid
        po = event.preferred_origin() or event.origins[0]
        pm = event.preferred_magnitude() or None
        ort = po.time
        lat = po.latitude
        lon = po.longitude
        dep = po.depth * 1e-3
        mag = pm.mag if pm else None
        gap = po.quality.azimuthal_gap or None
        rms = po.quality.standard_error or None
        nst = po.quality.used_station_count or None
        erh = (
            po.origin_uncertainty.horizontal_uncertainty
            if po.origin_uncertainty
            else None
        )
        erz = po.depth_errors.uncertainty * 1e-3 if po.depth_errors else None
        event_info = {
            "eid": eid,
            "ort": ort,
            "lat": lat,
            "lon": lon,
            "dep": dep,
            "mag": mag,
            "gap": gap,
            "rms": rms,
            "nst": nst,
            "erh": erh,
            "erz": erz,
        }
        phase_info = get_picks_info(event, eid, config.INPUTS.stations)
        phase_info.append(event_info)
        data.extend(phase_info)
    cat_df = DataFrame(data)
    cat_df.to_csv(cat_outfile, index=False, float_format="%.3f")
    return cat


def select(
    config, catalog, catalog_sum_file, select_cat_outfile, select_csv_outfile
):

    # process
    df = read_csv(catalog_sum_file)
    df.dropna(subset=["ort"], inplace=True)
    mask = df["lat"].between(
        config.STUDY_AREA.min_lat, config.STUDY_AREA.max_lat
    ) & df["lon"].between(config.STUDY_AREA.min_lon, config.STUDY_AREA.max_lon)

    limits = (
        ("erh", config.EVENT_SELECTION.max_horizontal_error_km, "le"),
        ("erz", config.EVENT_SELECTION.max_depth_error_km, "le"),
        ("rms", config.EVENT_SELECTION.max_rms_sec, "le"),
        ("gap", config.EVENT_SELECTION.max_gap_deg, "le"),
        ("nst", config.EVENT_SELECTION.min_used_stations, "ge"),
    )
    for column, limit, operator in limits:
        if limit is None:
            continue
        comparison = (
            df[column] <= limit if operator == "le" else df[column] >= limit
        )
        mask &= comparison | df[column].isna()

    filtered = df[mask].copy()
    selected = select_events_3d(filtered, config)
    selected.to_csv(select_csv_outfile, index=False, float_format="%.3f")
    cat_sel = Catalog()
    for eid in selected.eid:
        event = catalog[int(eid)-1]
        cat_sel.events.append(event)
    cat_sel.write(
        select_cat_outfile,
        format="NORDIC",
        high_accuracy=False,
        nordic_format="OLD",
    )
    return cat_sel


def make_origin(origin_data):
    origin = Origin(
        time=origin_data["time"],
        latitude=origin_data["latitude"],
        longitude=origin_data["longitude"],
        depth=(
            origin_data["depth_km"] * 1000.0
            if origin_data["depth_km"]
            else None
        ),
    )
    return origin


def make_magnitude(origin_data):
    magnitude = Magnitude(
        mag=origin_data["magnitude"],
        magnitude_type="ML",
    )
    return magnitude


def make_catalog(events_df):
    catalog = Catalog()
    events_df.replace({nan: None}, inplace=True)
    for eid, event_df in events_df.groupby(["eid"]):

        event_info = event_df.iloc[-1]
        event = Event()

        if not event_info.ort:
            origin_data = {
                "time": utc(1900, 1, 1),
                "latitude": 0,
                "longitude": 0,
                "depth_km": 0,
                "magnitude": 0,
            }
        else:
            origin_data = {
                "time": utc(event_info.ort),
                "latitude": event_info.lat,
                "longitude": event_info.lon,
                "depth_km": event_info.dep,
                "magnitude": event_info.mag,
            }

        origin = make_origin(origin_data)
        event.origins.append(origin)
        event.preferred_origin_id = origin.resource_id

        magnitude = make_magnitude(origin_data)
        event.magnitudes.append(magnitude)
        event.preferred_magnitude_id = magnitude.resource_id
        catalog.events.append(event)

        for r, row in event_df.iterrows():
            if row.ptt == None:
                continue
            pick = Pick(
                time=origin_data["time"] + row.ptt,
                waveform_id=WaveformStreamID(station_code=row.sta),
                phase_hint=row.pha,
                evaluation_mode="automatic",
            )
            pick.extra = {
                "nordic_pick_weight": {
                    "value": row.wet or 0,
                    "namespace": "velest",
                }
            }
            event.picks.append(pick)

            arrival = Arrival(
                pick_id=pick.resource_id,
                phase=row.pha,
                time_residual=row.res,
            )
            origin.arrivals.append(arrival)
    return catalog


def catalog_to_cnv(
    catalog,
    stations_inpfile,
    root_out,
    cnv_outfile,
    use_pick_time_errors=False,
):
    """
    Convert an ObsPy Catalog to VELEST CNV format.

    Parameters
    ----------
    catalog : obspy.core.event.Catalog
    csv_outfile : str or Path
    default_weight : int
        Weight assigned when no uncertainty exists.
    use_pick_time_errors : bool
        Convert pick.time_errors.uncertainty to VELEST weights.
    """

    # outputs
    root_out = Path(root_out)
    unused_st_file = root_out / "unused_st.csv"

    # processing
    stations_df = read_csv(stations_inpfile)
    station_codes = set(stations_df["code"].astype(str))
    default_weight = 0
    unused_st = []

    with open(cnv_outfile, "w") as f:

        for event in catalog:

            if event.preferred_origin():
                origin = event.preferred_origin()
            else:
                origin = event.origins[0]

            if event.preferred_magnitude():
                magnitude = event.preferred_magnitude().mag
            elif event.magnitudes:
                magnitude = event.magnitudes[0].mag
            else:
                magnitude = 0.0

            ot = origin.time.datetime.replace(tzinfo=timezone.utc)

            lat = abs(origin.latitude)
            lon = abs(origin.longitude)

            ns = "N" if origin.latitude >= 0 else "S"
            ew = "E" if origin.longitude >= 0 else "W"

            depth = origin.depth / 1000.0 if origin.depth else 0.0

            # -------------------------
            # Event header
            # -------------------------

            header = (
                f"{ot:%y%m%d} "
                f"{ot:%H%M} "
                f"{ot.second + ot.microsecond/1e6:05.2f} "
                f"{lat:7.4f}{ns} "
                f"{lon:8.4f}{ew} "
                f"{depth:7.2f}  "
                f"{magnitude:5.2f}"
            )

            f.write(header + "\n")

            observations = []
            repetitions = []

            for arrival in origin.arrivals:
                if arrival.pick_id is None:
                    continue
                pick = arrival.pick_id.get_referred_object()
                if pick is None or pick.waveform_id is None:
                    continue
                sta = pick.waveform_id.station_code[:4].ljust(4)
                if sta.strip() not in station_codes:
                    unused_st.append(sta.strip())
                    continue
                if not arrival.phase or not arrival.phase.strip():
                    continue
                phase = arrival.phase.upper()[0]
                tt = pick.time - origin.time
                weight = default_weight
                # 1. Nordic operator weight (highest priority)
                if (
                    hasattr(pick, "extra")
                    and "nordic_pick_weight" in pick.extra
                ):
                    try:
                        weight = int(pick.extra["nordic_pick_weight"]["value"])
                    except Exception:
                        pass

                # 2. Otherwise derive from pick uncertainty
                elif use_pick_time_errors and pick.time_errors:

                    err = pick.time_errors.uncertainty

                    if err is not None:
                        if err <= 0.05:
                            weight = 0
                        elif err <= 0.10:
                            weight = 1
                        elif err <= 0.20:
                            weight = 2
                        elif err <= 0.50:
                            weight = 3
                        else:
                            weight = 4
                obs = f"{sta}{phase}{weight:d}{tt:6.2f}"
                rep = f"{sta}{phase}"
                observations.append(obs) if rep not in repetitions else None
                repetitions.append(rep)

            # six observations per line
            for i in range(0, len(observations), 6):
                f.write("".join(observations[i : i + 6]) + "\n")

            f.write("\n")

        f.write("9999\n")

    unused_st = sorted(set(unused_st))
    unused_st_df = DataFrame({"code": unused_st})
    with open(unused_st_file, "w") as f:
        f.write(
            "# Stations in the catalog but outside the specified region are listed below:\n"
        )
        unused_st_df.to_csv(f, index=False)


def get_eid(line):
    args = line.split()
    return int(args[-3])


def add_origin(line, header):
    date = re.sub(r"\s+", "0", line[1:7])
    time = re.sub(r"\s+", "0", line[8:20])
    ort = dt.strptime(" ".join([date, time]), "%y%m%d %H:%M:%S.%f")
    lat = float(line[21:28])
    lon = float(line[30:37])
    dep = float(line[39:46])
    mag = float(line[47:51])
    nst = float(line[52:55])
    dmn = float(line[56:59])
    gap = float(line[60:63])
    rms = float(line[64:68])
    header.update(
        {
            "ort": ort,
            "lat": lat,
            "lon": lon,
            "dep": dep,
            "mag": mag,
            "nst": nst,
            "dmn": dmn,
            "gap": gap,
            "rms": rms,
        }
    )
    return header


def add_statistics(line, header):
    args = line.split()
    erx = float(args[0])
    ery = float(args[1])
    erz = float(args[2])
    erh = sqrt(erx**2 + ery**2)

    header.update({"erh": erh, "erz": erz})
    return header


def add_phase_info(line):
    sta = line[2:6]
    dis = float(line[7:11])
    azm = float(line[12:15])
    pha = line[21:23].strip()
    wet = line[23]
    ptt = float(line[37:43])
    res = float(line[58:65])
    data = {
        "sta": sta,
        "dis": dis,
        "azm": azm,
        "pha": pha,
        "wet": wet,
        "ptt": ptt,
        "res": res,
    }
    return data


def summarize_single_events(single_events_file, sum_outfile):
    events_df = []
    header = {}
    phases = []
    eid = None
    with open(single_events_file) as f:
        for line in f:
            if "E V E N T" in line:
                header = {}
                phases = []
                eid = get_eid(line)
            elif "ERROR" in line:
                header.update({"eid": eid})
                info = {"eid": eid}
                phases.append(info)
                phases.append(header)
                events_df.extend(phases)
            elif "DATE  ORIGIN" in line:
                line = next(f)
                header = add_origin(line, header)
                header.update({"eid": eid})
            elif "ERX  ERY  ERZ" in line:
                line = next(f)
                header = add_statistics(line, header)
                header.update({"eid": eid})
            elif "STN  DIST" in line:
                line = next(f)
                line = next(f)
                while line.strip() and "$$$" not in line:
                    info = add_phase_info(line)
                    info.update({"eid": eid})
                    phases.append(info)
                    line = next(f)
                    if "$$$" in line:
                        phases.append(header)
                        events_df.extend(phases)
    events_df = DataFrame(events_df)
    events_df.to_csv(sum_outfile, index=False, float_format="%.3f")
    return events_df


def convert_single_events(single_events_file, cat_outfile, sum_outfile):
    events_df = summarize_single_events(single_events_file, sum_outfile)
    catalog = make_catalog(events_df)
    catalog.write(
        cat_outfile, format="NORDIC", high_accuracy=False, nordic_format="OLD"
    )
    return catalog


def cnv2csv(cnv_inpfile, csv_outfile):
    data = []
    with open(cnv_inpfile) as fo:
        for line in fo:
            if len(line) > 35 and line[25] in "SN" and line[35] in "EW":
                ort = line[:17]
                ort = ort[:-5]+"59.99" if ort[-5] == "6" else ort
                ort = dt.strptime(ort, "%y%m%d %H%M %S.%f")
                lat = float(line[18:25])
                lon = float(line[27:35])
                dep = float(line[37:43])
                mag = float(line[46:50])
                info = {
                    "ort": ort,
                    "lat": lat,
                    "lon": lon,
                    "dep": dep,
                    "mag": mag,
                }
                data.append(info)
    cnv_df = DataFrame(data)
    cnv_df.to_csv(
        csv_outfile,
        index=False,
    )


def reselect_best(config, stage_n):

    rootpath = Path(f"outputs/stage_{stage_n:02d}")

    # inputs
    single_events_file = rootpath / "single_events.out"

    # outputs
    cat_rel_outfile = rootpath / "relocated.out"
    cat_sum_outfile = rootpath / "relocated.csv"

    select_cat_outfile = rootpath / "select.dat"
    select_cnv_outfile = rootpath / "select.cnv"
    select_csv_outfile = rootpath / "select.csv"

    cat_rel = convert_single_events(
        single_events_file, cat_rel_outfile, cat_sum_outfile
    )

    cat_sel = select(
        config,
        cat_rel,
        cat_sum_outfile,
        select_cat_outfile,
        select_csv_outfile,
    )

    stations_inpfile = rootpath / "stations_sel.csv"
    catalog_to_cnv(cat_sel, stations_inpfile, rootpath, select_cnv_outfile)


def perturb_catalog(
    config, run_n, cat_cnv_file, cat_csv_file, cnv_outfile, csv_outfile
):

    lat_std = k2d(config.STAGE_04.dx_km)
    lon_std = k2d(config.STAGE_04.dy_km)
    dep_std = config.STAGE_04.dz_km
    seed = config.STAGE_04.seed

    cat_csv = read_csv(cat_csv_file)

    rng = random.default_rng(seed + run_n)

    nEvt = len(cat_csv)
    cat_csv["lat"] += rng.normal(0, lat_std, size=nEvt)
    cat_csv["lon"] += rng.normal(0, lon_std, size=nEvt)
    cat_csv["dep"] += rng.normal(0, dep_std, size=nEvt)

    # Ensure all depths are positive
    cat_csv["dep"] = cat_csv["dep"].clip(lower=0)

    with open(cat_cnv_file) as fo, open(cnv_outfile, "w") as go:
        eid = 0
        for line in fo:
            if len(line) > 35 and line[25] in "NS" and line[35] in "EW":
                lat = cat_csv.iloc[eid].lat
                lon = cat_csv.iloc[eid].lon
                dep = cat_csv.iloc[eid].dep
                line = (
                    line[:18]
                    + f"{lat:7.4f}"
                    + line[25:27]
                    + f"{lon:8.4f}"
                    + line[35:38]
                    + f"{dep:6.2f}"
                    + line[44:]
                )
                eid += 1
            go.write(line)
    cat_csv.to_csv(csv_outfile, index=False, float_format="%.3f")


def prepare_catalog(config, cat_inp=None, stage_n=1, run_n=None):

    if stage_n == 1:

        root = Path("outputs/stage_01")

        # inputs
        catalog_csv_inpfile = Path("inputs/catalog.csv")
        stations_inpfile = root / "stations_sel.csv"

        # output
        select_cat_outfile = root / "select.dat"
        select_csv_outfile = root / "select.csv"
        select_cnv_outfile = root / "select.cnv"

        cat_sel = select(
            config,
            cat_inp,
            catalog_csv_inpfile,
            select_cat_outfile,
            select_csv_outfile,
        )

        run_dir = root / f"run_{run_n:02d}"
        for n_model in range(1, config.SYNTHETIC_MODELS.n_models + 1):
            model_dir = run_dir / f"model_{n_model:02d}"
            catalog_to_cnv(
                cat_sel, stations_inpfile, model_dir, select_cnv_outfile
            )

    elif stage_n == 2:
        # output
        root_out = Path("outputs/stage_02")
        stations_inpfile = root_out / "stations_sel.csv"
        original_cnv_outfile = root_out / "original.cnv"
        csv_outfile = root_out / "original.csv"
        catalog_to_cnv(
            cat_inp, stations_inpfile, root_out, original_cnv_outfile
        )
        cnv2csv(original_cnv_outfile, csv_outfile)

    elif stage_n == 3:

        root = Path("outputs/stage_03")

        cat_sel = Path("outputs/stage_02") / "select.cnv"

        # output
        run_dir = root / f"run_{run_n:02d}"

        for n_model in range(1, config.STAGE_03.n_models + 1):
            model_dir = run_dir / f"model_{n_model:02d}"
            copy(cat_sel, model_dir)

    elif stage_n == 4:

        root = Path("outputs/stage_04")

        # inputs
        cat_cnv_file = Path("outputs/stage_02") / "select.cnv"
        cat_csv_file = Path("outputs/stage_02") / "select.csv"

        # output
        run_dir = root / f"run_{run_n:02d}"
        cnv_outfile = root / f"run_{run_n:02d}" / "perturb.cnv"
        csv_outfile = root / f"run_{run_n:02d}" / "perturb.csv"

        perturb_catalog(
            config, run_n, cat_cnv_file, cat_csv_file, cnv_outfile, csv_outfile
        )

    elif stage_n == 5:
        # output
        root_out = Path("outputs/stage_05")
        stations_inpfile = root_out / "stations_sel.csv"
        original_cnv_outfile = root_out / "original.cnv"
        csv_outfile = root_out / "original.csv"
        catalog_to_cnv(
            cat_inp, stations_inpfile, root_out, original_cnv_outfile
        )
        cnv2csv(original_cnv_outfile, csv_outfile)
