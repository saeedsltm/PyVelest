from pandas import DataFrame
from numpy import linspace, clip, digitize


def select_events_3d(df: DataFrame, config) -> DataFrame:
    """
    Homogeneously select events from a 3D grid.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain:
            eid, lat, lon, dep
    config : object
        Must provide:
            max_events
            cell_count_lat
            cell_count_lon
            cell_count_depth

    Returns
    -------
    pandas.DataFrame
        Selected events.
    """

    df = df.copy()

    # -------------------------------------------------------------
    # build 3D grid
    # -------------------------------------------------------------
    lat_edges = linspace(df.lat.min(), df.lat.max(),
                            config.EVENT_SELECTION.cell_count_lat + 1)

    lon_edges = linspace(df.lon.min(), df.lon.max(),
                            config.EVENT_SELECTION.cell_count_lon + 1)

    dep_edges = linspace(df.dep.min(), df.dep.max(),
                            config.EVENT_SELECTION.cell_count_depth + 1)

    df["ilat"] = clip(
        digitize(df.lat, lat_edges) - 1,
        0, config.EVENT_SELECTION.cell_count_lat - 1
    )

    df["ilon"] = clip(
        digitize(df.lon, lon_edges) - 1,
        0, config.EVENT_SELECTION.cell_count_lon - 1
    )

    df["idep"] = clip(
        digitize(df.dep, dep_edges) - 1,
        0, config.EVENT_SELECTION.cell_count_depth - 1
    )

    # -------------------------------------------------------------
    # rank events inside each cell
    # best quality first
    # -------------------------------------------------------------
    sort_cols = ["rms", "erh", "erz", "gap", "mag"]
    ascending = [True, True, True, True, False]

    df = df.sort_values(sort_cols, ascending=ascending)

    groups = {
        key: grp.copy()
        for key, grp in df.groupby(["ilat", "ilon", "idep"])
    }

    selected = []

    # -------------------------------------------------------------
    # round-robin selection
    # -------------------------------------------------------------
    while len(selected) < config.EVENT_SELECTION.max_events:

        added = False

        for key in list(groups.keys()):

            grp = groups[key]

            if grp.empty:
                continue

            selected.append(grp.iloc[0])

            groups[key] = grp.iloc[1:]

            added = True

            if len(selected) >= config.EVENT_SELECTION.max_events:
                break

        if not added:
            break
    selected = (
        DataFrame(selected)
        .drop(columns=["ilat", "ilon", "idep"])
        .sort_values("eid")
        .reset_index(drop=True)
    )
    

    return selected