from pathlib import Path
import re
from numpy import abs, nan, min, max , std, sqrt, mean, where
from pandas import DataFrame


def read_velest_inversion(filename):
    """
    Parse a VELEST main inversion output file.

    The parser follows the output logic implemented in VELEST:

        RMSDATVAR()
        CHECKSOLUTION()
        NITTOUTPUT()
        BACKUP()
        STEPLENGTHCALC()

    Returns
    -------
    dict
        attempts
            One row for every normal iteration / BACKUP attempt.

        velocity_changes
            Actual velocity changes applied to the model.

        velocity_states
            Velocity model states printed by VELEST.

        velocity_summary
            Statistics of velocity changes.

        hypocenter_changes
            Actual hypocenter changes applied to events.

        hypocenter_summary
            Statistics of hypocenter changes.

    Notes
    -----
    For a normal iteration:

        vp_new = vp_old + dvp

    For BACKUP k:

        b = b / 2
        vp_new = vp_old - b

    Therefore the actual BACKUP change is:

        applied_change = -printed_dvp

    The final NITTOUTPUT after a BACKUP prints the remaining
    solution vector b. That b is NOT another applied change.
    """

    filename = Path(filename)

    text = filename.read_text(errors="replace")

    # =========================================================
    # Numeric pattern
    # =========================================================

    number = (
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)"
        r"(?:[Ee][-+]?\d+)?"
    )

    # =========================================================
    # 1. FIND ALL CHECKSOLUTION RESULTS
    # =========================================================
    #
    # Normal:
    #
    # Iteration nr  2 obtained:
    #
    # Backup:
    #
    # (Iteration nr  2)   BACKUP nr 1 obtained:
    #
    # =========================================================

    result_re = re.compile(
        rf"^\s*(?:"
        rf"Iteration nr\s+(?P<iteration>\d+)"
        rf"|"
        rf"\(Iteration nr\s+(?P<backup_iteration>\d+)\)"
        rf"\s+BACKUP nr\s+(?P<backup>\d+)"
        rf")\s+obtained:\s*$"
        rf"\n\s*DATVAR=\s*(?P<datavar>{number})"
        rf"\s+mean sqrd residual=\s*(?P<msres>{number})"
        rf"\s+RMS RESIDUAL=\s*(?P<rms>{number})",
        re.MULTILINE,
    )

    results = list(result_re.finditer(text))

    if not results:
        raise ValueError(
            "No VELEST iteration/backup results found "
            "in the inversion output."
        )

    # =========================================================
    # 2. ITERATION / BACKUP STATISTICS
    # =========================================================

    attempts = []

    # These are the values stored by CHECKSOLUTION().
    #
    # They are updated only when VELEST accepts a solution.
    previous_datavar = None
    previous_msres = None

    for i, match in enumerate(results):

        iteration = int(
            match.group("iteration")
            or match.group("backup_iteration")
        )

        backup = int(
            match.group("backup") or 0
        )

        datavar = float(
            match.group("datavar")
        )

        msres = float(
            match.group("msres")
        )

        rms = float(
            match.group("rms")
        )

        # -----------------------------------------------------
        # Find the step length immediately preceding this
        # CHECKSOLUTION output.
        #
        # This is important because VELEST prints:
        #
        #     (Applied) Step length = ...
        #
        # BEFORE the forward calculation and CHECKSOLUTION.
        # -----------------------------------------------------

        previous_end = (
            results[i - 1].end()
            if i > 0
            else 0
        )

        preceding_text = text[
            previous_end:match.start()
        ]

        step_matches = list(
            re.finditer(
                rf"\(Applied\) Step length\s*=\s*"
                rf"({number})",
                preceding_text,
            )
        )

        if step_matches:
            step_length = float(
                step_matches[-1].group(1)
            )
        else:
            step_length = nan

        # -----------------------------------------------------
        # Reproduce CHECKSOLUTION() exactly.
        #
        # VELEST accepts if:
        #
        #     old_datavar / new_datavar >= .99
        #
        # OR
        #
        #     old_msres / new_msres >= .99
        #
        # -----------------------------------------------------

        if iteration == 0:

            better = True

        else:

            if previous_datavar is None:
                raise ValueError(
                    "Cannot determine VELEST convergence status "
                    "because no previous accepted solution exists."
                )

            var_ratio = (
                previous_datavar / datavar
            )

            ms_ratio = (
                previous_msres / msres
            )

            better = (
                var_ratio >= 0.99
                or ms_ratio >= 0.99
            )

        attempts.append({
            "iteration": iteration,
            "backup": backup,
            "attempt": backup,

            "data_variance": datavar,
            "mean_squared_residual": msres,
            "rms_residual": rms,

            "step_length": step_length,

            "better": better,
        })

        # CHECKSOLUTION only changes its reference solution
        # when the current solution is better.
        if better:

            previous_datavar = datavar
            previous_msres = msres

    attempts_df = DataFrame(
        attempts
    )

    # ---------------------------------------------------------
    # Determine final result.
    # ---------------------------------------------------------

    attempts_df["accepted"] = (
        attempts_df["better"]
    )

    attempts_df["final"] = False

    # The last CHECKSOLUTION record preceding
    # "final solution reached" is the terminal attempt.
    final_marker = re.search(
        r"~~~ final solution reached",
        text,
        re.IGNORECASE,
    )

    if final_marker:

        final_result_index = None

        for i, result in enumerate(results):

            if result.start() < final_marker.start():
                final_result_index = i

        if final_result_index is not None:
            attempts_df.loc[
                final_result_index,
                "final"
            ] = True

    else:

        # For a truncated output, use the last result.
        attempts_df.loc[
            attempts_df.index[-1],
            "final"
        ] = True

    attempts_df["status"] = where(
        attempts_df["accepted"],
        "accepted",
        where(
            attempts_df["final"],
            "final_rejected",
            "rejected",
        ),
    )

    # =========================================================
    # 3. HYPOCENTER TABLE
    # =========================================================

    event_header_re = re.compile(
        r"^\s*eq\s+ot\s+x\s+y\s+z\s+rms\s+avres\s+"
        r"dot\s+dx\s+dy\s+dz\s*$",
        re.MULTILINE,
    )

    event_row_re = re.compile(
        rf"^\s*(\d+)\s+"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})\s*$",
        re.MULTILINE,
    )

    # event table belonging to each result
    event_tables = {}

    for i, result in enumerate(results):

        next_start = (
            results[i + 1].start()
            if i + 1 < len(results)
            else len(text)
        )

        block = text[
            result.end():next_start
        ]

        header = event_header_re.search(
            block
        )

        if header is None:
            continue

        rows = []

        for row in event_row_re.finditer(
            block[header.end():]
        ):

            rows.append({
                "event": int(row.group(1)),

                "ot": float(row.group(2)),
                "x": float(row.group(3)),
                "y": float(row.group(4)),
                "z": float(row.group(5)),

                "rms": float(row.group(6)),
                "avres": float(row.group(7)),

                "dot": float(row.group(8)),
                "dx": float(row.group(9)),
                "dy": float(row.group(10)),
                "dz": float(row.group(11)),
            })

        event_tables[i] = rows

    # =========================================================
    # 4. HYPOCENTER CHANGES
    # =========================================================

    hypocenter_changes = []

    # Store the original iteration adjustment.
    #
    # BACKUP() uses this recursively:
    #
    # backup 1 = -original / 2
    # backup 2 = -original / 4
    # backup 3 = -original / 8
    # backup 4 = -original / 16
    #
    original_adjustments = {}

    for i, attempt in attempts_df.iterrows():

        iteration = int(
            attempt["iteration"]
        )

        backup = int(
            attempt["backup"]
        )

        # -----------------------------------------------------
        # Normal iteration output
        # -----------------------------------------------------

        if backup == 0 and i in event_tables:

            rows = event_tables[i]

            original_adjustments[
                iteration
            ] = {}

            for row in rows:

                original_adjustments[
                    iteration
                ][row["event"]] = {
                    "dot": row["dot"],
                    "dx": row["dx"],
                    "dy": row["dy"],
                    "dz": row["dz"],
                }

                hypocenter_changes.append({
                    "iteration": iteration,
                    "backup": 0,
                    "event": row["event"],

                    "origin_time_change":
                        row["dot"],

                    "dx": row["dx"],
                    "dy": row["dy"],
                    "dz": row["dz"],

                    "change_type":
                        "iteration",

                    "applied": True,
                })

        # -----------------------------------------------------
        # BACKUP
        # -----------------------------------------------------
        #
        # BACKUP() does:
        #
        #     b = b / 2
        #     e = e - b
        #
        # Therefore the actual applied change is:
        #
        #     -original / 2^backup
        #
        # -----------------------------------------------------

        elif backup > 0:

            if iteration not in original_adjustments:

                raise ValueError(
                    "Cannot reconstruct hypocenter backup "
                    f"{backup} of iteration {iteration}: "
                    "the original iteration adjustment was "
                    "not found."
                )

            factor = -1.0 / (2 ** backup)

            for event, original in (
                original_adjustments[
                    iteration
                ].items()
            ):

                hypocenter_changes.append({
                    "iteration": iteration,
                    "backup": backup,
                    "event": event,

                    "origin_time_change":
                        original["dot"] * factor,

                    "dx":
                        original["dx"] * factor,

                    "dy":
                        original["dy"] * factor,

                    "dz":
                        original["dz"] * factor,

                    "change_type":
                        "backup",

                    "applied": True,
                })

    hypocenter_changes_df = DataFrame(
        hypocenter_changes
    )

    # =========================================================
    # 5. VELOCITY OUTPUT
    # =========================================================
    #
    # We must treat the two VELEST outputs differently.
    #
    # NITTOUTPUT:
    #
    #     Velocity adjustments:
    #
    # BACKUP:
    #
    #     Velocity readjustments:
    #
    # Also note that the BACKUP readjustment is printed BEFORE
    # the corresponding BACKUP CHECKSOLUTION line.
    #
    # Therefore:
    #
    #     readjustments -> next result
    #
    #     adjustments -> previous result
    #
    # =========================================================

    velocity_marker_re = re.compile(
        r"^\s*Velocity "
        r"(?P<kind>adjustments|readjustments):\s*$",
        re.MULTILINE,
    )

    velocity_model_re = re.compile(
        r"^\s*Velocity model\s+(\d+)\s*$",
        re.MULTILINE,
    )

    velocity_row_re = re.compile(
        rf"^\s*"
        rf"({number})\s+"
        rf"({number})\s+"
        rf"({number})"
        rf"(?:\s+\S)?\s*$",
        re.MULTILINE,
    )

    velocity_changes = []
    velocity_states = []

    markers = list(
        velocity_marker_re.finditer(text)
    )

    for marker_index, marker in enumerate(
        markers
    ):

        # -----------------------------------------------------
        # Find result before marker
        # -----------------------------------------------------

        previous_result = None

        for i, result in enumerate(results):

            if result.start() < marker.start():
                previous_result = i
            else:
                break

        # -----------------------------------------------------
        # Find result after marker
        # -----------------------------------------------------

        next_result = None

        for i, result in enumerate(results):

            if result.start() > marker.start():
                next_result = i
                break

        kind = marker.group("kind")

        if kind == "readjustments":

            if next_result is None:
                raise ValueError(
                    "Velocity readjustment without "
                    "a corresponding BACKUP result."
                )

            result_index = next_result
            change_type = "backup"

        else:

            if previous_result is None:
                raise ValueError(
                    "Velocity adjustment without "
                    "a corresponding iteration result."
                )

            result_index = previous_result

            # A normal NITTOUTPUT gives an actual applied
            # iteration adjustment.
            #
            # A NITTOUTPUT following a BACKUP is instead a
            # state report; its b vector is the remaining
            # solution vector, not a newly applied change.
            if attempts_df.iloc[
                result_index
            ]["backup"] == 0:

                change_type = "iteration"

            else:

                change_type = "state_output"

        attempt = attempts_df.iloc[
            result_index
        ]

        iteration = int(
            attempt["iteration"]
        )

        backup = int(
            attempt["backup"]
        )

        # -----------------------------------------------------
        # End of this marker section
        # -----------------------------------------------------

        if marker_index + 1 < len(markers):

            section_end = markers[
                marker_index + 1
            ].start()

        else:

            section_end = len(text)

        section = text[
            marker.end():section_end
        ]

        models = list(
            velocity_model_re.finditer(
                section
            )
        )

        for model_index, model_match in enumerate(
            models
        ):

            model = int(
                model_match.group(1)
            )

            if model_index + 1 < len(models):

                model_end = models[
                    model_index + 1
                ].start()

            else:

                model_end = len(section)

            model_section = section[
                model_match.end():model_end
            ]

            # -------------------------------------------------
            # IMPORTANT:
            #
            # Read until the exact VELEST average-velocity
            # section. This guarantees that all velocity rows
            # belonging to the model are captured.
            # -------------------------------------------------

            average_marker = re.search(
                r"^\s*Calculation of average velocity starts",
                model_section,
                re.MULTILINE,
            )

            if average_marker:

                model_section = model_section[
                    :average_marker.start()
                ]

            layer = 0

            for row in velocity_row_re.finditer(
                model_section
            ):

                layer += 1

                velocity = float(
                    row.group(1)
                )

                printed_change = float(
                    row.group(2)
                )

                depth = float(
                    row.group(3)
                )

                # ---------------------------------------------
                # Interpret DVP according to the source code.
                # ---------------------------------------------

                if change_type == "iteration":

                    applied_change = (
                        printed_change
                    )

                    semantics = (
                        "applied_iteration_adjustment"
                    )

                elif change_type == "backup":

                    # BACKUP():
                    #
                    #     b = b / 2
                    #     vp = vp - b
                    #
                    applied_change = (
                        -printed_change
                    )

                    semantics = (
                        "applied_backup_readjustment"
                    )

                else:

                    # NITTOUTPUT after a BACKUP.
                    #
                    # This is a state report. The printed b
                    # is the remaining inversion vector.
                    applied_change = nan

                    semantics = (
                        "remaining_solution_vector"
                    )

                velocity_changes.append({
                    "iteration": iteration,
                    "backup": backup,

                    "model": model,
                    "layer": layer,
                    "depth": depth,

                    "velocity":
                        velocity,

                    "printed_change":
                        printed_change,

                    "applied_change":
                        applied_change,

                    "change_type":
                        change_type,

                    "change_semantics":
                        semantics,

                    "applied":
                        change_type != "state_output",
                })

                velocity_states.append({
                    "iteration": iteration,
                    "backup": backup,

                    "model": model,
                    "layer": layer,
                    "depth": depth,

                    "velocity":
                        velocity,

                    "state_type":
                        (
                            "iteration"
                            if change_type == "iteration"
                            else
                            "backup"
                            if change_type == "backup"
                            else
                            "final_state"
                        ),
                })

    velocity_changes_df = DataFrame(
        velocity_changes
    )

    velocity_states_df = DataFrame(
        velocity_states
    )

    # =========================================================
    # 6. VELOCITY STATISTICS
    # =========================================================

    velocity_summary = []

    if not velocity_changes_df.empty:

        applied = velocity_changes_df[
            velocity_changes_df["applied"]
        ]

        for (
            iteration,
            backup,
            model,
        ), group in applied.groupby(
            [
                "iteration",
                "backup",
                "model",
            ]
        ):

            values = group[
                "applied_change"
            ].to_numpy(float)

            velocity_summary.append({
                "iteration":
                    int(iteration),

                "backup":
                    int(backup),

                "model":
                    int(model),

                "mean_change":
                    mean(values),

                "std_change":
                    std(
                        values,
                        ddof=0,
                    ),

                "rms_change":
                    sqrt(
                        mean(
                            values ** 2
                        )
                    ),

                "mean_abs_change":
                    mean(
                        abs(values)
                    ),

                "min_change":
                    min(values),

                "max_change":
                    max(values),

                "max_abs_change":
                    max(
                        abs(values)
                    ),

                "n_layers":
                    len(values),
            })

    velocity_summary_df = DataFrame(
        velocity_summary
    )

    # =========================================================
    # 7. HYPOCENTER STATISTICS
    # =========================================================

    hypocenter_summary = []

    if not hypocenter_changes_df.empty:

        for (
            iteration,
            backup,
        ), group in hypocenter_changes_df.groupby(
            [
                "iteration",
                "backup",
            ]
        ):

            for parameter, column in {
                "origin_time":
                    "origin_time_change",

                "x":
                    "dx",

                "y":
                    "dy",

                "z":
                    "dz",
            }.items():

                values = group[
                    column
                ].to_numpy(float)

                hypocenter_summary.append({
                    "iteration":
                        int(iteration),

                    "backup":
                        int(backup),

                    "parameter":
                        parameter,

                    "mean":
                        mean(values),

                    "std":
                        std(
                            values,
                            ddof=0,
                        ),

                    "rms":
                        sqrt(
                            mean(
                                values ** 2
                            )
                        ),

                    "mean_abs":
                        mean(
                            abs(values)
                        ),

                    "min":
                        min(values),

                    "max":
                        max(values),

                    "n_events":
                        len(values),
                })

    hypocenter_summary_df = DataFrame(
        hypocenter_summary
    )

    # =========================================================
    # RETURN
    # =========================================================

    return {
        "attempts":
            attempts_df,

        "velocity_changes":
            velocity_changes_df,

        "velocity_states":
            velocity_states_df,

        "velocity_summary":
            velocity_summary_df,

        "hypocenter_changes":
            hypocenter_changes_df,

        "hypocenter_summary":
            hypocenter_summary_df,
    }