from pathlib import Path
from subprocess import run
from shutil import copy
from src.io.velest import read_velest_inversion
from src.io.catalog import cnv2csv
from tqdm import tqdm


def run_velest(config, stage_n, run_n=None):
    here = Path()

    if stage_n == 1:

        root = Path("outputs/stage_01")

        # outputs

        run_dir = root / f"run_{run_n:02d}"

        # processing
        n_models = config.SYNTHETIC_MODELS.n_models
        desc = "+++ Running VELEST ..."
        for n_model in tqdm(range(1, n_models + 1), total=n_models, desc=desc):

            model_dir = run_dir / f"model_{n_model:02d}"

            # inputs
            cmn_inpfile = model_dir / "velest.cmn"
            inv_inpfile = model_dir / "velout.inv"
            velest_exec = Path("src/exc") / "velest"

            copy(cmn_inpfile, here)
            result = run([velest_exec], capture_output=True, text=True)

            if result.stderr == "":

                # outputs
                atm_outfile = model_dir / "attempts.csv"
                hyc_outfile = model_dir / "hypocenter_changes.csv"
                hps_outfile = model_dir / "hypocenter_summary.csv"
                vec_outfile = model_dir / "velocity_changes.csv"
                vlst_outfile = model_dir / "velocity_summary.csv"
                vlsu_outfile = model_dir / "velocity_states.csv"

                results = read_velest_inversion(inv_inpfile)

                results["attempts"].to_csv(atm_outfile, index=False)
                results["hypocenter_changes"].to_csv(hyc_outfile, index=False)
                results["hypocenter_summary"].to_csv(hps_outfile, index=False)
                results["velocity_changes"].to_csv(vec_outfile, index=False)
                results["velocity_summary"].to_csv(vlsu_outfile, index=False)
                results["velocity_states"].to_csv(vlst_outfile, index=False)

    if stage_n == 2:

        # outputs
        root = Path("outputs/stage_02")

        # processing
        cmn_inpfile = root / "velest.cmn"
        inv_inpfile = root / "velout.inv"
        velest_exec = Path("src/exc") / "velest"

        copy(cmn_inpfile, here)
        result = run([velest_exec], capture_output=True, text=True)

    if stage_n == 3:

        root = Path("outputs/stage_03")

        # outputs
        run_dir = root / f"run_{run_n:02d}"

        # processing
        n_models = config.STAGE_03.n_models
        desc = "+++ Running VELEST ..."
        for n_model in tqdm(range(1, n_models + 1), total=n_models, desc=desc):

            model_dir = run_dir / f"model_{n_model:02d}"

            # inputs
            cmn_inpfile = model_dir / "velest.cmn"
            inv_inpfile = model_dir / "velout.inv"
            velest_exec = Path("src/exc") / "velest"

            copy(cmn_inpfile, here)
            result = run([velest_exec], capture_output=True, text=True)

            if result.stderr == "":

                # outputs
                atm_outfile = model_dir / "attempts.csv"
                hyc_outfile = model_dir / "hypocenter_changes.csv"
                hps_outfile = model_dir / "hypocenter_summary.csv"
                vec_outfile = model_dir / "velocity_changes.csv"
                vlst_outfile = model_dir / "velocity_summary.csv"
                vlsu_outfile = model_dir / "velocity_states.csv"

                results = read_velest_inversion(inv_inpfile)

                results["attempts"].to_csv(atm_outfile, index=False)
                results["hypocenter_changes"].to_csv(hyc_outfile, index=False)
                results["hypocenter_summary"].to_csv(hps_outfile, index=False)
                results["velocity_changes"].to_csv(vec_outfile, index=False)
                results["velocity_summary"].to_csv(vlsu_outfile, index=False)
                results["velocity_states"].to_csv(vlst_outfile, index=False)


    if stage_n == 4:

        root = Path("outputs/stage_04")

        # outputs
        run_dir = root / f"run_{run_n:02d}"
        rel_cnv_outfile = root / f"run_{run_n:02d}" / "velout.cnv"
        rel_csv_outfile = root / f"run_{run_n:02d}" / "relocated.csv"

        # processing
        cmn_inpfile = run_dir / "velest.cmn"
        inv_inpfile = run_dir / "velout.inv"
        velest_exec = Path("src/exc") / "velest"

        copy(cmn_inpfile, here)
        result = run([velest_exec], capture_output=True, text=True)

        if result.stderr == "":

            # outputs
            atm_outfile = run_dir / "attempts.csv"
            hyc_outfile = run_dir / "hypocenter_changes.csv"
            hps_outfile = run_dir / "hypocenter_summary.csv"
            vec_outfile = run_dir / "velocity_changes.csv"
            vlst_outfile = run_dir / "velocity_summary.csv"
            vlsu_outfile = run_dir / "velocity_states.csv"

            results = read_velest_inversion(inv_inpfile)

            results["attempts"].to_csv(atm_outfile, index=False)
            results["hypocenter_changes"].to_csv(hyc_outfile, index=False)
            results["hypocenter_summary"].to_csv(hps_outfile, index=False)
            results["velocity_changes"].to_csv(vec_outfile, index=False)
            results["velocity_summary"].to_csv(vlsu_outfile, index=False)
            results["velocity_states"].to_csv(vlst_outfile, index=False)
            
            cnv2csv(rel_cnv_outfile, rel_csv_outfile)

    if stage_n == 5:

        # outputs
        root = Path("outputs/stage_05")

        # processing
        cmn_inpfile = root / "velest.cmn"
        inv_inpfile = root / "velout.inv"
        velest_exec = Path("src/exc") / "velest"

        copy(cmn_inpfile, here)
        result = run([velest_exec], capture_output=True, text=True)            