from pathlib import Path


def prepare_workspace(config):

    outputs = Path("outputs")
    outputs.mkdir(parents=True, exist_ok=True)

    for directory in [
        "stage_01",
        "stage_02",
        "stage_03",
        "stage_04",
        "stage_05",
    ]:

        target = Path("outputs") / directory
        target.mkdir(parents=True, exist_ok=True)

        if directory in ["stage_01", "stage_03"]:
            for n_run in range(1, config.STAGE_01.n_runs + 1):
                run_dir = target / f"run_{n_run:02d}"
                run_dir.mkdir(parents=True, exist_ok=True)
                sum_dir = run_dir / "summary"
                sum_dir.mkdir(parents=True, exist_ok=True)
                for n_model in range(1, config.SYNTHETIC_MODELS.n_models + 1):
                    model_dir = run_dir / f"model_{n_model:02d}"
                    model_dir.mkdir(parents=True, exist_ok=True)

        if directory in ["stage_04"]:
            for n_run in range(1, config.STAGE_04.n_runs + 1):
                run_dir = target / f"run_{n_run:02d}"
                run_dir.mkdir(parents=True, exist_ok=True)

                    
    figures = Path("figures")
    figures.mkdir(parents=True, exist_ok=True)

    if config.GLOBAL_SETTINGS.clean_workspace:
        for folder in [outputs, figures]:
            for item in folder.rglob("*"):
                if item.is_file() and item.exists():
                    item.unlink()
