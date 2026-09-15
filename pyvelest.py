from pathlib import Path
from src.io.config import read_config
from src.io.workspace import prepare_workspace
from src.io.station import prepare_station, finalized_inverted_stations
from src.io.model import (
    prepare_model,
    finalized_inverted_models,
)
from src.io.catalog import summarize_catalog, prepare_catalog, reselect_best
from src.io.cmn import prepare_cmn
from src.core.velest import run_velest
from src.core.visualize import (
    plot_events,
    plot_models,
    plot_statistics,
    plot_dislocations,
)
from src.core.final import (
    prepare_final_model,
    prepare_final_catalog,
    prepare_final_stations,
)


class Main:
    def __init__(self, config_path: str):
        self.config = read_config(config_path)
        prepare_workspace(self.config)

    def load_catalog(self):
        self.catalog = summarize_catalog(self.config)

    def stage_01(self):
        n_runs = self.config.STAGE_01.n_runs
        for n in range(1, n_runs + 1):
            print(f"+++ Execution number {n} from stage_01 ...")
            prepare_station(self.config, stage_n=1, run_n=n)
            prepare_model(self.config, stage_n=1, run_n=n)
            prepare_catalog(self.config, self.catalog, stage_n=1, run_n=n)
            prepare_cmn(self.config, stage_n=1, run_n=n)
            run_velest(self.config, stage_n=1, run_n=n)
            finalized_inverted_models(self.config, stage_n=1, run_n=n)
            finalized_inverted_stations(self.config, stage_n=1, run_n=n)
            plot_models(self.config, stage_n=1, run_n=n)
            plot_statistics(self.config, stage_n=1, run_n=n)
        plot_events(self.config)

    def stage_02(self):
        print("+++ Executing stage_02 ...")
        prepare_station(self.config, stage_n=2)
        prepare_model(self.config, stage_n=2)
        prepare_catalog(self.config, self.catalog, stage_n=2)
        prepare_cmn(self.config, stage_n=2)
        run_velest(self.config, stage_n=2)
        reselect_best(self.config, stage_n=2)

    def stage_03(self):
        n_runs = self.config.STAGE_03.n_runs
        for n in range(1, n_runs + 1):
            print(f"+++ Execution number {n} from stage_03 ...")
            prepare_station(self.config, stage_n=3, run_n=n)
            prepare_model(self.config, stage_n=3, run_n=n)
            prepare_catalog(self.config, stage_n=3, run_n=n)
            prepare_cmn(self.config, stage_n=3, run_n=n)
            run_velest(self.config, stage_n=3, run_n=n)
            finalized_inverted_models(self.config, stage_n=3, run_n=n)
            finalized_inverted_stations(self.config, stage_n=3, run_n=n)
            plot_models(self.config, stage_n=3, run_n=n)
            plot_statistics(self.config, stage_n=3, run_n=n)

    def stage_04(self):
        n_runs = self.config.STAGE_04.n_runs
        for n in range(1, n_runs + 1):
            print(f"+++ Execution number {n} from stage_04 ...")
            prepare_station(self.config, stage_n=4, run_n=n)
            prepare_model(self.config, stage_n=4, run_n=n)
            prepare_catalog(self.config, stage_n=4, run_n=n)
            prepare_cmn(self.config, stage_n=4, run_n=n)
            run_velest(self.config, stage_n=4, run_n=n)
            plot_dislocations(self.config, stage_n=4, run_n=n)

    def stage_05(self):
        # prepare_station(self.config, stage_n=5)
        # prepare_model(self.config, stage_n=5)
        # prepare_catalog(self.config, self.catalog, stage_n=5)
        # prepare_cmn(self.config, stage_n=5)
        # run_velest(self.config, stage_n=5)   
        reselect_best(self.config, stage_n=5)
        # prepare_final_model(self.config)
        # prepare_final_catalog(self.config)
        # prepare_final_stations(self.config)


if __name__ == "__main__":
    config = Path("configs") / "user.yaml"
    app = Main(config.as_posix())
    # app.load_catalog()
    # app.stage_01()
    # app.stage_02()
    # app.stage_03()
    # app.stage_04()
    app.stage_05()
