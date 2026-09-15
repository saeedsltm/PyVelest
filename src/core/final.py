from src.io.model import read_vmodel, write_vmodel
from pathlib import Path
from pandas import DataFrame


def prepare_final_model(config):

    # inputs
    run_n = config.STAGE_03.n_runs
    vmodel_inpfile = (
        Path("outputs/stage_03") / f"run_{run_n:02d}" / "summary" / "model.mod"
    )

    # outputs
    vmodel_outfile = Path("outputs/stage_05") / "model.mod"
    csv_outfile = Path("outputs/stage_05") / "model.csv"

    # process
    vmodel = read_vmodel(vmodel_inpfile)
    write_vmodel(vmodel, vmodel_outfile)

    vmodel_df = {
        "vp": [layer.vel for layer in vmodel.layers_vp],
        "depth": [layer.depth for layer in vmodel.layers_vp],
        "vs": (
            [layer.vel for layer in vmodel.layers_vs]
            if vmodel.layers_vs
            else None
        ),
    }
    vmodel_df = DataFrame(vmodel_df)
    vmodel_df.to_csv(csv_outfile, index=False)


def prepare_final_catalog(config):
    pass


def prepare_final_stations(config):
    pass
