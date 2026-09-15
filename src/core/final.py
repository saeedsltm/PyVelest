from src.io.model import read_vmodel, write_vmodel
from src.io.station import read_vstation
from pathlib import Path
from pandas import DataFrame
from shutil import copy

def prepare_final_model(config):

    # inputs
    vmodel_inpfile = (
        Path("outputs/stage_05") / "model.mod"
    )

    # outputs
    vmodel_outfile = Path("outputs/stage_fn") / "model.mod"
    csv_outfile = Path("outputs/stage_fn") / "model.csv"

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
    
    # inputs
    reloc_csv_file = Path("outputs/stage_05") / "relocated.csv"
    reloc_cat_file = Path("outputs/stage_05") / "relocated.out"
    
    # outputs
    reloc_csv_outfile = Path("outputs/stage_fn") / "relocated.csv"
    reloc_cat_outfile = Path("outputs/stage_fn") / "relocated.out"
    
    copy(reloc_csv_file, reloc_csv_outfile)
    copy(reloc_cat_file, reloc_cat_outfile)


def prepare_final_stations(config):
    
    # inputs
    vstations_file = Path("outputs/stage_05") / "stations.sta"
    
    # outputs
    vstations_outfile = Path("outputs/stage_fn") / "stations.sta"
    stations_csvfile = Path("outputs/stage_fn") / "stations.csv"
    
    vstations = read_vstation(config, vstations_file)
    
    copy(vstations_file, vstations_outfile)
    
    stations_df = DataFrame([i.dict() for i in vstations.stations])

    stations_df.drop(["mn", "icc"], axis=1, inplace=True)
    stations_df.to_csv(stations_csvfile, index=False, float_format="%.4f")
    
