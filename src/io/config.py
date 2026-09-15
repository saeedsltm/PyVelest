from src.model.Models import USER_CONFIG
from pathlib import Path
import yaml


def read_config(path: str | Path) -> USER_CONFIG:
    """Read and validate a PyVELEST YAML config file."""

    with Path(path).open("r", encoding="utf-8") as file:
        data = yaml.safe_load(file) or {}
    return USER_CONFIG.model_validate(data)
