from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import yaml


_ANALYSIS_PATH = Path("configs/analysis.yml")
_SOFTWARE_PATH = Path("configs/software.yml")


def _load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Missing config file: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config root must be a mapping: {path}")
    return data


def _resolve_path(data: dict[str, Any], path: Iterable[str]) -> Any:
    current: Any = data
    for key in path:
        if not isinstance(current, dict) or key not in current:
            joined = ".".join(path)
            raise KeyError(f"Missing config key: {joined}")
        current = current[key]
    return current


def get_analysis_value(path: list[str]) -> Any:
    analysis = _load_yaml(_ANALYSIS_PATH)
    return _resolve_path(analysis, path)


def get_results_dir() -> str:
    return str(get_analysis_value(["results_dir"]))


def get_software_module(tool: str) -> str:
    software = _load_yaml(_SOFTWARE_PATH)
    if tool not in software:
        raise KeyError(f"Unknown tool in software config: {tool}")
    module_name = software[tool].get("module")
    if not module_name:
        raise KeyError(f"Missing module for tool: {tool}")
    return str(module_name)


def get_software_value(tool: str, path: list[str], default: Any | None = None) -> Any:
    software = _load_yaml(_SOFTWARE_PATH)
    if tool not in software:
        if default is not None:
            return default
        raise KeyError(f"Unknown tool in software config: {tool}")

    current: Any = software[tool]
    for key in path:
        if not isinstance(current, dict) or key not in current:
            if default is not None:
                return default
            joined = ".".join(path)
            raise KeyError(f"Missing software key: {tool}.{joined}")
        current = current[key]
    return current
