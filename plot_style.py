"""Shared plot styling for INCI plotting modules and subprocesses."""

from __future__ import annotations

import json
import os
from typing import Any


DEFAULTS: dict[str, Any] = {
    "font_family": "DejaVu Sans",
    "font_size": 11.0,
    "title_size": 14.0,
    "line_width": 2.0,
    "grid_width": 0.7,
    "marker_size": 6.0,
    "dpi": 220,
    "figure_width": 12.0,
    "figure_height": 5.0,
}


def get_plot_settings() -> dict[str, Any]:
    settings = dict(DEFAULTS)
    try:
        incoming = json.loads(os.environ.get("INCI_PLOT_SETTINGS", "{}"))
    except json.JSONDecodeError:
        incoming = {}
    if isinstance(incoming, dict):
        for key in settings:
            if key in incoming:
                settings[key] = incoming[key]
    return settings


def apply_matplotlib_style() -> dict[str, Any]:
    import matplotlib

    settings = get_plot_settings()
    matplotlib.rcParams.update(
        {
            "font.family": settings["font_family"],
            "font.size": settings["font_size"],
            "axes.titlesize": settings["title_size"],
            "axes.labelsize": settings["font_size"],
            "lines.linewidth": settings["line_width"],
            "lines.markersize": settings["marker_size"],
            "grid.linewidth": settings["grid_width"],
            "legend.fontsize": max(6.0, float(settings["font_size"]) - 1.0),
            "figure.dpi": settings["dpi"],
            "savefig.dpi": settings["dpi"],
        }
    )
    return settings
