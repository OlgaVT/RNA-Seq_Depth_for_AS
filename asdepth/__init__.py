"""Helpers shared by the figure notebooks of the RNA-Seq depth guidelines paper.

Only file reading and plot formatting live here. Every statistical step stays
in the notebook that produces the figure, so a reader can follow the analysis
without opening this package.
"""

from .summaries import (
    DepthCounts,
    category_of,
    count_deep_sequenced,
    count_per_sample,
    read_depth_counts,
)
from .plotting import PLATFORM_COLORS, SINGLE_COLOR, STACK_COLORS, millions, thousands

__all__ = [
    "DepthCounts",
    "category_of",
    "count_deep_sequenced",
    "count_per_sample",
    "read_depth_counts",
    "PLATFORM_COLORS",
    "SINGLE_COLOR",
    "STACK_COLORS",
    "millions",
    "thousands",
]
