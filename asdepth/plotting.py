"""Axis formatters and the colours used across the figures."""

# Blue: detected in all 100 subsamples. Yellow: detected in fewer than 100.
STACK_COLORS = {"robust": "#1E88E5", "sporadic": "#FFC107"}
SINGLE_COLOR = "#004D40"

# Sequencing platforms priced by CCGA Kiel.
PLATFORM_COLORS = {"SP": "#E69F00", "S1": "#0072B2", "S4": "#D55E00"}


def thousands(value, _tick=None):
    """Tick label in thousands. Ticks that are not a multiple of 1000 stay blank,
    which is how the published axes were drawn."""
    if value % 1000 == 0 and value != 0:
        return "{:.0f}K".format(value / 1000)
    return ""


def millions(value, _tick=None):
    """Tick label in millions, blank unless the tick is a multiple of 10M."""
    if value % 10_000_000 == 0:
        return "{:.0f}M".format(value / 1_000_000)
    return ""
