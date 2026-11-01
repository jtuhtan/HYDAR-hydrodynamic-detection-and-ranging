"""HYDAR: Hydrodynamic Detection and Ranging.

Simulate objects dropped into water, record the surface waves they make with
pressure sensor arrays, and locate the drop by time of flight.
"""

from .solver import Grid, Solver, pml_profile
from .drop import DropObject, recommend_grid, cap_volume
from .sensors import SensorArray, depth_response, wavenumber
from .locate import (ArrayFix, BoatTrack, intersect_bearings, locate, locate_array,
                     pick_arrivals, track_boat, wake_arrival)
from .boat import Boat

__all__ = [
    "Grid", "Solver", "pml_profile",
    "DropObject", "recommend_grid", "cap_volume",
    "SensorArray", "depth_response", "wavenumber",
    "ArrayFix", "intersect_bearings", "locate", "locate_array", "pick_arrivals",
    "Boat", "BoatTrack", "track_boat", "wake_arrival",
]
__version__ = "0.1.0"
