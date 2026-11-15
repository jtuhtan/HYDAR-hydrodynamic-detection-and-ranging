"""HYDAR: Hydrodynamic Detection and Ranging.

Simulate objects in water (dropped spheres, boats, ships), record the surface
waves they make with pressure sensor arrays, and locate the source: drops by
time of flight, vessels from their wakes (hydar.locate, hydar.kelvin,
hydar.interferometry), and rocks in rivers from their Karman vortex
streets (hydar.karman).
"""

from .solver import Grid, Solver, pml_profile
from .drop import DropObject, recommend_grid, cap_volume
from .sensors import SensorArray, depth_response, wavenumber
from .locate import (ArrayFix, BoatTrack, intersect_bearings, locate, locate_array,
                     pick_arrivals, pick_wake, track_boat, wake_arrival,
                     KELVIN_HALF_ANGLE)
from .spectral import SpectralSolver
from .boat import Boat

__all__ = [
    "Grid", "Solver", "pml_profile",
    "DropObject", "recommend_grid", "cap_volume",
    "SensorArray", "depth_response", "wavenumber",
    "ArrayFix", "intersect_bearings", "locate", "locate_array", "pick_arrivals",
    "Boat", "BoatTrack", "track_boat", "wake_arrival", "pick_wake", "KELVIN_HALF_ANGLE",
    "SpectralSolver",
]
__version__ = "1.3.0"
