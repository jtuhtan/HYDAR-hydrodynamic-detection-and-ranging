"""Case study 4: design a 3 x 3 array for correlative interferometry and track a
ship without Kelvin wake theory.

1. design_array() finds the 3 x 3 grid (free column and row positions) with
   the smallest mean bearing error bound over the wake band, subject to a
   sidelobe limit (no ambiguous directions) and a footprint limit.
2. The ferry passage of case study 3 is recorded by the Rätsep et al. grid
   (2.5 m spacing) and by the designed grid, in the same wind sea.
3. Both records are analysed with the Kelvin method of case study 3
   (locate_vessel) and with correlative interferometry (track_vessel), which
   uses only the dispersion relation, wave stationarity and rays.

Run:  python examples/06_array_design.py   (about 2 minutes)
"""

import numpy as np

from hydar import SensorArray
from hydar.interferometry import design_array, evaluate_layout, grid_layout, track_vessel
from hydar.kelvin import KelvinShip, locate_vessel, record_passage

DEPTH = 3.0
BAND = (0.07, 0.3)            # wake band [Hz]; deep-water waves, as synthesized in hydar.kelvin

# 1. Array design ----------------------------------------------------------
base = grid_layout([-2.5, 0, 2.5], [-2.5, 0, 2.5])
crlb0, psl0 = evaluate_layout(base, BAND)
print(f"Rätsep grid, 2.5 m spacing: mean bearing bound {crlb0:.1f} deg at SNR 1, peak sidelobe {psl0:.2f}")
designs = {"same footprint (7.1 m)": design_array(BAND, max_aperture=7.08, min_spacing=0.5),
           "20 m footprint": design_array(BAND, max_aperture=20.0, min_spacing=0.5)}
for name, d in designs.items():
    print(f"design, {name}: {d}")
    print(f"   bearing bound {d.mean_crlb_deg / crlb0:.2f} x that of the Rätsep grid")

# 2. Ferry passage -----------------------------------------------------------
course = 90.0 - 159.0
e = np.array([np.cos(np.radians(course)), np.sin(np.radians(course))])
n_port = np.array([-e[1], e[0]])
d_true, U_true, t_abeam = 2500.0, 11.0, 300.0
p0 = -d_true * n_port - U_true * t_abeam * e
ship = KelvinShip(speed=U_true, course_deg=course, x0=p0[0], y0=p0[1], amplitude=0.25)

layouts = {"Rätsep grid": None, **{f"design, {k}": d.offsets for k, d in designs.items()}}
tt = np.linspace(0, 1800, 61)
for name, lay in layouts.items():
    kw = {"layout": lay} if lay is not None else {"spacing": 2.5}
    arr = SensorArray("H", 0.0, 0.0, z=-(DEPTH - 0.2), rotation_deg=22.5, pressure_model="airy", **kw)
    record_passage([arr], ship, t_end=1800.0, depth=DEPTH,
                   wind={"hs": 0.3, "f_peak": 0.4, "mean_dir_deg": 112.5})
    print(f"\n{name}  (wind sea Hs = 0.3 m)")
    print(f"  true : speed {U_true:.2f} m/s, course {course:.1f} deg")
    fix = locate_vessel(arr)
    ci = track_vessel([arr], BAND, depth=np.inf, window=120.0, hop=20.0)
    for label, spd, crs, pos in (("Kelvin (case 3)", fix.speed, fix.course_deg, fix.position),
                                 ("interferometry ", ci.speed, ci.heading_deg, ci.position)):
        rms = np.sqrt(np.mean(np.sum((pos(tt) - ship.position(tt)) ** 2, axis=1)))
        print(f"  {label}: speed {spd:.2f} m/s, course {crs:.1f} deg, "
              f"rms track error over 30 min {rms:.0f} m")
