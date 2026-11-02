"""Case study 3: locate a passing ferry from its Kelvin wake.

After Rätsep et al. (2020, 2021): a 3 x 3 grid of seabed pressure sensors
(2.5 m spacing, 0.2 m above the bed in 3 m of water, 100 Hz) records the wake
of a ferry sailing at 11 m/s (21 kn) on a course of 159 deg, 2.5 km offshore.
A wind sea (0.4 Hz) and sensor noise are added. From the pressures alone:

  spectrogram + 5-sensor stacking  -> detect the wake
  Kelvin chirp fit                 -> speed, distance from track, time abeam
  cross-correlation delays         -> wake direction -> course (cusp angle 35.26 deg)
  distance + course + side         -> ship track and position

Run:  python examples/05_vessel_localization.py   (about 1 minute)
"""

import numpy as np

from hydar import SensorArray
from hydar.kelvin import KelvinShip, locate_vessel, record_passage

DEPTH = 3.0
array = SensorArray("H", x=0.0, y=0.0, z=-(DEPTH - 0.2), spacing=2.5, fs=100.0,
                    rotation_deg=22.5, pressure_model="airy")

# Ship: nautical course 159 deg = 90 - 159 = -69 deg counter-clockwise from east.
course = 90.0 - 159.0
e = np.array([np.cos(np.radians(course)), np.sin(np.radians(course))])
n_port = np.array([-e[1], e[0]])
d_true, U_true, t_abeam = 2500.0, 11.0, 300.0
p_abeam = -d_true * n_port            # the array lies 2.5 km to port of the track
p0 = p_abeam - U_true * t_abeam * e
ship = KelvinShip(speed=U_true, course_deg=course, x0=p0[0], y0=p0[1], amplitude=0.25)

print("synthesizing wake, wind sea and noise ...")
record_passage([array], ship, t_end=1800.0, depth=DEPTH,
               wind={"hs": 0.15, "f_peak": 0.4, "mean_dir_deg": 112.5})

fix = locate_vessel(array)
print("HYDAR:", fix)
print(f"true : speed {U_true:.2f} m/s, distance {d_true:.0f} m, course {course:.1f} deg, "
      f"abeam at t = {t_abeam:.0f} s")
for t in (fix.t_abeam, fix.t_cusp):
    err = np.linalg.norm(fix.position(t) - ship.position(t))
    print(f"  ship position error at t = {t:.0f} s: {err:.0f} m")
print(f"  course candidates from the cusp angle: {fix.candidates_deg[0]:.1f} and "
      f"{fix.candidates_deg[1]:.1f} deg; transverse waves pick {fix.course_deg:.1f} deg")
