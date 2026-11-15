"""Case study 5: a fish locates a rock in a river from its Karman vortex street.

A gravel-bed river with a trapezoidal cross section (4 m bed, 1 : 2 banks,
0.25 m deep, slope 0.002, Manning n 0.035) flows past a hemispherical rock
that sticks out of the water by half its height. A 25 cm fish holds station
6 m downstream with its head upstream. Its lateral line has 10 pressure
sensors: 5 on the head third of the body, 5 on the middle third, none on the
tail third. The fish steps across the river (19 stations, 0.2 m apart, 60 s
at each) and, from the pressures alone, estimates the size of the rock and
where it is.

Run:  python examples/08_fish_rock.py   (about 15 s)
"""

import numpy as np

from hydar.karman import (Fish, HemisphereRock, KarmanStreet, RiverTurbulence,
                          TrapezoidalChannel, lateral_survey, locate_rock)

channel = TrapezoidalChannel(bed_width=4.0, side_slope=2.0, depth=0.25, slope=0.002, manning_n=0.035)
rock = HemisphereRock.from_protrusion(channel, x=0.0, y=0.4, protrusion=0.5)
street = KarmanStreet(channel, rock, phase=3.0)
fish = Fish(x=6.0 - 0.25 / 3, y=0.0)          # sensor line centred 6 m behind the rock

print(f"River: top width {channel.top_width:.1f} m, mean velocity {channel.mean_velocity:.2f} m/s, "
      f"discharge {channel.discharge:.2f} m3/s, Froude {channel.froude:.2f}")
print(f"Rock: radius {rock.radius:.2f} m, waterline diameter {rock.waterline_diameter(street.depth):.2f} m, "
      f"shedding diameter D {street.D:.2f} m")
print(f"Street: U {street.U:.2f} m/s, f {street.f:.3f} Hz, U_c {street.Uc:.2f} m/s, "
      f"vortex spacing {street.wavelength:.2f} m x {street.row_spacing:.2f} m, "
      f"Gamma {street.gamma0:.2f} m2/s, Re {street.reynolds:.1e}, shallow-wake S {street.stability_number:.3f}")
print("Fish sensors, distance from the snout (cm):", np.round(100 * fish.sensor_offsets, 1))

turbulence = RiverTurbulence(street.U, intensity=0.10, scale=channel.depth, seed=1)
stations = np.round(np.arange(-1.8, 1.81, 0.2), 2)
survey = lateral_survey(street, fish.x, stations, hold=60.0, turbulence=turbulence, seed=1)
print(f"Survey: {len(stations)} stations x 60 s, turbulence {turbulence.rms:.1f} Pa RMS, 2 Pa sensor noise")

fix = locate_rock(survey)
true_dist = survey.x_center - rock.x
print(f"\nDetected: {fix.detected}, street comes from upstream: {fix.upstream}")
print(f"             estimate   true")
print(f"f (Hz)       {fix.f:8.4f} {street.f:8.4f}")
print(f"U_c (m/s)    {fix.Uc:8.3f} {street.Uc:8.3f}")
print(f"D (m)        {fix.diameter:8.3f} {street.D:8.3f}")
print(f"distance (m) {fix.distance:8.2f} {true_dist:8.2f}")
print(f"rock x (m)   {fix.x:8.2f} {rock.x:8.2f}")
print(f"rock y (m)   {fix.y:8.2f} {rock.y:8.2f}")
