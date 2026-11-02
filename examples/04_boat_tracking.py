"""Case study 2: track a small boat from its Kelvin wake.

A 4.5 m Series 60 hull (beam 0.6 m, draft 0.24 m, 389 L displacement) runs
at 3 m/s (about 6 knots, length Froude number 0.45) through 3 m deep water.
Its wake is computed with the dispersive SpectralSolver (linear potential
flow, exact dispersion relation), so it has the real Kelvin pattern:
transverse waves behind the boat and diverging waves inside a 19.5 deg
wedge. Four 3x3 arrays record the pressure at 100 Hz at their own depth;
HYDAR fits the boat's track, speed and heading to the times the wake's cusp
line passes each sensor.

Run:  python examples/04_boat_tracking.py   ->  boat_timeseries.png, boat_map.png
Takes about 5 minutes.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from hydar import Boat, Grid, SensorArray, SpectralSolver, track_boat
from hydar.boat import Series60Hull
from hydar.plotting import plot_timeseries

grid = Grid(lx=120.0, ly=50.0, dx=0.125, n_pml=80)   # 10 m absorbing sponge on every side
sim = SpectralSolver(grid, depth=3.0, dt=0.01)
print(f"grid {grid.nx} x {grid.ny}, dt = {sim.dt * 1e3:.0f} ms")

# "solver": each sensor reads the dynamic pressure at its own depth directly
arrays = [SensorArray("A", 50.0, 14.0, z=-0.3, pressure_model="solver"),
          SensorArray("B", 75.0, 16.0, z=-0.5, pressure_model="solver"),
          SensorArray("C", 55.0, 36.0, z=-0.5, pressure_model="solver"),
          SensorArray("D", 80.0, 40.0, z=-0.8, pressure_model="solver")]
for a in arrays:
    sim.add(a)

# The boat starts inside the sponge on the left and fades in over 4 s, so the
# start-up waves are absorbed before they reach the arrays.
hull = Series60Hull(length=4.5)
print("hull coefficients:", {k: round(float(v), 3) for k, v in hull.coefficients().items()})
boat = Boat(hull=hull, speed=3.0, x=-8.0, y=22.0, heading_deg=4.0, ramp=4.0)
sim.add(boat)
print(f"boat: Fr_L = {boat.speed / np.sqrt(9.81 * boat.length):.2f}, "
      f"Fr_h = {boat.froude(sim.depth):.2f}, displaced volume {boat.volume * 1e3:.0f} L")

sim.run(30.0)
snapshot = sim.eta[grid.interior].copy()
sim.run(44.0)

track = track_boat(arrays)
print("HYDAR track:", track)
print(f"true: speed {boat.speed} m/s, heading {boat.heading_deg} deg")
t_check = np.array([20.0, 30.0, 40.0])
err = np.linalg.norm(track.position(t_check) - boat.position(t_check), axis=1)
for t, e in zip(t_check, err):
    print(f"  distance from the true hull centre at t = {t:.0f} s: {e:.2f} m")

plot_timeseries(arrays)
plt.tight_layout()
plt.savefig("boat_timeseries.png", dpi=120)

fig, ax = plt.subplots(figsize=(10, 5))
ax.imshow(snapshot.T, origin="lower", extent=grid.extent, cmap="RdBu_r", vmin=-0.04, vmax=0.04)
tt = np.linspace(0, 44, 50)
ax.plot(*boat.position(tt).T, "k-", lw=1, label="true track")
ax.plot(*track.position(tt).T, "r--", lw=1.5, label="HYDAR track")
for a in arrays:
    ax.plot(*a.positions.T, "k.", ms=3)
ax.set(xlim=(0, grid.lx), ylim=(0, grid.ly), xlabel="x [m]", ylabel="y [m]")
ax.legend()
plt.savefig("boat_map.png", dpi=120)
print("saved boat_timeseries.png and boat_map.png")
