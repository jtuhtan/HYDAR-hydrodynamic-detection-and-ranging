"""Case study 2: track a small boat crossing the domain from its wake.

A 4.5 m Series 60 hull (beam 0.6 m, draft 0.24 m, 389 L displacement) runs
at 5 m/s through 1 m deep water. The shallow-water wave speed
is sqrt(g h) = 3.13 m/s, so the boat is supercritical (Fr = 1.6) and drags a
V-shaped wake. Four 3x3 arrays record the wake at 100 Hz; HYDAR fits the
boat's track, speed and heading to the wake arrival times.

Run:  python examples/04_boat_tracking.py   ->  boat_timeseries.png, boat_map.png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from hydar import Boat, Grid, SensorArray, Solver, track_boat
from hydar.boat import Series60Hull
from hydar.plotting import plot_timeseries

grid = Grid(lx=40.0, ly=30.0, dx=0.05, n_pml=60)   # beam spans 12 cells
sim = Solver(grid, depth=1.0)
print(f"grid {grid.nx} x {grid.ny}, dt = {sim.dt * 1e3:.1f} ms, c = {sim.c:.2f} m/s")

arrays = [SensorArray("A", 10.0, 4.0, z=-0.3), SensorArray("B", 28.0, 5.0, z=-0.5),
          SensorArray("C", 12.0, 25.0, z=-0.5), SensorArray("D", 31.0, 25.0, z=-0.8)]
for a in arrays:
    sim.add(a)

# The boat starts inside the absorbing layer left of the domain, so it
# enters already up to speed.
hull = Series60Hull(length=4.5)
print("hull coefficients:", {k: round(float(v), 3) for k, v in hull.coefficients().items()})
boat = Boat(hull=hull, speed=5.0, x=-2.5, y=12.0, heading_deg=10.0)
sim.add(boat)
print(f"boat: Fr = {boat.froude(sim.depth):.2f}, displaced volume {boat.volume * 1e3:.0f} L")

sim.run(5.0)
snapshot = sim.eta[grid.interior].copy()
sim.run(12.0)

track = track_boat(arrays, c=sim.c)
print("HYDAR track:", track)
true_p0 = boat.position(0.0)
t_check = np.array([2.0, 5.0, 8.0])
bow = lambda t: boat.position(t) + 0.5 * boat.length * boat.e
err = np.linalg.norm(track.position(t_check) - bow(t_check), axis=1)
print(f"true: speed {boat.speed} m/s, heading {boat.heading_deg} deg, p(0) = {true_p0}")
for t, e in zip(t_check, err):
    print(f"  distance from the true bow at t = {t:.0f} s: {e:.2f} m")

plot_timeseries(arrays)
plt.tight_layout()
plt.savefig("boat_timeseries.png", dpi=120)

fig, ax = plt.subplots(figsize=(8, 6))
v = np.abs(snapshot).max()
ax.imshow(snapshot.T, origin="lower", extent=grid.extent, cmap="RdBu_r", vmin=-v, vmax=v)
tt = np.linspace(0, 10, 50)
ax.plot(*boat.position(tt).T, "k-", lw=1, label="true track")
ax.plot(*track.position(tt).T, "r--", lw=1.5, label="HYDAR track")
for a in arrays:
    ax.plot(*a.positions.T, "k.", ms=3)
ax.set(xlim=(0, grid.lx), ylim=(0, grid.ly), xlabel="x [m]", ylabel="y [m]")
ax.legend()
plt.savefig("boat_map.png", dpi=120)
print("saved boat_timeseries.png and boat_map.png")
