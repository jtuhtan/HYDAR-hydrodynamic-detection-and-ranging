"""Full HYDAR pipeline: drop a sphere, record with four 3x3 arrays, locate it.

Run:  python examples/03_hydar_localization.py
Writes hydar_timeseries.png and hydar_map.png in the current directory.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from hydar import DropObject, Grid, SensorArray, Solver, locate
from hydar.plotting import plot_map, plot_timeseries

# 1. Open-water domain: 30 m x 30 m, 1 m deep, absorbing walls.
grid = Grid(lx=30.0, ly=30.0, dx=0.05, n_pml=40)
sim = Solver(grid, depth=1.0)
print(f"grid {grid.nx} x {grid.ny}, dt = {sim.dt*1e3:.2f} ms, c = {sim.c:.2f} m/s")

# 2. Four 3x3 pressure arrays (1 m spacing, 100 Hz) at different places and depths.
arrays = [
    SensorArray("A", x=6.0, y=6.0, z=-0.2),
    SensorArray("B", x=24.0, y=6.0, z=-0.5),
    SensorArray("C", x=24.0, y=24.0, z=-0.8),
    SensorArray("D", x=6.0, y=24.0, z=-0.5),
]
for a in arrays:
    sim.add(a)

# 3. Drop a 0.5 m stone sphere (2500 kg/m^3) at (11, 17).
ball = DropObject(diameter=0.5, density=2500.0, x=11.0, y=17.0, t_release=0.0)
sim.add(ball)

# 4. Simulate 12 s (long enough for the waves to reach every array).
sim.run(3.0)
eta_snapshot = sim.eta[grid.interior].copy()  # wave field at t = 3 s, for the map
sim.run(12.0)
print(f"sphere came to rest on the bottom: {ball.on_bottom}")

# 5. HYDAR: bearing from TDOA, range from time of flight (drop time known).
c = sim.c
result = locate(arrays, c=c, t0=ball.t_release)
for f in result["fixes"]:
    print(" ", f)
err = np.linalg.norm(result["position"] - [ball.x, ball.y])
print(f"fused position ({result['position'][0]:.2f}, {result['position'][1]:.2f}) m,"
      f" true ({ball.x}, {ball.y}) m, error {err:.2f} m")

# Same, without knowing when the object was dropped (t0 estimated).
blind = locate(arrays, c=c)
err_b = np.linalg.norm(blind["position"] - [ball.x, ball.y])
print(f"unknown drop time: position ({blind['position'][0]:.2f}, {blind['position'][1]:.2f}) m,"
      f" t0 = {blind['t0']:.2f} s, error {err_b:.2f} m")

plot_timeseries(arrays, fixes=result["fixes"])
plt.tight_layout()
plt.savefig("hydar_timeseries.png", dpi=120)
ax = plot_map(sim, arrays, [ball], result, eta=eta_snapshot)
ax.set_title("surface at t = 3 s, bearings (dashed) and HYDAR fix")
plt.tight_layout()
plt.savefig("hydar_map.png", dpi=120)
print("saved hydar_timeseries.png and hydar_map.png")
