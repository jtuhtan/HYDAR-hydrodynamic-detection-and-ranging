# HYDAR: Hydrodynamic Detection and Ranging

<p align="center"><img src="docs/img/hydar_hero.gif" alt="A drone drops a sphere, the waves wash over four sensor arrays, and red lasers trace the HYDAR fix back to the drop point" width="720"></p>

**HYDAR** locates an object that enters the water using the surface waves it makes. Arrays of pressure sensors record the passing waves, and the **time of flight** of the wave front gives the **direction and distance** from each array's centre to the source. It works the way radar or sonar does, except that the signal is a hydrodynamic (water-wave) disturbance rather than an electromagnetic or acoustic echo.

This repository contains a self-contained Python simulator for developing and testing HYDAR:

1. a 2D shallow-water wave solver on a large rectangular domain whose walls do not reflect, plus a **dispersive** linear potential-flow solver (exact dispersion relation) for boat wakes,
2. a **drop object** module, where a sphere of given diameter and density is dropped at a chosen point and sinks or floats, making surface waves,
3. **pressure sensor arrays** (for example 3 × 3 sensors at 1 m spacing) placed at any x, y and depth z, sampling at 100 Hz,
4. the **HYDAR localization** step: bearing from time differences of arrival, range from time of flight, and fusion across arrays.

<p align="center"><img src="docs/img/fig04_hydar_map.png" alt="HYDAR localization map" width="420"></p>

*A 0.5 m stone sphere is dropped at (11, 17) m in 1 m deep water; the colours show the surface elevation η 3 s later. Four 3 × 3 arrays (A to D, 1 m spacing, at depths of 0.2 to 0.8 m) each estimate a bearing and a time-of-flight range (dashed lines). The fused HYDAR fix (red cross) lands 4 cm from the true drop point (open circle).*

**Created by Jeff Tuhtan** · GitHub [@jtuhtan](https://github.com/jtuhtan) · see [Author and contact](#author-and-contact)

---

## Quick start

```bash
git clone https://github.com/jtuhtan/HYDAR-hydrodynamic-detection-and-ranging.git
cd HYDAR-hydrodynamic-detection-and-ranging
pip install -e ".[test]"

python examples/03_hydar_localization.py   # full pipeline, about 40 s
pytest                                     # 12 tests
```

```python
from hydar import Grid, Solver, DropObject, SensorArray, locate
from hydar.plotting import plot_timeseries, plot_map

# Open-water domain: 30 m x 30 m, 1 m deep, absorbing walls
sim = Solver(Grid(lx=30, ly=30, dx=0.05), depth=1.0)

# Pressure arrays: 3x3 sensors, 1 m spacing, 100 Hz, placed at (x, y, z)
arrays = [SensorArray("A", x=6,  y=6,  z=-0.2),
          SensorArray("B", x=24, y=6,  z=-0.5),
          SensorArray("C", x=24, y=24, z=-0.8)]
for a in arrays:
    sim.add(a)

# Drop a 0.5 m sphere, density 2500 kg/m^3, at (11, 17)
ball = DropObject(diameter=0.5, density=2500, x=11, y=17)
sim.add(ball)
sim.run(12.0)

# HYDAR: direction and range from each array, plus the fused position
result = locate(arrays, c=sim.c, t0=ball.t_release)   # leave out t0 if the drop time is unknown
for fix in result["fixes"]:
    print(fix)                     # bearing, range, position per array
print(result["position"])          # fused estimate

plot_timeseries(arrays, fixes=result["fixes"])   # pressure time series of every sensor
```

## Examples

| Script | What it shows |
|---|---|
| `examples/01_open_boundary.py` | A hump of water spreads out and leaves the domain. After 8 s, 2 × 10⁻⁵ of the energy is left, so the walls act as open water. |
| `examples/02_drop_sphere.py` | A floating sphere (600 kg/m³) heaves and settles at its Archimedes draft. A sinking one (2500 kg/m³) falls to the bottom. Both radiate rings of waves. |
| `examples/03_hydar_localization.py` | The full HYDAR pipeline with four arrays, both with a known drop time and with an unknown one. |
| `examples/04_boat_tracking.py` | Case study 2: a 4.5 m Series 60 boat crosses a 120 m × 50 m domain at 3 m/s, its Kelvin wake is computed with the dispersive solver, and HYDAR tracks it from the wake (about 5 minutes). |
| `examples/05_vessel_localization.py` | Case study 3: a ferry 2.5 km offshore located from its Kelvin wake with one 3 × 3 seabed array, after Rätsep et al. (2020, 2021) (about 1 minute). |
| `examples/make_boat_animation.py` | Renders the case study 2 animation from a real simulation (needs ffmpeg, about 15 minutes). The blue lines mark the Kelvin cusp lines; the red lasers end where the fitted HYDAR track places the boat at the moment the cusp line first reaches each array. |
| `examples/make_hero_animation.py` | Renders the animation at the top of this page from a real simulation (needs ffmpeg, about 5 minutes). The waves, arrival times and the final fix are all computed by `hydar`; the drone and splash are decoration. |
| `docs/make_figures.py` | Regenerates the four figures below (PNG at 600 dpi, plus editable PDF and SVG) in the author's paper figure style (`docs/paperfig.py`), sized for a 170 mm text width. |

![Figure 1](docs/img/fig01_open_boundary.png)

**Figure 1.** Non-reflecting walls. **a–c**, a 1 cm Gaussian hump of water in a 20 m × 12 m, 1 m deep domain spreads out and leaves through the perfectly matched layer. **d**, wave energy left in the domain, relative to the start: it stays at 1 with rigid walls and falls below 10⁻⁴ with the PML.

![Figure 2](docs/img/fig02_drop_sphere.png)

**Figure 2.** Drop objects. A 0.2 m sphere is released 5 cm above 0.5 m deep water. **a**, a floating sphere (600 kg/m³) keeps radiating rings as it bobs; **b**, a sinking sphere (2500 kg/m³) makes a single ring. **c**, height of the sphere's lowest point: the floating sphere settles at its Archimedes draft (dashed), the sinking sphere stops on the bottom. **d**, surface elevation at the wave gauge 1 m from the drop (triangle in a and b).

![Figure 3](docs/img/fig03_hydar_timeseries.png)

**Figure 3.** Dynamic pressure recorded at 100 Hz by the nine sensors of each array (**a–d** = arrays A–D) for the drop in the map above. The shaded band spans the first to the last picked arrival across each array; this spread in arrival time gives the bearing, and the arrival time at the array centre gives the range.

---

## Case study 2: tracking a small boat from its Kelvin wake

<p align="center"><img src="docs/img/hydar_boat.gif" alt="A Series 60 boat crosses the domain, the cusp lines of its Kelvin wake sweep over the four arrays, and red lasers trace the HYDAR track" width="720"></p>

A 4.5 m boat with a **Series 60** hull (beam 0.60 m, draft 0.24 m, 389 L displacement) runs at 3 m/s (about 6 knots, length Froude number 0.45) on a heading of 4° through 3 m deep water. Real boat waves are dispersive: long waves travel faster than short ones, and that is what shapes the familiar **Kelvin wake**, with transverse waves behind the boat and diverging waves inside a wedge of half-angle 19.47°. The shallow-water solver cannot make this pattern, so the case study uses `SpectralSolver`, a linear potential-flow model with the exact dispersion relation (Theory, section 1b). Four 3 × 3 arrays, at depths of 0.3 to 0.8 m, read the dynamic pressure at their own depth at 100 Hz.

HYDAR picks the moment the wake's **cusp line** passes each of the 36 sensors and fits a straight, constant-speed track to those times.

| | True | HYDAR |
|---|---|---|
| Speed | 3.00 m/s | 3.06 m/s |
| Heading | 4.0° | 3.7° |
| Distance from the hull centre at t = 20, 30, 40 s | | 0.7, 1.3, 2.0 m |
| RMS misfit of the 36 picks | | 115 ms |

```python
from hydar import Boat, Grid, SensorArray, SpectralSolver, track_boat
from hydar.boat import Series60Hull

grid = Grid(lx=120.0, ly=50.0, dx=0.125, n_pml=80)        # 10 m absorbing sponge
sim = SpectralSolver(grid, depth=3.0, dt=0.01)
arrays = [SensorArray("A", 50.0, 14.0, z=-0.3, pressure_model="solver"), ...]
for a in arrays:
    sim.add(a)
sim.add(Boat(hull=Series60Hull(length=4.5), speed=3.0, x=-8.0, y=22.0, heading_deg=4.0, ramp=4.0))
sim.run(44.0)
track = track_boat(arrays)            # Kelvin wake: speed, heading, position at t = 0
print(track.position(30.0))           # where the boat was at t = 30 s
```

![Figure 4](docs/img/fig05_series60_hull.png)

**Figure 4.** The Series 60 hull used in the boat case study, scaled to L = 4.5 m. **a**, half-breadth plan with waterlines at 0, 25, 50, 75 and 95 % of the draft (dark to light). **b**, body plan: stations 1 to 4 of the aft body on the left and stations 6 to 9 of the fore body on the right, with the midship section (black) on both sides. As built, the hull has C_B = 0.600, C_P = 0.614, C_M = 0.977, C_WP = 0.706 and its LCB 1.5 % L aft of midships.

![Figure 5](docs/img/fig06_boat_tracking.png)

**Figure 5.** HYDAR tracking of the boat from its Kelvin wake. **a**, surface elevation 30 s after the start, computed with the dispersive solver: transverse waves behind the boat and diverging waves along the 19.47° cusp lines. The HYDAR track (dashed red) lies on the true track (black). **b**, dynamic pressure at the nine sensors of array A (0.3 m deep), with the picked cusp-line passages (red). **c**, picked cusp passages at all 36 sensors against those predicted by the fitted track; the RMS misfit is 115 ms.

**How the boat tracking works.** In deep water, a boat's wake is bounded by two cusp lines at the Kelvin angle μ = arcsin(1/3) = 19.47°, whatever the boat's speed. A point a distance s along the track and d to the side is crossed by the cusp line at t = (s + d cot μ) / U = (s + 2√2 d) / U. Across a cusp line the wave amplitude follows an Airy function: it grows from zero ahead of the line, reaches Ai(0) = 0.355 on it and peaks at 0.536 just behind it. `pick_wake` therefore band-passes each record, takes its envelope, and picks the time the envelope first reaches 0.355/0.536 = 66 % of its first maximum. A simple first-threshold pick does not work for a wake, because the oscillating signal crosses any threshold on a different crest at each sensor. `track_boat` then fits the boat's position at t = 0, its heading and its speed to the 36 picks by least squares from several starting guesses. The start-up of the boat sends a weak train of long waves ahead of it, so the boat starts inside the absorbing layer and the arrays sit at least 50 m down the track, where the wake is fully developed.

---

## Case study 3: locating a passing ferry from its Kelvin wake

This case study follows the field work of Rätsep et al. ([2020](https://doi.org/10.1175/JTECH-D-19-0192.1), [2021](https://doi.org/10.1016/j.oceaneng.2021.109156)) in Tallinn Bay. A 3 × 3 grid of seabed pressure sensors (2.5 m spacing, 0.2 m above the bed in 3 m of water, 100 Hz) records the wakes of fast ferries passing about 3 km offshore. These ships are slower than the long-wave speed of the deep fairway, so they make a **dispersive Kelvin wake**. Simulating a 3 km wide domain with the dispersive solver of case study 2 would be very costly, so the wake is synthesized from linear deep-water theory by integrating the Havelock wave spectrum numerically, so the cusp and both wave branches are exact. A random wind sea (H_s = 0.15 m, 0.4 Hz) and sensor noise are added, and the pressure at the sensor depth uses the Airy depth response.

In the example, a ferry sails at 11 m/s (21 kn) on a course of 159°, 2.5 km from the array. The analysis sees only the recorded pressures and follows the steps of the papers:

1. **Detection.** Spectrograms of five sensors (four corners and the centre) are averaged, which suppresses incoherent noise ([2020](https://doi.org/10.1175/JTECH-D-19-0192.1)). The wake is flagged where its band power rises above the background.
2. **Speed, distance and time abeam from the spectrogram.** A fixed observer a distance d from the track sees two branches: transverse waves whose frequency falls towards f_t = g/(2πU), and divergent waves whose frequency rises. Both start at the cusp, at √(3/2)·f_t. Fitting the two stationary-phase branches to the spectrogram ridge gives U, d and the time the ship passed abeam. The upper divergent branch is hidden under the wind sea, as noted in the papers, so the fit uses the band below it.
3. **Wake direction from time delays.** Cross-correlating all 36 sensor pairs gives the delay of the same waves between sensors. A plane-wave fit to those delays gives the direction the cusp waves travel ([2021](https://doi.org/10.1016/j.oceaneng.2021.109156)).
4. **Course and position.** At the cusp, the waves travel at the Kelvin angle 35°16′ to the ship's course, which leaves two candidate courses. In the field this ambiguity was resolved with the coastline. Here the direction of the late transverse waves, which travel along the course, picks the right one. The distance d and the side of the track the waves came from then place the ship.

| | True | HYDAR (case study 3) |
|---|---|---|
| Speed | 11.00 m/s | 10.99 m/s |
| Distance from track | 2500 m | 2390 m |
| Course | −69.0° (159° nautical) | −67.1° |
| Ship position when the cusp reaches the array | | 463 m off |

For comparison, the field study reported positions mostly within 0.5 km of the AIS positions, with about 400 m in its example case.

```python
from hydar import SensorArray
from hydar.kelvin import KelvinShip, record_passage, locate_vessel

array = SensorArray("H", 0, 0, z=-2.8, spacing=2.5, rotation_deg=22.5, pressure_model="airy")
ship = KelvinShip(speed=11.0, course_deg=-69.0, x0=..., y0=...)     # see examples/05
record_passage([array], ship, t_end=1800, depth=3.0, wind={"hs": 0.15, "f_peak": 0.4})
fix = locate_vessel(array)       # speed, distance, course, ship position at any time
```

![Figure 6](docs/img/fig07_vessel_localization.png)

**Figure 6.** Case study 3, after Rätsep et al. (2020, 2021). **a**, dynamic pressure at the centre sensor: wind-sea noise, then the ferry wake 15 min after the ship passed abeam. **b**, spectrogram averaged over five sensors. The fitted Kelvin branches (red) follow the ridge points (cyan); the upper divergent branch disappears under the wind sea. **c**, the true track (black) and the HYDAR track (dashed red) from the single 3 × 3 array (blue). The arrow shows the measured direction of the cusp waves; the markers show the ship when the cusp reached the array, true (circle) and estimated (cross).

**Assumptions.** Deep water is assumed along the whole wave path, so refraction and shoaling into the 3 m deep nearshore are ignored, and the hull is represented by a smooth amplitude spectrum rather than a specific ship. The spectrogram relations used here are derived in `hydar/kelvin.py` from stationary phase; the papers state the same physics in terms of the cusp and limiting transverse frequencies.

---

## Theory

### 1. Wave model

The domain is the horizontal free surface. The unknowns are the surface elevation η(x, y, t) and the depth-averaged velocity (u, v). For still water of uniform depth h, the linear shallow-water equations are

$$
\frac{\partial \eta}{\partial t} = -h\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right) + q(x,y,t),\qquad
\frac{\partial u}{\partial t} = -g\frac{\partial \eta}{\partial x},\qquad
\frac{\partial v}{\partial t} = -g\frac{\partial \eta}{\partial y}
$$

Waves travel at $c=\sqrt{gh}$, and $q$ is a volume source per unit area that objects use to push water aside. The equations are solved with second-order finite differences on a staggered (Arakawa C / Yee) grid with a time step $\Delta t = \mathrm{CFL}\,\Delta x / c$, using CFL = 0.5 by default.

### 1b. Dispersive solver (boat wakes)

`SpectralSolver` keeps the full linear free-surface physics for inviscid, irrotational flow. In Fourier space, every wavenumber k of the surface elevation $\hat\eta$ and of the surface velocity potential $\hat\phi$ is a harmonic oscillator,

$$
\frac{\partial \hat\eta}{\partial t} = k\tanh(kh)\,\hat\phi,\qquad
\frac{\partial \hat\phi}{\partial t} = -g\,\hat\eta - \frac{\hat P}{\rho},
$$

with the exact dispersion relation $\omega^2 = gk\tanh(kh)$, valid from shallow to deep water. Each oscillator is advanced exactly over a time step, with the applied surface pressure P held at its mid-step value, so the solver has no stability limit and no numerical dispersion. The boat is the classical **pressure patch** of linear ship-wave theory (Wehausen and Laitone, 1960), $P=\rho g D(x,y,t)$, where D is the hull draft. Sensors with `pressure_model="solver"` read the dynamic pressure at their own depth directly, $p_{dyn}/(\rho g)=\mathcal F^{-1}[(\hat\eta+\hat P/\rho g)\cosh k(h+z)/\cosh kh]$. The FFT makes the domain periodic, so the walls are a **sponge layer** in which η and φ are damped by $e^{-\sigma\Delta t}$, with σ growing quadratically to 3 s⁻¹ at the edge. Tests check the standing-wave frequency against the dispersion relation and that less than 10⁻³ of the wave energy is left once the waves have left the domain.

### 2. Non-reflecting walls

The physical domain is wrapped in a split-field **perfectly matched layer** (PML; Bérenger, 1994), 40 cells thick by default. Inside the layer, η is split into $\eta_x+\eta_y$, and each part is damped only in its own direction by $\sigma(d)=\sigma_{\max}(d/L)^2$, where $\sigma_{\max} = -3c\ln R/(2L)$ and R = 10⁻⁶. Waves enter the layer without reflection and decay there, so a finite box behaves like an unbounded lake. The test suite checks that less than 10⁻⁴ of the wave energy remains after the waves have crossed the domain twice. With rigid walls, over 30 % remains.

### 3. Drop object

A sphere of diameter D and density ρₛ is released from rest at (x, y), optionally from a height above the surface and at a chosen release time. Several drops can share one simulation. The sphere moves vertically under gravity, buoyancy, added mass and quadratic drag:

$$
(m + C_a\rho_w V_{sub})\,\dot w = \rho_w g V_{sub} - m g - \tfrac12 \rho_w C_d A\, w|w|
$$

$V_{sub}$ is the spherical-cap volume below the local water level, with $C_a = 0.5$ and $C_d = 0.5$ by default. As the sphere crosses the surface, the change in submerged volume, $\dot V_{sub}$, is injected as the source q over the waterline disc. The volume of water raised in the domain therefore equals the volume the sphere displaces, which the tests verify to within 2 %. The water level around the sphere is fed back into $V_{sub}$, so a floating sphere also responds to the waves it radiates (two-way coupling). Spheres lighter than water settle at the Archimedes draft. Heavier ones sink and stop on the bottom at $z = -h + D/2$.

**Resolution.** The sphere should span at least 8 to 10 grid cells (Δx ≤ D/10 is recommended), and the code warns if it spans fewer. `recommend_grid(D, h)` suggests Δx and a domain size.

### 4. Pressure sensors

A `SensorArray` is an n_x × n_y grid of sensors (3 × 3 at 1 m by default) centred at (x, y) and mounted at depth z < 0. η is interpolated bilinearly to each sensor at every solver step and resampled at `fs` = 100 Hz. The solver step must be shorter than 1/fs, and the code checks this. Gauge pressure is the hydrostatic head plus a dynamic part:

$$
p(t) = \rho g(-z) + p_{dyn}(t)
$$

* `pressure_model="hydrostatic"` (default): $p_{dyn}=\rho g\eta$. This matches the shallow-water solver and is used for HYDAR picking.
* `pressure_model="airy"`: each frequency component is attenuated with depth by linear (Airy) wave theory, $K_p=\cosh k(h+z)/\cosh kh$ with $\omega^2 = gk\tanh kh$. It shows how strongly a deep array filters short waves. Because the filter is zero-phase, it smears sharp fronts slightly earlier in time, so do not use it for arrival picking.

### 5. Boat

`Boat` moves a hull along a straight line at constant speed. In the dispersive `SpectralSolver` the hull's depth below the waterline, D(x, y, t), acts as the surface pressure head P/(ρg) (section 1b). In the shallow-water `Solver`, the rate of change q = ∂D/∂t is injected as the volume source instead. Water is pushed out ahead of the bow and flows back behind the stern; the source integrates to zero, and the water raised around the hull equals its displaced volume (tested to within 3 %). The boat can start inside the absorbing layer so that it enters the domain already at speed.

The default `Series60Hull` follows the Series 60 parent form (Todd, 1963): L/B = 7.5, B/T = 2.5, C_B = 0.60, C_P = 0.614, C_M = 0.977, C_WP = 0.706, LCB 1.5 % L aft. The published offset table is not reproduced. Instead, smooth curves are fitted to exactly these coefficients:

* the sectional-area curve $A(\xi)=A_M(1-|\xi|^{n})$, with separate fore and aft exponents that set C_P and the LCB;
* the waterline $b(\xi)=\tfrac{B}{2}(1-|\xi|^{m})$, with m set by C_WP;
* section depth profiles $D = T\,(1-|y/b|^{k})$, with k chosen to give each section its area.

Principal dimensions, displacement, prismatic distribution and LCB therefore match Series 60, while the detailed section shapes are approximate. `Series60Hull.coefficients()` reports the values achieved, and a test checks them to ±0.003. An `EllipseHull` is also available.

### 6. HYDAR localization

Each array acts as a small sonar receiver.

1. **Arrival picking.** At each sensor, the arrival time $t_i$ is the first time $|p_{dyn}|$ exceeds a fraction (10 % by default) of the array's peak, interpolated between samples.
2. **Direction from time differences of arrival.** Across a small array the front is nearly plane, so $t_i = \tau + \mathbf{s}\cdot(\mathbf r_i-\mathbf r_c)$. A linear least-squares fit gives the slowness vector **s**, and the source lies in the direction $-\mathbf s/|\mathbf s|$ from the array centre. The fitted apparent speed $1/|\mathbf s|$ should match $c$, which serves as a quality check.
3. **Range from time of flight.** When the drop time $t_0$ is known, $r = c\,(\tau - t_0)$, as in radar ranging. Bearing and range together give a position from a single array.
4. **Fusion across arrays.** All arrival times are fitted together to a circular front, $t_i = t_0 + |\mathbf r_i - \mathbf x_s|/c$, solving for the source position $\mathbf x_s$ and, if it is unknown, $t_0$. Bearing-line intersection provides the initial guess. Because $t_0$ can be fitted, a constant picking delay cancels out, and the method also works when nobody knows when the object entered the water.

In example 3, every bearing is within 0.4° of the truth, the single-array ranges are within 0.2 m over distances of 8 to 17 m, and the fused position is within 0.05 m, whether or not the drop time is known.

---

## Assumptions and limitations

* **Shallow-water physics (drop case).** The drop examples use the shallow-water `Solver`, which is linear and non-dispersive, so it is valid for waves much longer than the depth and for amplitudes much smaller than the depth. Real waves from objects smaller than the water depth are dispersive: longer waves outrun shorter ones, and the front spreads out. The dispersive `SpectralSolver` (used for the boat) could carry the drop case too, but arrival picking would then have to account for the group velocity.
* Uniform depth, flat bottom, no current, no wind waves or sensor noise.
* The sphere moves vertically only and is a source only while it crosses the surface. A fully submerged moving body does not radiate in this model.
* Sensors are ideal point gauges.
* **Boats.** Case study 2 uses the dispersive solver, which gives the Kelvin wake, but it is still linear: the hull is a pressure patch rather than a body boundary, and wave breaking and nonlinear steepening near the hull are left out. The tracking fit assumes a straight track at constant speed and the deep-water Kelvin angle, which holds while the depth Froude number U/√(gh) stays below about 0.6 (it is 0.55 here). In the shallow-water `Solver`, only supercritical boats (U > √(gh)) make a wake, a Mach cone, and `track_boat(..., wake="mach")` tracks that case.

## Repository layout

```
hydar/
  solver.py     2D shallow-water solver with PML walls
  spectral.py   dispersive linear potential-flow solver with sponge walls (boat wakes)
  drop.py       DropObject (sphere sinks or floats, two-way coupled source)
  boat.py       Boat with Series 60 (or elliptical) hull moving on a straight track
  kelvin.py     Kelvin wake synthesis, wind sea, spectrogram chirp fit, wake direction,
                vessel localization (case study 3)
  sensors.py    SensorArray (3x3 at 1 m, 100 Hz, x/y/z placement, depth response)
  locate.py     HYDAR: arrival picking, TDOA bearing, TOF range, multi-array fusion,
                and boat tracking from the Kelvin or Mach wake (pick_wake, track_boat)
  plotting.py   time-series and map plots
examples/       runnable scripts (01 to 05) and the two animations
tests/          pytest suite
docs/img/       figures used in this README
```

## Author and contact

HYDAR was conceived and created by **Jeff Tuhtan**.

* GitHub: [@jtuhtan](https://github.com/jtuhtan)
* Questions, ideas and bug reports: please [open an issue](https://github.com/jtuhtan/HYDAR-hydrodynamic-detection-and-ranging/issues)

If you use HYDAR in your work, please credit the author and link to this repository.

## References

Bérenger, J.-P. (1994). A perfectly matched layer for the absorption of electromagnetic waves. *Journal of Computational Physics*, 114(2), 185–200. https://doi.org/10.1006/jcph.1994.1159

Rätsep, M., Parnell, K. E., Soomere, T., Kruusmaa, M., Ristolainen, A., & Tuhtan, J. A. (2020). Using spectrograms from underwater total pressure sensors to detect passing vessels in a coastal environment. *Journal of Atmospheric and Oceanic Technology*, 37(8), 1353–1363. https://doi.org/10.1175/JTECH-D-19-0192.1

Rätsep, M., Parnell, K. E., Soomere, T., Kruusmaa, M., Ristolainen, A., & Tuhtan, J. A. (2021). Surface vessel localization from wake measurements using an array of pressure sensors in the littoral zone. *Ocean Engineering*, 233, 109156. https://doi.org/10.1016/j.oceaneng.2021.109156

Todd, F. H. (1963). *Series 60: Methodical experiments with models of single-screw merchant ships*. David Taylor Model Basin Report 1712, Washington, DC. DTIC AD0419990, https://apps.dtic.mil/sti/html/tr/AD0419990/index.html

Wehausen, J. V., & Laitone, E. V. (1960). Surface waves. In S. Flügge (Ed.), *Handbuch der Physik*, Vol. 9, 446–778. Springer, Berlin. https://doi.org/10.1007/978-3-642-45944-3_6
