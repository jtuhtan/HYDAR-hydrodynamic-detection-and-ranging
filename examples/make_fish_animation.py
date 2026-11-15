"""Render the animation for case study 5 (a fish finds a rock from its vortex street).

The camera starts close to the fish, whose 10 lateral-line sensors (5 on the
head third, 5 on the middle third) glow purple as the vortices of the Karman
street pass over them, then pulls back over the river. The fish steps across
the river; an inset builds up the amplitude of the shedding frequency at each
station. A red laser then marks where HYDAR places the rock. The pressure
field, the sensor amplitudes and the fix all come from a real hydar.karman
run (the scenario of example 08 and Figure 8).

Time: the street runs at SPEED times real time. The survey is shortened for
the video: in the simulation the fish holds each of the 19 stations for 60 s.
The water surface shows the street's pressure, exaggerated (the real surface
moves by about 2 mm). The water is drawn translucent over the gravel bed, so the
whole rock shows. Dye released at the rock's two flanks is carried by the
river and rolled up by the street's vortices, as in a flume with dye
injection; faint streaks across the river move with the current. The dye and
the streaks are advected by the same kinematic street (velocity U plus the
velocity the vortices induce) that gives the pressure.

Run:   python examples/make_fish_animation.py [out_dir]
Needs: ffmpeg on the PATH. Takes about 15 minutes.
"""

import os
import shutil
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from hydar.karman import (Fish, HemisphereRock, KarmanStreet, RiverTurbulence, TrapezoidalChannel,
                          lateral_survey, locate_rock)

FPS = 15
SPEED = 3.0                       # street: simulated seconds per video second
OUT = sys.argv[1] if len(sys.argv) > 1 else "docs/img"
W_PX, H_PX, DPI = 960, 540, 120

# ---------------------------------------------------------------- simulation
channel = TrapezoidalChannel()
rock = HemisphereRock.from_protrusion(channel, 0.0, 0.4)
street = KarmanStreet(channel, rock, phase=3.0)
L = 0.25
X_SNOUT = 6.0 - L / 3
STATIONS = np.round(np.arange(-1.8, 1.81, 0.2), 2)
print("recording the survey ...")
survey = lateral_survey(street, X_SNOUT, STATIONS, seed=1,
                        turbulence=RiverTurbulence(street.U, intensity=0.10, scale=channel.depth, seed=1))
fix = locate_rock(survey)
print(f"fix: x {fix.x:.2f} m, y {fix.y:.2f} m, D {fix.diameter:.2f} m (true 0, {rock.y}, {street.D:.2f})")
yf = np.linspace(-2, 2, 201)
M1, _ = fix.model.harmonic_profile(fix.distance, yf - fix.y)
M1 = np.sqrt(M1 ** 2 + fix.noise1 ** 2)
SENS = Fish(X_SNOUT, 0.0).sensor_offsets

# ---------------------------------------------------------------- timeline (video seconds)
T_ZOOM = (3.5, 7.0)
Y_INTRO = float(STATIONS[np.argmax(fix.A1)])     # the fish starts under a vortex row
T_SURVEY = (7.5, 17.5)
T_LASER = 18.0
T_END = 23.0
N_FRAMES = int(T_END * FPS)
W = channel.top_width

BG = "#02040c"
NAVY = np.array([0.02, 0.06, 0.16])
DEEP = np.array([0.01, 0.02, 0.07])
CYAN = np.array([0.25, 0.85, 1.00])
PURPLE = "#b84dff"
RED = "#ff2a3c"
P_REF = 15.0                      # Pa, pressure that saturates the glow


# ---------------------------------------------------------------- dye and flow streaks
DT_SIM = 0.1                       # s, advection step
WARM = 45.0                        # s of street before the video starts, so the dye fills the view
X_DROP = 17.0                      # particles leave the simulation here


def velocity(x, y, ts):
    """Velocity (u, v) [m/s] of the river plus the street at points (x, y), time ts."""
    xv, yv, gv, rv = street.vortices(ts, rock.x + X_DROP + 4)
    u = np.full(x.shape, street.U)
    v = np.zeros(x.shape)
    if xv.size:
        yi, gi = street._with_images(yv, gv)
        n = len(yi) // xv.size
        xi, ri = np.tile(xv, n), np.tile(rv, n)
        dx = x[:, None] - xi[None, :]
        dy = y[:, None] - yi[None, :]
        r2 = np.maximum(dx ** 2 + dy ** 2, 1e-12)
        fac = gi / (2 * np.pi * r2) * (-np.expm1(-r2 / ri ** 2))
        u = u - (fac * dy).sum(1)
        v = v + (fac * dx).sum(1)
    return u, v


def advect_particles():
    """Positions of dye (from the rock's flanks) and streak particles at every frame.

    Returns dict kind -> (x, y, age) arrays of shape (n_frames, n_particles), NaN where absent.
    """
    rng = np.random.default_rng(5)
    ts0 = survey.t[0, 0]
    steps_per_frame = int(round(SPEED / FPS / DT_SIM))
    n_warm = int(WARM / DT_SIM)
    n_steps = n_warm + N_FRAMES * steps_per_frame
    r_flank = 0.5 * rock.waterline_diameter(street.depth)
    out = {}
    for kind in ("dye", "streak"):
        xs, ys, born = [], [], []
        X = np.empty(0); Y = np.empty(0); B = np.empty(0); ID = np.empty(0, int)
        nid = 0
        rec = []
        for k in range(n_steps):
            ts = ts0 - WARM + k * DT_SIM
            if kind == "dye":             # two dye ports on each flank, a little behind the widest point
                nx = np.full(4, rock.x + 0.05) + rng.normal(0, 0.01, 4)
                ny = rock.y + np.array([1, 1, -1, -1]) * (r_flank + np.array([0.03, 0.10, 0.03, 0.10]))
                ny = ny + rng.normal(0, 0.008, 4)
            elif k % 4 == 0:              # streaks: random points across the river, upstream
                m = 5
                nx = rng.uniform(-3.5, -3.0, m)
                ny = rng.uniform(-0.5 * W + 0.05, 0.5 * W - 0.05, m)
            else:
                nx = ny = np.empty(0)
            X = np.r_[X, nx]; Y = np.r_[Y, ny]; B = np.r_[B, np.full(nx.size, ts)]
            ID = np.r_[ID, np.arange(nid, nid + nx.size)]; nid += nx.size
            # midpoint (RK2) step, plus a little turbulent diffusion
            u1, v1 = velocity(X, Y, ts)
            u2, v2 = velocity(X + 0.5 * DT_SIM * u1, Y + 0.5 * DT_SIM * v1, ts + 0.5 * DT_SIM)
            X = X + DT_SIM * u2 + rng.normal(0, 0.004, X.size)
            Y = np.clip(Y + DT_SIM * v2 + rng.normal(0, 0.004, Y.size), -0.5 * W + 0.02, 0.5 * W - 0.02)
            if kind == "streak":          # the rock is solid: streaks that hit it are dropped
                keep = np.hypot(X - rock.x, Y - rock.y) > r_flank
            else:
                keep = np.ones(X.size, bool)
            keep &= X < rock.x + X_DROP
            X, Y, B, ID = X[keep], Y[keep], B[keep], ID[keep]
            if k >= n_warm and (k - n_warm) % steps_per_frame == steps_per_frame - 1:
                rec.append((ID.copy(), X.copy(), Y.copy(), ts + DT_SIM - B))
        Xa = np.full((N_FRAMES, nid), np.nan); Ya = Xa.copy(); Aa = Xa.copy()
        for i, (ids, x, y, a) in enumerate(rec[:N_FRAMES]):
            Xa[i, ids], Ya[i, ids], Aa[i, ids] = x, y, a
        used = ~np.all(np.isnan(Xa), axis=0)
        out[kind] = (Xa[:, used], Ya[:, used], Aa[:, used])
        print(f"  {kind}: {used.sum()} particles")
    return out


def ease(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def lerp(a, b, s):
    return a + (b - a) * s


def fish_y(t):
    """Lateral position of the fish at video time t (steps between stations)."""
    if t < T_SURVEY[0]:          # hold under a vortex row, then move to the first station
        return lerp(Y_INTRO, STATIONS[0], ease((t - T_ZOOM[0]) / (T_SURVEY[0] - T_ZOOM[0])))
    u = (t - T_SURVEY[0]) / (T_SURVEY[1] - T_SURVEY[0]) * (len(STATIONS) - 1)
    k = int(np.clip(np.floor(u), 0, len(STATIONS) - 2))
    s = ease((u - k - 0.55) / 0.45)          # hold, then move to the next station
    return lerp(STATIONS[k], STATIONS[k + 1], s) if u < len(STATIONS) - 1 else STATIONS[-1]


def stations_done(t):
    if t < T_SURVEY[0]:
        return 0
    u = (t - T_SURVEY[0]) / (T_SURVEY[1] - T_SURVEY[0]) * (len(STATIONS) - 1)
    return int(np.clip(np.floor(u + 0.55), 0, len(STATIONS)))


def camera(t):
    """Window centre, half-widths (x, y), elevation, azimuth at video time t."""
    s = ease((t - T_ZOOM[0]) / (T_ZOOM[1] - T_ZOOM[0]))
    near = np.array([X_SNOUT + 0.1, fish_y(t)])
    c = lerp(near, np.array([5.9, 0.0]), s)
    hx = lerp(1.1, 7.9, s)
    hy = lerp(1.1, 0.5 * W + 0.6, s)
    elev = lerp(28.0, 40.0, s)
    azim = lerp(222.0, 240.0, s) + 5.0 * t / T_END
    return c, hx, hy, elev, azim


def field(t, c, hx, hy, n=200):
    xs = np.linspace(c[0] - hx, c[0] + hx, n)
    ys = np.linspace(max(c[1] - hy, -W / 2), min(c[1] + hy, W / 2), max(int(n * hy / hx), 50))
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    ts = survey.t[0, 0] + SPEED * t
    P = street.pressure(X.ravel(), Y.ravel(), [ts]).reshape(X.shape)
    P[np.hypot(X - rock.x, Y - rock.y) < 0.5 * rock.waterline_diameter(street.depth)] = 0.0
    return X, Y, P


def surface_colors(X, Y, P, hx, hy, ts):
    gx, gy = np.gradient(P, X[1, 0] - X[0, 0], Y[0, 1] - Y[0, 0])
    light = np.clip(0.5 + 0.02 * (-0.6 * gx - 0.8 * gy), 0, 1)
    base = lerp(DEEP, NAVY, light[..., None]) * 1.25
    glow = np.clip(-P / P_REF, 0, 1)[..., None] ** 1.1        # vortex cores: low pressure glows
    high = np.clip(P / P_REF, 0, 1)[..., None]
    xa = X - street.U * ts                                       # ripples drift with the current
    shimmer = 0.016 * (np.sin(17 * xa + 6 * Y) * np.sin(9 * Y - 13 * xa))[..., None]
    col = base * (1 - 0.4 * high) + glow * (CYAN - base) + shimmer
    r = np.maximum(np.abs(X - X.mean()) / hx, np.abs(Y - Y.mean()) / hy)
    fade = np.clip((1.0 - r) / 0.10, 0, 1)[..., None]
    col = lerp(np.array([0.008, 0.016, 0.047]), col, fade)
    alpha = np.clip(0.34 + 0.55 * glow[..., 0] + 0.15 * high[..., 0], 0, 0.92) * fade[..., 0] \
        + (1 - fade[..., 0])                                     # translucent water, opaque glow
    return np.concatenate([np.clip(col, 0, 1), alpha[..., None]], axis=-1)


def glow_line(ax, xs, ys, zs, color, width, alpha=1.0, layers=((6, .10), (3, .22), (1, 1.0)), zorder=6):
    for mult, a in layers:
        ax.plot(xs, ys, zs, color=color, lw=width * mult, alpha=a * alpha,
                solid_capstyle="round", zorder=zorder)


def glow_points(ax, x, y, z, color, size, alpha=1.0, zorder=7):
    for mult, a in ((9, .06), (4, .14), (1.8, .35), (1, 1.0)):
        ax.scatter(x, y, z, s=size * mult, color=color, alpha=min(a * alpha, 1),
                   edgecolors="none", depthshade=False, zorder=zorder)


def draw_rock(ax, zscale, above):
    """The part of the hemisphere above (above=True) or below the water surface."""
    R = rock.radius
    h = street.depth
    v0 = np.arcsin(h / R)                      # elevation angle of the waterline
    v = np.linspace(v0, np.pi / 2, 10) if above else np.linspace(0, v0, 8)
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 40), v)
    zc = R * np.sin(v) - h                     # bed at -h, water surface at 0
    X = rock.x + R * np.cos(v) * np.cos(u)
    Y = rock.y + R * np.cos(v) * np.sin(u)
    shade = 0.35 + 0.45 * np.clip(np.cos(u - 2.3) * np.cos(v) * 0.5 + np.sin(v) * 0.6, 0, 1)
    tint = np.array([0.55, 0.57, 0.62]) if above else np.array([0.30, 0.42, 0.52])
    cols = np.concatenate([shade[..., None] * tint, np.ones(shade.shape + (1,))], -1)
    ax.plot_surface(X, Y, zc * zscale, facecolors=cols, shade=False, linewidth=0, antialiased=False,
                    zorder=8 if above else 0.5)


def draw_bed(ax, c, hx, hy, zscale):
    """Gravel bed and the submerged bank slopes, seen through the water."""
    h = street.depth
    xs = np.linspace(c[0] - hx, c[0] + hx, 90)
    ys = np.linspace(max(c[1] - hy, -0.5 * W), min(c[1] + hy, 0.5 * W), 40)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    Z = -np.asarray(channel.local_depth(Y.ravel())).reshape(Y.shape) * zscale
    g = np.random.default_rng(2).random(X.shape)
    lum = 0.55 + 0.45 * g
    col = np.stack([0.09 * lum, 0.12 * lum, 0.15 * lum, np.ones_like(lum)], -1)
    r = np.maximum(np.abs(X - c[0]) / hx, np.abs(Y - c[1]) / hy)
    fade = np.clip((1.0 - r) / 0.30, 0, 1)[..., None]          # bed fades into the background
    col[..., :3] = lerp(np.array([0.008, 0.016, 0.047]), col[..., :3], fade)
    ax.plot_surface(X, Y, Z, facecolors=col, shade=False, linewidth=0, antialiased=False, zorder=0)


def draw_banks(ax, c, hx, zscale):
    x0, x1 = c[0] - hx, c[0] + hx
    for sgn in (1, -1):
        y0, y1 = sgn * W / 2, sgn * (W / 2 + 0.7)
        quad = [[(x0, y0, 0), (x1, y0, 0), (x1, y1, 0.35 * zscale), (x0, y1, 0.35 * zscale)]]
        ax.add_collection3d(Poly3DCollection(quad, facecolor="#1b2233", edgecolor="#2c3a55", lw=0.4,
                                             zorder=2))


def draw_fish(ax, y, z, amp, zoom):
    """Fish outline (head upstream, -x) with its 10 sensors; amp = sensor glow 0..1."""
    s = np.linspace(0, 1, 60)
    half = 0.11 * L * np.sqrt(s) * (1 - s) ** 0.8 * 2.2
    xs = X_SNOUT + s * L
    outline_x = np.r_[xs, xs[::-1], X_SNOUT]
    outline_y = y + np.r_[half, -half[::-1], 0]
    tail_x = [X_SNOUT + L, X_SNOUT + L + 0.03, X_SNOUT + L + 0.03, X_SNOUT + L]
    tail_y = [y, y + 0.022, y - 0.022, y]
    lw = 0.9 if zoom else 0.6
    glow_line(ax, outline_x, outline_y, np.full(outline_x.size, z), "#cfe9ff", lw, alpha=0.9, zorder=9)
    glow_line(ax, tail_x, tail_y, np.full(4, z), "#cfe9ff", lw, alpha=0.9, zorder=9)
    sx = X_SNOUT + SENS
    size = 10 if zoom else 2.5
    ax.scatter(sx, np.full(10, y), np.full(10, z), s=size, color="#5b3b7a", depthshade=False, zorder=10)
    for xk, ak in zip(sx, amp):
        glow_points(ax, [xk], [y], [z], PURPLE, size * (0.6 + 1.6 * ak), alpha=0.4 + 0.6 * ak, zorder=11)


def draw_particles(ax, i, c, hx, hy, zs, zoom):
    """Dye from the rock (pale cyan) and flow streaks (short white dashes) at frame i."""
    inside = lambda x, y: (np.abs(x - c[0]) < 0.97 * hx) & (np.abs(y - c[1]) < 0.97 * hy)
    # streaks: a dash from where the particle was 3 frames ago
    xs, ys, _ = PARTICLES["streak"]
    j = max(i - (1 if zoom else 2), 0)
    ok = ~np.isnan(xs[i]) & ~np.isnan(xs[j]) & inside(xs[i], ys[i]) & inside(xs[j], ys[j])
    for x0, y0, x1, y1 in zip(xs[j, ok], ys[j, ok], xs[i, ok], ys[i, ok]):
        ax.plot([x0, x1], [y0, y1], [0.6 * zs * 0.2] * 2, color="#dff4ff", lw=0.5, alpha=0.2, zorder=3)
    # dye: fades as it spreads downstream
    xd, yd, ad = PARTICLES["dye"]
    ok = ~np.isnan(xd[i]) & inside(xd[i], yd[i])
    if ok.any():
        fade = np.clip(1.15 - ad[i, ok] / 40.0, 0.15, 1.0)
        size = 1.6 if not zoom else 5.0
        cols = np.zeros((ok.sum(), 4))
        cols[:, :3] = [0.80, 0.97, 1.0]
        ax.scatter(xd[i, ok], yd[i, ok], np.full(ok.sum(), 0.12 * zs), s=size * 5,
                   color=np.c_[cols[:, :3], 0.05 * fade], edgecolors="none", depthshade=False, zorder=3)
        ax.scatter(xd[i, ok], yd[i, ok], np.full(ok.sum(), 0.12 * zs), s=size,
                   color=np.c_[cols[:, :3], 0.55 * fade], edgecolors="none", depthshade=False, zorder=4)


def render(i, path):
    t = i / FPS
    c, hx, hy, elev, azim = camera(t)
    zoom = hx < 2.5
    zscale = lerp(0.6, 1.0, ease((t - T_ZOOM[0]) / (T_ZOOM[1] - T_ZOOM[0])))
    zs = 0.035 * hy                         # height of the pressure "surface" at P_REF
    fig = plt.figure(figsize=(W_PX / DPI, H_PX / DPI), dpi=DPI, facecolor=BG)
    ax = fig.add_axes([-0.12, -0.2, 1.24, 1.4], projection="3d", computed_zorder=False)
    ax.set_facecolor(BG)
    ax.set_axis_off()
    ax.set_xlim(c[0] - hx, c[0] + hx)
    ax.set_ylim(c[1] - hy, c[1] + hy)
    ax.set_zlim(-0.35 * hy, 0.45 * hy)
    ax.set_box_aspect((hx / hy, 1, 0.7), zoom=1.45 if hx / hy < 1.3 else 1.3)
    ax.view_init(elev=elev, azim=azim)

    ts = survey.t[0, 0] + SPEED * t
    rock_in_view = abs(rock.x - c[0]) < hx + 0.6
    draw_bed(ax, c, hx, hy, zscale)
    if rock_in_view:
        draw_rock(ax, zscale, above=False)
    X, Y, P = field(t, c, hx, hy)
    Z = zs * np.tanh(P / P_REF)
    ax.plot_surface(X, Y, Z, facecolors=surface_colors(X, Y, P, hx, hy, ts), shade=False,
                    rstride=1, cstride=1, linewidth=0, antialiased=False, zorder=1)
    draw_particles(ax, i, c, hx, hy, zs, zoom)
    if not zoom:
        draw_banks(ax, c, hx, zscale * 1.0)
    if rock_in_view:
        draw_rock(ax, zscale, above=True)

    # fish and its sensors, glowing with the pressure they feel
    y = fish_y(t)
    pos = Fish(X_SNOUT, y).sensor_positions
    p_now = street.pressure(pos[:, 0], pos[:, 1], [ts])[:, 0]
    amp = np.clip(-p_now / P_REF, 0, 1)
    draw_fish(ax, y, -0.02 * hy, amp, zoom)

    # stations already surveyed, coloured by the measured amplitude of f
    nd = stations_done(t)
    if nd and not zoom:
        a1 = fix.A1[:nd] / fix.A1.max()
        for yk, ak in zip(STATIONS[:nd], a1):
            glow_points(ax, [survey.x_center], [yk], [0.0], PURPLE, 2 + 10 * ak, alpha=0.3 + 0.7 * ak, zorder=5)

    # HYDAR: laser up the street's centreline to the rock
    if t > T_LASER:
        s = ease((t - T_LASER) / 1.4)
        start = np.array([survey.x_center, fix.y])
        q = start + s * (np.array([fix.x, fix.y]) - start)
        zl = 0.02 * hy
        glow_line(ax, [start[0], q[0]], [start[1], q[1]], [zl, zl], RED, 1.4, zorder=12)
        if s >= 1:
            glow_points(ax, [fix.x], [fix.y], [0.3], RED, 30, zorder=13)

    fig.text(0.035, 0.92, "HYDAR", color="white", fontsize=22, weight="bold", family="DejaVu Sans")
    fig.text(0.035, 0.875, "Hydrodynamic Detection and Ranging", color="#9fb3d9", fontsize=10)
    if zoom:
        fig.text(0.035, 0.82, "25 cm fish, 10 lateral-line pressure sensors\n(5 head third, 5 middle third, none on the tail)",
                 color="#d9b3ff", fontsize=9, va="top")
    else:
        fig.text(0.035, 0.82, "Kármán vortex street behind a rock", color="#7fdfff", fontsize=9, va="top")

    # inset: amplitude of the shedding frequency across the river
    if t > T_SURVEY[0] - 0.3:
        ia = fig.add_axes([0.70, 0.08, 0.26, 0.24], facecolor=(0.02, 0.04, 0.1, 0.85))
        for sp in ia.spines.values():
            sp.set_color("#3a4a6a")
        ia.tick_params(colors="#9fb3d9", labelsize=7, length=2)
        ia.set_xlim(-2, 2)
        ia.set_ylim(0, 1.25 * fix.A1.max())
        ia.set_xlabel("y (m)", color="#9fb3d9", fontsize=7, labelpad=1)
        ia.set_ylabel("Amplitude of f (Pa)", color="#9fb3d9", fontsize=7, labelpad=1)
        if nd:
            ia.plot(STATIONS[:nd], fix.A1[:nd], "o", color=PURPLE, ms=3, mew=0)
        if t > T_LASER:
            s = ease((t - T_LASER) / 1.0)
            ia.plot(yf, M1, color="#3fd0ff", lw=1.0, alpha=s)
            ia.axvline(fix.y, color=RED, lw=1.0, alpha=s)
    if t > T_LASER + 1.6:
        a = ease((t - T_LASER - 1.6) / 0.8)
        fig.text(0.035, 0.75, f"Rock found {fix.distance:.2f} m upstream, y = {fix.y:.2f} m, D = {fix.diameter:.2f} m\n"
                 f"(true: {survey.x_center - rock.x:.2f} m, {rock.y:.2f} m, {street.D:.2f} m)",
                 color="#ffb3ba", fontsize=9.5, alpha=a, va="top")
    fig.savefig(path, facecolor=BG)
    plt.close(fig)


PARTICLES = None


def setup():
    global PARTICLES
    print("advecting dye and flow streaks ...")
    PARTICLES = advect_particles()


def main():
    setup()
    os.makedirs(OUT, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="hydar_fish_frames_")
    print(f"rendering {N_FRAMES} frames ...")
    for i in range(N_FRAMES):
        render(i, os.path.join(tmp, f"f{i:04d}.png"))
        if i % 20 == 0:
            print(f"  {i}/{N_FRAMES}")
    pat = os.path.join(tmp, "f%04d.png")
    mp4 = os.path.join(OUT, "hydar_fish.mp4")
    gif = os.path.join(OUT, "hydar_fish.gif")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pat,
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "22", mp4], check=True)
    pal = os.path.join(tmp, "pal.png")
    scale = "scale=720:-1:flags=lanczos"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pat,
                    "-vf", f"{scale},palettegen=max_colors=160:stats_mode=diff", pal], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pat, "-i", pal,
                    "-lavfi", f"{scale}[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle",
                    "-loop", "0", gif], check=True)
    shutil.rmtree(tmp)
    for f in (gif, mp4):
        print(f"wrote {f} ({os.path.getsize(f) / 1e6:.1f} MB)")


if __name__ == "__main__":
    if os.environ.get("HYDAR_PREVIEW"):
        setup()
        os.makedirs(OUT, exist_ok=True)
        for t in map(float, os.environ["HYDAR_PREVIEW"].split(",")):
            render(int(round(t * FPS)), os.path.join(OUT, f"preview_{t:05.2f}.png"))
    else:
        main()
