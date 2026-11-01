"""Render the animation for case study 2 (boat tracking).

A Series 60 boat crosses the domain while the camera chases it, then the
camera pulls back. Each 3x3 array glows purple as the blue wake washes over
it. Red lasers then mark where HYDAR placed the boat when each array heard
it, and the fitted HYDAR track lights up across the domain. Waves, arrival
times and the track all come from a real hydar simulation.

Run:   python examples/make_boat_animation.py [out_dir]
Needs: ffmpeg on the PATH. Takes several minutes.
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
from scipy.ndimage import gaussian_filter, map_coordinates

from hydar import Boat, Grid, SensorArray, Solver, track_boat
from hydar.boat import Series60Hull

FPS = 15
OUT = sys.argv[1] if len(sys.argv) > 1 else "docs/img"
W_PX, H_PX, DPI = 960, 540, 120

# ---------------------------------------------------------------- simulation
grid = Grid(lx=40.0, ly=30.0, dx=0.05, n_pml=60)
sim = Solver(grid, depth=1.0)
arrays = [SensorArray("A", 10.0, 4.0, z=-0.3), SensorArray("B", 28.0, 5.0, z=-0.5),
          SensorArray("C", 12.0, 25.0, z=-0.5), SensorArray("D", 31.0, 25.0, z=-0.8)]
for a in arrays:
    sim.add(a)
hull = Series60Hull(length=4.5)
boat = Boat(hull=hull, speed=5.0, x=-2.5, y=12.0, heading_deg=10.0)
sim.add(boat)

T_SIM = 10.0
snaps, next_t = [], 0.0
si, sj = grid.interior
print("simulating ...")
while sim.t < T_SIM:
    sim.step()
    if sim.t >= next_t:
        snaps.append(sim.eta[si, sj][::2, ::2].astype(np.float32))  # 0.1 m grid
        next_t += 1.0 / FPS
snaps = np.array(snaps)
SNAP_DX = 0.1
tr = track_boat(arrays, c=sim.c)
print("HYDAR track:", tr)
detect = [float(np.nanmin(p)) for p in tr.arrivals]   # first wake arrival per array
fix_pts = [tr.position(td) for td in detect]          # HYDAR boat position at that moment

# ---------------------------------------------------------------- timeline
T_ZOOM = (3.2, 6.2)
T_LASER = max(detect) + 0.5
T_END = T_LASER + 4.5
N_FRAMES = int(T_END * FPS)

BG = "#02040c"
NAVY = np.array([0.02, 0.06, 0.16])
DEEP = np.array([0.01, 0.02, 0.07])
CYAN = np.array([0.25, 0.85, 1.00])
PURPLE = "#b84dff"
RED = "#ff2a3c"
ETA_REF = 0.09


def ease(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def lerp(a, b, s):
    return a + (b - a) * s


def camera(t):
    """Window centre, half-width, elevation, azimuth at time t."""
    s = ease((t - T_ZOOM[0]) / (T_ZOOM[1] - T_ZOOM[0]))
    chase = boat.position(min(t, T_ZOOM[1])) - 1.0 * boat.e
    chase = np.clip(chase, [5.0, 5.0], [35.0, 25.0])
    c = lerp(chase, np.array([20.0, 15.0]), s)
    hw = lerp(5.5, 21.5, s)
    elev = lerp(20.0, 42.0, s)
    azim = lerp(boat.heading_deg + 180 + 38, 300.0, s) + 4.0 * t / T_END
    return c, hw, elev, azim


def field(t, c, hw, n=150):
    xs = np.linspace(c[0] - hw, c[0] + hw, n)
    ys = np.linspace(c[1] - hw, c[1] + hw, n)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    k = min(t * FPS, len(snaps) - 1)
    k0, w = int(k), k - int(k)
    k1 = min(k0 + 1, len(snaps) - 1)
    coords = np.array([X.ravel() / SNAP_DX - 0.5, Y.ravel() / SNAP_DX - 0.5])
    e0 = map_coordinates(snaps[k0], coords, order=1, mode="constant")
    e1 = map_coordinates(snaps[k1], coords, order=1, mode="constant")
    eta = ((1 - w) * e0 + w * e1).reshape(X.shape)
    return X, Y, gaussian_filter(eta, 1.4)


def surface_colors(X, Y, eta, hw):
    gx, gy = np.gradient(eta, X[1, 0] - X[0, 0])
    light = np.clip(0.5 + 6 * (-0.6 * gx - 0.8 * gy), 0, 1)
    base = lerp(DEEP, NAVY, light[..., None])
    glow = np.clip(eta / (0.7 * ETA_REF), 0, 1)[..., None] ** 0.9
    trough = np.clip(-eta / ETA_REF, 0, 1)[..., None]
    shimmer = 0.012 * (np.sin(1.7 * X + 0.6 * Y) * np.sin(0.9 * Y - 1.3 * X))[..., None]
    col = base * (1 - 0.6 * trough) + glow * (CYAN - base) + shimmer
    r = np.maximum(np.abs(X - X.mean()), np.abs(Y - Y.mean())) / hw
    fade = np.clip((1.0 - r) / 0.12, 0, 1)[..., None]
    col = lerp(np.array([0.008, 0.016, 0.047]), col, fade)
    return np.concatenate([np.clip(col, 0, 1), np.ones(eta.shape + (1,))], axis=-1)


def glow_line(ax, xs, ys, zs, color, width, alpha=1.0, layers=((6, .10), (3, .22), (1, 1.0))):
    for mult, a in layers:
        ax.plot(xs, ys, zs, color=color, lw=width * mult, alpha=a * alpha,
                solid_capstyle="round", zorder=6)


def glow_points(ax, x, y, z, color, size, alpha=1.0, zorder=7):
    for mult, a in ((9, .06), (4, .14), (1.8, .35), (1, 1.0)):
        ax.scatter(x, y, z, s=size * mult, color=color, alpha=min(a * alpha, 1),
                   edgecolors="none", depthshade=False, zorder=zorder)


# Series 60 waterline outline in hull coordinates (bow at +L/2)
XI = np.linspace(-1, 1, 41)
HB = 0.5 * hull.beam * (1 - np.abs(XI) ** hull.m)
OUTLINE = np.r_[np.column_stack([XI * hull.length / 2, HB]),
                np.column_stack([XI[::-1] * hull.length / 2, -HB[::-1]])]


def draw_boat(ax, t):
    p = boat.position(t)
    if not (-3 < p[0] < 43):
        return
    e, n = boat.e, boat.n
    world = p + OUTLINE[:, :1] * e + OUTLINE[:, 1:] * n
    fb = 0.32                                  # freeboard (m)
    deck = np.column_stack([world, np.full(len(world), fb)])
    water = np.column_stack([world, np.zeros(len(world))])
    sides = [[water[i], water[i + 1], deck[i + 1], deck[i]] for i in range(len(world) - 1)]
    ax.add_collection3d(Poly3DCollection(sides, facecolor="#9fb0c8", edgecolor="none", zorder=8))
    ax.add_collection3d(Poly3DCollection([deck], facecolor="#eef3fb", edgecolor="#ffffff",
                                         lw=0.5, zorder=9))
    # cabin
    cx = p - 0.35 * e
    hl, hb, h = 0.55, 0.2, 0.32
    corners = [cx + sx * hl * e + sy * hb * n for sx, sy in ((1, 1), (1, -1), (-1, -1), (-1, 1))]
    top = [np.r_[c_, fb + h] for c_ in corners]
    bot = [np.r_[c_, fb] for c_ in corners]
    faces = [top] + [[bot[i], bot[(i + 1) % 4], top[(i + 1) % 4], top[i]] for i in range(4)]
    ax.add_collection3d(Poly3DCollection(faces, facecolor="#c9d4e6", edgecolor="#ffffff",
                                         lw=0.4, zorder=10))
    # navigation lights: red port, green starboard, white stern
    bow = p + 0.42 * hull.length * e
    glow_points(ax, [bow[0] + 0.2 * n[0]], [bow[1] + 0.2 * n[1]], [fb + 0.05], "#ff3b4e", 5, zorder=11)
    glow_points(ax, [bow[0] - 0.2 * n[0]], [bow[1] - 0.2 * n[1]], [fb + 0.05], "#39ff88", 5, zorder=11)
    stern = p - 0.5 * hull.length * e
    glow_points(ax, [stern[0]], [stern[1]], [fb + 0.1], "#ffffff", 4, zorder=11)
    # bow spray
    k = np.arange(24)
    ph = (t * 3.0 + k / 24) % 1.0
    side = np.where(k % 2, 1, -1)
    sp = bow[None, :] - ph[:, None] * 1.2 * e + side[:, None] * (0.25 + ph[:, None] * 0.9) * n
    ax.scatter(sp[:, 0], sp[:, 1], 0.05 + 0.35 * ph * (1 - ph), s=6, color="#e6faff",
               alpha=0.7, edgecolors="none", depthshade=False, zorder=7)


def render(i, path):
    t = i / FPS
    c, hw, elev, azim = camera(t)
    zs = 0.014 * hw
    fig = plt.figure(figsize=(W_PX / DPI, H_PX / DPI), dpi=DPI, facecolor=BG)
    ax = fig.add_axes([-0.12, -0.2, 1.24, 1.4], projection="3d", computed_zorder=False)
    ax.set_facecolor(BG)
    ax.set_axis_off()
    ax.set_xlim(c[0] - hw, c[0] + hw)
    ax.set_ylim(c[1] - hw, c[1] + hw)
    ax.set_zlim(-0.25 * hw, 0.45 * hw)
    ax.set_box_aspect((1, 1, 0.7), zoom=1.45)
    ax.view_init(elev=elev, azim=azim)

    X, Y, eta = field(t, c, hw)
    Z = zs * np.tanh(eta / ETA_REF)
    ax.plot_surface(X, Y, Z, facecolors=surface_colors(X, Y, eta, hw), shade=False,
                    rstride=1, cstride=1, linewidth=0, antialiased=False, zorder=1)
    inwin = lambda px, py: abs(px - c[0]) < hw * 0.97 and abs(py - c[1]) < hw * 0.97

    # glowing blue Mach-cone front trailing the bow
    bow = boat.position(t) + 0.5 * hull.length * boat.e
    mu = np.arcsin(sim.c / boat.speed)
    if t > 0.8:
        for sgn in (1, -1):
            d = -np.cos(mu) * boat.e + sgn * np.sin(mu) * boat.n
            L = np.linspace(0, 40, 300)
            fx, fy = bow[0] + L * d[0], bow[1] + L * d[1]
            keep = (np.abs(fx - c[0]) < hw * 0.95) & (np.abs(fy - c[1]) < hw * 0.95) & \
                   (fx > 0.3) & (fy > 0.3) & (fy < 29.7)
            glow_line(ax, np.where(keep, fx, np.nan), np.where(keep, fy, np.nan),
                      np.full_like(fx, zs * 0.3), "#3fd0ff", 1.1, alpha=0.7)

    # arrays glow purple when the wake reaches them
    for a, td in zip(arrays, detect):
        if not inwin(a.x, a.y):
            continue
        dt = t - td
        lit = 0.0 if dt < -0.1 else (1.0 + 2.0 * np.exp(-max(dt, 0) / 0.35))
        px, py = a.positions[:, 0], a.positions[:, 1]
        ax.scatter(px, py, np.full(9, zs * 0.3), s=6, color="#5b3b7a", depthshade=False, zorder=5)
        if lit > 0:
            glow_points(ax, px, py, np.full(9, zs * 0.3), PURPLE, 7 * lit, alpha=min(lit / 1.5, 1))
        if hw > 9:
            ax.text(a.x + 1.4, a.y - 1.4, zs * 0.4, a.name, color="#d9b3ff" if lit else "#6a5a80",
                    fontsize=9, zorder=9)

    # HYDAR: lasers to where each array placed the boat, then the fitted track
    if t > T_LASER:
        zl = zs * 0.35
        for n_, (a, fp) in enumerate(zip(arrays, fix_pts)):
            s = ease((t - T_LASER - 0.3 * n_) / 0.8)
            if s <= 0:
                continue
            q = a.center + s * (fp - a.center)
            glow_line(ax, [a.x, q[0]], [a.y, q[1]], [zl, zl], RED, 1.4)
            if s >= 1:
                glow_points(ax, [fp[0]], [fp[1]], [zl], RED, 18)
        t_track = T_LASER + 0.3 * 3 + 0.8
        if t > t_track:
            s = ease((t - t_track) / 1.2)
            t0 = (0.0 - tr.p0[0]) / (tr.speed * tr.direction[0])
            t1 = (40.0 - tr.p0[0]) / (tr.speed * tr.direction[0])
            tt = np.linspace(t0, lerp(t0, t1, s), 120)
            P = tr.position(tt)
            glow_line(ax, P[:, 0], P[:, 1], np.full(len(tt), zl), RED, 1.8)

    draw_boat(ax, t)

    fig.text(0.035, 0.92, "HYDAR", color="white", fontsize=22, weight="bold", family="DejaVu Sans")
    fig.text(0.035, 0.875, "Hydrodynamic Detection and Ranging", color="#9fb3d9", fontsize=10)
    fig.savefig(path, facecolor=BG)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="hydar_boat_frames_")
    print(f"rendering {N_FRAMES} frames ...")
    for i in range(N_FRAMES):
        render(i, os.path.join(tmp, f"f{i:04d}.png"))
        if i % 20 == 0:
            print(f"  {i}/{N_FRAMES}")
    pat = os.path.join(tmp, "f%04d.png")
    mp4 = os.path.join(OUT, "hydar_boat.mp4")
    gif = os.path.join(OUT, "hydar_boat.gif")
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
        os.makedirs(OUT, exist_ok=True)
        for t in map(float, os.environ["HYDAR_PREVIEW"].split(",")):
            render(int(t * FPS), os.path.join(OUT, f"preview_{t:05.2f}.png"))
    else:
        main()
