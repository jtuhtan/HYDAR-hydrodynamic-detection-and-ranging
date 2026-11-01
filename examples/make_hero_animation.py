"""Render the HYDAR hero animation shown at the top of the README.

A drone flies in, hovers and drops a sphere. The splash makes waves (taken
from a real hydar simulation), the camera pulls back, each 3x3 array glows
purple as the blue wave front washes over it, and red lasers trace the HYDAR
time-of-flight fix back to the drop point.

Run:   python examples/make_hero_animation.py [out_dir]
Needs: ffmpeg on the PATH for the GIF/MP4 encoding. Takes a few minutes.
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

from hydar import DropObject, Grid, SensorArray, Solver, locate

FPS = 15
OUT = sys.argv[1] if len(sys.argv) > 1 else "docs/img"
W_PX, H_PX, DPI = 960, 540, 120

# ---------------------------------------------------------------- simulation
DROP = np.array([11.0, 17.0])
grid = Grid(lx=30.0, ly=30.0, dx=0.05, n_pml=40)
sim = Solver(grid, depth=1.0)
arrays = [SensorArray("A", 6.0, 6.0, z=-0.2), SensorArray("B", 24.0, 6.0, z=-0.5),
          SensorArray("C", 24.0, 24.0, z=-0.8), SensorArray("D", 6.0, 24.0, z=-0.5)]
for a in arrays:
    sim.add(a)
ball = DropObject(diameter=0.5, density=2500.0, x=DROP[0], y=DROP[1])
sim.add(ball)

T_SIM = 9.5
snaps, next_t = [], 0.0
si, sj = grid.interior
print("simulating ...")
while sim.t < T_SIM:
    sim.step()
    if sim.t >= next_t:
        snaps.append(sim.eta[si, sj][::2, ::2].astype(np.float32))  # 0.1 m grid
        next_t += 1.0 / FPS
snaps = np.array(snaps)
fix = locate(arrays, c=sim.c, t0=0.0)
est = fix["position"]
arrivals = [np.nanmin(f.arrivals) for f in fix["fixes"]]
err_cm = 100 * np.linalg.norm(est - DROP)
print(f"HYDAR fix ({est[0]:.2f}, {est[1]:.2f}), error {err_cm:.1f} cm")
SNAP_DX = 0.1

# ---------------------------------------------------------------- timeline
T_IMPACT = 3.2            # film time when the sphere hits the water (sim t = 0)
T_RELEASE = 2.6
T_ZOOM = (T_IMPACT + 0.7, T_IMPACT + 3.6)
T_LASER = T_IMPACT + max(arrivals) + 0.6
T_END = T_LASER + 4.0
N_FRAMES = int(T_END * FPS)

# colours
BG = "#02040c"
NAVY = np.array([0.02, 0.06, 0.16])
DEEP = np.array([0.01, 0.02, 0.07])
CYAN = np.array([0.25, 0.85, 1.00])
PURPLE = "#b84dff"
RED = "#ff2a3c"


def ease(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def lerp(a, b, s):
    return a + (b - a) * s


def camera(t):
    """Window centre, half-width, elevation, azimuth at film time t."""
    s = ease((t - T_ZOOM[0]) / (T_ZOOM[1] - T_ZOOM[0]))
    c = lerp(DROP, np.array([15.0, 15.0]), s)
    hw = lerp(4.0, 16.5, s)
    elev = lerp(24.0, 42.0, s)
    azim = -62.0 + 6.0 * t / T_END * 2 + 10 * s
    return c, hw, elev, azim


def field(t, c, hw, n=150):
    """Surface elevation on an n x n grid covering the camera window."""
    ts = t - T_IMPACT
    xs = np.linspace(c[0] - hw, c[0] + hw, n)
    ys = np.linspace(c[1] - hw, c[1] + hw, n)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    if ts <= 0:
        return X, Y, np.zeros_like(X)
    k = min(ts * FPS, len(snaps) - 1)
    k0, w = int(k), k - int(k)
    k1 = min(k0 + 1, len(snaps) - 1)
    coords = np.array([X.ravel() / SNAP_DX - 0.5, Y.ravel() / SNAP_DX - 0.5])
    e0 = map_coordinates(snaps[k0], coords, order=1, mode="constant")
    e1 = map_coordinates(snaps[k1], coords, order=1, mode="constant")
    eta = ((1 - w) * e0 + w * e1).reshape(X.shape)
    return X, Y, gaussian_filter(eta, 0.8)


def surface_colors(X, Y, eta, hw):
    """Dark water with fake lighting and a blue glow on wave crests."""
    gx, gy = np.gradient(eta, X[1, 0] - X[0, 0])
    light = np.clip(0.5 + 40 * (-0.6 * gx - 0.8 * gy), 0, 1)
    base = lerp(DEEP, NAVY, light[..., None])
    glow = np.clip(eta / 0.006, 0, 1)[..., None] ** 0.8
    trough = np.clip(-eta / 0.008, 0, 1)[..., None]
    # faint shimmer so calm water does not look flat (colour only, no height)
    shimmer = 0.012 * (np.sin(1.7 * X + 0.6 * Y) * np.sin(0.9 * Y - 1.3 * X))[..., None]
    col = base * (1 - 0.6 * trough) + glow * (CYAN - base) + shimmer
    # fade the edge of the window into the background for a soft vignette
    r = np.maximum(np.abs(X - X.mean()), np.abs(Y - Y.mean())) / hw
    fade = np.clip((1.0 - r) / 0.12, 0, 1)[..., None]
    col = lerp(np.array([0.008, 0.016, 0.047]), col, fade)
    return np.concatenate([np.clip(col, 0, 1), np.ones(eta.shape + (1,))], axis=-1)


def glow_line(ax, xs, ys, zs, color, width, alpha=1.0, layers=((6, .10), (3, .22), (1, 1.0))):
    for mult, a in layers:
        ax.plot(xs, ys, zs, color=color, lw=width * mult, alpha=a * alpha,
                solid_capstyle="round", zorder=6)


def glow_points(ax, x, y, z, color, size, alpha=1.0):
    for mult, a in ((9, .06), (4, .14), (1.8, .35), (1, 1.0)):
        ax.scatter(x, y, z, s=size * mult, color=color, alpha=min(a * alpha, 1),
                   edgecolors="none", depthshade=False, zorder=7)


def draw_drone(ax, p, scale, t, carrying):
    """Quadcopter at position p (x, y, z)."""
    x, y, z = p
    a = 0.55 * scale
    body = np.array([[x - .18 * scale, y - .18 * scale, z], [x + .18 * scale, y - .18 * scale, z],
                     [x + .18 * scale, y + .18 * scale, z], [x - .18 * scale, y + .18 * scale, z]])
    ax.add_collection3d(Poly3DCollection([body], facecolor="#cfd6e6", edgecolor="#ffffff",
                                         lw=0.6, zorder=8))
    for k, (dx, dy) in enumerate(((1, 1), (-1, 1), (-1, -1), (1, -1))):
        rx, ry = x + dx * a, y + dy * a
        ax.plot([x, rx], [y, ry], [z, z], color="#9aa6bf", lw=2.2, zorder=8)
        th = np.linspace(0, 2 * np.pi, 24)
        rr = 0.32 * scale
        ax.plot(rx + rr * np.cos(th), ry + rr * np.sin(th), np.full_like(th, z + 0.03 * scale),
                color="#e8f4ff", lw=1.0, alpha=0.55, zorder=8)
        blade = th[0] + 9 * t + k
        ax.plot([rx - rr * np.cos(blade), rx + rr * np.cos(blade)],
                [ry - rr * np.sin(blade), ry + rr * np.sin(blade)],
                [z + 0.03 * scale] * 2, color="#ffffff", lw=1.4, alpha=0.8, zorder=8)
    led = "#39ff88" if int(t * 4) % 2 else "#ff3b4e"
    glow_points(ax, [x + a], [y + a], [z + 0.02], led, 6)
    if carrying:
        ax.plot([x, x], [y, y], [z, z - 0.45 * scale], color="#c0c8d8", lw=0.8, zorder=8)
        glow_points(ax, [x], [y], [z - 0.55 * scale], "#f1f4fa", 70 * scale ** 2, alpha=0.9)


rng = np.random.default_rng(3)
N_SPRAY = 160
spray_v = np.column_stack([rng.normal(0, 1, (N_SPRAY, 2)) * 0.9, rng.uniform(1.5, 4.0, N_SPRAY)])


def render(i, path):
    t = i / FPS
    c, hw, elev, azim = camera(t)
    zs = 0.022 * hw  # vertical display scale of the waves

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
    Z = zs * np.tanh(eta / 0.012)
    ax.plot_surface(X, Y, Z, facecolors=surface_colors(X, Y, eta, hw), shade=False,
                    rstride=1, cstride=1, linewidth=0, antialiased=False, zorder=1)

    ts = t - T_IMPACT
    inwin = lambda px, py: abs(px - c[0]) < hw * 0.97 and abs(py - c[1]) < hw * 0.97

    # blue glowing wave front
    if ts > 0.15:
        r = sim.c * ts - 0.15
        th = np.linspace(0, 2 * np.pi, 400)
        fx, fy = DROP[0] + r * np.cos(th), DROP[1] + r * np.sin(th)
        keep = (np.abs(fx - c[0]) < hw * 0.95) & (np.abs(fy - c[1]) < hw * 0.95)
        fx, fy = np.where(keep, fx, np.nan), np.where(keep, fy, np.nan)
        fade = np.clip(1.4 - ts / 9.0, 0.25, 1)
        glow_line(ax, fx, fy, np.full_like(fx, zs * 0.25), "#3fd0ff", 1.3, alpha=fade)

    # HYDAR arrays: purple glow when the wave reaches them
    for a, ta in zip(arrays, arrivals):
        dt = ts - ta
        lit = 0.0 if dt < -0.15 else (1.0 + 2.0 * np.exp(-max(dt, 0) / 0.35))
        if not inwin(a.x, a.y):
            continue
        px, py = a.positions[:, 0], a.positions[:, 1]
        ax.scatter(px, py, np.full(9, zs * 0.3), s=6, color="#5b3b7a", depthshade=False, zorder=5)
        if lit > 0:
            glow_points(ax, px, py, np.full(9, zs * 0.3), PURPLE, 7 * lit, alpha=min(lit / 1.5, 1))
        if hw > 9 and inwin(a.x, a.y):
            ax.text(a.x + 1.3, a.y - 1.3, zs * 0.4, a.name, color="#d9b3ff" if lit else "#6a5a80",
                    fontsize=9, zorder=9)

    # red lasers from each array to the HYDAR fix
    if t > T_LASER:
        for n, f in enumerate(fix["fixes"]):
            s = ease((t - T_LASER - 0.25 * n) / 0.9)
            if s <= 0:
                continue
            p = f.center + s * (est - f.center)
            zl = zs * 0.35
            glow_line(ax, [f.center[0], p[0]], [f.center[1], p[1]], [zl, zl], RED, 1.6)
        t_hit = T_LASER + 0.25 * 3 + 0.9
        if t > t_hit:
            k = t - t_hit
            th = np.linspace(0, 2 * np.pi, 80)
            for ph in (0.0, 0.5, 1.0):
                rr = 0.4 + 2.2 * ((k + ph) % 1.5) / 1.5
                al = 1 - ((k + ph) % 1.5) / 1.5
                glow_line(ax, est[0] + rr * np.cos(th), est[1] + rr * np.sin(th),
                          np.full(80, zs * 0.35), RED, 1.0, alpha=al)
            glow_points(ax, [est[0]], [est[1]], [zs * 0.4], RED, 30)

    # drone and sphere
    hover = np.array([DROP[0], DROP[1], 1.0])
    if t < T_IMPACT + 3.0:
        if t < 2.0:
            s = ease(t / 2.0)
            pos = hover + np.array([-9.0, 3.0, 1.0]) * (1 - s)
        elif t < T_IMPACT + 0.4:
            pos = hover.copy()
        else:
            s = ease((t - T_IMPACT - 0.4) / 2.6)
            pos = hover + np.array([7.0, 6.0, 5.0]) * s
        pos = pos + np.array([0, 0, 0.08 * np.sin(5 * t)])
        draw_drone(ax, pos, 1.0, t, carrying=t < T_RELEASE)
        if T_RELEASE <= t < T_IMPACT:
            s = ((t - T_RELEASE) / (T_IMPACT - T_RELEASE)) ** 2
            zb = lerp(hover[2] - 0.55, 0.2, s)
            glow_points(ax, [DROP[0]], [DROP[1]], [zb], "#f1f4fa", 70, alpha=0.9)

    # splash
    if 0 <= ts < 1.4:
        tt = ts
        P = np.column_stack([DROP[0] + spray_v[:, 0] * tt, DROP[1] + spray_v[:, 1] * tt,
                             (spray_v[:, 2] * tt - 4.9 * tt ** 2)])
        up = P[:, 2] > 0
        al = max(1 - tt / 1.4, 0)
        ax.scatter(P[up, 0], P[up, 1], P[up, 2], s=7, color="#e6faff", alpha=al,
                   edgecolors="none", depthshade=False, zorder=7)
        ax.scatter(P[up, 0], P[up, 1], P[up, 2], s=40, color="#7fe3ff", alpha=0.15 * al,
                   edgecolors="none", depthshade=False, zorder=7)

    # captions
    fig.text(0.035, 0.92, "HYDAR", color="white", fontsize=22, weight="bold", family="DejaVu Sans")
    fig.text(0.035, 0.875, "Hydrodynamic Detection and Ranging", color="#9fb3d9", fontsize=10)
    if t < T_IMPACT:
        cap = "drone drops a sphere"
    elif t < T_IMPACT + max(arrivals) + 0.3:
        cap = "surface waves reach the sensor arrays"
    elif t < T_LASER + 1.7:
        cap = "time of flight → bearing and range"
    else:
        cap = f"HYDAR fix: {err_cm:.0f} cm from the true drop point"
    col = RED if t > T_LASER + 1.7 else "#c9d6f0"
    fig.text(0.035, 0.06, cap, color=col, fontsize=11)
    fig.savefig(path, facecolor=BG)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="hydar_frames_")
    print(f"rendering {N_FRAMES} frames ...")
    for i in range(N_FRAMES):
        render(i, os.path.join(tmp, f"f{i:04d}.png"))
        if i % 20 == 0:
            print(f"  {i}/{N_FRAMES}")
    pat = os.path.join(tmp, "f%04d.png")
    mp4 = os.path.join(OUT, "hydar_hero.mp4")
    gif = os.path.join(OUT, "hydar_hero.gif")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pat,
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "22", mp4], check=True)
    pal = os.path.join(tmp, "pal.png")
    scale = "scale=720:-1:flags=lanczos"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pat,
                    "-vf", f"{scale},palettegen=max_colors=160:stats_mode=diff", pal], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pat, "-i", pal,
                    "-lavfi", f"{scale}[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle",
                    "-loop", "0", gif], check=True)
    shutil.copy(os.path.join(tmp, f"f{N_FRAMES - 1:04d}.png"), os.path.join(OUT, "hydar_hero_last.png"))
    shutil.rmtree(tmp)
    for f in (gif, mp4):
        print(f"wrote {f} ({os.path.getsize(f) / 1e6:.1f} MB)")


if __name__ == "__main__":
    if os.environ.get("HYDAR_PREVIEW"):  # render a few stills for checking
        os.makedirs(OUT, exist_ok=True)
        for t in map(float, os.environ["HYDAR_PREVIEW"].split(",")):
            render(int(t * FPS), os.path.join(OUT, f"preview_{t:05.2f}.png"))
    else:
        main()
