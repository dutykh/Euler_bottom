#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figures of the README, drawn from the data that example_data.m and frequency_data.m write.

Run first, from this folder, in MATLAB or GNU Octave:
    example_data                                        (the example, version 2.0.0)
    frequency_data                                      (steady waves, version 2.0.0)
and, for the two comparisons, the same scripts on a copy of version 1.0 beside this one
(git worktree add ../Euler_bottom-v1.0 v1.0, with the semicolon of line 33 of its
code/methods/RHS.m restored, so that it does not print the bed at every evaluation):
    codedir = '../../../Euler_bottom-v1.0/code'; tag = 'v1'; example_data
    codedir = '../../../Euler_bottom-v1.0/code'; tag = 'v1'; frequency_data
Then
    python3 plot_figures.py
writes example.png, corrections-frequency.png and corrections-example.png. The style is the
stamped copy manuscript_style.py beside this script.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

G, K0 = 9.8, 0.05            # gravity and carrier wavenumber of the example
hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 264
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER


def group_speed(k, h):
    omega = np.sqrt(G * k * np.tanh(k * h))
    return 0.5 * omega / k * (1.0 + 2.0 * k * h / np.sinh(2.0 * k * h))


def read(name):
    return np.genfromtxt(HERE / name, delimiter=",", names=True)


def line(index):
    style = hs.series(index)
    style["marker"] = None
    return style


def dots(index):
    style = hs.series(index)
    style["linestyle"] = "none"
    return style


# ---------------------------------------------------------------- the corrected example
bed = read("example_v2_bed.csv")
prof = read("example_v2_profiles.csv")
times = (0.0, 100.0, 200.0)

fig = plt.figure(figsize=(hs.WIDTH, 5.0), constrained_layout=True)
grid = fig.add_gridspec(4, 1, height_ratios=(1.5, 1.0, 1.0, 1.0))
top = fig.add_subplot(grid[0])
order = np.argsort(bed["x"])
xb, yb = bed["x"][order] / 1000.0, bed["y"][order]
top.fill_between(xb, -60.0, yb, color=hs.SAND, lw=0)
first = prof[prof["t"] == 0.0]
o = np.argsort(first["x"])
xs, es = first["x"][o] / 1000.0, first["eta"][o]
top.fill_between(xs, np.interp(xs, xb, yb), es, color=hs.NAVY, alpha=0.10, lw=0)
top.plot(xb, yb, color=hs.INK, lw=1.0)
top.plot(xs, es, color=hs.NAVY, lw=0.7)
top.axvline(10.0, color=hs.SLATE, lw=0.7, linestyle=(0, (1.6, 1.6)))
top.annotate("", xy=(4.6, -12.0), xytext=(2.5, -12.0),
             arrowprops=dict(arrowstyle="->", color=hs.SLATE, lw=0.8))
top.annotate("", xy=(15.4, -12.0), xytext=(17.5, -12.0),
             arrowprops=dict(arrowstyle="->", color=hs.SLATE, lw=0.8))
top.text(10.25, -44.0, "mirror line", color=hs.SLATE, fontsize=7.5, ha="left", va="center")
top.set_xlim(0.0, 20.0)
top.set_ylim(-56.0, 6.0)
top.set_ylabel(r"$y$ (m)")
top.tick_params(labelbottom=False)
hs.panel(top, "(a) the bed and the free surface at $t = 0$: two mirror-image packets over 50 m of water")

for row, t in enumerate(times, start=1):
    ax = fig.add_subplot(grid[row], sharex=top)
    snap = prof[prof["t"] == t]
    o = np.argsort(snap["x"])
    for xstep in (5.0, 15.0):
        ax.axvline(xstep, color=hs.MIST, lw=0.8)
    ax.axvline(10.0, color=hs.SLATE, lw=0.7, linestyle=(0, (1.6, 1.6)))
    ax.plot(snap["x"][o] / 1000.0, snap["eta"][o], color=hs.NAVY, lw=0.6)
    ax.set_ylim(-1.4, 1.4)
    ax.set_ylabel(r"$\eta$ (m)")
    hs.time_stamp(ax, t, label=rf"$t = {t:.0f}$ s", corner="upper right")
    if row < len(times):
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel(r"$x$ (km)")
    hs.panel(ax, f"({'bcd'[row - 1]})")
hs.save(fig, HERE / "example.png")
print("written", HERE / "example.png")

# ---------------------------------------------------------------- frequency of steady waves
f1, f2 = read("frequency_v1.csv"), read("frequency_v2.csv")
fig, (left, right) = plt.subplots(1, 2, figsize=(hs.WIDTH, 2.7), constrained_layout=True)
for i, kh in enumerate((0.5, 1.0, 2.0)):
    a, b = f1[np.isclose(f1["kh"], kh)], f2[np.isclose(f2["kh"], kh)]
    left.loglog(a["kA"], a["predicted_v1"], lw=0.9, **line(i))
    left.loglog(a["kA"], a["error"], ms=4.5, **dots(i))
    right.loglog(b["kA"], np.abs(b["error"]), ms=4.5, **dots(i))
for ax in (left, right):
    ax.set_xlabel(r"steepness $kA$")
    ax.set_xlim(7e-4, 0.14)
left.set_ylim(1e-8, 0.1)
right.set_ylim(1e-15, 1e-9)
left.set_ylabel(r"$\omega/(kc) - 1$")
right.set_ylabel(r"$|\omega/(kc) - 1|$")
hs.panel(left, "(a) version 1.0: markers measured, lines $-3q/(2c)$")
hs.panel(right, "(b) version 2.0.0")
handles = [Line2D([], [], ms=4.5, lw=0.9, **hs.series(i)) for i in range(3)]
hs.legend_outside(fig, where="lower", ncol=3, handles=handles,
                  labels=[r"$kh = 0.5$", r"$kh = 1$", r"$kh = 2$"])
hs.save(fig, HERE / "corrections-frequency.png")
print("written", HERE / "corrections-frequency.png")

# ---------------------------------------------------------------- the example, both versions
b1, b2 = read("example_v1_bed.csv"), read("example_v2_bed.csv")
s1, s2 = read("example_v1_series.csv"), read("example_v2_series.csv")
p1, p2 = read("example_v1_profiles.csv"), read("example_v2_profiles.csv")
depth_v1 = -b1["y"].max()                    # depth under the packet in version 1.0

fig = plt.figure(figsize=(hs.WIDTH, 4.3), constrained_layout=True)
grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.15))
top = fig.add_subplot(grid[0, :])
for bed_, index in ((b2, 0), (b1, 1)):
    o = np.argsort(bed_["x"])
    top.plot(bed_["x"][o] / 1000.0, bed_["y"][o], lw=1.0, **line(index))
top.set_xlim(0.0, 20.0)
top.set_ylim(-60.0, -20.0)
top.set_xlabel(r"$x$ (km)")
top.set_ylabel(r"bed $y$ (m)")
hs.panel(top, "(a) the bed described by the conformal map of each version")

mid = fig.add_subplot(grid[1, 0])
for s, index in ((s2, 0), (s1, 1)):
    mid.plot(s["t"], s["centroid"] - s["centroid"][0], ms=3.5, **dots(index))
t = np.linspace(0.0, 200.0, 2)
for h, index in ((50.0, 2), (depth_v1, 3)):
    mid.plot(t, group_speed(K0, h) * t, lw=0.8, **line(index))
mid.set_xlim(0.0, 205.0)
mid.set_ylim(0.0, 1800.0)
mid.set_xlabel(r"$t$ (s)")
mid.set_ylabel("displacement (m)")
hs.panel(mid, "(b) travel of the packet")

low = fig.add_subplot(grid[1, 1])
for p, index in ((p2, 0), (p1, 1)):
    snap = p[(p["t"] == 200.0) & (p["x"] < 10000.0)]
    o = np.argsort(snap["x"])
    low.plot(snap["x"][o] / 1000.0, snap["eta"][o], lw=0.7, **line(index))
low.set_xlim(3.0, 5.2)
low.set_ylim(-1.3, 1.3)
low.set_xlabel(r"$x$ (km)")
low.set_ylabel(r"$\eta$ (m)")
hs.time_stamp(low, 200.0, label=r"$t = 200$ s", corner="upper left")
hs.panel(low, "(c) the packet")
handles = [Line2D([], [], ms=3.5, lw=1.0, **hs.series(0)), Line2D([], [], ms=3.5, lw=1.0, **hs.series(1)),
           Line2D([], [], lw=0.8, **line(2)), Line2D([], [], lw=0.8, **line(3))]
labels = ["version 2.0.0", "version 1.0", "linear group speed over 50 m",
          f"linear group speed over {depth_v1:.2f} m"]
hs.legend_outside(fig, where="lower", ncol=2, handles=handles, labels=labels)
hs.save(fig, HERE / "corrections-example.png")
print("written", HERE / "corrections-example.png")
