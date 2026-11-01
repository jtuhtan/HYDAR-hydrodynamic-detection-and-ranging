"""paperfig: Jeff Tuhtan's paper figure style for matplotlib.

A blend of
  * Tufte      - high data-ink ratio, no chart junk, direct labels, range frames
  * gramm      - grouped colour/marker mapping, raw data shown next to summaries
  * ggplot2    - consistent theme, faint grid as reading aid, small multiples
  * Jeff's own - Arial, black/grey bars with thin black edges, blue/red group
                 pairs, cyan accent for events and call-outs, bold lowercase
                 panel letters, value labels on bars, units in parentheses.

Usage
-----
    import paperfig as pf
    pf.use()                                  # apply the style once
    pf.page(170)                              # page text width (ask the author)
    fig, ax = pf.figure("half")               # 79 mm: two figures side by side
    ax.plot(t, y, color=pf.C["blue"])
    pf.save(fig, "fig01_timeseries")          # -> pdf + png + svg
"""
from __future__ import annotations

import string
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

# ---------------------------------------------------------------- palette ---
# Named colours. Groups come in blue/red pairs first (as in the BDS papers),
# greys carry "neutral" data, cyan is reserved for highlights and events.
C = {
    "ink": "#1A1A1A",        # text, axes, black edges
    "dark": "#595959",       # dark grey bars / secondary text
    "grey": "#7F7F7F",       # primary bar fill (group 1)
    "light": "#E3E3E3",      # secondary bar fill (group 2)
    "grid": "#E6E6E6",       # grid lines
    "blue": "#1F6FB4",       # group A (e.g. smolt, Fall 1, left axis)
    "red": "#C0262D",        # group B (e.g. kelt, Fall 2, right axis)
    "cyan": "#00A3D9",       # accent: events, call-outs, water
    "orange": "#E08A1E",
    "green": "#3E8E41",
    "purple": "#6A4C93",
}

# Default colour cycle: blue, red, dark grey, cyan, orange, green, purple.
CYCLE = [C["blue"], C["red"], C["dark"], C["cyan"], C["orange"], C["green"], C["purple"]]
# Grey ramp for bars / pies where colour is not needed (light last).
GREYS = ["#404040", C["grey"], "#A6A6A6", "#CCCCCC", C["light"]]
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]   # gramm-style marker mapping

# Spacing rules, in units of the font size (em). Measured from Jeff's published
# figures (see references/placement_rules.md, section 6):
TEXT_GAP = 1.0    # min clear space between two independent labels
GROUP_GAP = 0.3   # min clear space between lines of one group (legend rows, stacked labels)
DATA_GAP = 0.3    # min clear space between text and any data mark (line, point, bar)

# ------------------------------------------------------------ page sizes ----
# Figures are sized as a fraction of the page TEXT width (the space between the
# margins, or one column of a two-column journal page counts as "half"). Ask the
# author for the text width before making figures and set it with page().
# Measured in Jeff's papers: printed figures fill 97-99 % of their column/text
# width; his own drafts leave more room. Rule: every figure keeps WIGGLE_MM of
# clear space on its left and right, so n figures side by side each get
#     width = text_width / n - 2 * WIGGLE_MM
FRACTIONS = {"full": 1, "half": 2, "third": 3}
ALIASES = {"double": "full", "single": "half", "column": "half"}
# default height/width per size, medians of single-row figures in his papers
ASPECT = {"full": 0.42, "half": 0.70, "third": 0.85}
TEXT_WIDTH_MM = None      # set with page(); 170 mm (A4, 20 mm margins) if unset
WIGGLE_MM = 3.0
COMMON_TEXT_WIDTHS = {    # mm, for reference when the author is unsure
    "A4, 20 mm margins (Word default for drafts, Sci Rep, reports)": 170,
    "A4, 25 mm margins": 160,
    "Letter, 1 in margins": 165,
    "Elsevier / IEEE / MDPI two-column page": 180,
    "Springer, Wiley two-column page": 174,
    "IOP single column": 160,
    "PLOS ONE (figures max)": 132,
}
MM = 1 / 25.4

FONT_STACK = ["Arial", "Liberation Sans", "Helvetica", "Nimbus Sans", "DejaVu Sans"]


def _font_available(name: str) -> bool:
    return any(f.name == name for f in font_manager.fontManager.ttflist)


def use(base_size: float = 8.0, grid: bool = True) -> str:
    """Apply the paper style globally. Returns the font family actually used.

    base_size: tick/legend size in pt (8 pt at every figure size);
    axis labels use base_size + 1, panel letters base_size + 2.
    """
    font = next((f for f in FONT_STACK if _font_available(f)), "sans-serif")
    small = base_size
    mpl.rcParams.update({
        # fonts
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "font.size": small,
        "mathtext.fontset": "custom",
        "mathtext.rm": font, "mathtext.cal": font, "mathtext.sf": font, "mathtext.it": f"{font}:italic", "mathtext.bf": f"{font}:bold",
        "axes.labelsize": small + 1,
        "axes.titlesize": small + 1,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "xtick.labelsize": small,
        "ytick.labelsize": small,
        "legend.fontsize": small - 0.5,
        # colours / lines
        "axes.prop_cycle": mpl.cycler(color=CYCLE),
        "text.color": C["ink"],
        "axes.labelcolor": C["ink"],
        "axes.edgecolor": C["ink"],
        "xtick.color": C["ink"],
        "ytick.color": C["ink"],
        "lines.linewidth": 1.2,
        "lines.markersize": 4,
        "patch.linewidth": 0.6,
        "patch.edgecolor": C["ink"],
        # axes: Tufte - only left/bottom spines, ticks out
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.facecolor": "white",
        "axes.axisbelow": True,
        "axes.grid": grid,
        "axes.grid.axis": "both",
        "grid.color": C["grid"],
        "grid.linewidth": 0.6,
        "grid.linestyle": "-",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.visible": False,
        "ytick.minor.visible": False,
        # legend: thin grey box, as in the BDS figures
        "legend.frameon": True,
        "legend.framealpha": 1.0,
        "legend.edgecolor": "#A6A6A6",
        "legend.fancybox": False,
        "legend.borderpad": 0.4,
        "legend.handlelength": 1.6,
        "legend.labelspacing": 0.3,
        # bars / errorbars
        "errorbar.capsize": 2.5,
        # boxplots (gramm/ggplot: thin, black, open median)
        "boxplot.flierprops.marker": "o",
        "boxplot.flierprops.markersize": 3,
        "boxplot.flierprops.markeredgecolor": C["dark"],
        "boxplot.medianprops.color": C["ink"],
        "boxplot.medianprops.linewidth": 1.2,
        # output
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": None,          # keep the exact page-fraction size
        "savefig.pad_inches": 0.0,
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,   # embed TrueType: text stays editable in Illustrator/Inkscape
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "figure.constrained_layout.use": True,
        # >= 1 em of clear space between neighbouring panels (TEXT_GAP)
        "figure.constrained_layout.w_pad": TEXT_GAP * base_size / 72,
        "figure.constrained_layout.h_pad": TEXT_GAP * base_size / 72 / 2,
    })
    return font


# ---------------------------------------------------------------- layout ----
def page(text_width_mm: float, wiggle_mm: float = 3.0):
    """Set the page text width (mm) that figure widths are fractions of.

    text_width_mm: page width minus left and right margins (e.g. A4 with 20 mm
    margins = 170). For a two-column journal give the full text block (both
    columns + gutter); 'half' then gives a one-column figure.
    wiggle_mm: clear space kept on each side of every figure.
    """
    global TEXT_WIDTH_MM, WIGGLE_MM
    TEXT_WIDTH_MM, WIGGLE_MM = float(text_width_mm), float(wiggle_mm)
    return {k: width_mm(k) for k in FRACTIONS}


def width_mm(width="half") -> float:
    """Figure width in mm for 'full', 'half', 'third' (or a number in mm)."""
    if not isinstance(width, str):
        return float(width)
    key = ALIASES.get(width, width)
    if key not in FRACTIONS:
        raise ValueError(f"width must be one of {list(FRACTIONS)} or mm, not {width!r}")
    global TEXT_WIDTH_MM
    if TEXT_WIDTH_MM is None:
        print("[size] page text width not set: using 170 mm (A4, 20 mm margins). "
              "Ask the author and call pf.page(text_width_mm).")
        TEXT_WIDTH_MM = 170.0
    return TEXT_WIDTH_MM / FRACTIONS[key] - 2 * WIGGLE_MM


def figure(width="half", aspect: float | None = None, nrows: int = 1, ncols: int = 1,
           height_mm: float | None = None, **subplot_kw):
    """Create a figure sized as a fraction of the page text width.

    width: 'full', 'half' (2 side by side, or one journal column), 'third'
    (3 side by side) or a number in mm. See page().
    aspect: height/width of the whole figure; default from ASPECT per size,
    multiplied by the number of rows (ignored if height_mm is given).
    The saved file has exactly this size: insert it at 100 % and never rescale.
    """
    w_mm = width_mm(width)
    if aspect is None:
        key = ALIASES.get(width, width) if isinstance(width, str) else "half"
        aspect = ASPECT.get(key, 0.7) * (1 + 0.8 * (nrows - 1))
    h_mm = height_mm if height_mm is not None else w_mm * aspect
    fig, axes = plt.subplots(nrows, ncols, figsize=(w_mm * MM, h_mm * MM), **subplot_kw)
    fig._pf_size = (width, w_mm, h_mm)
    return fig, axes


PANEL_FORMATS = {"a": "{}", "(a)": "({})", "a)": "{})", "A": "{}"}


def panel_labels(axes, style: str = "a", loc: str = "outside", start: int = 0,
                 bold=None, size=None, pad: float = 5, **kw):
    """Number the panels of a figure.

    style: 'a' (default, bold lowercase - current house style), '(a)' (IOP / Elsevier
           style, regular weight), 'a)' or 'A' (bold capitals). Follow the journal.
    loc:   'outside' - top-left, aligned with the left edge of the y-axis label and
           with the top of the axes (default; never collides with data);
           'inside'  - inside the top-left corner of the axes (use for maps, images,
           and tightly packed grids with no space above the axes).
    The letter is anchored to the y-label at draw time, so it stays aligned across
    stacked panels even when tick labels have different widths.
    """
    axes = np.atleast_1d(axes).ravel()
    letters = string.ascii_uppercase if style == "A" else string.ascii_lowercase
    fmt = PANEL_FORMATS.get(style, "{}")
    if bold is None:
        bold = style != "(a)"
    fs = size or mpl.rcParams["font.size"] + (2 if bold else 1)
    out = []
    for ax, letter in zip(axes, letters[start:]):
        text = fmt.format(letter)
        common = dict(fontweight="bold" if bold else "normal", fontsize=fs, **kw)
        if loc == "inside":
            t = ax.annotate(text, xy=(0, 1), xycoords="axes fraction", xytext=(4, -3),
                            textcoords="offset points", ha="left", va="top", zorder=10,
                            bbox=dict(fc="white", ec="none", pad=1, alpha=0.8), **common)
        else:
            has_ylab = bool(ax.yaxis.label.get_text()) and ax.axison

            def anchor(renderer, ax=ax, has_ylab=has_ylab):
                bb = ax.get_window_extent(renderer)
                if has_ylab:
                    bb = mpl.transforms.Bbox.union([bb, ax.yaxis.label.get_window_extent(renderer)])
                # letters of panels in the same row share one baseline
                row = [o for o in axes if o.get_subplotspec() is not None and ax.get_subplotspec() is not None
                       and o.get_subplotspec().rowspan.start == ax.get_subplotspec().rowspan.start]
                top = max([bb.y1] + [o.get_window_extent(renderer).y1 for o in row])
                return mpl.transforms.Bbox([[bb.x0, bb.y0], [bb.x1, top]])
            t = ax.annotate(text, xy=(0, 1), xycoords=anchor,
                            xytext=(0 if has_ylab else -14, pad), textcoords="offset points",
                            ha="left", va="bottom", annotation_clip=False, **common)
        t.set_gid("panel")
        out.append(t)
    return out


def stack(nrows: int = 3, width="half", height_mm=None, row_mm: float = 28,
          sharey: bool = False, **kw):
    """Vertically stacked panels sharing the x-axis (time series, profiles).

    Only the bottom panel carries x tick labels and the x label; the y labels are
    aligned. Returns fig, axes (array).
    """
    h = height_mm or row_mm * nrows + 10
    fig, axes = figure(width, nrows=nrows, ncols=1, height_mm=h, sharex=True, sharey=sharey, **kw)
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    fig.align_ylabels(axes)
    return fig, axes


def legend(ax, loc: str = "best", last=(), frame: bool = True, **kw):
    """House legend: inside the axes, in the emptiest corner, thin grey frame.

    Data series come first; entries whose labels are listed in `last` (reference
    and fit lines) are moved to the end. Use frame=False on plots without a grid.
    """
    handles, labels = ax.get_legend_handles_labels()
    order = [i for i, l in enumerate(labels) if l not in last] + \
            [i for i, l in enumerate(labels) if l in last]
    return ax.legend([handles[i] for i in order], [labels[i] for i in order],
                     loc=loc, frameon=frame, **kw)


def figure_legend(fig, axes, ncol=None, **kw):
    """One shared legend for a multi-panel figure, centred above the panels."""
    handles, labels = np.atleast_1d(axes).ravel()[0].get_legend_handles_labels()
    return fig.legend(handles, labels, loc="outside upper center",
                      ncol=ncol or len(labels), frameon=False, **kw)


def despine(ax, left=True, bottom=True, offset: float = 0):
    """Keep only the requested spines; optionally offset them (Tufte 'floating' axes)."""
    for side, keep in (("left", left), ("bottom", bottom), ("top", False), ("right", False)):
        ax.spines[side].set_visible(keep)
        if keep and offset:
            ax.spines[side].set_position(("outward", offset))


def range_frame(ax, x=None, y=None, offset: float = 4):
    """Tufte range frame: spines span only the data range, offset from the plot."""
    ax.grid(False)
    despine(ax, offset=offset)
    if x is not None:
        ax.spines["bottom"].set_bounds(np.nanmin(x), np.nanmax(x))
    if y is not None:
        ax.spines["left"].set_bounds(np.nanmin(y), np.nanmax(y))


def grid_y_only(ax):
    """Horizontal grid only (bar charts, box plots)."""
    ax.grid(False, axis="x")
    ax.grid(True, axis="y")


# ------------------------------------------------------------ annotation ----
def label_bars(ax, bars, fmt: str = "{:.2f}", pad: float = 3, tops=None,
               skip_zero: bool = True, **kw):
    """Write the value above each bar (direct labelling instead of reading the axis).

    tops: optional y positions for the labels (e.g. upper CI) so they clear error bars.
    skip_zero: leave empty bars unlabelled (less clutter).
    """
    style = dict(fontsize=mpl.rcParams["font.size"] - 1, color=C["ink"], ha="center", va="bottom")
    style.update(kw)
    for i, b in enumerate(bars):
        h = b.get_height()
        if skip_zero and h == 0:
            continue
        x = b.get_x() + b.get_width() / 2
        y = h if tops is None else tops[i]
        ax.annotate(fmt.format(h), xy=(x, y), xytext=(0, pad), textcoords="offset points", **style)


def label_lines(ax, lines=None, labels=None, dy=None, pad: float = 4, min_gap=None, **kw):
    """Direct labels at the right end of each line, coloured like the line (no legend).

    Labels are kept in the vertical order of the line ends and pushed apart to at
    least `min_gap` (axes fraction; default = 1.3 text lines) so they never overlap.
    dy: extra nudges in points.
    """
    if min_gap is None:
        h_pt = ax.get_position().height * ax.figure.get_figheight() * 72 * 0.8  # after layout
        min_gap = 1.3 * mpl.rcParams["font.size"] / h_pt
    lines = lines or ax.get_lines()
    labels = labels or [l.get_label() for l in lines]
    dy = dy or [0] * len(lines)
    lo, hi = ax.get_ylim()
    ends = []
    for line in lines:
        x, y = (np.asarray(v, float) for v in line.get_data())
        ok = np.isfinite(y)
        ends.append((x[ok][-1], (y[ok][-1] - lo) / (hi - lo)))
    order = np.argsort([e[1] for e in ends])
    pos = np.array([ends[i][1] for i in order])
    for k in range(1, len(pos)):                      # push up
        pos[k] = max(pos[k], pos[k - 1] + min_gap)
    shift = max(0, pos[-1] - 0.98) if len(pos) else 0  # keep inside the axes
    pos -= shift
    for k in range(len(pos) - 2, -1, -1):             # push down if needed
        pos[k] = min(pos[k], pos[k + 1] - min_gap)
    for k, i in enumerate(order):
        line = lines[i]
        t = ax.annotate(labels[i], xy=(1, pos[k]), xycoords="axes fraction", xytext=(pad, dy[i]),
                        textcoords="offset points", color=line.get_color(), va="center", ha="left",
                        fontsize=mpl.rcParams["font.size"], annotation_clip=False, **kw)
        t.set_gid(f"group-lines-{id(ax)}")       # one group: rows may sit GROUP_GAP apart


def callout(ax, text, xy, xytext, color=None, **kw):
    """Rounded call-out box with an arrow (cyan edge by default)."""
    color = color or C["cyan"]
    return ax.annotate(
        text, xy=xy, xytext=xytext, fontsize=mpl.rcParams["font.size"] - 0.5,
        ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.3,rounding_size=0.25", fc="white", ec=color, lw=0.9),
        arrowprops=dict(arrowstyle="-|>", color=C["ink"], lw=0.8, shrinkA=0, shrinkB=2),
        **kw)


def event_line(ax, x, label=None, color=None, **kw):
    """Vertical marker for an event (e.g. release, impact).

    The label sits above the axes frame, centred on the line, so it can never
    cover data (text is not allowed to occlude lines or points).
    """
    color = color or C["cyan"]
    ax.axvline(x, color=color, lw=1.0, zorder=1, **kw)
    if label:
        ax.annotate(label, xy=(x, 1), xycoords=("data", "axes fraction"), xytext=(0, 3),
                    textcoords="offset points", ha="center", va="bottom", color=color,
                    fontsize=mpl.rcParams["font.size"] - 0.5, annotation_clip=False)


def shade_span(ax, x0, x1, color=None, alpha=0.12, label=None):
    """Light shaded interval (e.g. treatment period); the label sits above the frame,
    centred on the span, so it never covers the data."""
    color = color or C["cyan"]
    ax.axvspan(x0, x1, color=color, alpha=alpha, lw=0, zorder=0)
    if label:
        xm = x0 + (x1 - x0) / 2
        ax.annotate(label, xy=(xm, 1), xycoords=("data", "axes fraction"), xytext=(0, 3),
                    textcoords="offset points", ha="center", va="bottom", color=color,
                    fontsize=mpl.rcParams["font.size"] - 0.5, annotation_clip=False)


def ref_line(ax, y=None, x=None, text=None, color=None, ls="--", where="left", **kw):
    """Reference / threshold line labelled on the line itself, in the line colour.

    e.g. ref_line(ax, y=0.5, text="KLD = 0.5", color=pf.C["red"]). No legend entry.
    Choose `where` so the label sits over empty space; check_layout() flags it if not.
    """
    color = color or C["ink"]
    fs = mpl.rcParams["font.size"] - 0.5
    gap = DATA_GAP * fs + 1
    if y is not None:
        ax.axhline(y, color=color, ls=ls, lw=0.8, zorder=1, **kw)
        if text:
            xf = 0.02 if where == "left" else 0.98
            ax.annotate(text, xy=(xf, y), xycoords=("axes fraction", "data"), xytext=(0, gap),
                        textcoords="offset points", color=color, fontsize=fs,
                        ha="left" if where == "left" else "right", va="bottom")
    else:
        ax.axvline(x, color=color, ls=ls, lw=0.8, zorder=1, **kw)
        if text:
            ax.annotate(text, xy=(x, 1), xycoords=("data", "axes fraction"), xytext=(0, 2),
                        textcoords="offset points", color=color, fontsize=fs,
                        ha="center", va="bottom", annotation_clip=False)


def one_to_one(ax, lo=None, hi=None, color=None, label="1:1"):
    """1:1 line for observed-vs-estimated plots, equal aspect and equal limits."""
    if lo is None or hi is None:
        lims = np.r_[ax.get_xlim(), ax.get_ylim()]
        lo, hi = lims.min() if lo is None else lo, lims.max() if hi is None else hi
    ax.plot([lo, hi], [lo, hi], color=color or C["red"], lw=0.9, zorder=2, label=label)
    ax.set(xlim=(lo, hi), ylim=(lo, hi))
    ax.set_aspect("equal", adjustable="box")


def marker_label(ax, x, y, text, color=None, offset=(0, 0), size=None, **kw):
    """Circled number/letter identifying a point, peak or part (diagrams, peaks)."""
    color = color or C["ink"]
    fs = size or mpl.rcParams["font.size"] - 1
    return ax.annotate(text, xy=(x, y), xytext=offset, textcoords="offset points",
                       ha="center", va="center", fontsize=fs, color=color, zorder=6,
                       bbox=dict(boxstyle="circle,pad=0.2", fc="white", ec=color, lw=0.8), **kw)


def highlight_box(ax, x0, x1, y0=None, y1=None, color=None, lw=1.2, **kw):
    """Outline a region of interest with a coloured rectangle (no fill)."""
    from matplotlib.patches import Rectangle
    color = color or C["red"]
    if y0 is None:
        y0, y1 = ax.get_ylim()
    r = Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=color, lw=lw, zorder=5,
                  clip_on=False, **kw)
    ax.add_patch(r)
    return r


def flow_arrow(ax, xy, dxdy=(0.08, 0), text="Flow", color=None, coords="axes fraction"):
    """Block arrow indicating flow direction, labelled above (cyan by default)."""
    color = color or C["cyan"]
    x, y = xy
    ax.annotate("", xy=(x + dxdy[0], y + dxdy[1]), xytext=(x, y), xycoords=coords,
                arrowprops=dict(arrowstyle="simple,head_width=0.8,head_length=0.8,tail_width=0.35",
                                fc=color, ec="none"))
    if text:
        ax.annotate(text, xy=(x + dxdy[0] / 2, y), xycoords=coords, xytext=(0, 5),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=mpl.rcParams["font.size"] - 0.5, color=color)


def colorbar(fig, mappable, ax, label, orientation="vertical", **kw):
    """Colour bar: vertical bars carry the label (with units) horizontally on top;
    horizontal bars sit below the panel with the label underneath."""
    if orientation == "vertical":
        cb = fig.colorbar(mappable, ax=ax, orientation="vertical", fraction=0.05, pad=0.02,
                          aspect=25, **kw)
        cb.ax.set_title(label, fontsize=mpl.rcParams["font.size"], loc="left", pad=4)
    else:
        cb = fig.colorbar(mappable, ax=ax, orientation="horizontal", fraction=0.06, pad=0.12,
                          aspect=30, **kw)
        cb.set_label(label)
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(length=2)
    return cb


# ---------------------------------------------------------------- recipes ---
def twin_y(ax, left_color=None, right_color=None):
    """Second y-axis whose spine, ticks and label match the right-hand series colour.

    Returns the right axes. Colour the left axis to match its series too, so the
    reader never needs a legend to pair line and axis.
    """
    left_color = left_color or C["blue"]
    right_color = right_color or C["red"]
    ax2 = ax.twinx()
    ax2.grid(False)
    ax2.spines["right"].set_visible(True)
    for a, side, col in ((ax, "left", left_color), (ax2, "right", right_color)):
        a.spines[side].set_color(col)
        a.tick_params(axis="y", colors=col)
        a.yaxis.label.set_color(col)
    return ax2


def grouped_bars(ax, categories, groups: dict, colors=None, width=0.8,
                 errors: dict | None = None, labels: bool = True, fmt="{:.2f}"):
    """Grouped bar chart: greys with thin black edges (Jeff's default)."""
    n = len(groups)
    colors = colors or ([C["grey"], C["light"]] if n == 2 else GREYS[:n])
    x = np.arange(len(categories))
    w = width / n
    out = []
    for i, (name, vals) in enumerate(groups.items()):
        err = None if errors is None else errors.get(name)
        bars = ax.bar(x + (i - (n - 1) / 2) * w, vals, w, label=name, color=colors[i],
                      edgecolor=C["ink"], linewidth=0.6, yerr=err,
                      error_kw=dict(ecolor=C["ink"], elinewidth=0.8))
        if labels:
            label_bars(ax, bars, fmt=fmt)
        out.append(bars)
    ax.set_xticks(x, categories)
    ax.tick_params(axis="x", length=0)
    grid_y_only(ax)
    return out


def box_with_points(ax, data, labels, colors=None, width=0.5, jitter=0.08,
                    point_size=8, style="tufte", seed=0):
    """Distribution per group with every observation shown.

    style='tufte' (default): Tufte's quartile plot drawn to the LEFT of the raw
        points - thin line for each whisker (to the furthest point within 1.5 IQR),
        thick line for the interquartile range, median as a short gap with a white
        dot. The jittered observations sit to the right of it.
    style='box': gramm/ggplot2 open box (IQR), whiskers, thick median line, with
        the jittered raw points on top.
    """
    rng = np.random.default_rng(seed)
    colors = colors or CYCLE
    pos = np.arange(1, len(data) + 1)
    if style == "box":
        bp = ax.boxplot(data, positions=pos, widths=width, patch_artist=True, showfliers=False)
        for i, box in enumerate(bp["boxes"]):
            col = colors[i % len(colors)]
            box.set(facecolor=mpl.colors.to_rgba(col, 0.15), edgecolor=col, linewidth=1.0)
            for part in ("whiskers", "caps"):
                for art in bp[part][2 * i: 2 * i + 2]:
                    art.set(color=col, linewidth=1.0)
            bp["medians"][i].set(color=col, linewidth=1.8)
    span = np.nanmax([np.max(d) for d in data]) - np.nanmin([np.min(d) for d in data])
    for i, d in enumerate(data):
        col = colors[i % len(colors)]
        y = np.asarray(d, float)
        if style == "tufte":
            q1, med, q3 = np.percentile(y, [25, 50, 75])
            lo = y[y >= q1 - 1.5 * (q3 - q1)].min()
            hi = y[y <= q3 + 1.5 * (q3 - q1)].max()
            xq = pos[i] - width * 0.3                       # quartile plot on the left
            gap = 0.012 * span
            seg = dict(color=col, solid_capstyle="butt")
            ax.plot([xq, xq], [lo, q1], lw=1.0, **seg)
            ax.plot([xq, xq], [q3, hi], lw=1.0, **seg)
            ax.plot([xq, xq], [q1, med - gap], lw=3.0, **seg)
            ax.plot([xq, xq], [med + gap, q3], lw=3.0, **seg)
            ax.plot([xq], [med], "o", ms=3, mfc="white", mec=col, mew=0.9, zorder=4)
            xj = pos[i] + rng.uniform(-jitter, jitter, size=y.size) * width * 2 + width * 0.05
        else:
            xj = pos[i] + rng.uniform(-jitter, jitter, size=y.size) * width * 2
        ax.scatter(xj, y, s=point_size, color=col, alpha=0.6, lw=0, zorder=3)
    ax.set_xticks(pos - (width * 0.12 if style == "tufte" else 0), labels)
    ax.tick_params(axis="x", length=0)
    grid_y_only(ax)


def _bbox_perimeter(bb, n=12):
    xs = np.linspace(bb.x0, bb.x1, n)
    ys = np.linspace(bb.y0, bb.y1, n)
    return np.vstack([np.column_stack([xs, np.full(n, bb.y0)]), np.column_stack([xs, np.full(n, bb.y1)]),
                      np.column_stack([np.full(n, bb.x0), ys]), np.column_stack([np.full(n, bb.x1), ys])])


def _inside_patch(text, patch, renderer, margin=0.0):
    """True if the text box (grown by margin px) lies wholly inside the patch."""
    bb = mpl.text.Text.get_window_extent(text, renderer).padded(margin)
    return all(patch.contains_point(pt) for pt in _bbox_perimeter(bb))


def _spread(pos, min_gap, lo=-np.inf, hi=np.inf):
    """Push sorted 1-D positions apart to at least min_gap, staying within [lo, hi]."""
    pos = np.array(pos, float)
    for k in range(1, len(pos)):
        pos[k] = max(pos[k], pos[k - 1] + min_gap)
    if len(pos) and pos[-1] > hi:
        pos -= pos[-1] - hi
    for k in range(len(pos) - 2, -1, -1):
        pos[k] = min(pos[k], pos[k + 1] - min_gap)
    if len(pos) and pos[0] < lo:
        pos += lo - pos[0]
    return pos


def pie(ax, values, labels, colors=None, highlight: int | None = None,
        fmt="{:.0f}%", donut: bool = False, startangle=90, inside_min=15):
    """Pie with direct labels (no legend), clockwise from 12 o'clock, largest first.

    Slices >= inside_min % are labelled inside; smaller ones outside, in a column
    left or right of the pie, connected by elbow leader lines. Outside labels are
    spread so independent labels keep >= TEXT_GAP (1 em) of clear space.
    Use greys and colour at most one highlighted slice. Keep to <= 5 slices;
    otherwise use a sorted horizontal bar chart instead (Tufte).
    """
    values = np.asarray(values, float)
    if colors is None:
        colors = list(GREYS[: len(values)])
        if highlight is not None:
            colors[highlight] = C["cyan"]
    pct = 100 * values / values.sum()
    width = 0.42 if donut else 1.0
    if donut:
        inside_min = 101          # ring is too thin for inside labels
    wedges, _ = ax.pie(values, colors=colors, startangle=startangle, counterclock=False,
                       radius=1.0, wedgeprops=dict(edgecolor="white", linewidth=1.2, width=width))
    ax.set_aspect("equal")
    ax.grid(False)
    fs = mpl.rcParams["font.size"]
    inside, outside = [], {1: [], -1: []}
    for w, lab, p, col in zip(wedges, labels, pct, colors):
        ang = np.deg2rad((w.theta2 + w.theta1) / 2)
        x, y = np.cos(ang), np.sin(ang)
        if p >= inside_min:
            inside.append((x, y, lab, p, col))
        else:
            # slices near 12 or 6 o'clock go to the right-hand column, keeping the
            # top-left corner free for the panel letter
            side = 1 if (x >= 0 or abs(x) < 0.35) else -1
            outside[side].append((y, x, f"{lab} {fmt.format(p)}"))
    ax.set(xlim=(-1.1, 1.1), ylim=(-1.25, 1.25))
    fig = ax.figure
    placed = []
    for (x, y, lab, p, col), w in zip(inside, [w for w, pp in zip(wedges, pct) if pp >= inside_min]):
        r = 1 - width / 2 if donut else 0.62
        lum = np.dot(mpl.colors.to_rgb(col), [0.299, 0.587, 0.114])
        t = ax.text(r * x, r * y, f"{lab}\n{fmt.format(p)}", ha="center", va="center",
                    fontsize=fs, color="white" if lum < 0.55 else C["ink"])
        placed.append((t, w, x, y, lab, p))
    for _ in range(4):                                  # inside label must fit its slice
        ax.set(xlim=(-2.1 if outside[-1] else -1.1, 2.1 if outside[1] else 1.1))
        fig.canvas.draw()                               # real scale after layout
        moved = False
        for item in list(placed):
            t, w, x, y, lab, p = item
            if not _inside_patch(t, w, fig.canvas.get_renderer(), margin=DATA_GAP * fs * fig.dpi / 72):
                t.remove()
                placed.remove(item)
                side = 1 if (x >= 0 or abs(x) < 0.35) else -1
                outside[side].append((y, x, f"{lab} {fmt.format(p)}"))
                moved = True
        if not moved:
            break
    pt_per_du = (ax.transData.transform((0, 1))[1] - ax.transData.transform((0, 0))[1]) * 72 / fig.dpi
    step = (1.2 + TEXT_GAP) * fs / pt_per_du            # line height + 1 em clear space
    texts = []
    for side, items in outside.items():
        if not items:
            continue
        items.sort()
        ys = _spread([1.12 * it[0] for it in items], step, -1.2, 1.2)
        for (y, x, text), yl in zip(items, ys):
            xl = side * 1.12
            ax.annotate(text, xy=(x, y), xytext=(xl, yl), ha="left" if side > 0 else "right",
                        va="center", fontsize=fs,
                        arrowprops=dict(arrowstyle="-", color=C["dark"], lw=0.6, shrinkA=3, shrinkB=0,
                                        connectionstyle="angle,angleA=0,angleB=90,rad=0"))
            texts.append(ax.texts[-1])
    if texts:                                           # widen the axes to hold every label
        fig.canvas.draw()
        inv = ax.transData.inverted()
        xs = [inv.transform(mpl.text.Text.get_window_extent(t))[:, 0] for t in texts]
        lo = min(-1.1, min(x.min() for x in xs) - 0.05)
        hi = max(1.1, max(x.max() for x in xs) + 0.05)
        ax.set_xlim(lo, hi)
    return wedges


# ---------------------------------------------------------------- checking ---
def _text_boxes(fig, renderer):
    """Visible text artists with their display bbox, font size (px) and group id."""
    out = []
    px = fig.dpi / 72
    for ax in fig.axes:
        def shown(axis):
            lo, hi = sorted(axis.get_view_interval())
            eps = (hi - lo) * 1e-6
            return [tk.label1 for tk in axis.get_major_ticks()
                    if tk.label1.get_visible() and lo - eps <= tk.get_loc() <= hi + eps]
        groups = [("xtick", shown(ax.xaxis)), ("ytick", shown(ax.yaxis))]
        for g, texts in groups:
            for t in texts:
                out.append((t, f"{id(ax)}{g}", "axis"))
        out.append((ax.xaxis.label, f"{id(ax)}xl", "axis"))
        out.append((ax.yaxis.label, f"{id(ax)}yl", "axis"))
        out.append((ax.title, f"{id(ax)}t", "axis"))
        for t in ax.texts:
            gid = t.get_gid() or ""
            kind = "panel" if gid == "panel" else "data"
            out.append((t, gid if gid.startswith("group") else id(t), kind))
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                out.append((t, f"leg{id(leg)}", "legend"))
    for t in fig.texts:
        out.append((t, id(t), "fig"))
    for leg in fig.legends:
        for t in leg.get_texts():
            out.append((t, f"leg{id(leg)}", "legend"))
    boxes = []
    for t, g, kind in out:
        if not t.get_visible() or not t.get_text().strip():
            continue
        try:   # the text itself plus its box, never the arrow of an annotation
            bb = mpl.text.Text.get_window_extent(t, renderer)
            if t.get_bbox_patch() is not None:
                bb = mpl.transforms.Bbox.union([bb, t.get_bbox_patch().get_window_extent(renderer)])
        except Exception:
            continue
        if bb.width < 1:
            continue
        boxes.append((t, bb, t.get_fontsize() * px, g, kind))
    return boxes


def _data_points(ax, step_px=1.5):
    """Display-coordinate points sampling every data mark of an axes."""
    pts = []
    for line in ax.get_lines():
        if not line.get_visible():
            continue
        xy = line.get_transform().transform(line.get_path().vertices)
        xy = xy[np.all(np.isfinite(xy), axis=1)]
        if line.get_linestyle() not in ("None", "none", "", " "):
            for a, b in zip(xy[:-1], xy[1:]):
                n = max(2, int(np.hypot(*(b - a)) / step_px))
                pts.append(np.linspace(a, b, n))
        if len(xy) == 1 or line.get_marker() not in ("None", None, "", " "):
            pts.append(xy)
    for coll in ax.collections:
        try:
            off = coll.get_offsets()
            if len(off):
                pts.append(coll.get_offset_transform().transform(off))
        except Exception:
            pass
    for patch in ax.patches:
        if isinstance(patch, mpl.patches.Rectangle) and patch.get_fill() and patch.get_visible():
            bb = patch.get_window_extent()
            gx, gy = np.meshgrid(np.arange(bb.x0, bb.x1 + 1, 3), [bb.y1])
            pts.append(np.column_stack([gx.ravel(), gy.ravel()]))  # top edge of bars
    return np.vstack(pts) if pts else np.empty((0, 2))


def check_layout(fig, text_gap=TEXT_GAP, group_gap=GROUP_GAP, data_gap=DATA_GAP, verbose=True):
    """Check the spacing rules; return a list of problems (empty list = clean).

    * text over data: annotation/legend text must keep >= data_gap em from any
      line, point or bar top (text must never occlude data);
    * text on text: independent labels keep >= text_gap em of clear space,
      lines within one group (legend, tick labels) >= group_gap em.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    boxes = _text_boxes(fig, r)
    problems = []
    # text cut off at the figure edge (the figure is saved at its exact size)
    W, H = fig.bbox.width, fig.bbox.height
    for t, b, f, g, k in boxes:
        if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W + 0.5 or b.y1 > H + 0.5:
            problems.append(f"text '{t.get_text()[:20]}' is cut off at the figure edge: "
                            "make the figure taller/wider or move the text in")
    # text vs text
    for i in range(len(boxes)):
        t1, b1, f1, g1, k1 = boxes[i]
        for j in range(i + 1, len(boxes)):
            t2, b2, f2, g2, k2 = boxes[j]
            same = g1 == g2
            if same and k1 == "axis":
                continue                      # neighbouring tick labels: matplotlib spaces them
            need = (group_gap if same else text_gap) * min(f1, f2)
            if k1 == "axis" and k2 == "axis":
                need = 0                      # axis furniture: only true overlaps matter
            elif "panel" in (k1, k2) and "axis" in (k1, k2):
                need = 0.5 * min(f1, f2)      # panel letter next to tick/axis labels
            dx = max(b2.x0 - b1.x1, b1.x0 - b2.x1, 0)
            dy = max(b2.y0 - b1.y1, b1.y0 - b2.y1, 0)
            gap = np.hypot(dx, dy) if (dx and dy) else max(dx, dy)
            overlap = dx == 0 and dy == 0
            if overlap or gap < need:
                problems.append(f"text '{t1.get_text()[:20]}' and '{t2.get_text()[:20]}': "
                                f"{'overlap' if overlap else f'gap {gap / min(f1, f2):.2f} em < {need / min(f1, f2):.2f} em'}")
    # text vs data
    for ax in fig.axes:
        pts = _data_points(ax)
        if not len(pts):
            continue
        for t, bb, f, g, kind in boxes:
            if kind not in ("data", "legend") or getattr(t, "axes", None) is not ax and kind == "data":
                continue
            if kind == "legend":
                leg = ax.get_legend()
                if leg is None or t not in leg.get_texts():
                    continue
                bb = leg.get_window_extent(r)
            m = data_gap * f
            inside = ((pts[:, 0] > bb.x0 - m) & (pts[:, 0] < bb.x1 + m) &
                      (pts[:, 1] > bb.y0 - m) & (pts[:, 1] < bb.y1 + m))
            if inside.any():
                problems.append(f"{'legend' if kind == 'legend' else 'text'} '{t.get_text()[:20]}' "
                                f"is within {data_gap} em of data ({inside.sum()} samples)")
    # text crossing the edge of a filled shape (pie slice, bar): partly in, partly out
    for ax in fig.axes:
        shapes = [p_ for p_ in ax.patches if p_.get_fill() and p_.get_visible()]
        for t, bb, f, g, kind in boxes:
            if kind != "data" or t.axes is not ax:
                continue
            per = _bbox_perimeter(bb)
            for sh in shapes:
                hits = [sh.contains_point(pt) for pt in per]
                if any(hits) and not all(hits):
                    problems.append(f"text '{t.get_text()[:20]}' crosses the edge of a shape")
                    break
    problems = list(dict.fromkeys(problems))
    if verbose:
        for p_ in problems:
            print("  [layout]", p_)
    return problems


# ------------------------------------------------------------------ output ---
def save(fig, name, outdir=".", formats=("pdf", "png", "svg"), dpi=600, check=True, **kw):
    """Save in vector (pdf/svg, editable text) and raster (png, 600 dpi) formats.

    check=True runs check_layout() first and prints any spacing problems; fix them
    before using the figure.
    """
    if check:
        probs = check_layout(fig, verbose=False)
        for p_ in probs:
            print(f"  [layout] {name}: {p_}")
    w_in, h_in = fig.get_size_inches()
    size = getattr(fig, "_pf_size", (None, w_in / MM, h_in / MM))
    tw = f" ({size[0]} of {TEXT_WIDTH_MM:.0f} mm text width)" if isinstance(size[0], str) else ""
    print(f"  [size] {name}: {w_in / MM:.1f} x {h_in / MM:.1f} mm{tw}")
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in formats:
        p = out / f"{name}.{ext}"
        fig.savefig(p, dpi=dpi, **kw)
        paths.append(p)
    return paths
