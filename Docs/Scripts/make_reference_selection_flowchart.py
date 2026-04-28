"""
Generate a high-level organism-specific reference selection flowchart.
Output: reference_selection_flowchart.png

Aesthetic direction: editorial scientific infographic.
- Warm off-white paper, dark ink typography, restrained accents.
- Three vertical lanes (Fly / Human / Mouse) for candidate collection.
- Single horizontal pipeline below where the lanes converge.
- Source databases shown as discrete pill chips beneath each lane,
  so multi-source lists never overflow a box.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle


# ── Canvas ───────────────────────────────────────────────────────────────────
FIG_W, FIG_H = 18.5, 13.0
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis("off")

PAPER = "#F4F1EA"
fig.patch.set_facecolor(PAPER)
ax.set_facecolor(PAPER)


# ── Palette ──────────────────────────────────────────────────────────────────
INK = "#0F1B2D"
INK_SOFT = "#3A4A63"
RULE = "#1A2A40"
MUTED = "#7A8499"
HAIRLINE = "#C9C2B2"

# Organism accents — chosen to read as a curated trio, not a rainbow.
FLY = "#1B4F72"     # deep teal navy
HUMAN = "#8C2F3F"   # oxblood
MOUSE = "#B5651D"   # warm ochre

# Pipeline accent
PIPE = "#2F5D3A"    # forest green

CARD = "#FFFFFF"


# ── Typography helpers ───────────────────────────────────────────────────────
TITLE_FONT = {"family": "serif", "weight": "bold"}
EYEBROW_FONT = {"family": "sans-serif", "weight": "bold"}
BODY_FONT = {"family": "sans-serif"}


def text(x, y, s, *, size=10, color=INK, weight="normal", family="sans-serif",
         ha="center", va="center", italic=False, alpha=1.0, zorder=6):
    style = "italic" if italic else "normal"
    ax.text(
        x, y, s,
        ha=ha, va=va,
        fontsize=size, color=color, weight=weight,
        family=family, style=style, alpha=alpha, zorder=zorder,
    )


# ── Shape primitives ─────────────────────────────────────────────────────────
def card(cx, cy, w, h, *, edge=INK, fill=CARD, lw=1.4, radius=0.10, zorder=3):
    patch = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={radius}",
        linewidth=lw, edgecolor=edge, facecolor=fill, zorder=zorder,
    )
    ax.add_patch(patch)


def solid_card(cx, cy, w, h, *, color, lw=0, radius=0.10, zorder=3):
    patch = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={radius}",
        linewidth=lw, edgecolor=color, facecolor=color, zorder=zorder,
    )
    ax.add_patch(patch)


def lane_band(x, y, w, h, color):
    """Soft tinted vertical band that contains a lane."""
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.18",
        linewidth=0, facecolor=color, alpha=0.07, zorder=1,
    )
    ax.add_patch(patch)


def pill(cx, cy, label, *, color):
    """Source-database chip rendered as a thin outlined pill."""
    pad_x = 0.18
    pad_y = 0.16
    # Estimate chip width from label length so chips are tightly fitted.
    char_w = 0.085
    w = max(1.20, len(label) * char_w + pad_x * 2)
    h = 0.42
    patch = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.01,rounding_size=0.18",
        linewidth=1.0, edgecolor=color, facecolor=CARD, zorder=4,
    )
    ax.add_patch(patch)
    text(cx, cy, label, size=8.6, color=color, weight="bold", zorder=5)
    return w


def stacked_pills(cx, cy_top, items, *, color, gap=0.46):
    """Render a column of pill chips, top-down, centered on cx."""
    cy = cy_top
    for label in items:
        pill(cx, cy, label, color=color)
        cy -= gap
    return cy + gap  # y of last chip


def arrow(x1, y1, x2, y2, *, color=INK_SOFT, lw=1.6, head=10, alpha=1.0,
          shrink=4, zorder=5):
    ax.annotate(
        "",
        xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle=f"-|>,head_width=0.28,head_length=0.40",
            color=color, lw=lw, alpha=alpha,
            shrinkA=shrink, shrinkB=shrink,
        ),
        zorder=zorder,
    )


def hairline(x1, y, x2, color=HAIRLINE, lw=0.8):
    ax.plot([x1, x2], [y, y], color=color, lw=lw, zorder=2)


# ── Layout constants ─────────────────────────────────────────────────────────
CX = FIG_W / 2

# Vertical anchors (top → bottom)
Y_EYEBROW   = 12.45
Y_TITLE     = 11.95
Y_SUBTITLE  = 11.30
Y_RULE      = 10.95

Y_INPUT     = 10.20

# Lane region
LANE_TOP    = 9.65
LANE_BOTTOM = 4.30
LANE_HEADER = 9.32
Y_IDENTITY  = 8.55
Y_COLLECT   = 7.60   # header for the source list
Y_PILLS_TOP = 7.05   # first pill top center
Y_TAG       = 4.65

# Convergence + pipeline
Y_MERGE     = 3.55
Y_PIPELINE  = 2.20
Y_OUTPUT    = 0.70

# Lane horizontal positions
LANE_W = 5.55
GAP    = 0.45
LANE_TOTAL = LANE_W * 3 + GAP * 2
LANE_LEFT_X0 = (FIG_W - LANE_TOTAL) / 2
lane_centers = [
    LANE_LEFT_X0 + LANE_W / 2,
    LANE_LEFT_X0 + LANE_W + GAP + LANE_W / 2,
    LANE_LEFT_X0 + LANE_W * 2 + GAP * 2 + LANE_W / 2,
]


# ── Header ───────────────────────────────────────────────────────────────────
text(CX, Y_EYEBROW, "fl.AI  /  Pipeline Overview",
     size=10.5, color=MUTED, weight="bold",
     family="sans-serif")
text(CX, Y_TITLE, "Reference Selection Across Organisms",
     size=26, color=INK, weight="bold", family="serif")
text(
    CX, Y_SUBTITLE,
    "How candidate literature is gathered, filtered, and carried forward "
    "for fly, human-ortholog, and mouse-ortholog gene sets",
    size=11.2, color=INK_SOFT, family="sans-serif", italic=True,
)

# Horizontal rule under the title.
ax.plot([1.6, FIG_W - 1.6], [Y_RULE, Y_RULE], color=RULE, lw=1.6, zorder=2)
# Tiny center tick to give the rule character.
ax.plot([CX - 0.18, CX + 0.18], [Y_RULE, Y_RULE], color=PAPER, lw=4, zorder=3)
ax.plot([CX], [Y_RULE], marker="o", color=RULE, markersize=4.5, zorder=4)


# ── Input node ───────────────────────────────────────────────────────────────
INPUT_W, INPUT_H = 5.4, 0.78
solid_card(CX, Y_INPUT, INPUT_W, INPUT_H, color=INK, radius=0.18)
text(CX, Y_INPUT + 0.12, "INPUT  ·  Fly gene set",
     size=12.5, color=PAPER, weight="bold", family="serif")
text(CX, Y_INPUT - 0.18, "Symbols normalized before reference selection",
     size=9.6, color="#D7D9DE", family="sans-serif")


# ── Lane bands ───────────────────────────────────────────────────────────────
lane_colors = [FLY, HUMAN, MOUSE]
lane_labels = ["Fly", "Human ortholog", "Mouse ortholog"]
lane_left_xs = [c - LANE_W / 2 for c in lane_centers]

for x_left, color in zip(lane_left_xs, lane_colors):
    lane_band(x_left, LANE_BOTTOM, LANE_W, LANE_TOP - LANE_BOTTOM, color)


# ── Lane headers (eyebrow + title) ───────────────────────────────────────────
for cx_lane, color, label in zip(lane_centers, lane_colors, lane_labels):
    # Small vertical accent bar on the left of the lane header.
    bar = Rectangle(
        (cx_lane - LANE_W / 2 + 0.30, LANE_HEADER - 0.18),
        0.10, 0.36, facecolor=color, edgecolor="none", zorder=4,
    )
    ax.add_patch(bar)
    text(cx_lane - LANE_W / 2 + 0.55, LANE_HEADER + 0.08,
         f"LANE  ·  {label.upper()}",
         size=9.0, color=color, weight="bold", ha="left",
         family="sans-serif")
    text(cx_lane - LANE_W / 2 + 0.55, LANE_HEADER - 0.18,
         "Candidate collection", size=9.0, color=MUTED, ha="left",
         family="sans-serif", italic=True)


# ── Identity step (top of each lane) ─────────────────────────────────────────
identity_subtitles = {
    "Fly": "No ortholog step.\nReference search uses fly identity.",
    "Human ortholog": "Each fly gene mapped to its human ortholog.\nFly to human link preserved.",
    "Mouse ortholog": "Each fly gene mapped to its mouse ortholog.\nFly to mouse link preserved.",
}
ID_W, ID_H = LANE_W - 0.80, 1.20
for cx_lane, color, label in zip(lane_centers, lane_colors, lane_labels):
    card(cx_lane, Y_IDENTITY, ID_W, ID_H, edge=color, lw=1.4)
    if label == "Fly":
        title_line = "Use fly gene identity"
    else:
        title_line = f"Map fly  →  {label.split()[0].lower()}"
    text(cx_lane, Y_IDENTITY + 0.32, title_line,
         size=12.0, color=color, weight="bold", family="serif")
    text(cx_lane, Y_IDENTITY - 0.22, identity_subtitles[label],
         size=9.2, color=INK_SOFT, family="sans-serif")


# ── "Collect candidate references" header for each lane ──────────────────────
for cx_lane, color in zip(lane_centers, lane_colors):
    text(cx_lane, Y_COLLECT + 0.10, "Collect candidate references",
         size=10.6, color=INK, weight="bold", family="serif")
    text(cx_lane, Y_COLLECT - 0.18,
         "Sources queried in parallel:",
         size=8.8, color=MUTED, family="sans-serif", italic=True)


# ── Source pill chips per lane ───────────────────────────────────────────────
fly_sources = [
    "Curated fly references",
    "PubMed search",
    "Europe PMC search",
]
mammal_sources = [
    "Gene-linked PubMed records",
    "GeneRIF notes",
    "UniProt literature",
    "PubMed search",
    "Europe PMC search",
]
lane_sources = [fly_sources, mammal_sources, mammal_sources]

last_pill_y = []
for cx_lane, color, sources in zip(lane_centers, lane_colors, lane_sources):
    last_y = stacked_pills(cx_lane, Y_PILLS_TOP, sources, color=color)
    last_pill_y.append(last_y)


# ── "Tag references with their source" footer per lane ───────────────────────
TAG_W, TAG_H = LANE_W - 0.70, 0.65
for cx_lane, color in zip(lane_centers, lane_colors):
    card(cx_lane, Y_TAG, TAG_W, TAG_H, edge=color, lw=1.2,
         fill=CARD, radius=0.14)
    text(cx_lane, Y_TAG + 0.08, "Tag each reference with its source",
         size=9.6, color=color, weight="bold", family="sans-serif")
    text(cx_lane, Y_TAG - 0.16,
         "Source labels carry through downstream ranking",
         size=8.4, color=MUTED, family="sans-serif", italic=True)


# ── Vertical arrows within lanes (only between aligned lane elements) ────────
for cx_lane, color in zip(lane_centers, lane_colors):
    # Identity card → "Collect" header
    arrow(cx_lane, Y_IDENTITY - ID_H / 2,
          cx_lane, Y_COLLECT + 0.30,
          color=color, lw=1.3, head=8, alpha=0.85)
    # Last pill → Tag footer
    arrow(cx_lane, last_pill_y[lane_centers.index(cx_lane)] - 0.26,
          cx_lane, Y_TAG + TAG_H / 2,
          color=color, lw=1.3, head=8, alpha=0.85)


# ── Input → each lane band top (single fan-out, lands above lane header) ────
for cx_lane, color in zip(lane_centers, lane_colors):
    arrow(CX, Y_INPUT - INPUT_H / 2,
          cx_lane, LANE_TOP - 0.04,
          color=color, lw=1.4, head=10, alpha=0.9)


# ── Convergence row: "Merge duplicate papers per organism" ───────────────────
MERGE_W, MERGE_H = 13.0, 0.95
solid_card(CX, Y_MERGE, MERGE_W, MERGE_H, color=PIPE, radius=0.16)
text(CX, Y_MERGE + 0.16, "Merge duplicate papers within each organism",
     size=13.0, color=PAPER, weight="bold", family="serif")
text(CX, Y_MERGE - 0.20,
     "A paper found by multiple sources is kept once, with source labels combined",
     size=10.0, color="#D9E5DC", family="sans-serif", italic=True)


# ── Lane → merge convergence arrows (no crossing other lanes) ────────────────
# Each lane's tag box drops straight down to a y above the merge bar, then
# routes diagonally to the merge bar's top edge in its own column.
for cx_lane, color in zip(lane_centers, lane_colors):
    # Drop arrow head lands on the merge bar's top edge, in-column.
    arrow(cx_lane, Y_TAG - TAG_H / 2,
          cx_lane, Y_MERGE + MERGE_H / 2,
          color=color, lw=1.5, head=10, alpha=0.95)


# ── Horizontal pipeline below the merge bar ──────────────────────────────────
pipeline_steps = [
    ("Rank candidates", "Source agreement, then recency"),
    ("Apply keyword filter", "Match requested biology terms"),
    ("Retrieve evidence text", "Full text or title and abstract"),
    ("Summarize and classify", "Carry forward best evidence"),
]
N = len(pipeline_steps)
PIPE_W = 4.00
PIPE_H = 1.05
PIPE_GAP = 0.40
total_pipe_w = PIPE_W * N + PIPE_GAP * (N - 1)
pipe_x0 = (FIG_W - total_pipe_w) / 2
pipe_centers = [pipe_x0 + PIPE_W / 2 + i * (PIPE_W + PIPE_GAP) for i in range(N)]

# Center connector from merge bar down to the pipeline row.
arrow(CX, Y_MERGE - MERGE_H / 2,
      CX, Y_PIPELINE + PIPE_H / 2 + 0.05,
      color=PIPE, lw=1.6, head=11)

for i, (cx_p, (title, subtitle)) in enumerate(zip(pipe_centers, pipeline_steps)):
    card(cx_p, Y_PIPELINE, PIPE_W, PIPE_H, edge=PIPE, lw=1.3, radius=0.12)
    # Step number label.
    text(cx_p - PIPE_W / 2 + 0.34, Y_PIPELINE + PIPE_H / 2 - 0.22,
         f"0{i + 1}", size=9.2, color=PIPE, weight="bold",
         ha="left", family="sans-serif")
    text(cx_p, Y_PIPELINE + 0.18, title,
         size=11.4, color=INK, weight="bold", family="serif")
    text(cx_p, Y_PIPELINE - 0.22, subtitle,
         size=9.0, color=INK_SOFT, family="sans-serif")
    if i > 0:
        # Horizontal arrow between adjacent pipeline cards.
        arrow(pipe_centers[i - 1] + PIPE_W / 2,
              Y_PIPELINE,
              cx_p - PIPE_W / 2,
              Y_PIPELINE,
              color=PIPE, lw=1.4, head=10)


# ── Final output node ────────────────────────────────────────────────────────
OUT_W, OUT_H = 13.0, 0.95
solid_card(CX, Y_OUTPUT, OUT_W, OUT_H, color=INK, radius=0.18)
text(CX, Y_OUTPUT + 0.16, "OUTPUT  ·  Per-organism classification workbook",
     size=12.8, color=PAPER, weight="bold", family="serif")
text(CX, Y_OUTPUT - 0.20,
     "Selected references with source attribution, summaries, and gene-level classifications",
     size=10.0, color="#D7D9DE", family="sans-serif", italic=True)

# Pipeline → output connector (single straight arrow from the row's center).
arrow(CX, Y_PIPELINE - PIPE_H / 2,
      CX, Y_OUTPUT + OUT_H / 2,
      color=PIPE, lw=1.6, head=11)


# ── Save ─────────────────────────────────────────────────────────────────────
out_path = "reference_selection_flowchart.png"
fig.savefig(out_path, dpi=240, bbox_inches="tight",
            facecolor=fig.get_facecolor())
print(f"Saved {out_path}")
