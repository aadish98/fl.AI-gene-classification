"""
Generate a high-level algorithm flowchart for fl.AI gene classification.
Output: fl.AI_algorithm_flowchart.png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# ── Canvas ───────────────────────────────────────────────────────────────────
FIG_W, FIG_H = 12, 26
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis("off")
fig.patch.set_facecolor("#F7F9FC")
ax.set_facecolor("#F7F9FC")

# ── Palette ──────────────────────────────────────────────────────────────────
NAVY   = "#0D2B55"
BLUE   = "#1A5FA8"
TEAL   = "#00796B"
PURPLE = "#6A1B9A"
GREEN  = "#2E7D32"
ORANGE = "#D95B00"
LGRAY  = "#DDE3EA"
MGRAY  = "#7A90A0"
WHITE  = "#FFFFFF"

# ── Layout constants ──────────────────────────────────────────────────────────
CX  = 6.0    # horizontal centre
BW  = 6.8    # standard box width
BH  = 0.70   # standard box height
SBW = 2.75   # search box width
SBH = 0.85   # search box height (2-line)
RADIUS = 0.14

# Centre-Y of every node (high → low)
Y_INPUT      = 24.0
Y_GCONV      = 22.65
Y_SEARCH     = 20.9
Y_MERGE      = 19.2
Y_RANK       = 17.95
Y_DIAMOND    = 16.45
Y_FULLTEXT   = 14.9
Y_FALLBACK   = 14.05   # thin annotation row
Y_CHUNK      = 13.2
Y_AIEXTRACT  = 11.65
Y_MERGESUM   = 10.15
Y_AGG        =  8.95
Y_CLASSIFY   =  7.35
Y_OUTPUT     =  5.75
Y_PILLS      =  4.8
Y_LEGEND     =  1.4
Y_FOOTER     =  0.5

# ── Helper functions ──────────────────────────────────────────────────────────

def rbox(cx, cy, w, h, fill, text,
         fc=WHITE, fs=11, bold=False, alpha=1.0,
         sub=None, sub_fc=None):
    """Rounded rectangle centred at (cx, cy)."""
    patch = FancyBboxPatch(
        (cx - w/2, cy - h/2), w, h,
        boxstyle=f"round,pad=0,rounding_size={RADIUS}",
        linewidth=0, facecolor=fill, alpha=alpha, zorder=3,
    )
    ax.add_patch(patch)
    ty = cy if sub is None else cy + h * 0.16
    ax.text(cx, ty, text, ha="center", va="center",
            fontsize=fs, fontweight="bold" if bold else "normal",
            color=fc, zorder=4, multialignment="center")
    if sub:
        ax.text(cx, cy - h * 0.22, sub,
                ha="center", va="center",
                fontsize=fs - 1.5, color=sub_fc or fc,
                alpha=0.85, zorder=4, multialignment="center")


def dmd(cx, cy, w, h, fill, text, fc=WHITE, fs=10):
    """Decision diamond."""
    pts = [[cx, cy+h/2], [cx+w/2, cy], [cx, cy-h/2], [cx-w/2, cy]]
    ax.add_patch(plt.Polygon(pts, closed=True, facecolor=fill,
                             edgecolor="none", zorder=3))
    ax.text(cx, cy, text, ha="center", va="center",
            fontsize=fs, fontweight="bold", color=fc, zorder=4,
            multialignment="center")


def arr(x1, y1, x2, y2, col=MGRAY, lw=1.8, hw=0.16):
    """Arrow from (x1,y1) to (x2,y2)."""
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(
                    arrowstyle=f"-|>,head_width={hw},head_length={hw*0.65}",
                    color=col, lw=lw),
                zorder=5)


def pill(cx, cy, w, h, fill, text, fs=8.5):
    """Small rounded pill."""
    patch = FancyBboxPatch(
        (cx - w/2, cy - h/2), w, h,
        boxstyle=f"round,pad=0,rounding_size=0.10",
        linewidth=0, facecolor=fill, alpha=0.92, zorder=4,
    )
    ax.add_patch(patch)
    ax.text(cx, cy, text, ha="center", va="center",
            fontsize=fs, fontweight="bold", color=WHITE, zorder=5,
            multialignment="center")


def phase_bg(y_top, y_bot, col, label):
    """Lightly shaded phase band."""
    rect = FancyBboxPatch(
        (0.3, y_bot), FIG_W - 0.6, y_top - y_bot,
        boxstyle="round,pad=0,rounding_size=0.2",
        linewidth=1.0, edgecolor=col,
        facecolor=col, alpha=0.06, zorder=1,
    )
    ax.add_patch(rect)
    # Label on left outside
    ax.text(0.08, (y_top + y_bot) / 2, label,
            ha="center", va="center", fontsize=7.5,
            fontweight="bold", color=col, alpha=0.75,
            rotation=90, zorder=2)


# ── Phase bands ───────────────────────────────────────────────────────────────
phase_bg(25.1, 21.9,  BLUE,   "INPUT")
phase_bg(21.9, 16.9,  BLUE,   "LITERATURE\nCOLLECTION")
phase_bg(16.9, 12.75, ORANGE, "FILTER &\nRETRIEVAL")
phase_bg(12.75, 8.35, PURPLE, "AI EXTRACTION")
phase_bg( 8.35, 3.9,  GREEN,  "CLASSIFY &\nOUTPUT")

# ── Title ─────────────────────────────────────────────────────────────────────
ax.text(CX, 25.55, "fl.AI Gene Classification — Core Algorithm",
        ha="center", va="center",
        fontsize=16, fontweight="bold", color=NAVY, zorder=6)
ax.text(CX, 25.12, "High-Level Pipeline Overview  ·  Drosophila melanogaster",
        ha="center", va="center", fontsize=10, color=MGRAY, zorder=6)
ax.plot([1.2, FIG_W-1.2], [24.85, 24.85], color=TEAL, lw=2, zorder=6)

# ── 1. Input ──────────────────────────────────────────────────────────────────
rbox(CX, Y_INPUT, BW, BH, NAVY, "INPUT:  Gene List  (CSV)",
     fs=13, bold=True)
ax.text(CX + BW/2 + 0.22, Y_INPUT + 0.05,
        "One row per gene\ne.g.  period, timeless, sifamide",
        ha="left", va="center", fontsize=8, color=MGRAY, zorder=6)

# ── 2. Gene conversion ────────────────────────────────────────────────────────
arr(CX, Y_INPUT - BH/2, CX, Y_GCONV + BH/2, col=BLUE)
rbox(CX, Y_GCONV, BW, BH, BLUE,
     "Convert gene symbols  →  FlyBase IDs  (FBgn)",
     sub="Using FlyBase synonym tables  ·  handles Greek letters & aliases",
     sub_fc="#B3D4F5", fs=11)

# ── 3. Three parallel searches ────────────────────────────────────────────────
search_xs = [CX - 3.7, CX, CX + 3.7]
search_labels = ["FlyBase\nReferences", "PubMed\nSearch", "Europe PMC\nSearch"]
search_subs   = ["Curated gene→paper\nmappings", "Gene name variants\n+ keyword queries",
                 "International OA\ncoverage"]

# fan-out arrows
for sx in search_xs:
    arr(CX, Y_GCONV - BH/2, sx, Y_SEARCH + SBH/2, col=BLUE, lw=1.4, hw=0.13)

for sx, lbl, sub in zip(search_xs, search_labels, search_subs):
    rbox(sx, Y_SEARCH, SBW, SBH, BLUE, lbl,
         sub=sub, sub_fc="#B3D4F5", fs=10.5)

# fan-in arrows
for sx in search_xs:
    arr(sx, Y_SEARCH - SBH/2, CX, Y_MERGE + BH/2, col=BLUE, lw=1.4, hw=0.13)

# ── 4. Merge ──────────────────────────────────────────────────────────────────
rbox(CX, Y_MERGE, BW, BH, BLUE,
     "Merge & Deduplicate References",
     sub="Combine all PMCIDs  ·  remove duplicates across sources",
     sub_fc="#B3D4F5", fs=11)

# ── 5. Rank ───────────────────────────────────────────────────────────────────
arr(CX, Y_MERGE - BH/2, CX, Y_RANK + BH/2, col=BLUE)
rbox(CX, Y_RANK, BW, BH, BLUE,
     "Rank by Source Consensus  +  Recency",
     sub="Papers found by multiple databases rank higher  ·  capped at user limit",
     sub_fc="#B3D4F5", fs=10.5)

# ── 6. Keyword filter diamond ─────────────────────────────────────────────────
DW, DH = 4.4, 0.95
arr(CX, Y_RANK - BH/2, CX, Y_DIAMOND + DH/2, col=ORANGE)
dmd(CX, Y_DIAMOND, DW, DH, ORANGE,
    "Does title / abstract\ncontain keywords?", fs=10.5)

# "No" → Skip paper  (right branch)
SKIP_X = CX + 4.6
SKIP_W = 2.2
arr(CX + DW/2, Y_DIAMOND, SKIP_X - SKIP_W/2, Y_DIAMOND,
    col=MGRAY, lw=1.4, hw=0.13)
ax.text((CX + DW/2 + SKIP_X - SKIP_W/2) / 2, Y_DIAMOND + 0.22,
        "No", ha="center", va="bottom", fontsize=9.5,
        fontweight="bold", color=MGRAY)
rbox(SKIP_X, Y_DIAMOND, SKIP_W, BH * 0.85, MGRAY, "Skip paper", fs=10)

# "Yes" label below diamond
ax.text(CX + 0.18, Y_DIAMOND - DH/2 - 0.18,
        "Yes", ha="left", va="top", fontsize=9.5,
        fontweight="bold", color=ORANGE)

# ── 7. Fetch full text ────────────────────────────────────────────────────────
arr(CX, Y_DIAMOND - DH/2, CX, Y_FULLTEXT + BH/2, col=ORANGE)
rbox(CX, Y_FULLTEXT, BW, BH, ORANGE,
     "Fetch Full Text  (per matching paper)",
     sub="Attempted in order: PMC OA PDF  →  PMC XML  →  Europe PMC API  →  Unpaywall  →  Abstract fallback",
     sub_fc="#FFD4A8", fs=10.5)

# ── 8. Chunk text ─────────────────────────────────────────────────────────────
arr(CX, Y_FULLTEXT - BH/2, CX, Y_CHUNK + BH/2, col=ORANGE)
rbox(CX, Y_CHUNK, BW, BH, ORANGE,
     "Split paper text into chunks  (≤ 5 000 words each)",
     sub="Each chunk fits within the AI model's context window",
     sub_fc="#FFD4A8", fs=10.5)

# ── 9. AI evidence extraction ────────────────────────────────────────────────
AI_BH = 1.15
arr(CX, Y_CHUNK - BH/2, CX, Y_AIEXTRACT + AI_BH/2, col=PURPLE)
rbox(CX, Y_AIEXTRACT, BW, AI_BH, PURPLE,
     "AI Evidence Extraction  (GPT — per chunk)",
     fs=12, bold=True)
items = [
    "Function:   What does the gene do biologically?",
    "Phenotype:  What happens when the gene is mutated or knocked down?",
    "Reagents:   Which fly stocks / RNAi lines / tools were used?",
]
for i, t in enumerate(items):
    ax.text(CX - BW/2 + 0.45, Y_AIEXTRACT + 0.28 - i * 0.34,
            f"•  {t}", ha="left", va="center",
            fontsize=9, color="#DDB3FF", zorder=6)

# ── 10. Merge chunk summaries ─────────────────────────────────────────────────
arr(CX, Y_AIEXTRACT - AI_BH/2, CX, Y_MERGESUM + BH/2, col=PURPLE)
rbox(CX, Y_MERGESUM, BW, BH, PURPLE,
     "Merge chunk summaries  →  Paper-level evidence summary",
     sub="One coherent function + phenotype summary per reviewed paper",
     sub_fc="#DDB3FF", fs=10.5)

# ── 11. Aggregate per gene ────────────────────────────────────────────────────
arr(CX, Y_MERGESUM - BH/2, CX, Y_AGG + BH/2, col=PURPLE)
rbox(CX, Y_AGG, BW, BH, PURPLE,
     "Aggregate all paper summaries  →  Gene evidence block",
     sub="All reviewed papers combined into one evidence block per gene",
     sub_fc="#DDB3FF", fs=10.5)

# ── 12. AI classification ─────────────────────────────────────────────────────
CL_BH = 1.15
arr(CX, Y_AGG - BH/2, CX, Y_CLASSIFY + CL_BH/2, col=GREEN)
rbox(CX, Y_CLASSIFY, BW, CL_BH, GREEN,
     "AI Gene Classification  (GPT)",
     fs=12, bold=True)
cl_items = [
    ("Category:  ", "e.g.  sleep  ·  circadian  ·  synapse  (zero or more)"),
    ("Confidence:", "0–100 score based on strength & consistency of evidence"),
    ("Rationale:  ", "Plain-language explanation citing specific papers"),
]
for i, (lbl, val) in enumerate(cl_items):
    y = Y_CLASSIFY + 0.28 - i * 0.35
    ax.text(CX - BW/2 + 0.45, y, f"•  {lbl}",
            ha="left", va="center", fontsize=9,
            fontweight="bold", color="#ADEFB0", zorder=6)
    ax.text(CX - BW/2 + 1.55, y, val,
            ha="left", va="center", fontsize=9,
            color="#ADEFB0", zorder=6)

# ── 13. Output ────────────────────────────────────────────────────────────────
arr(CX, Y_CLASSIFY - CL_BH/2, CX, Y_OUTPUT + BH/2, col=GREEN)
rbox(CX, Y_OUTPUT, BW, BH, NAVY,
     "OUTPUT:  Excel Workbook  (.xlsx)", fs=13, bold=True)

# Sheet pills
sheet_data = [
    ("Gene Set",        NAVY),
    ("Classification",  GREEN),
    ("Ref. Summaries",  TEAL),
    ("Reagents",        BLUE),
]
pill_w, pill_h = 2.35, 0.46
pill_xs = [CX - 3.6, CX - 1.2, CX + 1.2, CX + 3.6]
ax.text(CX, Y_PILLS + 0.42, "4 output sheets:",
        ha="center", va="center", fontsize=8.5, color=MGRAY)
for px, (lbl, col) in zip(pill_xs, sheet_data):
    pill(px, Y_PILLS - 0.05, pill_w, pill_h, col, lbl, fs=8.5)

# ── Legend ────────────────────────────────────────────────────────────────────
leg_items = [
    (BLUE,   "Data collection"),
    (ORANGE, "Filtering & retrieval"),
    (PURPLE, "AI analysis"),
    (GREEN,  "Classify & output"),
]
ax.text(0.4, Y_LEGEND + 0.45, "Phase colour key:",
        ha="left", va="center", fontsize=8.5,
        fontweight="bold", color=MGRAY)
leg_xs = [1.5, 4.2, 6.9, 9.6]
for lx, (col, lbl) in zip(leg_xs, leg_items):
    pill(lx, Y_LEGEND, 2.4, 0.42, col, lbl, fs=8.5)

# ── Footer ────────────────────────────────────────────────────────────────────
ax.plot([0.8, FIG_W - 0.8], [Y_FOOTER + 0.32, Y_FOOTER + 0.32],
        color=LGRAY, lw=1, zorder=6)
ax.text(CX, Y_FOOTER,
        "fl.AI-gene-classification  ·  Allada Lab, University of Michigan  ·  2026",
        ha="center", va="center", fontsize=8, color=MGRAY)

# ── Save ──────────────────────────────────────────────────────────────────────
OUT = "/Users/aadishms/Desktop/Projects.nosync/fl.AI-gene-classification/fl.AI_algorithm_flowchart.png"
plt.tight_layout(pad=0)
plt.savefig(OUT, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
print(f"Saved: {OUT}")
