"""
Generate a strategy deck for validating the fl.AI gene classification algorithm.

Run: python3 Docs/Scripts/make_validation_deck.py
Output: Docs/Decks/fl.AI_Algorithm_Validation_Strategy.pptx
"""

from __future__ import annotations

from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Emu, Inches, Pt


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "Docs" / "Decks" / "fl.AI_Algorithm_Validation_Strategy.pptx"

SLIDE_W = 13.333
SLIDE_H = 7.5

# Topic-specific palette: lab notebook cream, evidence ink, assay teal, alert coral.
INK = RGBColor(0x0D, 0x1B, 0x2A)
INK_2 = RGBColor(0x17, 0x2A, 0x3A)
CREAM = RGBColor(0xF7, 0xF1, 0xE4)
SAND = RGBColor(0xEC, 0xDE, 0xC7)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
TEAL = RGBColor(0x0B, 0x8F, 0x8C)
TEAL_DARK = RGBColor(0x06, 0x66, 0x64)
CORAL = RGBColor(0xE8, 0x5D, 0x4A)
AMBER = RGBColor(0xD7, 0xA6, 0x27)
MOSS = RGBColor(0x69, 0x8F, 0x3A)
SLATE = RGBColor(0x31, 0x3C, 0x4A)
MUTED = RGBColor(0x6D, 0x78, 0x86)
PALE_TEAL = RGBColor(0xD9, 0xEF, 0xEC)
PALE_CORAL = RGBColor(0xF7, 0xD8, 0xD1)

HEADER_FONT = "Georgia"
BODY_FONT = "Calibri"


def add_rect(slide, x, y, w, h, fill, line=None, width=Pt(0.5), shape=MSO_SHAPE.RECTANGLE):
    obj = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    obj.fill.solid()
    obj.fill.fore_color.rgb = fill
    if line is None:
        obj.line.fill.background()
    else:
        obj.line.color.rgb = line
        obj.line.width = width
    obj.shadow.inherit = False
    return obj


def add_round(slide, x, y, w, h, fill, line=None, radius=0.12):
    obj = add_rect(slide, x, y, w, h, fill, line=line, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    obj.adjustments[0] = radius
    return obj


def add_oval(slide, x, y, w, h, fill, line=None, width=Pt(0.5)):
    return add_rect(slide, x, y, w, h, fill, line=line, width=width, shape=MSO_SHAPE.OVAL)


def add_text(
    slide,
    text,
    x,
    y,
    w,
    h,
    *,
    size=14,
    bold=False,
    italic=False,
    color=SLATE,
    font=BODY_FONT,
    align=PP_ALIGN.LEFT,
    valign=MSO_ANCHOR.TOP,
    char_spacing=None,
):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Emu(0)
    tf.margin_top = tf.margin_bottom = Emu(0)
    tf.vertical_anchor = valign
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.name = font
    run.font.color.rgb = color
    if char_spacing is not None:
        run._r.get_or_add_rPr().set("spc", str(char_spacing))
    return box


def add_paragraphs(slide, items, x, y, w, h, *, size=13, color=SLATE, gap=5, font=BODY_FONT):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Emu(0)
    tf.margin_top = tf.margin_bottom = Emu(0)
    for i, item in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(gap)
        p.alignment = PP_ALIGN.LEFT
        run = p.add_run()
        run.text = item
        run.font.size = Pt(size)
        run.font.name = font
        run.font.color.rgb = color
    return box


def add_header(slide, title, eyebrow, slide_no, total):
    add_text(slide, eyebrow.upper(), 0.65, 0.46, 5.2, 0.25, size=9, bold=True, color=TEAL_DARK, char_spacing=250)
    add_text(slide, title, 0.62, 0.78, 8.8, 0.55, size=28, bold=True, color=INK, font=HEADER_FONT)
    add_text(slide, f"{slide_no:02d} / {total:02d}", 12.0, 0.55, 0.72, 0.2, size=8, color=MUTED, align=PP_ALIGN.RIGHT)


def add_footer(slide):
    add_text(slide, "fl.AI validation strategy", 0.65, 7.08, 3.2, 0.22, size=8, color=MUTED, char_spacing=120)


def add_chip(slide, text, x, y, w, h, *, fill=PALE_TEAL, color=INK, size=10):
    add_round(slide, x, y, w, h, fill, radius=0.45)
    add_text(slide, text, x + 0.06, y + 0.06, w - 0.12, h - 0.08, size=size, bold=True, color=color, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE)


def add_card(slide, x, y, w, h, title, body, *, accent=TEAL, fill=WHITE, title_size=15, body_size=11):
    add_round(slide, x, y, w, h, fill, line=SAND, radius=0.08)
    add_rect(slide, x, y, 0.08, h, accent)
    add_text(slide, title, x + 0.22, y + 0.18, w - 0.38, 0.28, size=title_size, bold=True, color=INK)
    if isinstance(body, str):
        add_text(slide, body, x + 0.22, y + 0.6, w - 0.38, h - 0.74, size=body_size, color=SLATE)
    else:
        add_paragraphs(slide, body, x + 0.22, y + 0.58, w - 0.38, h - 0.68, size=body_size, color=SLATE, gap=3)


def add_metric(
    slide,
    x,
    y,
    w,
    h,
    value,
    label,
    *,
    fill=INK,
    accent=TEAL,
    value_size=31,
    value_color=WHITE,
    label_color=WHITE,
):
    add_round(slide, x, y, w, h, fill, radius=0.08)
    add_oval(slide, x + w - 0.58, y + 0.18, 0.34, 0.34, accent)
    add_text(slide, value, x + 0.18, y + 0.15, w - 0.56, 0.46, size=value_size, bold=True, color=value_color, font=HEADER_FONT)
    add_text(slide, label, x + 0.2, y + 0.76, w - 0.34, 0.42, size=10.5, color=label_color)


def add_flow_step(slide, num, title, body, x, y, w, *, fill=WHITE, accent=TEAL):
    add_round(slide, x, y, w, 1.05, fill, line=SAND, radius=0.1)
    add_oval(slide, x + 0.18, y + 0.24, 0.48, 0.48, accent)
    add_text(slide, str(num), x + 0.18, y + 0.28, 0.48, 0.25, size=13, bold=True, color=WHITE, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE)
    add_text(slide, title, x + 0.82, y + 0.18, w - 1.0, 0.28, size=14, bold=True, color=INK)
    add_text(slide, body, x + 0.82, y + 0.52, w - 1.0, 0.36, size=10.5, color=SLATE)


def add_background(slide, dark=False):
    if dark:
        add_rect(slide, 0, 0, SLIDE_W, SLIDE_H, INK)
        add_oval(slide, 11.25, -1.55, 2.85, 2.85, INK_2, line=INK_2)
        add_oval(slide, -1.25, 6.05, 1.95, 1.95, TEAL_DARK, line=TEAL_DARK)
    else:
        add_rect(slide, 0, 0, SLIDE_W, SLIDE_H, CREAM)
        add_rect(slide, 0, 0, 0.22, SLIDE_H, TEAL_DARK)
        add_oval(slide, 12.25, 6.38, 1.0, 1.0, SAND, line=SAND)


def make_deck():
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)
    blank = prs.slide_layouts[6]
    total = 10

    # Slide 1
    slide = prs.slides.add_slide(blank)
    add_background(slide, dark=True)
    add_text(slide, "Validation Strategy", 0.72, 1.05, 7.4, 0.7, size=44, bold=True, color=WHITE, font=HEADER_FONT)
    add_text(slide, "for the fl.AI Gene Classification Algorithm", 0.78, 1.88, 6.9, 0.38, size=18, italic=True, color=SAND, font=HEADER_FONT)
    add_text(
        slide,
        "A practical framework for proving that gene labels are accurate, evidence-backed, and stable under reasonable perturbations.",
        0.8,
        2.58,
        6.7,
        0.8,
        size=17,
        color=WHITE,
    )
    for i, (label, fill) in enumerate([("Classification", TEAL), ("References", CORAL), ("Robustness", AMBER)]):
        add_metric(slide, 0.78 + i * 2.12, 4.4, 1.85, 1.22, str(i + 1), label, fill=INK_2, accent=fill, value_size=30)
    add_round(slide, 8.3, 1.35, 3.85, 4.65, WHITE, radius=0.1)
    add_text(slide, "Validation target", 8.65, 1.72, 2.4, 0.32, size=16, bold=True, color=INK)
    add_paragraphs(
        slide,
        [
            "Every accepted classification should have measurable category accuracy.",
            "Every cited paper should contain gene-specific evidence.",
            "The answer should survive source, keyword, and reference-limit stress tests.",
        ],
        8.65,
        2.28,
        3.1,
        2.35,
        size=13,
        color=SLATE,
        gap=8,
    )
    add_chip(slide, "Outcome: defensible hit lists", 8.68, 5.1, 2.8, 0.45, fill=PALE_TEAL, color=TEAL_DARK, size=11)

    # Slide 2
    slide = prs.slides.add_slide(blank)
    add_background(slide)
    add_header(slide, "Validate the pipeline in layers", "Why a single score is not enough", 2, total)
    add_footer(slide)
    layers = [
        ("1", "Identity", "Fly symbol, FBgn, ortholog, and synonym resolution", TEAL),
        ("2", "Reference retrieval", "Candidate papers from FlyBase, PubMed, Europe PMC, GeneRIF, UniProt", AMBER),
        ("3", "Evidence extraction", "Gene-specific function and phenotype claims from title, abstract, or full text", CORAL),
        ("4", "Classification", "Allowed category labels plus confidence and rationale", MOSS),
        ("5", "Robustness", "Sensitivity to limits, sources, keywords, and repeated runs", TEAL_DARK),
    ]
    for i, (num, title, body, color) in enumerate(layers):
        x = 0.85 + i * 2.45
        add_oval(slide, x + 0.68, 2.0, 0.8, 0.8, color)
        add_text(slide, num, x + 0.68, 2.12, 0.8, 0.35, size=18, bold=True, color=WHITE, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE)
        add_card(slide, x, 3.05, 2.15, 2.0, title, body, accent=color, body_size=10.3)
    add_text(slide, "Measure each layer independently, then combine them into release gates.", 1.0, 6.1, 10.8, 0.35, size=18, bold=True, color=INK, font=HEADER_FONT, align=PP_ALIGN.CENTER)

    # Slide 3
    slide = prs.slides.add_slide(blank)
    add_background(slide)
    add_header(slide, "Build a small, trusted gold standard", "The test set is the anchor for every metric", 3, total)
    add_footer(slide)
    add_card(
        slide,
        0.8,
        1.65,
        3.55,
        4.55,
        "Gold file schema",
        ["species, gene_id, symbol", "expected_categories", "gold_pmids", "excluded_categories", "reviewer notes"],
        accent=TEAL,
        body_size=12,
    )
    add_card(
        slide,
        4.75,
        1.65,
        3.55,
        4.55,
        "Panel composition",
        ["Canonical positives per category", "Hard negatives that should be None", "Ambiguous genes for adjudication", "High- and low-literature genes"],
        accent=CORAL,
        body_size=12,
    )
    add_card(
        slide,
        8.7,
        1.65,
        3.55,
        4.55,
        "Review process",
        ["Two reviewers label independently", "Resolve disagreements into a final key", "Keep provenance for each accepted PMID", "Refresh quarterly as categories evolve"],
        accent=AMBER,
        body_size=12,
    )
    add_chip(slide, "Start with 100-200 genes, not the whole universe", 4.3, 6.55, 4.8, 0.42, fill=INK, color=WHITE, size=12)

    # Slide 4
    slide = prs.slides.add_slide(blank)
    add_background(slide)
    add_header(slide, "Classification metrics answer: did we label it right?", "Score category quality before trusting confidence", 4, total)
    add_footer(slide)
    add_metric(slide, 0.82, 1.65, 2.28, 1.25, "F1", "per category plus macro average", fill=TEAL_DARK, accent=WHITE)
    add_metric(slide, 3.35, 1.65, 2.28, 1.25, "J", "Jaccard for multi-label outputs", fill=INK, accent=CORAL)
    add_metric(slide, 5.88, 1.65, 2.28, 1.25, "None", "false positive rate on negatives", fill=CORAL, accent=WHITE)
    add_metric(slide, 8.41, 1.65, 2.28, 1.25, "Cal", "confidence calibration by bucket", fill=AMBER, accent=INK, value_color=INK, label_color=INK)
    add_card(
        slide,
        0.95,
        3.55,
        5.3,
        2.35,
        "Acceptance view",
        "Report all predictions, then show selective accuracy after accepting only classifications above confidence thresholds such as 70, 80, and 90.",
        accent=TEAL,
        body_size=13,
    )
    add_card(
        slide,
        6.9,
        3.55,
        5.3,
        2.35,
        "Failure view",
        "Keep a confusion table by category and a per-gene error list. The error list is where prompt, retrieval, synonym, and category-definition issues become visible.",
        accent=CORAL,
        body_size=13,
    )

    # Slide 5
    slide = prs.slides.add_slide(blank)
    add_background(slide)
    add_header(slide, "Reference metrics answer: is the evidence robust?", "A correct label with weak papers should not pass", 5, total)
    add_footer(slide)
    metrics = [
        ("Recall@K", "Do top-ranked references recover known gold PMIDs?", TEAL),
        ("Precision@support", "How many accepted supporting refs are truly relevant?", MOSS),
        ("Source consensus", "How often do independent sources recover the same paper?", AMBER),
        ("Full-text rate", "How often do accepted refs use full paper text versus metadata only?", CORAL),
        ("Evidence specificity", "Does the text make a direct claim about the target gene?", TEAL_DARK),
        ("Citation stability", "Do the same refs survive repeated runs and perturbations?", INK_2),
    ]
    for i, (title, body, color) in enumerate(metrics):
        x = 0.78 + (i % 3) * 4.05
        y = 1.7 + (i // 3) * 2.0
        add_card(slide, x, y, 3.55, 1.43, title, body, accent=color, body_size=10.8)
    add_text(slide, "Recommended minimum: a high-confidence hit should cite at least two gene-specific supporting references, or one canonical reference with reviewer-approved strength.", 1.05, 6.16, 11.0, 0.42, size=15, bold=True, color=INK, align=PP_ALIGN.CENTER)

    # Slide 6
    slide = prs.slides.add_slide(blank)
    add_background(slide)
    add_header(slide, "Stress tests reveal fragile classifications", "Perturb the inputs and watch what changes", 6, total)
    add_footer(slide)
    add_flow_step(slide, 1, "Reference-limit sweep", "Compare outputs at 25, 50, 100, 250, and 500 refs.", 0.8, 1.55, 5.55, accent=TEAL)
    add_flow_step(slide, 2, "Source ablation", "Drop one source family at a time and re-score.", 6.95, 1.55, 5.55, accent=AMBER)
    add_flow_step(slide, 3, "Keyword perturbation", "Test synonyms, broader terms, narrower terms, and reordered categories.", 0.8, 2.95, 5.55, accent=CORAL)
    add_flow_step(slide, 4, "Leave-one-ref-out", "Remove each supporting paper and quantify label sensitivity.", 6.95, 2.95, 5.55, accent=MOSS)
    add_flow_step(slide, 5, "Repeated-run stability", "Hold references constant and repeat model classification.", 0.8, 4.35, 5.55, accent=TEAL_DARK)
    add_flow_step(slide, 6, "Negative controls", "Run genes and categories expected to have no match.", 6.95, 4.35, 5.55, accent=INK_2)
    add_chip(slide, "Output: stability score per gene and per category", 4.4, 6.35, 4.55, 0.42, fill=INK, color=WHITE, size=12)

    # Slide 7
    slide = prs.slides.add_slide(blank)
    add_background(slide)
    add_header(slide, "Quality gates turn metrics into decisions", "Do not make users interpret every raw score", 7, total)
    add_footer(slide)
    bands = [
        ("Auto-accept", "High confidence, calibrated category, >=2 strong refs, stable under stress", MOSS, "Use in final hit list"),
        ("Review", "One weak dimension: low source diversity, single ref, or medium confidence", AMBER, "Route to scientist review"),
        ("Hold out", "No direct evidence, unstable label, or fails negative control behavior", CORAL, "Exclude from automated claims"),
    ]
    for i, (title, body, color, action) in enumerate(bands):
        x = 0.95 + i * 4.05
        add_round(slide, x, 1.75, 3.45, 3.85, WHITE, line=SAND, radius=0.08)
        add_rect(slide, x, 1.75, 3.45, 0.55, color)
        add_text(slide, title, x + 0.22, 1.89, 3.0, 0.28, size=16, bold=True, color=WHITE)
        add_text(slide, body, x + 0.28, 2.58, 2.9, 1.1, size=12.2, color=SLATE)
        add_chip(slide, action, x + 0.32, 4.62, 2.8, 0.42, fill=PALE_TEAL if color != CORAL else PALE_CORAL, color=INK, size=10.5)
    add_text(slide, "The goal is not to maximize labels. The goal is to maximize reliable labels.", 1.1, 6.35, 11.0, 0.38, size=18, bold=True, color=INK, font=HEADER_FONT, align=PP_ALIGN.CENTER)

    # Slide 8
    slide = prs.slides.add_slide(blank)
    add_background(slide)
    add_header(slide, "Evaluation harness for current outputs", "Use workbook sheets and run-store telemetry", 8, total)
    add_footer(slide)
    add_card(slide, 0.85, 1.9, 3.1, 4.15, "Inputs", ["Output workbook", "run_store JSON", "gold labels CSV", "stress-run manifests"], accent=TEAL, body_size=12)
    add_card(slide, 5.05, 1.9, 3.1, 4.15, "Evaluator", ["Normalize labels", "Match PMIDs/PMCIDs", "Score classification", "Score references", "Compare perturbations"], accent=CORAL, body_size=12)
    add_card(slide, 9.25, 1.9, 3.1, 4.15, "Outputs", ["Metrics JSON", "failure workbook", "calibration plots", "category dashboard"], accent=AMBER, body_size=12)
    add_rect(slide, 4.18, 3.82, 0.52, 0.1, TEAL_DARK)
    add_rect(slide, 4.72, 3.76, 0.14, 0.22, TEAL_DARK, shape=MSO_SHAPE.RIGHT_TRIANGLE)
    add_rect(slide, 8.38, 3.82, 0.52, 0.1, TEAL_DARK)
    add_rect(slide, 8.92, 3.76, 0.14, 0.22, TEAL_DARK, shape=MSO_SHAPE.RIGHT_TRIANGLE)
    add_text(slide, "Same pipeline, plus an evidence-aware scoring layer.", 2.3, 6.42, 8.7, 0.34, size=16, bold=True, color=INK, align=PP_ALIGN.CENTER)

    # Slide 9
    slide = prs.slides.add_slide(blank)
    add_background(slide)
    add_header(slide, "Dashboard view for a validation run", "Make quality visible at a glance", 9, total)
    add_footer(slide)
    add_metric(slide, 0.85, 1.55, 2.2, 1.2, "0.82", "macro F1 target >=0.80", fill=TEAL_DARK, accent=WHITE, value_size=29)
    add_metric(slide, 3.35, 1.55, 2.2, 1.2, "0.91", "supporting-ref precision", fill=INK, accent=CORAL, value_size=29)
    add_metric(slide, 5.85, 1.55, 2.2, 1.2, "86%", "stable classifications", fill=AMBER, accent=INK, value_size=29, value_color=INK, label_color=INK)
    add_metric(slide, 8.35, 1.55, 2.2, 1.2, "14", "review-queue genes", fill=CORAL, accent=WHITE, value_size=29)
    chart_x, chart_y = 1.05, 3.55
    add_text(slide, "Quality by category", chart_x, 3.05, 3.2, 0.3, size=15, bold=True, color=INK)
    bars = [("sleep", 0.86, TEAL), ("circadian", 0.81, MOSS), ("synapse", 0.74, AMBER), ("metabolism", 0.67, CORAL)]
    for i, (label, val, color) in enumerate(bars):
        y = chart_y + i * 0.52
        add_text(slide, label, chart_x, y, 1.15, 0.22, size=10.5, color=SLATE)
        add_rect(slide, chart_x + 1.28, y + 0.04, 3.5, 0.18, SAND)
        add_rect(slide, chart_x + 1.28, y + 0.04, 3.5 * val, 0.18, color)
        add_text(slide, f"{val:.2f}", chart_x + 4.94, y - 0.01, 0.5, 0.2, size=9.5, color=MUTED)
    add_card(slide, 7.0, 3.05, 4.85, 2.75, "Failure drilldown", ["Top error genes by confidence", "Papers accepted by one source only", "Categories with poor calibration", "Predictions sensitive to one reference"], accent=TEAL, body_size=12)
    add_chip(slide, "A scientist should be able to decide what to trust in under five minutes.", 3.05, 6.45, 7.35, 0.42, fill=INK, color=WHITE, size=12)

    # Slide 10
    slide = prs.slides.add_slide(blank)
    add_background(slide, dark=True)
    add_text(slide, "Implementation Roadmap", 0.72, 0.72, 7.4, 0.56, size=34, bold=True, color=WHITE, font=HEADER_FONT)
    add_text(slide, "Move from useful prototype to auditable algorithm.", 0.78, 1.37, 6.7, 0.35, size=16, italic=True, color=SAND, font=HEADER_FONT)
    roadmap = [
        ("Phase 1", "Instrument", "Expose candidate, queried, supporting, source-count, full-text, and QC fields."),
        ("Phase 2", "Gold set", "Curate labels and gold PMIDs for 100-200 genes across categories."),
        ("Phase 3", "Evaluator", "Score classification, references, calibration, and abstention quality."),
        ("Phase 4", "Stress suite", "Automate source ablation, limit sweeps, keyword perturbations, and ref dropout."),
        ("Phase 5", "Quality gates", "Publish auto-accept, review, and hold-out rules with every run."),
    ]
    for i, (phase, title, body) in enumerate(roadmap):
        y = 2.12 + i * 0.82
        add_oval(slide, 0.92, y + 0.04, 0.42, 0.42, [TEAL, CORAL, AMBER, MOSS, WHITE][i])
        add_text(slide, phase, 1.55, y, 1.3, 0.22, size=9.5, bold=True, color=SAND, char_spacing=120)
        add_text(slide, title, 2.8, y - 0.02, 1.2, 0.28, size=15, bold=True, color=WHITE)
        add_text(slide, body, 4.15, y - 0.01, 7.6, 0.28, size=12.5, color=SAND)
    add_round(slide, 8.35, 0.82, 3.7, 0.84, WHITE, radius=0.08)
    add_text(slide, "Deliverable", 8.63, 0.98, 1.2, 0.22, size=10, bold=True, color=TEAL_DARK, char_spacing=130)
    add_text(slide, "a validation report with every run", 9.72, 0.98, 2.0, 0.26, size=12, bold=True, color=INK)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    prs.save(OUTPUT)
    print(OUTPUT)


if __name__ == "__main__":
    make_deck()
