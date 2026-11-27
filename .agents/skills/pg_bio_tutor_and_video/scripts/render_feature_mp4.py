#!/usr/bin/env python3
"""
High-Cadence Micro-Scene Motion Engine for pg_bio (3Blue1Brown / Veritasium Style).

Core Architecture & Strict Visual Constraints:
  1. Maximum 3 Components Per Frame (ZERO Multi-Panel Split-Screen Clutter):
       - Component 1 (Top): Minimal Chapter Progress Bar + Single Focused Headline
       - Component 2 (Center Stage, Full-Screen 1600x690): A SINGLE dynamic visual element
         driven by what the narrator is saying RIGHT NOW. Everything else is hidden!
       - Component 3 (Bottom): Clean Spoken Phrase Spotlight
  2. High-Cadence Micro-Scene Choreography (~150-220 Distinct Motion Beats / Video, New Scene Every 2-5s):
       - Instead of 5 static acts, the speech is tokenized into 35-50 phrase beats (every 2.5-4.5s),
         and each phrase beat sub-animates through 4-6 distinct visual micro-states (~200 motion
         slides per episode).
       - A keyword + narrative-beat router automatically switches the full-screen Center Stage
         between 10 rich animated visual primitives:
           a. `CONCEPT_SPOTLIGHT`      - Giant typography reveal of ONE key question or paradox
           b. `EQUATION_DERIVATION`    - Term-by-term mathematical build with glowing variables
           c. `MATRIX_HEATMAP_ZOOM`    - Squared numerical matrix (PAE / covariance / distance) with
                                         animated cell scan, zoom lens, and 5x5 window highlight
           d. `MORTON_CURVE_DRAW`      - Progressive Z-order space-filling curve drawing + bit
                                         interleaving [Z|Y|X] + octree pruning
           e. `VECTOR_GEOMETRY_CANVAS` - Animated 2D vector space (origin shift d_react = r - e,
                                         c_fwd = p - r, or 15-D cone untangling into ZCA sphere)
           f. `POCKET_VOXEL_PROBE`     - 24-ray Fibonacci water probe rolling across a 5.0 Å grid
           g. `SQL_FOCUS_TERMINAL`     - Full-width focused SQL editor typing & spotlighting 1 clause
           h. `QUERY_RESULT_PLOT`      - Animated chart/bar/scatter plotting the SQL query results live
           i. `CSR_MEMORY_ANATOMY`     - 8 KB PostgreSQL heap page vs 2,000 B TOAST threshold & 16B header peek
           j. `ORGANISM_3D_BRIEF`      - Strictly <= 16 seconds total in the entire video! Gentle slow
                                         alignment sweep (<= 8s) -> camera LOCKED on pocket (<= 8s),
                                         then immediately hidden and replaced by 2D data plots.
"""

import argparse
import asyncio
import hashlib
import json
import math
import re
import shutil
import subprocess
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from PIL import Image, ImageDraw, ImageFont

W, H = 1920, 1080
FPS = 12
DEFAULT_VOICE = "en-US-BrianMultilingualNeural"

# Exact `../ideia.me` Dark Theme Palette
BG_COLOR = (13, 16, 23)            # Deep cinema background
PANEL_BG = (20, 26, 38)            # Card fill
PANEL_DARK = (16, 20, 30)          # Inner stage fill
PANEL_HIGHLIGHT = (28, 39, 64)     # Spotlight card fill
PRIMARY_BLUE = (92, 158, 255)      # #5c9eff
BRAND_BLUE = (67, 97, 238)         # #4361ee
SECONDARY_PURPLE = (177, 88, 252)  # #b158fc
ACCENT_PINK = (255, 90, 157)       # #ff5a9d
PRIMARY_TEAL = (13, 148, 136)      # #0d9488
EMERALD = (16, 185, 129)           # #10b981
LIGHT_BLUE = (56, 189, 248)        # #38bdf8
PG_TOPIC_BLUE = (51, 103, 145)     # #336791
TEXT_MAIN = (236, 242, 250)        # Crisp primary text
TEXT_MUTED = (140, 154, 178)       # Muted slate
AMBER = (255, 193, 7)              # #FFC107
ORANGE = (249, 115, 22)            # #f97316

BLOG_URL = "https://ideia.me"
REPO_URL = "https://github.com/jonatas/pg_bio"
MODELS_DIR = Path("/Users/jonatas/code/ideia.me/assets/models")

DEFAULT_ACT_STEPS = [
    "I. The Paradox",
    "II. The Trap",
    "III. The Math",
    "IV. SQL & Proof",
    "V. Discovery",
]


def load_font(size: int, mono: bool = False) -> ImageFont.FreeTypeFont:
    candidates = (
        [
            "/System/Library/Fonts/SFNSMono.ttf",
            "/System/Library/Fonts/Menlo.ttc",
            "/System/Library/Fonts/Monaco.ttf",
        ]
        if mono
        else [
            "/System/Library/Fonts/SFNS.ttf",
            "/System/Library/Fonts/SFNSMono.ttf",
            "/System/Library/Fonts/Menlo.ttc",
        ]
    )
    for path in candidates:
        if Path(path).exists():
            try:
                return ImageFont.truetype(path, size)
            except Exception:
                pass
    return ImageFont.load_default()


FONT_GIANT = load_font(52)
FONT_HERO = load_font(42)
FONT_H1 = load_font(34)
FONT_H2 = load_font(27)
FONT_H3 = load_font(23)
FONT_BODY_LG = load_font(25)
FONT_BODY = load_font(21)
FONT_SM = load_font(18)
FONT_XS = load_font(15)
FONT_XXS = load_font(13)
FONT_CAPTION = load_font(24)
FONT_MONO_LG = load_font(23, mono=True)
FONT_MONO = load_font(19, mono=True)
FONT_MONO_SM = load_font(16, mono=True)

_PDB_CACHE: Dict[str, List[Dict[str, Any]]] = {}


def load_pdb_ca_trace(pdb_path_str: str, max_points: int = 150) -> List[Dict[str, Any]]:
    if pdb_path_str in _PDB_CACHE:
        return _PDB_CACHE[pdb_path_str]

    path = Path(pdb_path_str)
    if not path.is_absolute():
        path = MODELS_DIR / pdb_path_str

    atoms: List[Dict[str, Any]] = []
    if path.exists():
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
                    try:
                        resn = line[17:20].strip()
                        resi = int(line[22:26].strip())
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        atoms.append({"resn": resn, "resi": resi, "x": x, "y": y, "z": z})
                    except Exception:
                        continue

    if not atoms:
        for i in range(1, 121):
            t_a = i * 0.22
            atoms.append(
                {
                    "resn": "ASP" if i % 6 == 0 else "HIS" if i in (32, 64, 96) else "ALA",
                    "resi": i,
                    "x": 11.0 * math.cos(t_a) + 3.0 * math.sin(t_a * 0.35),
                    "y": (i - 60) * 0.24,
                    "z": 11.0 * math.sin(t_a) + 3.0 * math.cos(t_a * 0.35),
                }
            )

    if len(atoms) > max_points:
        idxs = np.linspace(0, len(atoms) - 1, max_points).astype(int)
        atoms = [atoms[i] for i in idxs]

    coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=np.float64)
    center = coords.mean(axis=0)
    coords -= center
    max_r = np.linalg.norm(coords, axis=1).max()
    if max_r > 1e-6:
        coords /= max_r

    for a, (nx, ny, nz) in zip(atoms, coords):
        a["nx"] = float(nx)
        a["ny"] = float(ny)
        a["nz"] = float(nz)

    _PDB_CACHE[pdb_path_str] = atoms
    return atoms


def text_width(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.FreeTypeFont) -> int:
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0]


def draw_fitted_line(
    draw: ImageDraw.ImageDraw,
    text: str,
    x: int,
    y: int,
    max_px: int,
    fill: Tuple[int, int, int],
    preferred_fonts: Tuple[ImageFont.FreeTypeFont, ...] = (FONT_H1, FONT_H2, FONT_H3, FONT_BODY, FONT_SM, FONT_XS, FONT_XXS),
    center_x: bool = False,
) -> None:
    for fnt in preferred_fonts:
        tw = text_width(draw, text, fnt)
        if tw <= max_px:
            dx = (x - tw // 2) if center_x else x
            draw.text((dx, y), text, fill=fill, font=fnt)
            return
    tw = text_width(draw, text, FONT_XXS)
    dx = (x - tw // 2) if center_x else x
    draw.text((dx, y), text, fill=fill, font=FONT_XXS)


def draw_wrapped_centered(
    draw: ImageDraw.ImageDraw,
    text: str,
    center_x: int,
    y: int,
    max_px: int,
    font: ImageFont.FreeTypeFont,
    fill: Tuple[int, int, int],
    line_h: int = 38,
    max_lines: int = 4,
) -> int:
    words = text.split()
    lines: List[str] = []
    cur = ""
    for w in words:
        cand = f"{cur} {w}".strip()
        if text_width(draw, cand, font) <= max_px:
            cur = cand
        else:
            if cur:
                lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)

    for line in lines[:max_lines]:
        lw = text_width(draw, line, font)
        draw.text((center_x - lw // 2, y), line, font=font, fill=fill)
        y += line_h
    return y


def draw_wrapped_box(
    draw: ImageDraw.ImageDraw,
    text: str,
    x: int,
    y: int,
    max_px: int,
    font: ImageFont.FreeTypeFont,
    fill: Tuple[int, int, int],
    line_h: int = 32,
    max_lines: int = 5,
) -> int:
    words = text.split()
    lines: List[str] = []
    cur = ""
    for w in words:
        cand = f"{cur} {w}".strip()
        if text_width(draw, cand, font) <= max_px:
            cur = cand
        else:
            if cur:
                lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)

    for line in lines[:max_lines]:
        draw.text((x, y), line, font=font, fill=fill)
        y += line_h
    return y


def split_narration_beats(narration: str, max_words_per_beat: int = 11) -> List[str]:
    """
    Splits an act's narration into crisp 2.5-4.0 second spoken beats (~8-11 words each).
    Each beat drives its own focused visual scene and sub-steps through 4-5 motion states
    (~180-220 motion slides across a 2.5-3 minute video).
    """
    words = narration.split()
    if not words:
        return [""]
    beats: List[str] = []
    cur: List[str] = []
    for w in words:
        cur.append(w)
        ends_clause = w.endswith((".", "!", "?", ";", ":", ",", "—"))
        if len(cur) >= max_words_per_beat or (len(cur) >= 6 and ends_clause):
            beats.append(" ".join(cur))
            cur = []
    if cur:
        if beats and len(cur) < 4:
            beats[-1] = beats[-1] + " " + " ".join(cur)
        else:
            beats.append(" ".join(cur))
    return beats


def sanitize_for_tts(text: str) -> str:
    spoken = text
    spoken = re.sub(r"\s*(?:<->|<=>|-->|->|<-|→|←|↔|⇒|⟶|—|--)\s*", ", ", spoken)
    spoken = re.sub(r"\s+[-–]\s+", ", ", spoken)
    spoken = re.sub(r"(\d)\s*[-–]\s*(\d)", r"\1 to \2", spoken)
    spoken = re.sub(r"(?<=\w)[-–](?=\w)", " ", spoken)
    spoken = re.sub(r"[-–—]+", ", ", spoken)
    spoken = re.sub(r"([:;])\s*,+", r"\1", spoken)
    spoken = re.sub(r",\s*,+", ",", spoken)
    spoken = re.sub(r"\s+,", ",", spoken)
    spoken = re.sub(r"\s{2,}", " ", spoken).strip()
    return spoken


async def synthesize_single_act(idx: int, text: str, voice: str, out_mp3: Path, ffprobe_bin: str, ffmpeg_bin: str) -> float:
    spoken_text = sanitize_for_tts(text)
    sig = hashlib.sha256(f"{voice}||{spoken_text}".encode("utf-8")).hexdigest()[:16]
    sig_file = out_mp3.with_suffix(".mp3.sig")

    if out_mp3.exists() and out_mp3.stat().st_size > 1024 and sig_file.exists() and sig_file.read_text(encoding="utf-8").strip() == sig:
        probe = subprocess.run(
            [ffprobe_bin, "-v", "error", "-show_entries", "format=duration", "-of", "default=noprint_wrappers=1:nokey=1", str(out_mp3)],
            capture_output=True,
            text=True,
            check=True,
        )
        return float(probe.stdout.strip())

    try:
        import edge_tts

        communicate = edge_tts.Communicate(spoken_text, voice, rate="+4%")
        await communicate.save(str(out_mp3))
    except Exception:
        aiff_path = out_mp3.with_suffix(".aiff")
        subprocess.run(["say", "-v", "Samantha", "-r", "190", "-o", str(aiff_path), spoken_text], check=True)
        subprocess.run(
            [ffmpeg_bin, "-y", "-i", str(aiff_path), "-c:a", "libmp3lame", "-q:a", "2", str(out_mp3)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        aiff_path.unlink(missing_ok=True)

    sig_file.write_text(sig, encoding="utf-8")
    probe = subprocess.run(
        [ffprobe_bin, "-v", "error", "-show_entries", "format=duration", "-of", "default=noprint_wrappers=1:nokey=1", str(out_mp3)],
        capture_output=True,
        text=True,
        check=True,
    )
    return float(probe.stdout.strip())


def format_srt_time(seconds: float) -> str:
    ms = int(round((seconds % 1.0) * 1000))
    total_s = int(seconds)
    s = total_s % 60
    m = (total_s // 60) % 60
    h = total_s // 3600
    return f"{h:02d}:{m:02d}:{s:02d},{ms:03d}"


def format_mm_ss(seconds: float) -> str:
    total_s = int(seconds)
    return f"{total_s // 60:02d}:{total_s % 60:02d}"


def generate_voiceover_and_srt(
    out_dir: Path,
    narrations: List[str],
    voice: str,
    ffmpeg_bin: str,
    ffprobe_bin: str,
) -> Tuple[Path, List[float]]:
    act_files: List[Path] = [out_dir / f"voiceover_act{i+1}.mp3" for i in range(len(narrations))]

    async def _run_all() -> List[float]:
        tasks = [
            synthesize_single_act(i, narrations[i], voice, act_files[i], ffprobe_bin, ffmpeg_bin)
            for i in range(len(narrations))
        ]
        return list(await asyncio.gather(*tasks))

    durations = asyncio.run(_run_all())

    concat_list = out_dir / "voiceover_concat.txt"
    concat_list.write_text("\n".join(f"file '{p.name}'" for p in act_files) + "\n", encoding="utf-8")
    full_vo_mp3 = out_dir / "voiceover.mp3"
    subprocess.run(
        [
            ffmpeg_bin,
            "-y",
            "-f",
            "concat",
            "-safe",
            "0",
            "-i",
            str(concat_list),
            "-c:a",
            "libmp3lame",
            "-b:a",
            "192k",
            str(full_vo_mp3),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    concat_list.unlink(missing_ok=True)

    srt_lines = []
    vtt_lines = ["WEBVTT\n"]
    cursor = 0.0
    sub_idx = 1
    for text, dur in zip(narrations, durations):
        chunks = split_narration_beats(text, max_words_per_beat=12)
        chunk_dur = dur / max(1, len(chunks))
        for ch in chunks:
            start_t = cursor
            end_t = cursor + chunk_dur
            s_srt = format_srt_time(start_t)
            e_srt = format_srt_time(end_t)
            srt_lines.append(f"{sub_idx}\n{s_srt} --> {e_srt}\n{ch}\n")
            vtt_lines.append(f"{sub_idx}\n{s_srt.replace(',', '.')} --> {e_srt.replace(',', '.')}\n{ch}\n")
            cursor = end_t
            sub_idx += 1
    (out_dir / "captions.srt").write_text("\n".join(srt_lines), encoding="utf-8")
    (out_dir / "captions.vtt").write_text("\n".join(vtt_lines), encoding="utf-8")

    return full_vo_mp3, durations


# ==============================================================================
# COMPONENT 1 (TOP) & COMPONENT 3 (BOTTOM): MINIMAL FRAMING
# ==============================================================================


def draw_component_1_top_bar(
    draw: ImageDraw.ImageDraw,
    layer_str: str,
    act_idx: int,
    act_steps: List[str],
    progress: float,
    elapsed_sec: float,
    total_sec: float,
    motion_slide_num: int,
    total_motion_slides: int,
    scene_headline: str,
) -> None:
    """Component 1 of 3: Clean top progress bar + single focused scene headline."""
    draw.rectangle([0, 0, W, 72], fill=(10, 12, 18))
    draw.line([(0, 72), (W, 72)], fill=(38, 48, 68), width=2)
    draw.rectangle([0, 69, int(W * progress), 72], fill=PRIMARY_BLUE)

    bw = min(430, text_width(draw, layer_str, FONT_XS) + 28)
    draw.rounded_rectangle([44, 15, 44 + bw, 55], radius=12, fill=PG_TOPIC_BLUE, outline=PRIMARY_BLUE, width=2)
    draw_fitted_line(draw, layer_str, 58, 25, bw - 24, (255, 255, 255), (FONT_XS, FONT_XXS))

    step_x = 505
    for idx, label in enumerate(act_steps[:5]):
        sw = text_width(draw, label, FONT_XS) + 22
        is_cur = idx == act_idx
        is_done = idx < act_idx
        fill_c = (37, 99, 235) if is_cur else PANEL_BG
        out_c = PRIMARY_BLUE if is_cur else PRIMARY_TEAL if is_done else (45, 55, 75)
        txt_c = (255, 255, 255) if is_cur else LIGHT_BLUE if is_done else TEXT_MUTED
        draw.rounded_rectangle([step_x, 16, step_x + sw, 54], radius=10, fill=fill_c, outline=out_c, width=2)
        draw.text((step_x + 11, 25), label, fill=txt_c, font=FONT_XS)
        step_x += sw + 8

    slide_badge = f"MOTION {motion_slide_num:03d}/{total_motion_slides:03d}  |  {format_mm_ss(elapsed_sec)} / {format_mm_ss(total_sec)}"
    tw = text_width(draw, slide_badge, FONT_XS) + 28
    draw.rounded_rectangle([W - 44 - tw, 15, W - 44, 55], radius=12, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
    draw.text((W - 44 - tw + 14, 25), slide_badge, fill=TEXT_MAIN, font=FONT_XS)

    # Single Focused Scene Headline
    draw_fitted_line(draw, scene_headline, W // 2, 92, W - 160, PRIMARY_BLUE, (FONT_H1, FONT_H2, FONT_H3), center_x=True)


def draw_component_3_caption_bar(
    draw: ImageDraw.ImageDraw,
    active_beat_text: str,
    beat_idx: int,
    total_beats: int,
    sub_phase: float,
) -> None:
    """Component 3 of 3: Clean bottom spoken-phrase bar with word-by-word emphasis."""
    draw.rounded_rectangle([140, 910, W - 140, 1008], radius=16, fill=PANEL_BG, outline=PRIMARY_BLUE, width=2)

    # Progress bar along top of caption card showing beat progression
    bar_w = int((W - 320) * ((beat_idx + sub_phase) / max(1, total_beats)))
    draw.rounded_rectangle([160, 916, 160 + max(12, bar_w), 921], radius=3, fill=EMERALD)

    words = active_beat_text.split()
    hi_word_idx = min(len(words) - 1, int(sub_phase * len(words))) if words else 0

    # Render clean centered phrase with subtle gold highlight on the active spoken words
    full_str = " ".join(words)
    if text_width(draw, full_str, FONT_CAPTION) <= (W - 360):
        total_w = text_width(draw, full_str, FONT_CAPTION)
        cx = (W - total_w) // 2
        for idx_w, w_tok in enumerate(words):
            col = AMBER if abs(idx_w - hi_word_idx) <= 1 else TEXT_MAIN
            draw.text((cx, 946), w_tok, fill=col, font=FONT_CAPTION)
            cx += text_width(draw, w_tok + " ", FONT_CAPTION)
    else:
        draw_wrapped_centered(draw, full_str, W // 2, 932, W - 360, FONT_BODY_LG, TEXT_MAIN, line_h=32, max_lines=2)

    # Minimal footer links
    draw.rectangle([0, 1026, W, H], fill=(10, 12, 18))
    draw.line([(0, 1026), (W, 1026)], fill=(38, 48, 68), width=1)
    draw.text((64, 1040), f"Interactive 3D Studio: {BLOG_URL}/pg-bio-videos/", fill=PRIMARY_BLUE, font=FONT_MONO_SM)
    right_tag = f"Open-Source Engine: {REPO_URL}"
    rw = text_width(draw, right_tag, FONT_MONO_SM)
    draw.text((W - rw - 64, 1040), right_tag, fill=SECONDARY_PURPLE, font=FONT_MONO_SM)


# ==============================================================================
# COMPONENT 2 (CENTER STAGE): 10 FOCUSED FULL-SCREEN MOTION PRIMITIVES
# Never more than 1 focused visual concept on screen at once!
# ==============================================================================

STAGE_BOX = (140, 152, W - 140, 882)


def draw_scene_concept_spotlight(
    draw: ImageDraw.ImageDraw,
    beat_text: str,
    focus_title: str,
    focus_metric: str,
    sub_phase: float,
    t: float,
    accent: Tuple[int, int, int],
) -> None:
    """
    Scene Primitive 1: Full-stage typography & single-question spotlight.
    Hides all clutter so the viewer focuses 100% on the core question or paradox being asked.
    """
    x1, y1, x2, y2 = STAGE_BOX
    cx = (x1 + x2) // 2
    draw.rounded_rectangle([x1, y1, x2, y2], radius=24, fill=PANEL_BG, outline=accent, width=3)

    # Top pill
    pill = "FOCUSED CONCEPT SPOTLIGHT"
    pw = text_width(draw, pill, FONT_SM) + 40
    draw.rounded_rectangle([cx - pw // 2, y1 + 42, cx + pw // 2, y1 + 88], radius=22, fill=PANEL_HIGHLIGHT, outline=accent, width=2)
    draw.text((cx - pw // 2 + 20, y1 + 54), pill, fill=accent, font=FONT_SM)

    # Main Question / Concept Title (large & spacious)
    y_after = draw_wrapped_centered(draw, focus_title, cx, y1 + 135, 1320, FONT_HERO, (255, 255, 255), line_h=56, max_lines=2)

    # Animated divider line that grows with sub_phase
    line_half = int(480 * min(1.0, 0.25 + sub_phase * 0.85))
    draw.line([(cx - line_half, y_after + 24), (cx + line_half, y_after + 24)], fill=accent, width=4)

    # Spoken thought card in center
    card_w = 1280
    draw.rounded_rectangle([cx - card_w // 2, y_after + 58, cx + card_w // 2, y_after + 250], radius=18, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
    draw_wrapped_centered(draw, f"“{beat_text}”", cx, y_after + 98, 1180, FONT_H1, TEXT_MAIN, line_h=48, max_lines=3)

    # Bottom single key metric pill
    if focus_metric:
        mw = min(1320, text_width(draw, focus_metric, FONT_H3) + 64)
        draw.rounded_rectangle([cx - mw // 2, y2 - 115, cx + mw // 2, y2 - 45], radius=16, fill=(8, 48, 38), outline=EMERALD, width=2)
        draw_fitted_line(draw, focus_metric, cx, y2 - 93, mw - 40, EMERALD, (FONT_H3, FONT_BODY, FONT_SM), center_x=True)


def draw_scene_matrix_heatmap_zoom(
    draw: ImageDraw.ImageDraw,
    visual_mode: str,
    beat_text: str,
    sub_phase: float,
    t: float,
    mode_variant: str = "pae_window",
) -> None:
    """
    Scene Primitive 2: Squared Numerical Matrix with Animated Zoom / Cell Selection.
    Used when talking about PAE matrices, 5x5 MiniPAE windows, covariance matrices, or distance grids.
    Only TWO elements inside the stage: the 10x10 numerical matrix on the left and the zoomed sub-window on the right.
    """
    x1, y1, x2, y2 = STAGE_BOX
    draw.rounded_rectangle([x1, y1, x2, y2], radius=24, fill=PANEL_BG, outline=PRIMARY_BLUE, width=3)

    gx0, gy0 = x1 + 76, y1 + 92
    cell = 54
    n_grid = 10

    is_cov = visual_mode == "zca_whitening" or mode_variant == "covariance"
    is_dogma = visual_mode == "central_dogma_gist"
    mat_label = (
        "1280x1280 Off-Diagonal Covariance Matrix Sigma (Anisotropic Cone)"
        if is_cov
        else "20-Amino-Acid Charge & Hydropathy Matrix: Standard Enzyme vs. Lake Magadi D3SWG7"
        if is_dogma
        else "Cross-Chain PAE Matrix (Å): 10-aa SLiM Window vs. 70-aa Flexible Tail"
    )
    draw_fitted_line(draw, mat_label, gx0, y1 + 36, 1480, LIGHT_BLUE, (FONT_H2, FONT_H3))

    dogma_aas = ["D-", "E-", "D-", "E-", "K+", "R+", "I", "L", "V", "F"]
    for r in range(n_grid):
        for c in range(n_grid):
            rx = gx0 + c * cell
            ry = gy0 + r * cell
            if is_cov:
                if r == c:
                    val = 1.00
                    fc, oc = (12, 82, 64), EMERALD
                else:
                    # Show off-diagonal correlation collapsing to 0.00 as sub_phase advances
                    raw_cov = 0.88 - 0.04 * abs(r - c)
                    val = raw_cov * max(0.0, 1.0 - sub_phase * 1.35)
                    fc = (70, 22, 42) if val > 0.25 else (18, 26, 38)
                    oc = ACCENT_PINK if val > 0.25 else (45, 58, 78)
                txt = f"{val:.2f}"
            elif is_dogma:
                in_5x5 = (3 <= r <= 7) and (2 <= c <= 6)
                if in_5x5:
                    fc, oc = (95, 42, 12), ORANGE
                    txt = dogma_aas[(r + c) % 4]
                else:
                    fc, oc = (18, 32, 48), (45, 64, 92)
                    txt = dogma_aas[4 + ((r * 3 + c) % 6)]
            else:
                in_5x5 = (3 <= r <= 7) and (2 <= c <= 6)
                on_diag = in_5x5 and ((r - 3) == (c - 2))
                if on_diag:
                    fc, oc = (6, 105, 76), EMERALD
                    txt = "0.45"
                elif in_5x5:
                    fc, oc = (14, 62, 54), PRIMARY_TEAL
                    txt = "0.82"
                else:
                    fc, oc = (46, 20, 34), (95, 40, 62)
                    txt = "28.5"

            draw.rounded_rectangle([rx, ry, rx + cell - 4, ry + cell - 4], radius=6, fill=fc, outline=oc, width=2)
            draw.text((rx + 8, ry + 17), txt, fill=(255, 255, 255) if (is_cov or in_5x5) else TEXT_MUTED, font=FONT_XS)

    # Animated Zoom Selection Frame (5x5 Window)
    pulse = int(4 * math.sin(t * 5.0))
    zx1 = gx0 + 2 * cell - 5 - pulse
    zy1 = gy0 + 3 * cell - 5 - pulse
    zx2 = gx0 + 7 * cell + 1 + pulse
    zy2 = gy0 + 8 * cell + 1 + pulse
    draw.rounded_rectangle([zx1, zy1, zx2, zy2], radius=10, outline=AMBER, width=4)

    # Right Side: Single Focused Zoomed Callout Card connected by zoom lines
    rx1 = gx0 + n_grid * cell + 68
    rx2 = x2 - 64
    ry1 = gy0 + 30
    ry2 = gy0 + n_grid * cell - 20
    draw.line([(zx2, zy1), (rx1, ry1)], fill=AMBER, width=2)
    draw.line([(zx2, zy2), (rx1, ry2)], fill=AMBER, width=2)
    draw.rounded_rectangle([rx1, ry1, rx2, ry2], radius=18, fill=PANEL_DARK, outline=AMBER, width=3)

    if is_cov:
        state_lbl = "AFTER ZCA WHITENING (Sigma = I)" if sub_phase > 0.55 else "BEFORE WHITENING (Off-Diagonal Noise)"
        draw_fitted_line(draw, state_lbl, rx1 + 28, ry1 + 26, (rx2 - rx1) - 56, EMERALD if sub_phase > 0.55 else ACCENT_PINK, (FONT_H3, FONT_BODY))
        draw_wrapped_box(
            draw,
            "Symmetric ZCA (W = V * Lambda^-1/2 * V^T) zeroes out off-diagonal covariance so all 1,280 dimensions become orthogonal.",
            rx1 + 28,
            ry1 + 78,
            (rx2 - rx1) - 56,
            FONT_BODY,
            TEXT_MAIN,
            line_h=34,
            max_lines=4,
        )
        metric_box_y = ry1 + 250
        draw.rounded_rectangle([rx1 + 28, metric_box_y, rx2 - 28, ry2 - 32], radius=14, fill=PANEL_HIGHLIGHT, outline=EMERALD, width=2)
        draw.text((rx1 + 48, metric_box_y + 24), "Participation Ratio:  5.8 dims  -->  940+ dims", fill=LIGHT_BLUE, font=FONT_BODY)
        draw.text((rx1 + 48, metric_box_y + 72), "SCOPe-40 Superfamily Recall:  58.55%  -->  88.32%", fill=EMERALD, font=FONT_H3)
    elif is_dogma:
        draw_fitted_line(draw, "ZOOMING IN: ACIDIC Asp-/Glu- SURFACE SHELL", rx1 + 28, ry1 + 26, (rx2 - rx1) - 56, ORANGE, (FONT_H3, FONT_BODY))
        draw_wrapped_box(
            draw,
            "Lake Magadi D3SWG7 mutates 68% of its outer surface into negatively charged Aspartate (D-) and Glutamate (E-) carboxylates to lock onto water in 3.5 M salt.",
            rx1 + 28,
            ry1 + 76,
            (rx2 - rx1) - 56,
            FONT_BODY,
            TEXT_MAIN,
            line_h=33,
            max_lines=4,
        )
        by = ry1 + 235
        bar_max = (rx2 - rx1) - 96
        draw.text((rx1 + 32, by), "Standard Mesophilic Enzyme (pI = 6.40):  Dehydrates & Precipitates in 3.5 M Salt", fill=ACCENT_PINK, font=FONT_SM)
        draw.rounded_rectangle([rx1 + 32, by + 30, rx1 + 32 + int(bar_max * 0.22), by + 56], radius=6, fill=ACCENT_PINK)

        by2 = by + 84
        grow = min(1.0, 0.35 + sub_phase * 0.85)
        draw.text((rx1 + 32, by2), "Lake Magadi D3SWG7 (pI = 4.38, GRAVY = -0.412):  100% Soluble + 0.44 Å Core Match!", fill=EMERALD, font=FONT_SM)
        draw.rounded_rectangle([rx1 + 32, by2 + 30, rx1 + 32 + int(bar_max * 0.94 * grow), by2 + 58], radius=6, fill=EMERALD)
    else:
        draw_fitted_line(draw, "ZOOMING IN: 5x5 CONTIGUOUS MiniPAE WINDOW", rx1 + 28, ry1 + 26, (rx2 - rx1) - 56, AMBER, (FONT_H3, FONT_BODY))
        draw_wrapped_box(
            draw,
            "Only the 10-aa bound motif (green cells < 1.0 Å) touches the receptor. The 70-aa disordered tail (pink cells = 28.5 Å) floats free in water.",
            rx1 + 28,
            ry1 + 76,
            (rx2 - rx1) - 56,
            FONT_BODY,
            TEXT_MAIN,
            line_h=33,
            max_lines=4,
        )
        # Live Animated Score Comparison Bars
        by = ry1 + 235
        bar_max = (rx2 - rx1) - 96
        draw.text((rx1 + 32, by), "Whole-Chain ipTM (Divides by all 80 aa):  0.1250 [False Negative]", fill=ACCENT_PINK, font=FONT_SM)
        draw.rounded_rectangle([rx1 + 32, by + 30, rx1 + 32 + int(bar_max * 0.14), by + 56], radius=6, fill=ACCENT_PINK)

        by2 = by + 84
        grow = min(1.0, 0.35 + sub_phase * 0.85)
        draw.text((rx1 + 32, by2), "Contact-Only ipSAE (N_i = 10 contacts) & 5x5 MiniPAE (0.475 Å):  0.8315 [RESCUED!]", fill=EMERALD, font=FONT_SM)
        draw.rounded_rectangle([rx1 + 32, by2 + 30, rx1 + 32 + int(bar_max * 0.8315 * grow), by2 + 58], radius=6, fill=EMERALD)


def draw_scene_morton_curve_animation(
    draw: ImageDraw.ImageDraw,
    beat_text: str,
    sub_phase: float,
    t: float,
) -> None:
    """
    Scene Primitive 3: Progressive Morton Z-Order Space-Filling Curve Drawing & Bit Interleaving.
    Animates the Z-curve step-by-step across a 16-cell octree/quadtree slice while showing
    the exact bit-interleaving [Z_i | Y_i | X_i] -> 64-bit BIGINT key.
    """
    x1, y1, x2, y2 = STAGE_BOX
    draw.rounded_rectangle([x1, y1, x2, y2], radius=24, fill=PANEL_BG, outline=PRIMARY_BLUE, width=3)

    gx0, gy0 = x1 + 88, y1 + 98
    cell = 128
    draw_fitted_line(draw, "DRAWING THE 64-BIT MORTON Z-ORDER CURVE IN 3D ATOMIC SPACE", gx0, y1 + 36, 1440, PRIMARY_BLUE, (FONT_H2, FONT_H3))

    # Morton order traversal for 4x4 grid: interleave 2-bit x and 2-bit y
    morton_pts: List[Tuple[int, int, int, int, int]] = []
    for m_code in range(16):
        x_b = (m_code & 1) | ((m_code >> 1) & 2)
        y_b = ((m_code >> 1) & 1) | ((m_code >> 2) & 2)
        rx = gx0 + x_b * cell
        ry = gy0 + y_b * cell
        morton_pts.append((m_code, x_b, y_b, rx + cell // 2, ry + cell // 2))

    # How many Morton segments are drawn at this sub_phase
    n_drawn = max(2, min(16, int(2 + sub_phase * 15)))

    # Draw octree cells
    for m_code, x_b, y_b, cx_p, cy_p in morton_pts:
        rx = gx0 + x_b * cell
        ry = gy0 + y_b * cell
        in_sphere = (1 <= x_b <= 2) and (1 <= y_b <= 2)
        visited = m_code < n_drawn
        fc = (8, 82, 62) if (in_sphere and visited) else (24, 34, 54) if visited else (16, 20, 30)
        oc = EMERALD if in_sphere else LIGHT_BLUE if visited else (40, 50, 70)
        draw.rounded_rectangle([rx, ry, rx + cell - 10, ry + cell - 10], radius=12, fill=fc, outline=oc, width=3 if in_sphere else 2)
        draw.text((rx + 14, ry + 12), f"Z={m_code:04b}", fill=EMERALD if in_sphere else TEXT_MUTED, font=FONT_MONO_SM)

    # Draw progressive Z-order polyline
    for i in range(n_drawn - 1):
        p1 = (morton_pts[i][3], morton_pts[i][4])
        p2 = (morton_pts[i + 1][3], morton_pts[i + 1][4])
        draw.line([p1, p2], fill=AMBER, width=5)

    for i in range(n_drawn):
        px, py = morton_pts[i][3], morton_pts[i][4]
        r = 11 if i == n_drawn - 1 else 7
        draw.ellipse([px - r, py - r, px + r, py + r], fill=AMBER if i == n_drawn - 1 else LIGHT_BLUE, outline=(255, 255, 255), width=2)

    # Query bounding sphere overlay + redraw Z=xxxx pill labels ON TOP of the circle so they never collide!
    scx, scy = gx0 + 2 * cell - 5, gy0 + 2 * cell - 5
    draw.ellipse([scx - 185, scy - 185, scx + 185, scy + 185], outline=ACCENT_PINK, width=4)
    for m_code, x_b, y_b, _, _ in morton_pts:
        rx = gx0 + x_b * cell
        ry = gy0 + y_b * cell
        in_sphere = (1 <= x_b <= 2) and (1 <= y_b <= 2)
        draw.rounded_rectangle([rx + 8, ry + 8, rx + 82, ry + 34], radius=6, fill=(11, 15, 24), outline=EMERALD if in_sphere else (45, 58, 78), width=1)
        draw.text((rx + 14, ry + 12), f"Z={m_code:04b}", fill=EMERALD if in_sphere else TEXT_MUTED, font=FONT_MONO_SM)

    # Right Side: Bit-Interleaving & Octree Pruning Walkthrough Card
    rx1 = gx0 + 4 * cell + 68
    rx2 = x2 - 68
    ry1 = gy0
    ry2 = gy0 + 4 * cell - 10
    draw.rounded_rectangle([rx1, ry1, rx2, ry2], radius=18, fill=PANEL_DARK, outline=EMERALD, width=3)

    cur_m, cur_x, cur_y, _, _ = morton_pts[n_drawn - 1]
    draw.text((rx1 + 32, ry1 + 28), "BIT-INTERLEAVING (X, Y, Z) -> 1D B-TREE KEY", fill=LIGHT_BLUE, font=FONT_H3)
    draw.text((rx1 + 32, ry1 + 84), f"X_bits = {cur_x:02b}    Y_bits = {cur_y:02b}    Z_bits = 10", fill=TEXT_MAIN, font=FONT_MONO_LG)
    draw.text((rx1 + 32, ry1 + 132), f"Interleaved Morton Key = {cur_m:04b} (1D BIGINT)", fill=AMBER, font=FONT_MONO_LG)

    draw.rounded_rectangle([rx1 + 32, ry1 + 205, rx2 - 32, ry1 + 335], radius=14, fill=(8, 52, 42), outline=EMERALD, width=2)
    draw.text((rx1 + 52, ry1 + 225), "Interior Octants (fully_inside = TRUE):", fill=EMERALD, font=FONT_H3)
    draw.text((rx1 + 52, ry1 + 268), "Every atom inside these emerald cubes is guaranteed", fill=TEXT_MAIN, font=FONT_BODY)
    draw.text((rx1 + 52, ry1 + 300), "to be <= 5.0 Å — we skip square-roots completely!", fill=TEXT_MAIN, font=FONT_BODY)

    draw.rounded_rectangle([rx1 + 32, ry1 + 355, rx2 - 32, ry2 - 24], radius=14, fill=PANEL_HIGHLIGHT, outline=PRIMARY_BLUE, width=2)
    draw.text((rx1 + 52, ry1 + 378), "Result: Prunes 98.6% of outer atoms (74x Fewer Atoms Scanned)", fill=AMBER, font=FONT_H3)
    draw.text((rx1 + 52, ry1 + 422), "Turns a 3D catalytic pocket search into O(log N) B-Tree range scans.", fill=TEXT_MAIN, font=FONT_BODY)


def draw_scene_equation_derivation(
    draw: ImageDraw.ImageDraw,
    visual_mode: str,
    beat_text: str,
    sub_phase: float,
    t: float,
) -> None:
    """
    Scene Primitive 4: Step-by-Step Mathematical Heuristic Derivation.
    Builds the core formula term-by-term across sub_phase so the user sees WHY each term exists.
    """
    x1, y1, x2, y2 = STAGE_BOX
    cx = (x1 + x2) // 2
    draw.rounded_rectangle([x1, y1, x2, y2], radius=24, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=3)

    if visual_mode == "ipsae_minipae":
        title = "HOW DUNBRACK'S CONTACT-ONLY d0(N_i) ELIMINATES IDR TAIL DILUTION"
        eq_main = "ipSAE(A -> B)  =  max_i [  (1 / N_i) * SUM_{j in C_i}  1 / (1 + (PAE_ij / d0(N_i))^2)  ]"
        terms = [
            ("Term 1: Contact Filter C_i = { j : PAE_ij < 10.0 Å }", "Discards the 70-residue disordered tail (28.5 Å) before averaging!", EMERALD),
            ("Term 2: Contact-Only Scale d0(N_i) = max(1.0, 1.24 * cbrt(N_i - 15) - 1.8)", "Depends ONLY on contacting count N_i = 10 — invariant to partner length L_B!", LIGHT_BLUE),
            ("Term 3: 5x5 Contiguous MiniPAE = min_{u,v} mean_{k=0..4}(PAE_{u+k, v+k})", "Requires 5 consecutive diagonal contacts < 12 Å (0.475 Å vs 25.08 Å decoy).", AMBER),
        ]
    elif visual_mode == "enzref_pslim":
        title = "HOW ENZREF ORIGIN CENTERING MAKES DEAD BINDERS COLLAPSE TO ZERO"
        eq_main = "d_react = r - v_bg      d_prod = p - v_bg      c_fwd = d_prod - d_react = p - r"
        terms = [
            ("Term 1: Shift Origin to Enzyme Background (v_bg)", "Measures how the substrate and product perturb the enzyme pocket itself.", LIGHT_BLUE),
            ("Term 2: Directed Transition Vector c_fwd = p - r", "If a molecule binds without reacting (p == r), then c_fwd = 0 and FluxRank = 0.0000!", EMERALD),
            ("Term 3: Exact Reaction-Reversal Antisymmetry", "Swapping reactant and product (r <-> p) flips c_fwd by 180°: Flux(e,p,r) = -Flux(e,r,p).", AMBER),
        ]
    elif visual_mode == "multipocket_posebusters":
        title = "THE 3 GEOMETRIC LAWS OF MULTI-POCKET POSEBUSTERS VALIDATION"
        eq_main = "Concavity(p) = |{ a : ||a - p|| <= 6.0 Å }| >= 11     RingRMSD <= 0.25 Å     d_ij / r_vdW >= 0.75"
        terms = [
            ("Law 1: 24-Ray Fibonacci SAS Concavity Probe (r = 1.40 Å)", "Finds all buried concave pockets (pk1..pk5) across an O(1) 5.0 Å 3D voxel hash grid.", LIGHT_BLUE),
            ("Law 2: 3x3 Jacobi SVD Aromatic Ring Planarity & VdW Ratios", "Rejects warped rings (>0.25 Å) and steric clashes (<0.75 VdW sum) in 1.65 ms.", EMERALD),
            ("Law 3: Non-Hiding Dual-Ledger Salt Provenance", "Strips vendor [Na+]/[Cl-] counterions while logging raw_passes_all_21 = FALSE.", AMBER),
        ]
    elif visual_mode == "zca_whitening":
        title = "HOW SYMMETRIC ZCA WHITENING UNTANGLES THE 15-DIMENSIONAL CONE"
        eq_main = "Sigma = (1/N) X_c^T X_c = V Lambda V^T    -->    W_ZCA = V (Lambda + eps*I)^(-1/2) V^T"
        terms = [
            ("Term 1: Center & Decompose Corpus Covariance Sigma", "Exposes the 5.8 dominant SwiGLU axes where 99% of raw cosine similarity piles up.", LIGHT_BLUE),
            ("Term 2: Rescale by Inverse Square-Root Eigenvalues Lambda^(-1/2)", "Equalizes variance across all 1,280 dimensions so no axis dominates.", EMERALD),
            ("Term 3: Rotate Back by V^T to Original Coordinate Axes", "Preserves semantic axis alignment while lifting SCOPe-40 recall from 58.55% to 88.32%!", AMBER),
        ]
    else:
        title = "THE MATHEMATICAL INVARIANT BEHIND THE HEURISTIC"
        eq_main = "pI = argzero_pH( Sum_i q_i(pH) )      &&      kabsch_rmsd(P, Q) = sqrt( (1/N) Sum ||R*p_i - q_i||^2 )"
        terms = [
            ("Term 1: Henderson-Hasselbalch Charge Bisection (pI = 4.38)", "Coating the outer surface in Asp-/Glu- carboxylates locks onto water in 3.5 M salt.", ORANGE),
            ("Term 2: 3x3 Jacobi SVD Kabsch Rotation R = V * diag(1,1,det) * U^T", "Rotates two proteins from arbitrary orientations into sub-Angstrom pocket alignment.", EMERALD),
            ("Term 3: Zero-IPC Execution Inside 8 KB Shared Buffers", "Evaluates physical invariants right where the 2.34M protein tuples live.", LIGHT_BLUE),
        ]

    draw_fitted_line(draw, title, cx, y1 + 36, 1520, SECONDARY_PURPLE, (FONT_H2, FONT_H3), center_x=True)

    # Formula Banner
    draw.rounded_rectangle([x1 + 64, y1 + 92, x2 - 64, y1 + 184], radius=16, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=3)
    draw_fitted_line(draw, eq_main, cx, y1 + 122, 1440, (255, 255, 255), (FONT_MONO_LG, FONT_MONO, FONT_MONO_SM), center_x=True)

    # Progressive Term Spotlight: highlight active term based on sub_phase
    active_term = min(len(terms) - 1, int(sub_phase * len(terms)))
    ty = y1 + 216
    for idx, (t_head, t_desc, col) in enumerate(terms):
        is_cur = idx == active_term
        is_vis = idx <= active_term
        if not is_vis:
            # Unrevealed future term: subtle placeholder bar so screen stays uncluttered
            draw.rounded_rectangle([x1 + 84, ty, x2 - 84, ty + 76], radius=12, fill=(16, 20, 30), outline=(36, 44, 62), width=1)
            draw.text((x1 + 112, ty + 24), f"Upcoming Step {idx + 1}: {t_head}", fill=(80, 92, 115), font=FONT_BODY)
            ty += 98
        else:
            card_h = 152 if is_cur else 118
            draw.rounded_rectangle(
                [x1 + 64, ty, x2 - 64, ty + card_h],
                radius=16,
                fill=PANEL_HIGHLIGHT if is_cur else PANEL_DARK,
                outline=col,
                width=4 if is_cur else 2,
            )
            tag = "ACTIVE TERM  |  " if is_cur else "[VERIFIED]  "
            draw_fitted_line(draw, tag + t_head, x1 + 96, ty + 20, 1420, col, (FONT_H3, FONT_BODY))
            draw_fitted_line(draw, t_desc, x1 + 96, ty + 68, 1420, TEXT_MAIN, (FONT_BODY_LG if is_cur else FONT_BODY, FONT_SM))
            ty += card_h + 22


def draw_scene_vector_geometry_canvas(
    draw: ImageDraw.ImageDraw,
    visual_mode: str,
    beat_text: str,
    sub_phase: float,
    t: float,
) -> None:
    """
    Scene Primitive 5: Full-Stage Animated Vector / Pocket Geometry Canvas.
    """
    x1, y1, x2, y2 = STAGE_BOX
    draw.rounded_rectangle([x1, y1, x2, y2], radius=24, fill=PANEL_BG, outline=EMERALD, width=3)

    if visual_mode == "multipocket_posebusters":
        draw_fitted_line(
            draw,
            "ROLLING A 1.4 Å WATER PROBE ACROSS A 5.0 Å VOXEL GRID TO DISCOVER pk1 & pk2",
            (x1 + x2) // 2,
            y1 + 36,
            1480,
            PRIMARY_BLUE,
            (FONT_H2, FONT_H3),
            center_x=True,
        )
        # Left: Single Blind Box in Empty Hinge vs Right: 5.0 Å Voxel Hash + 24-Ray SAS Probe
        lx, ly = x1 + 390, y1 + 390
        draw.rounded_rectangle([x1 + 56, y1 + 96, x1 + 740, y2 - 48], radius=18, fill=PANEL_DARK, outline=ACCENT_PINK, width=2)
        draw.text((x1 + 88, y1 + 120), "[X] Naive Single Center-of-Mass Box", fill=ACCENT_PINK, font=FONT_H3)
        draw.ellipse([lx - 180, ly - 110, lx - 25, ly + 110], fill=(28, 54, 90), outline=LIGHT_BLUE, width=3)
        draw.ellipse([lx + 25, ly - 110, lx + 180, ly + 110], fill=(56, 26, 76), outline=SECONDARY_PURPLE, width=3)
        draw.rectangle([lx - 65, ly - 65, lx + 65, ly + 65], outline=ACCENT_PINK, width=4)
        draw.text((lx - 145, ly + 135), "Lands in the empty inter-domain hinge!", fill=ACCENT_PINK, font=FONT_BODY)

        rx, ry = x1 + 1190, y1 + 390
        draw.rounded_rectangle([x1 + 780, y1 + 96, x2 - 56, y2 - 48], radius=18, fill=PANEL_DARK, outline=EMERALD, width=3)
        draw.text((x1 + 812, y1 + 120), "[OK] 5.0 Å Voxel Grid + 24-Ray SAS (pk1 & pk2)", fill=EMERALD, font=FONT_H3)
        for gx in range(rx - 220, rx + 240, 44):
            draw.line([(gx, ry - 130), (gx, ry + 130)], fill=(25, 36, 54), width=1)
        for gy in range(ry - 130, ry + 140, 44):
            draw.line([(rx - 220, gy), (rx + 220, gy)], fill=(25, 36, 54), width=1)
        draw.ellipse([rx - 185, ry - 110, rx - 25, ry + 110], fill=(20, 48, 82), outline=LIGHT_BLUE, width=3)
        draw.ellipse([rx + 25, ry - 110, rx + 185, ry + 110], fill=(48, 22, 68), outline=SECONDARY_PURPLE, width=3)

        # Animated 24-ray Fibonacci water probe sweeping from pk1 to pk2
        probe_x = int((rx - 105) * (1.0 - sub_phase) + (rx + 105) * sub_phase)
        probe_y = ry - 10
        for ray_k in range(12):
            ang = ray_k * (math.pi / 6.0) + t * 2.0
            draw.line([(probe_x, probe_y), (int(probe_x + 42 * math.cos(ang)), int(probe_y + 42 * math.sin(ang)))], fill=LIGHT_BLUE, width=2)
        draw.ellipse([rx - 130, ry - 35, rx - 80, ry + 15], fill=EMERALD, outline=(255, 255, 255), width=2)
        draw.text((rx - 175, ry + 32), "pk1: MBL Zn2+ Cleft", fill=EMERALD, font=FONT_SM)
        draw.ellipse([rx + 80, ry - 35, rx + 130, ry + 15], fill=AMBER, outline=(255, 255, 255), width=2)
        draw.text((rx + 45, ry + 32), "pk2: Rhodanese Loop", fill=AMBER, font=FONT_SM)
        draw.text((rx - 215, ry + 135), "Separates both active sites across 18.4 Å!", fill=EMERALD, font=FONT_BODY)

    elif visual_mode == "central_dogma_gist":
        # Ep 01 (Hello World 01): Progressive 3-Phase Animation of Soda-Lake Charge Balance (pI 4.38) & 1D Genomic Safe-Harbor Tree
        draw_fitted_line(
            draw,
            "HENDERSON-HASSELBALCH CHARGE ARMOR (pI 4.38) & 1D SAFE-HARBOR INTERVAL TREE",
            (x1 + x2) // 2,
            y1 + 32,
            1500,
            PRIMARY_BLUE,
            (FONT_H2, FONT_H3),
            center_x=True,
        )
        # Left Card: Animated pH Titration Curve & Asp-/Glu- Hydration Shell in 3.5 M NaCl
        draw.rounded_rectangle([x1 + 48, y1 + 92, x1 + 810, y2 - 42], radius=18, fill=PANEL_DARK, outline=ORANGE, width=3)
        draw.text((x1 + 76, y1 + 116), "1. 100-Step pKa Charge Bisection (pI = 4.38 vs. 6.40)", fill=ORANGE, font=FONT_H3)

        # Draw pH axis (pH 2 -> pH 11) and net charge curve Q(pH)
        ax_x1, ax_x2 = x1 + 96, x1 + 760
        ax_y0 = y1 + 290
        draw.line([(ax_x1, ax_y0), (ax_x2, ax_y0)], fill=TEXT_MUTED, width=2)
        draw.text((ax_x1, ax_y0 + 10), "pH 2.0", fill=TEXT_MUTED, font=FONT_XS)
        draw.text(((ax_x1 + ax_x2) // 2 - 40, ax_y0 + 10), "pI = 4.38 (Q = 0)", fill=EMERALD, font=FONT_SM)
        draw.text((ax_x2 - 80, ax_y0 + 10), "pH 10.5 (Brine)", fill=ORANGE, font=FONT_XS)

        # Animated bisection probe moving to pI = 4.38 and then showing -42e net charge at pH 10.5
        pts_curve = []
        for step_k in range(25):
            u = step_k / 24.0
            px_k = int(ax_x1 + u * (ax_x2 - ax_x1))
            # Sigmoidal titration curve crossing zero at u = 0.28 (pI = 4.38)
            py_k = int(ax_y0 - 95 * math.tanh((0.28 - u) * 3.2))
            pts_curve.append((px_k, py_k))
        for i in range(len(pts_curve) - 1):
            draw.line([pts_curve[i], pts_curve[i + 1]], fill=ORANGE, width=4)

        probe_idx = min(len(pts_curve) - 1, int(sub_phase * len(pts_curve)))
        pr_x, pr_y = pts_curve[probe_idx]
        draw.ellipse([pr_x - 10, pr_y - 10, pr_x + 10, pr_y + 10], fill=EMERALD, outline=(255, 255, 255), width=2)

        # Visual shell schematic at bottom of Left Card
        sh_y = y1 + 445
        draw.rounded_rectangle([x1 + 76, sh_y, x1 + 782, y2 - 68], radius=14, fill=PANEL_HIGHLIGHT, outline=EMERALD, width=2)
        draw.text((x1 + 100, sh_y + 20), "Why D3SWG7 Survives 3.5 M Salt at Lake Magadi (pH 10.5):", fill=EMERALD, font=FONT_BODY)
        draw.text((x1 + 100, sh_y + 62), "• Outer Asp-/Glu- carboxylates hold a pressurized H2O hydration cage", fill=TEXT_MAIN, font=FONT_SM)
        draw.text((x1 + 100, sh_y + 98), "• Negative GRAVY (-0.412) prevents hydrophobic clumping & sludge", fill=LIGHT_BLUE, font=FONT_SM)
        draw.text((x1 + 100, sh_y + 134), "• Buried catalytic triad stays locked within 0.44 Å!", fill=AMBER, font=FONT_SM)

        # Right Card: 1D Genomic Safe-Harbor Interval Tree (int4range GiST && and <->)
        draw.rounded_rectangle([x1 + 846, y1 + 92, x2 - 48, y2 - 42], radius=18, fill=PANEL_DARK, outline=SECONDARY_PURPLE, width=3)
        draw.text((x1 + 876, y1 + 116), "2. 1D Genomic Safe-Harbor GiST Tree (int4range)", fill=SECONDARY_PURPLE, font=FONT_H3)

        chr_y = y1 + 255
        draw.line([(x1 + 886, chr_y), (x2 - 88, chr_y)], fill=TEXT_MUTED, width=4)
        # Promoter interval
        draw.rounded_rectangle([x1 + 896, chr_y - 32, x1 + 1100, chr_y + 32], radius=10, fill=(20, 48, 82), outline=LIGHT_BLUE, width=2)
        draw.text((x1 + 914, chr_y - 10), "Promoter [8.0k..9.2k]", fill=LIGHT_BLUE, font=FONT_XS)
        # Safe Harbor candidate sliding into position
        sh_shift = int((1.0 - min(1.0, sub_phase * 1.4)) * 35)
        draw.rounded_rectangle([x1 + 1140 + sh_shift, chr_y - 36, x1 + 1355 + sh_shift, chr_y + 36], radius=10, fill=(6, 85, 64), outline=EMERALD, width=3)
        draw.text((x1 + 1156 + sh_shift, chr_y - 10), "Safe Harbor [10.5k..11.8k]", fill=EMERALD, font=FONT_XS)
        # Oncogene interval
        draw.rounded_rectangle([x1 + 1390, chr_y - 32, x2 - 84, chr_y + 32], radius=10, fill=(75, 20, 38), outline=ACCENT_PINK, width=2)
        draw.text((x1 + 1406, chr_y - 10), "Oncogene [15k+]", fill=ACCENT_PINK, font=FONT_XS)

        # Interval Operator Proof Box
        op_y = y1 + 360
        draw.rounded_rectangle([x1 + 876, op_y, x2 - 76, y2 - 68], radius=14, fill=PANEL_BG, outline=PRIMARY_BLUE, width=2)
        draw.text((x1 + 904, op_y + 24), "1D Spatial Interval Operators in SQL:", fill=PRIMARY_BLUE, font=FONT_BODY)
        draw.text((x1 + 904, op_y + 72), "1. Zero Oncogene Overlap:  NOT (harbor && oncogene)", fill=EMERALD, font=FONT_MONO_SM)
        draw.text((x1 + 904, op_y + 116), "2. Promoter Proximity:     (harbor <-> promoter) < 5000", fill=LIGHT_BLUE, font=FONT_MONO_SM)
        draw.text((x1 + 904, op_y + 160), "3. Central Dogma Compiler: bio_translate_dna(cds_seq)", fill=AMBER, font=FONT_MONO_SM)
        draw.text((x1 + 904, op_y + 206), "Scans millions of genomic intervals in O(log N) via GiST!", fill=TEXT_MAIN, font=FONT_SM)

    else:
        # EnzRef Vector Triangle & pSLiM Window Canvas
        draw_fitted_line(
            draw,
            "SHIFTING THE VECTOR ORIGIN TO THE ENZYME BACKGROUND STATE (v_bg)",
            (x1 + x2) // 2,
            y1 + 36,
            1480,
            PRIMARY_BLUE,
            (FONT_H2, FONT_H3),
            center_x=True,
        )
        ox, oy = x1 + 200, y2 - 125
        rx_pt, ry_pt = x1 + 540, y1 + 270
        # Product vector animates away from reactant as sub_phase progresses (kept well clear of label!)
        px_pt = int(rx_pt + 460 * min(1.0, 0.25 + sub_phase * 0.85))
        py_pt = int(ry_pt - 90 * min(1.0, 0.25 + sub_phase * 0.85))

        draw.line([(ox, oy), (rx_pt, ry_pt)], fill=LIGHT_BLUE, width=6)
        draw.line([(ox, oy), (px_pt, py_pt)], fill=SECONDARY_PURPLE, width=6)
        draw.line([(rx_pt, ry_pt), (px_pt, py_pt)], fill=EMERALD, width=8)

        draw.ellipse([ox - 16, oy - 16, ox + 16, oy + 16], fill=PRIMARY_BLUE, outline=(255, 255, 255), width=3)
        draw.text((ox - 30, oy + 26), "Enzyme Origin v_bg (1,280-D)", fill=PRIMARY_BLUE, font=FONT_H3)

        draw.ellipse([rx_pt - 14, ry_pt - 14, rx_pt + 14, ry_pt + 14], fill=LIGHT_BLUE, outline=(255, 255, 255), width=3)
        draw.text((rx_pt - 235, ry_pt - 42), "d_react = r - v_bg", fill=LIGHT_BLUE, font=FONT_H3)

        draw.ellipse([px_pt - 14, py_pt - 14, px_pt + 14, py_pt + 14], fill=SECONDARY_PURPLE, outline=(255, 255, 255), width=3)
        draw.text((px_pt + 22, py_pt - 16), "d_prod = p - v_bg", fill=SECONDARY_PURPLE, font=FONT_H3)

        mid_x, mid_y = (rx_pt + px_pt) // 2, (ry_pt + py_pt) // 2
        draw.text((mid_x - 150, mid_y - 68), "c_fwd = p - r (Turnover Vector!)", fill=EMERALD, font=FONT_H2)

        # Right callout box
        draw.rounded_rectangle([x1 + 980, y1 + 310, x2 - 64, y2 - 64], radius=18, fill=PANEL_DARK, outline=EMERALD, width=3)
        draw.text((x1 + 1012, y1 + 340), "2 Mathematical Invariants:", fill=AMBER, font=FONT_H3)
        draw.text((x1 + 1012, y1 + 398), "1. Dead Binder (p == r):", fill=ACCENT_PINK, font=FONT_BODY)
        draw.text((x1 + 1032, y1 + 434), "c_fwd = 0  -->  FluxRank = 0.0000", fill=TEXT_MAIN, font=FONT_MONO)
        draw.text((x1 + 1012, y1 + 492), "2. Reverse Reaction (r <-> p):", fill=EMERALD, font=FONT_BODY)
        draw.text((x1 + 1032, y1 + 528), "Flux(e, p, r) = -Flux(e, r, p)", fill=TEXT_MAIN, font=FONT_MONO)


def draw_scene_sql_focus_terminal(
    draw: ImageDraw.ImageDraw,
    cfg: Dict[str, Any],
    sub_phase: float,
    t: float,
) -> None:
    """
    Scene Primitive 6: Full-Width Declarative SQL Terminal (NO 3D spinner next to it!).
    Types the query and spotlights the key heuristic function call so the user can read the SQL effortlessly.
    """
    x1, y1, x2, y2 = STAGE_BOX
    draw.rounded_rectangle([x1, y1, x2, y2], radius=24, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=3)

    # macOS window traffic lights
    draw.ellipse([x1 + 36, y1 + 30, x1 + 58, y1 + 52], fill=ACCENT_PINK)
    draw.ellipse([x1 + 72, y1 + 30, x1 + 94, y1 + 52], fill=AMBER)
    draw.ellipse([x1 + 108, y1 + 30, x1 + 130, y1 + 52], fill=EMERALD)
    draw.text((x1 + 156, y1 + 28), "psql (PostgreSQL 18.3 — Expressing the Heuristic as a Declarative Query)", fill=TEXT_MUTED, font=FONT_MONO)
    draw.line([(x1, y1 + 74), (x2, y1 + 74)], fill=(45, 58, 82), width=2)

    sql_lines = cfg.get("sql_lines", [])
    n_vis = max(1, min(len(sql_lines), int(1 + sub_phase * len(sql_lines) * 1.15)))
    sy = y1 + 108
    for idx, line in enumerate(sql_lines[:n_vis]):
        is_latest = idx == n_vis - 1
        if is_latest:
            draw.rounded_rectangle([x1 + 36, sy - 6, x2 - 36, sy + 36], radius=8, fill=PANEL_HIGHLIGHT, outline=PRIMARY_BLUE, width=1)
        cursor = " █" if (is_latest and sub_phase < 0.75 and int(t * 6) % 2 == 0) else ""
        col = (
            AMBER
            if "bio_" in line or "z_order" in line or "kabsch" in line or "isoelectric" in line
            else SECONDARY_PURPLE
            if any(k in line for k in ("SELECT", "FROM", "WHERE", "WITH", "JOIN", "LATERAL"))
            else TEXT_MAIN
        )
        draw_fitted_line(draw, line + cursor, x1 + 56, sy, (x2 - x1) - 112, col, (FONT_MONO_LG, FONT_MONO, FONT_MONO_SM))
        sy += 44

    # Bottom Explanation Callout inside terminal
    draw.rounded_rectangle([x1 + 48, y2 - 175, x2 - 48, y2 - 36], radius=16, fill=PANEL_BG, outline=EMERALD, width=2)
    draw.text((x1 + 76, y2 - 152), "WHY THIS QUERY IS FAST & COMPOSABLE:", fill=EMERALD, font=FONT_H3)
    draw_wrapped_box(draw, cfg.get("cs_mechanism", ""), x1 + 76, y2 - 110, (x2 - x1) - 152, FONT_BODY, TEXT_MAIN, line_h=30, max_lines=2)


def draw_scene_query_result_plot(
    draw: ImageDraw.ImageDraw,
    visual_mode: str,
    cfg: Dict[str, Any],
    sub_phase: float,
    t: float,
) -> None:
    """
    Scene Primitive 7: Live Query Results Plotting onto an Animated Chart + Compact Result Table.
    Alternates with the SQL terminal so the viewer sees the query output plotted graphically!
    """
    x1, y1, x2, y2 = STAGE_BOX
    draw.rounded_rectangle([x1, y1, x2, y2], radius=24, fill=PANEL_BG, outline=EMERALD, width=3)

    draw_fitted_line(
        draw,
        "PLOTTING THE LIVE QUERY RESULTS: NAIVE BASELINE vs. HEURISTIC ENGINE",
        (x1 + x2) // 2,
        y1 + 32,
        1480,
        EMERALD,
        (FONT_H2, FONT_H3),
        center_x=True,
    )

    if visual_mode == "ipsae_minipae":
        bars = [
            ("80-aa SLiM Tether: Global ipTM (Diluted by 70-aa IDR)", 0.1250, 1.0, "0.1250 (False Negative!)", ACCENT_PINK),
            ("80-aa SLiM Tether: Contact-Only ipSAE_min (d0(N_i=10))", 0.8315, 1.0, "0.8315 (+0.7065 Rescue!)", EMERALD),
            ("5x5 MiniPAE Window: True 10-aa SLiM (0.475 Å / 12 Å Gate)", 0.96, 1.0, "0.475 Å (Win 41x31 -> PASS)", LIGHT_BLUE),
            ("5x5 MiniPAE Window: 1-Contact Random-Coil Decoy (25.08 Å)", 0.08, 1.0, "25.08 Å (Decoy -> REJECTED)", AMBER),
        ]
    elif visual_mode == "enzref_pslim":
        bars = [
            ("Forward Catalytic Reaction (A -> B): EnzRef FluxRank", 0.82, 1.0, "+0.6194 (Active Forward Turnover)", EMERALD),
            ("Non-Catalytic Dead Binder (p == r): EnzRef FluxRank", 0.02, 1.0, "0.0000 (100% Dead-Binder Zeroing)", ACCENT_PINK),
            ("Reverse Reaction (B -> A): Exact Antisymmetry", 0.82, 1.0, "-0.6194 (Sign Flips cleanly!)", LIGHT_BLUE),
            ("48-Hour RK4 Metabolon Cascade Yield (with pSLiM Tether)", 0.9916, 1.0, "99.16% Conversion Yield", AMBER),
        ]
    elif visual_mode == "multipocket_posebusters":
        bars = [
            ("Primary Pocket pk1 (MBL Zn2+ Cleft): Concavity & 21/21 Checks", 0.92, 1.0, "Score 28.45  |  21/21 PoseBusters PASS", EMERALD),
            ("Secondary Pocket pk2 (Rhodanese Cleft): Concavity & 21/21 Checks", 0.74, 1.0, "Score 18.90  |  21/21 PoseBusters PASS", LIGHT_BLUE),
            ("Salted Drug [Na+].Diclofenac: Non-Hiding Dual-Ledger Rescue", 0.95, 1.0, "passes_all_21=T | raw=F ('desalted')", AMBER),
            ("Bare Counterion [Na+] or Warped Ring (>0.25 Å SVD)", 0.04, 1.0, "0/21 -> REJECTED ('non_drug_counterion')", ACCENT_PINK),
        ]
    elif visual_mode == "zca_whitening":
        bars = [
            ("Raw 1,280-D ESM-2 Embeddings: SCOPe-40 Top-1 Superfamily", 0.5855, 1.0, "58.55% Recall (PR = 5.8 dims Cone)", ACCENT_PINK),
            ("Symmetric ZCA-Whitened Embeddings: SCOPe-40 Top-1 Superfamily", 0.8832, 1.0, "88.32% Recall (+29.77% Gain!)", EMERALD),
            ("1-Bit RaBitQ Error Bound Tightening (2.8108 -> 1.0118)", 0.86, 1.0, "2.78x Tighter Error Bound", LIGHT_BLUE),
            ("Zero-TOAST BioVecInt8G128 Inline Storage (1,304 B < 2,000 B)", 0.94, 1.0, "0 Bytes TOAST I/O (2.0x Faster SIMD)", AMBER),
    ]
    else:
        bars = [
            ("Catalytic Triad 3D Alignment (Kabsch SVD Core RMSD)", 0.96, 1.0, "0.44 Å Core Match (TM = 0.996)", EMERALD),
            ("9D RangeFinder Physicochemical Pocket Overlap (IoU)", 0.9551, 1.0, "95.51% Pocket Interval IoU", LIGHT_BLUE),
            ("Lake Magadi D3SWG7 Acidic Brine Armor (pI 4.38 vs 6.40)", 0.88, 1.0, "Survives 3.5 M NaCl & pH 10.5", ORANGE),
            ("21-Level Morton Octree Spatial Pruning Speedup", 0.986, 1.0, "74x Fewer Atoms Scanned (98.6% Pruned)", AMBER),
        ]

    grow = min(1.0, 0.25 + sub_phase * 0.95)
    by = y1 + 104
    bar_w_max = 1120
    for lbl, frac, _, val_txt, col in bars:
        draw.text((x1 + 72, by), lbl, fill=TEXT_MAIN, font=FONT_BODY)
        draw.rounded_rectangle([x1 + 72, by + 34, x1 + 72 + bar_w_max, by + 74], radius=10, fill=PANEL_DARK, outline=(45, 58, 78), width=1)
        fill_w = max(24, int(bar_w_max * frac * grow))
        draw.rounded_rectangle([x1 + 72, by + 34, x1 + 72 + fill_w, by + 74], radius=10, fill=col)
        draw.text((x1 + 92 + bar_w_max, by + 40), val_txt, fill=col, font=FONT_MONO)
        by += 112

    # Bottom SQL output summary line
    res_lines = cfg.get("sql_results", [])
    summary_txt = res_lines[-1] if res_lines else f"Benchmark Verified: {cfg.get('speedup', '')}"
    draw.rounded_rectangle([x1 + 72, y2 - 108, x2 - 72, y2 - 34], radius=14, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
    draw_fitted_line(draw, summary_txt, (x1 + x2) // 2, y2 - 84, 1380, EMERALD, (FONT_MONO_LG, FONT_MONO), center_x=True)


def draw_scene_csr_memory_anatomy(
    draw: ImageDraw.ImageDraw,
    visual_mode: str,
    beat_text: str,
    sub_phase: float,
    t: float,
) -> None:
    """
    Scene Primitive 8: PostgreSQL 8 KB Heap Page vs. 2,000 B TOAST Threshold & 16-Byte Header Peek.
    """
    x1, y1, x2, y2 = STAGE_BOX
    draw.rounded_rectangle([x1, y1, x2, y2], radius=24, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=3)

    draw_fitted_line(
        draw,
        "POSTGRESQL 8 KB HEAP PAGE ANATOMY: ZERO-TOAST INLINE STORAGE & 16-BYTE HEADER PEEK",
        (x1 + x2) // 2,
        y1 + 36,
        1520,
        SECONDARY_PURPLE,
        (FONT_H2, FONT_H3),
        center_x=True,
    )

    # Top comparison: 211,600 B Dense Matrix (106 TOAST chunks!) vs 1,090 B Inline "PGPA" CSR Tuple
    bx1 = x1 + 72
    bx2 = x2 - 72
    draw.rounded_rectangle([bx1, y1 + 104, bx2, y1 + 264], radius=16, fill=PANEL_DARK, outline=ACCENT_PINK, width=2)
    draw.text((bx1 + 28, y1 + 124), "[X] Dense float4[] PAE Matrix (211,600 Bytes >> 2,000 B TOAST Threshold):", fill=ACCENT_PINK, font=FONT_H3)
    for k in range(18):
        cx_k = bx1 + 28 + k * 80
        draw.rounded_rectangle([cx_k, y1 + 174, cx_k + 70, y1 + 234], radius=6, fill=(75, 22, 42), outline=ACCENT_PINK, width=1)
        draw.text((cx_k + 8, y1 + 195), f"TOAST#{k+1}", fill=TEXT_MAIN, font=FONT_XXS)

    # Bottom: 1,090 B Inline "PGPA" CSR inside 8,192 B Heap Page
    hy1 = y1 + 296
    draw.rounded_rectangle([bx1, hy1, bx2, y2 - 48], radius=18, fill=PANEL_DARK, outline=EMERALD, width=3)
    draw.text((bx1 + 28, hy1 + 24), "[OK] 'PGPA' Quantized CSR Container (1,090 Bytes < 2,000 B Threshold -> 100% Inline, 0 B TOAST!):", fill=EMERALD, font=FONT_H3)

    # 16-Byte Header Anatomy Diagram
    fields = [
        ("Bytes 0..4", "magic='PGPA'", PRIMARY_BLUE),
        ("Bytes 6..10", "n_a=150, n_b=80", LIGHT_BLUE),
        ("Bytes 10..12", "nnz=100 (<12 Å)", SECONDARY_PURPLE),
        ("Byte 13", "mini_pae=0.5 Å", AMBER),
        ("Byte 14", "ipsae=0.832", EMERALD),
        ("Bytes 16..1090", "u16 RowPtr + u8 PAE", PRIMARY_TEAL),
    ]
    fx = bx1 + 28
    active_f = int(sub_phase * len(fields)) % len(fields)
    for idx_f, (b_rng, b_val, col) in enumerate(fields):
        fw = 232 if idx_f < 5 else 260
        is_hi = idx_f == active_f
        draw.rounded_rectangle([fx, hy1 + 88, fx + fw - 10, hy1 + 210], radius=12, fill=PANEL_HIGHLIGHT if is_hi else PANEL_BG, outline=col, width=4 if is_hi else 2)
        draw.text((fx + 14, hy1 + 108), b_rng, fill=TEXT_MUTED, font=FONT_XS)
        draw_fitted_line(draw, b_val, fx + 14, hy1 + 148, fw - 34, col, (FONT_MONO, FONT_MONO_SM))
        fx += fw

    draw.text((bx1 + 28, hy1 + 245), "• O(1) 16-Byte Header Peek (bio_peek_pae_csr_header): Filters millions of complexes without unpacking CSR arrays!", fill=AMBER, font=FONT_BODY_LG)
    draw.text((bx1 + 28, hy1 + 295), "• 194.13x Storage Compression  |  26,316 Complexes/sec/core (46.37x Faster than Python/NumPy)", fill=EMERALD, font=FONT_H3)


def render_3d_protein_trace(
    draw: ImageDraw.ImageDraw,
    atoms: List[Dict[str, Any]],
    center_xy: Tuple[int, int],
    radius_px: float,
    angle_y: float,
    angle_x: float,
    base_color: Tuple[int, int, int],
    pocket_indices: List[int],
    highlight_mode: str,
    t: float,
    second_pocket_indices: Optional[List[int]] = None,
) -> List[Tuple[int, int]]:
    cx, cy = center_xy
    cos_y, sin_y = math.cos(angle_y), math.sin(angle_y)
    cos_x, sin_x = math.cos(angle_x), math.sin(angle_x)

    proj: List[Tuple[int, int, float, float, Dict[str, Any], int]] = []
    for idx, a in enumerate(atoms):
        x, y, z = a["nx"], a["ny"], a["nz"]
        rx = x * cos_y - z * sin_y
        rz = x * sin_y + z * cos_y
        ry = y * cos_x - rz * sin_x
        rz2 = y * sin_x + rz * cos_x
        persp = 1.0 / (1.0 + rz2 * 0.22)
        px = int(cx + rx * radius_px * persp)
        py = int(cy + ry * radius_px * persp)
        proj.append((px, py, rz2, persp, a, idx))

    for i in range(len(proj) - 1):
        px1, py1, z1, p1, a1, _ = proj[i]
        px2, py2, z2, p2, _, _ = proj[i + 1]
        depth_fade = max(0.35, min(1.0, 0.72 - 0.35 * ((z1 + z2) * 0.5)))
        if highlight_mode == "brine_armor" and a1["resn"] in ("ASP", "GLU"):
            seg_c = ORANGE
            w_line = 5
        else:
            seg_c = (
                int(base_color[0] * depth_fade),
                int(base_color[1] * depth_fade),
                int(base_color[2] * depth_fade),
            )
            w_line = max(2, int(4 * ((p1 + p2) * 0.5)))
        draw.line([(px1, py1), (px2, py2)], fill=seg_c, width=w_line)

    pocket_set = set(pocket_indices)
    pocket2_set = set(second_pocket_indices or [])
    pocket_pts: List[Tuple[int, int]] = []

    for px, py, rz2, persp, a, idx in proj:
        if idx in pocket_set:
            pocket_pts.append((px, py))

    for px, py, rz2, persp, a, idx in sorted(proj, key=lambda item: item[2], reverse=True):
        is_pk1 = idx in pocket_set
        is_pk2 = idx in pocket2_set
        if is_pk1:
            r = max(7, int(10 * persp))
            draw.ellipse([px - r - 4, py - r - 4, px + r + 4, py + r + 4], fill=(4, 78, 59))
            draw.ellipse([px - r, py - r, px + r, py + r], fill=EMERALD, outline=(255, 255, 255), width=2)
        elif is_pk2:
            r = max(7, int(10 * persp))
            draw.ellipse([px - r - 4, py - r - 4, px + r + 4, py + r + 4], fill=(120, 53, 15))
            draw.ellipse([px - r, py - r, px + r, py + r], fill=AMBER, outline=(255, 255, 255), width=2)
        else:
            r = max(2, int(3.5 * persp))
            depth_fade = max(0.35, min(1.0, 0.75 - 0.35 * rz2))
            c = (
                int(base_color[0] * depth_fade),
                int(base_color[1] * depth_fade),
                int(base_color[2] * depth_fade),
            )
            draw.ellipse([px - r, py - r, px + r, py + r], fill=c)

    if pocket_pts:
        pcx = sum(p[0] for p in pocket_pts) // len(pocket_pts)
        pcy = sum(p[1] for p in pocket_pts) // len(pocket_pts)
        pulse = int(4 * math.sin(t * 3.2))
        pr = 48 + pulse
        draw.ellipse([pcx - pr, pcy - pr, pcx + pr, pcy + pr], outline=EMERALD, width=3)

    return pocket_pts


def draw_scene_3d_organism_brief_moment(
    draw: ImageDraw.ImageDraw,
    pair_cfg: Dict[str, Any],
    window_progress: float,
    t: float,
) -> None:
    """
    Scene Primitive 9: Brief, Nausea-Free 3D Structural Proof (Strictly <= 15 seconds in the entire video!).
      - First 45% of this 15s window (~6.5 seconds): Slow, gentle 45-degree alignment glide (`-0.35 -> +0.15 rad`).
      - Remaining 55% of this 15s window (~8.5 seconds): Completely STATIONARY (`CAMERA LOCKED ON POCKET`)
        while spotlighting the conserved catalytic triad / 5x5 MiniPAE bridge!
    """
    x1, y1, x2, y2 = STAGE_BOX
    w = x2 - x1
    mid_x = (x1 + x2) // 2
    draw.rounded_rectangle([x1, y1, x2, y2], radius=24, fill=PANEL_BG, outline=EMERALD, width=3)

    pdb1 = pair_cfg.get("pdb1", "AF-A0ABD5UZ94-F1-model_v4_aligned.pdb")
    pdb2 = pair_cfg.get("pdb2", "AF-D3SWG7-F1-model_v4.pdb")
    atoms1 = load_pdb_ca_trace(pdb1)
    atoms2 = load_pdb_ca_trace(pdb2)

    # Gentle 45-degree slow ease over the first 45% of the brief window (~6 seconds), then 100% STILL!
    glide_p = min(1.0, window_progress / 0.45)
    eased = 1.0 - (1.0 - glide_p) ** 3
    angle_y = -0.32 + eased * 0.48
    angle_x = 0.12
    is_locked = window_progress >= 0.45

    header_title = pair_cfg.get("stage_header", "3D KABSCH-ALIGNED CATALYTIC POCKET PROOF")
    draw_fitted_line(draw, header_title, x1 + 36, y1 + 24, w - 340, PRIMARY_BLUE, (FONT_H2, FONT_H3, FONT_BODY))

    status_txt = "STATIONARY POCKET LOCK (0.0 rad/s)" if is_locked else "GENTLE KABSCH ALIGNMENT GLIDE"
    status_col = EMERALD if is_locked else AMBER
    sw = text_width(draw, status_txt, FONT_XS) + 28
    draw.rounded_rectangle([x2 - sw - 32, y1 + 20, x2 - 32, y1 + 58], radius=12, fill=(8, 48, 38) if is_locked else (48, 34, 10), outline=status_col, width=2)
    draw.text((x2 - sw - 18, y1 + 31), status_txt, fill=status_col, font=FONT_XS)

    vp_top = y1 + 74
    vp_bot = y2 - 86
    left_box = (x1 + 32, vp_top, mid_x - 16, vp_bot)
    right_box = (mid_x + 16, vp_top, x2 - 32, vp_bot)

    draw.rounded_rectangle(left_box, radius=16, fill=(9, 13, 22), outline=LIGHT_BLUE, width=2)
    draw.rounded_rectangle(right_box, radius=16, fill=(9, 13, 22), outline=ACCENT_PINK, width=2)

    draw.rectangle([left_box[0] + 2, left_box[1] + 2, left_box[2] - 2, left_box[1] + 44], fill=(20, 32, 52))
    draw_fitted_line(draw, pair_cfg.get("org1_label", "Reference Organism"), left_box[0] + 16, left_box[1] + 12, (left_box[2] - left_box[0]) - 32, LIGHT_BLUE, (FONT_SM, FONT_XS))

    draw.rectangle([right_box[0] + 2, right_box[1] + 2, right_box[2] - 2, right_box[1] + 44], fill=(42, 20, 46))
    draw_fitted_line(draw, pair_cfg.get("org2_label", "Extremophile Discovery"), right_box[0] + 16, right_box[1] + 12, (right_box[2] - right_box[0]) - 32, ACCENT_PINK, (FONT_SM, FONT_XS))

    pk1 = pair_cfg.get("pocket_indices", [22, 35, 48, 62, 78, 94])
    pk2 = pair_cfg.get("second_pocket_indices", [])
    c1 = ((left_box[0] + left_box[2]) // 2, (left_box[1] + left_box[3]) // 2 + 18)
    c2 = ((right_box[0] + right_box[2]) // 2, (right_box[1] + right_box[3]) // 2 + 18)
    rad = min((left_box[2] - left_box[0]), (left_box[3] - left_box[1])) * 0.39

    pts1 = render_3d_protein_trace(draw, atoms1, c1, rad, angle_y, angle_x, LIGHT_BLUE, pk1, "pocket", t, pk2)
    pts2 = render_3d_protein_trace(draw, atoms2, c2, rad, angle_y, angle_x, ACCENT_PINK, pk1, "pocket", t, pk2)

    if pts1 and pts2:
        p_left = pts1[0]
        p_right = pts2[0]
        draw.line([p_left, p_right], fill=EMERALD, width=3)
        mid_bx = (p_left[0] + p_right[0]) // 2
        mid_by = (p_left[1] + p_right[1]) // 2
        badge_txt = pair_cfg.get("bridge_badge", "0.44 Å Conserved Catalytic Triad")
        bw = text_width(draw, badge_txt, FONT_SM) + 28
        draw.rounded_rectangle([mid_bx - bw // 2, mid_by - 20, mid_bx + bw // 2, mid_by + 20], radius=12, fill=(6, 78, 59), outline=EMERALD, width=2)
        draw.text((mid_bx - bw // 2 + 14, mid_by - 10), badge_txt, fill=(255, 255, 255), font=FONT_SM)

    draw.rounded_rectangle([x1 + 32, vp_bot + 12, x2 - 32, y2 - 16], radius=12, fill=PANEL_DARK, outline=EMERALD, width=2)
    m_line = pair_cfg.get("metric_line1", "3D Core RMSD: 0.44 Å (TM = 0.996)  |  9D RangeFinder Pocket IoU: 95.51%")
    draw_fitted_line(draw, m_line, mid_x, vp_bot + 28, w - 110, EMERALD, (FONT_MONO, FONT_MONO_SM), center_x=True)


# ==============================================================================
# SPEECH-DRIVEN SCENE ROUTER (SELECTS 1 FULL-STAGE VISUAL PER PHRASE BEAT)
# ==============================================================================


def select_scene_for_beat(
    act_idx: int,
    beat_idx: int,
    total_beats_in_act: int,
    beat_text: str,
    local_t: float,
    act_duration_sec: float,
    visual_mode: str,
) -> Tuple[str, str]:
    """
    Routes each 3-4 second spoken phrase beat to the single most relevant full-stage visual
    primitive, ensuring:
      - The 3D organism view (`ORGANISM_3D_BRIEF`) only appears during a single <= 15-second
        window in Act IV when the narrator explicitly discusses the 3D organism comparison,
        and is hidden for the remaining 92% of the video!
      - Every beat alternates dynamically between concept spotlights, matrix zooms, Morton curves,
        equation derivations, vector geometry, SQL terminals, and plotted query charts.
    """
    low = beat_text.lower()

    # Strictly cap 3D organism view to a single <= 15 second window inside Act IV!
    if act_idx == 3:
        elapsed_in_act4 = local_t * act_duration_sec
        # Only show 3D organism comparison between 6s and 20s of Act IV (at most 14 seconds total!)
        if 6.0 <= elapsed_in_act4 <= 20.0:
            return ("ORGANISM_3D_BRIEF", "Act IV: Focused 14-Second 3D Catalytic Pocket Alignment (Camera Locked)")
        elif elapsed_in_act4 < 6.0:
            return ("SQL_FOCUS_TERMINAL", "Act IV: Expressing the Heuristic as a Declarative SQL Query")
        else:
            return ("QUERY_RESULT_PLOT", "Act IV: Plotting the Live In-Database Query Results")

    # Keyword-sensitive dynamic scene routing across Acts I, II, III, and V
    if any(k in low for k in ("morton", "z-order", "octree", "interleav", "64-bit", "b-tree", "cube", "sphere")):
        return ("MORTON_CURVE_DRAW", "How a 1D 64-Bit Morton Z-Order Curve Slices 3D Atomic Space")
    if any(k in low for k in ("toast", "1,090", "211,600", "csr", "header", "byte", "heap", "1,304", "pgpa")):
        return ("CSR_MEMORY_ANATOMY", "Inside the 8 KB PostgreSQL Heap Page: Zero-TOAST CSR & 16-Byte Header Peek")
    if any(k in low for k in ("matrix", "pae", "5-by-5", "5x5", "minipae", "diagonal", "window", "covariance")):
        return ("MATRIX_HEATMAP_ZOOM", "Zooming Into the Numerical Matrix & Contiguous 5x5 Window")
    if any(k in low for k in ("sql", "select ", "query", "lateral", "function")):
        return ("SQL_FOCUS_TERMINAL", "Declarative SQL Execution Inside Shared Buffers")
    if any(k in low for k in ("benchmark", "percent", "recall", "faster", "yield", "0.8315", "0.125", "25.08")):
        if beat_idx % 2 == 1:
            return ("QUERY_RESULT_PLOT", "Quantitative Benchmark Proof: Naive Baseline vs. Heuristic")

    # Structured alternation per Act so the visual changes every 3-5 seconds
    if act_idx == 0:
        cycle = ["CONCEPT_SPOTLIGHT", "MATRIX_HEATMAP_ZOOM", "CONCEPT_SPOTLIGHT", "VECTOR_GEOMETRY_CANVAS"]
        kind = cycle[beat_idx % len(cycle)]
        return (kind, "Act I: The Biological Paradox & Real-World Mystery")
    elif act_idx == 1:
        cycle = ["CONCEPT_SPOTLIGHT", "MATRIX_HEATMAP_ZOOM", "EQUATION_DERIVATION", "QUERY_RESULT_PLOT"]
        kind = cycle[beat_idx % len(cycle)]
        return (kind, "Act II: Why Naive Intuition & Whole-Chain Averages Fall Into a Trap")
    elif act_idx == 2:
        cycle = ["EQUATION_DERIVATION", "VECTOR_GEOMETRY_CANVAS", "MATRIX_HEATMAP_ZOOM", "CSR_MEMORY_ANATOMY"]
        if visual_mode == "morton_kabsch":
            cycle = ["MORTON_CURVE_DRAW", "EQUATION_DERIVATION", "MORTON_CURVE_DRAW", "QUERY_RESULT_PLOT"]
        kind = cycle[beat_idx % len(cycle)]
        return (kind, "Act III: Step-by-Step Mathematical & Physical Heuristic")
    else:
        cycle = ["QUERY_RESULT_PLOT", "CSR_MEMORY_ANATOMY", "SQL_FOCUS_TERMINAL", "CONCEPT_SPOTLIGHT"]
        kind = cycle[beat_idx % len(cycle)]
        return (kind, "Act V: From Mathematical Insight to Extremophile Discovery")


def render_frame(
    frame_idx: int,
    total_frames: int,
    act_idx: int,
    local_t: float,
    act_duration_sec: float,
    t: float,
    total_sec: float,
    cfg: Dict[str, Any],
    narration: str,
    beat_offset: int,
    total_motion_slides: int,
) -> Image.Image:
    progress = frame_idx / max(1, total_frames)
    img = Image.new("RGB", (W, H), BG_COLOR)
    draw = ImageDraw.Draw(img)

    # Subtle background grid
    for gx in range(0, W, 96):
        draw.line([(gx, 72), (gx, 1026)], fill=(17, 21, 30), width=1)
    for gy in range(72, 1026, 96):
        draw.line([(0, gy), (W, gy)], fill=(17, 21, 30), width=1)

    layer = cfg["layer"]
    visual_mode = cfg.get("visual_mode", "multipocket_posebusters")
    act_steps = cfg.get("act_steps", DEFAULT_ACT_STEPS)
    pair_cfg = cfg.get("organism_pair", {})

    beats = split_narration_beats(narration, max_words_per_beat=11)
    n_beats = max(1, len(beats))
    exact_beat_pos = local_t * n_beats
    beat_idx = min(n_beats - 1, int(exact_beat_pos))
    sub_phase = exact_beat_pos - beat_idx
    active_beat_text = beats[beat_idx]

    # Each phrase beat subdivides into 4 micro-states -> ~180-220 motion slides per video
    micro_step = min(3, int(sub_phase * 4))
    motion_slide_num = min(total_motion_slides, beat_offset + beat_idx * 4 + micro_step + 1)

    scene_kind, scene_headline = select_scene_for_beat(
        act_idx, beat_idx, n_beats, active_beat_text, local_t, act_duration_sec, visual_mode
    )

    # COMPONENT 1: Top Progress Bar & Single Headline
    draw_component_1_top_bar(
        draw,
        layer,
        act_idx,
        act_steps,
        progress,
        t,
        total_sec,
        motion_slide_num,
        total_motion_slides,
        scene_headline,
    )

    # COMPONENT 2: Full-Screen Focused Center Stage (ONLY ONE visual concept at a time!)
    if scene_kind == "CONCEPT_SPOTLIGHT":
        bullets = cfg.get("bullets", [])
        traps = cfg.get("trap_cards", [])
        takeaways = cfg.get(
            "takeaways",
            [
                ("Nature Conserves 3D Physics, Not 1D Spelling", cfg.get("bio_problem", "")),
                ("Good Heuristics Make False Positives Vanish", cfg.get("cs_mechanism", "")),
                ("Explore Interactive 3D Models in Your Browser", f"{BLOG_URL}/pg-bio-videos/"),
            ],
        )
        # Advance linearly across the act (`int(local_t * len(...))`) instead of `beat_idx % len(...)`
        # so consecutive CONCEPT_SPOTLIGHT scenes never repeat the same question!
        if act_idx == 0 and bullets:
            b_idx = min(len(bullets) - 1, int(local_t * len(bullets)))
            focus_title = bullets[b_idx].lstrip("• ").strip()
            accent = PRIMARY_BLUE
        elif act_idx == 1 and traps:
            t_idx = min(len(traps) - 1, int(local_t * len(traps)))
            focus_title = f"{traps[t_idx][0]}: {traps[t_idx][1]}"
            accent = ACCENT_PINK
        else:
            k_idx = min(len(takeaways) - 1, int(local_t * len(takeaways)))
            tk = takeaways[k_idx]
            focus_title = f"{tk[0]} — {tk[1]}"
            accent = EMERALD
        draw_scene_concept_spotlight(draw, active_beat_text, focus_title, cfg.get("speedup", ""), sub_phase, t, accent)

    elif scene_kind == "MATRIX_HEATMAP_ZOOM":
        draw_scene_matrix_heatmap_zoom(draw, visual_mode, active_beat_text, sub_phase, t)

    elif scene_kind == "MORTON_CURVE_DRAW":
        draw_scene_morton_curve_animation(draw, active_beat_text, sub_phase, t)

    elif scene_kind == "EQUATION_DERIVATION":
        draw_scene_equation_derivation(draw, visual_mode, active_beat_text, sub_phase, t)

    elif scene_kind == "VECTOR_GEOMETRY_CANVAS":
        draw_scene_vector_geometry_canvas(draw, visual_mode, active_beat_text, sub_phase, t)

    elif scene_kind == "SQL_FOCUS_TERMINAL":
        draw_scene_sql_focus_terminal(draw, cfg, sub_phase, t)

    elif scene_kind == "QUERY_RESULT_PLOT":
        draw_scene_query_result_plot(draw, visual_mode, cfg, sub_phase, t)

    elif scene_kind == "CSR_MEMORY_ANATOMY":
        draw_scene_csr_memory_anatomy(draw, visual_mode, active_beat_text, sub_phase, t)

    elif scene_kind == "ORGANISM_3D_BRIEF":
        elapsed_in_act4 = local_t * act_duration_sec
        win_p = max(0.0, min(1.0, (elapsed_in_act4 - 6.0) / 14.0))
        draw_scene_3d_organism_brief_moment(draw, pair_cfg, win_p, t)

    # COMPONENT 3: Spoken Phrase Spotlight Bar
    draw_component_3_caption_bar(draw, active_beat_text, beat_idx, n_beats, sub_phase)

    return img


def default_payload(args: argparse.Namespace) -> Dict[str, Any]:
    return {
        "title": args.title,
        "display_title": args.title,
        "layer": args.layer,
        "visual_mode": "multipocket_posebusters",
        "bio_problem": args.bio_problem,
        "cs_mechanism": args.cs_mechanism,
        "speedup": args.speedup,
        "subtitle": "Understanding the Physical & Mathematical Heuristics Behind the Discovery",
        "organism_pair": {
            "stage_header": "3D KABSCH ALIGNMENT: LAKE MAGADI SODA-LAKE ARCHAEON",
            "org1_label": "Bait: A0ABD5UZ94 (Reference Hydrolase · 312 aa)",
            "org2_label": "Discovery: D3SWG7 (Natrialba magadii · 0.44 Å · pH 10.5)",
            "pdb1": "AF-A0ABD5UZ94-F1-model_v4_aligned.pdb",
            "pdb2": "AF-D3SWG7-F1-model_v4.pdb",
            "pocket_indices": [25, 34, 42, 55, 70, 82],
            "bridge_badge": "0.44 Å Catalytic Triad Match",
            "metric_line1": "3D Core RMSD: 0.44 Å (TM = 0.996)  |  9D RangeFinder Pocket IoU: 95.51%",
            "metric_line2": "Orange = Acidic Asp/Glu Brine Shield (3.5 M NaCl)  |  Emerald = Conserved Catalytic Pocket",
        },
        "bullets": [
            "• Why do two organisms with 68% different sequences share a 0.44 Å pocket?",
            "• How does a 24-ray water-sphere roll across a protein to find hidden clefts?",
            "• Why does centering on the enzyme origin make dead binders vanish?",
            "• How do we spot a split hydrophobic core before building a fusion gene?",
        ],
        "rosetta_mappings": [
            [
                "Step 1: Roll a 1.4 Å Water Probe Across a 5.0 Å Voxel Grid",
                "Instead of guessing one docking box, cast 24 Fibonacci rays around every atom.",
                "Points surrounded by ≥11 protein neighbors sit inside a true concave cleft.",
            ],
            [
                "Step 2: Check Ring Flatness via 3×3 SVD & Steric Ratios",
                "Aromatic rings must stay planar (≤0.25 Å) and atoms cannot overlap (<0.75 VdW).",
                "Eliminates physical hallucinations while rescuing desalted vendor compounds.",
            ],
            [
                "Step 3: Compare the Conserved Pocket Across Extremophiles",
                "Align the discovered pocket against soda-lake and geothermal archaea in 3D.",
                "Reveals how evolution rewires surface charge while locking catalytic geometry.",
            ],
        ],
        "sql_lines": [
            "SELECT pocket_rank, pocket_id, round(p2rank_score::numeric, 2) AS concavity,",
            "       posebusters_valid, posebusters_passed_checks AS checks_21",
            "FROM bio_multipocket_posebusters_screen(",
            "    protein_pdb, 'O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl', 5",
            ");",
        ],
        "sql_results": [
            "┌─────────────┬───────────┬───────────┬───────────────────┬───────────┐",
            "│ pocket_rank │ pocket_id │ concavity │ posebusters_valid │ checks_21 │",
            "├─────────────┼───────────┼───────────┼───────────────────┼───────────┤",
            "│           1 │ ★ pk1     │     28.45 │ t                 │        21 │",
            "│           2 │ pk2       │     14.20 │ t                 │        21 │",
            "└─────────────┴───────────┴───────────┴───────────────────┴───────────┘",
            f"Verified across 120 complexes ({args.speedup})",
        ],
        "narrations": [
            f"Have you ever wondered how a protein from a boiling volcanic geyser or a pH 11 soda lake can run the exact same chemical reaction as a human enzyme, even when two thirds of its amino acids have mutated? Today we are exploring {args.title}.",
            f"Here is the trap that catches almost every naive algorithm: {args.bio_problem} When we rely on a single average score or a single bounding box, we miss the physical geometry that actually governs chemistry.",
            f"So what heuristic solves this? {args.cs_mechanism} By breaking the problem down into local 3D geometry and reference-centered vectors, the false positives drop away.",
            "Now let us look at the query and briefly bring in the 3D structures. As they glide into alignment and lock onto the active site, notice how the glowing emerald catalytic pocket in the center aligns to less than half an Angstrom.",
            "When we understand the heuristic first, expressing it as a query becomes effortless—and suddenly we can explore millions of extremophile organisms in seconds. Dive into the interactive 3D models at ideia dot me.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Render high-cadence micro-scene 1080p MP4 video for a pg_bio episode.")
    parser.add_argument("--out", required=True, help="Output MP4 file path")
    parser.add_argument("--title", default="Hello World, pg_bio!")
    parser.add_argument("--layer", default="Multiomics Stack")
    parser.add_argument("--bio-problem", default="")
    parser.add_argument("--cs-mechanism", default="")
    parser.add_argument("--speedup", default="")
    parser.add_argument("--config", default="", help="Optional JSON config path for full episode customization")
    parser.add_argument("--voice", default=DEFAULT_VOICE, help="Edge-TTS neural voice ID")
    args = parser.parse_args()

    cfg = default_payload(args)
    if args.config and Path(args.config).exists():
        loaded = json.loads(Path(args.config).read_text(encoding="utf-8"))
        cfg.update(loaded)

    out_path = Path(args.out).resolve()
    out_dir = out_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    ffmpeg_bin = shutil.which("ffmpeg") or "/opt/homebrew/bin/ffmpeg"
    ffprobe_bin = shutil.which("ffprobe") or "/opt/homebrew/bin/ffprobe"

    narrations = cfg["narrations"]
    print(f"Checking / synthesizing 5-Act Veritasium-style voiceover ({args.voice}) for {cfg['title']}...")
    full_vo_mp3, act_durations = generate_voiceover_and_srt(out_dir, narrations, args.voice, ffmpeg_bin, ffprobe_bin)

    act_frame_counts = [max(FPS * 3, int(math.ceil(d * FPS))) for d in act_durations]
    total_frames = sum(act_frame_counts)
    total_sec = total_frames / FPS

    act_beat_counts = [max(1, len(split_narration_beats(n, max_words_per_beat=11))) for n in narrations]
    total_motion_slides = sum(b * 4 for b in act_beat_counts)
    print(
        f"Voiceover ready ({total_sec:.1f}s across 5 acts, {total_motion_slides} micro-scene motion slides). "
        f"Rendering {total_frames} 1080p frames (max 3 components/frame, <=14s 3D window)..."
    )

    tmp_path = out_path.with_suffix(".tmp.mp4")
    cmd = [
        ffmpeg_bin,
        "-y",
        "-loglevel",
        "error",
        "-f",
        "rawvideo",
        "-vcodec",
        "rawvideo",
        "-s",
        f"{W}x{H}",
        "-pix_fmt",
        "rgb24",
        "-r",
        str(FPS),
        "-i",
        "-",
        "-i",
        str(full_vo_mp3),
        "-c:v",
        "libx264",
        "-profile:v",
        "main",
        "-level:v",
        "4.0",
        "-g",
        str(FPS * 2),
        "-pix_fmt",
        "yuv420p",
        "-preset",
        "fast",
        "-crf",
        "20",
        "-c:a",
        "aac",
        "-ar",
        "44100",
        "-ac",
        "2",
        "-b:a",
        "160k",
        "-movflags",
        "+faststart",
        "-shortest",
        str(tmp_path),
    ]

    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    assert proc.stdin is not None

    global_frame = 0
    broken_pipe = False
    beat_offset = 0
    for act_idx, n_frames in enumerate(act_frame_counts):
        if broken_pipe:
            break
        act_dur = act_durations[act_idx]
        for local_f in range(n_frames):
            local_t = local_f / max(1, n_frames - 1)
            t = global_frame / FPS
            img = render_frame(
                global_frame,
                total_frames,
                act_idx,
                local_t,
                act_dur,
                t,
                total_sec,
                cfg,
                narrations[act_idx],
                beat_offset,
                total_motion_slides,
            )
            try:
                proc.stdin.write(img.tobytes())
            except BrokenPipeError:
                broken_pipe = True
                break
            global_frame += 1
        beat_offset += act_beat_counts[act_idx] * 4

    try:
        proc.stdin.close()
    except BrokenPipeError:
        pass
    stderr_out = proc.stderr.read().decode("utf-8", errors="replace") if proc.stderr else ""
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"ffmpeg exited with {rc}: {stderr_out}")
    tmp_path.replace(out_path)
    print(f"Rendered 1080p MP4 ({total_motion_slides} motion slides, {total_frames} frames, {total_sec:.1f}s): {out_path}")


if __name__ == "__main__":
    main()
