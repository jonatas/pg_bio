#!/usr/bin/env python3
"""
Veritasium-Style 1080p MP4 Video + Neural Voiceover Renderer for pg_bio Episodes.

Key Pedagogical & Visual Design Principles:
  1. Multi-Slide, Low-Density, Speech-Synchronized Progression (18–20 Sub-Slides per Video):
       - Never cram 4 dense boxes onto one static screen!
       - Every Act automatically splits into 3–4 progressive sub-slides synchronized
         to the spoken voiceover sentences (`chunk_idx` / `total_chunks`).
       - On each sub-slide, ONLY the currently spoken bullet, trap, heuristic step,
         or takeaway is spotlighted in a spacious, high-contrast card with a live
         step indicator ("Slide 2 of 4 in Act II") and progressive checkmarks.
  2. Focused 3D Organism Motion ONLY When Explaining the 3D Alignment (Act IV):
       - Acts I, II, III, and V use clean, topic-specific 2D animated scientific diagrams
         (PAE matrices & 5x5 MiniPAE windows, reference-centered EnzRef vectors,
         24-ray SAS concavity probes, 64-bit Morton octrees, ZCA whitening cones,
         allosteric graphs, and pH charge balance) — NO distracting constant 3D spinning!
       - Act IV is the dedicated 3D Structural Proof moment: the Kabsch-aligned
         organisms rotate smoothly to show the 3D fold (`local_t < 0.55`) and then
         LOCK onto the aligned catalytic pocket / SLiM interface (`local_t >= 0.55`)
         while pulsing the exact matching residues.
  3. Zero Truncated Titles or Labels:
       - Uses multi-line wrapping or automatic font scaling so headers and labels
         never clip with `...`.
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
BG_COLOR = (15, 17, 23)            # --bg-color: #0f1117
PANEL_BG = (22, 27, 39)            # --card-bg-color: #161b27
PANEL_DARK = (26, 31, 46)          # --code-bg-color: #1a1f2e
PANEL_HIGHLIGHT = (30, 41, 68)     # Active spotlight card background
PRIMARY_BLUE = (92, 158, 255)      # --primary-color: #5c9eff
BRAND_BLUE = (67, 97, 238)         # #4361ee
SECONDARY_PURPLE = (177, 88, 252)  # --secondary-color: #b158fc
ACCENT_PINK = (255, 90, 157)       # --accent-color: #ff5a9d
PRIMARY_TEAL = (13, 148, 136)      # --primary-teal: #0d9488
EMERALD = (16, 185, 129)           # #10b981
LIGHT_BLUE = (56, 189, 248)        # #38bdf8
PG_TOPIC_BLUE = (51, 103, 145)     # #336791
TEXT_MAIN = (226, 232, 240)        # High-contrast readable text
TEXT_MUTED = (148, 163, 184)       # Muted slate
AMBER = (255, 193, 7)              # #FFC107
ORANGE = (249, 115, 22)            # #f97316

BLOG_URL = "https://ideia.me"
REPO_URL = "https://github.com/jonatas/pg_bio"
MODELS_DIR = Path("/Users/jonatas/code/ideia.me/assets/models")

DEFAULT_ACT_STEPS = [
    "I. The Paradox",
    "II. The Hidden Trap",
    "III. The Heuristic",
    "IV. 3D Pocket Focus",
    "V. Discovery Impact",
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


FONT_H1 = load_font(36)
FONT_H1_MD = load_font(30)
FONT_H2 = load_font(25)
FONT_H3 = load_font(22)
FONT_BODY_LG = load_font(24)
FONT_BODY = load_font(21)
FONT_SM = load_font(18)
FONT_XS = load_font(15)
FONT_XXS = load_font(13)
FONT_CAPTION = load_font(22)
FONT_MONO = load_font(18, mono=True)
FONT_MONO_SM = load_font(15, mono=True)

_PDB_CACHE: Dict[str, List[Dict[str, Any]]] = {}


def load_pdb_ca_trace(pdb_path_str: str, max_points: int = 160) -> List[Dict[str, Any]]:
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
                        plddt = float(line[60:66]) if len(line) >= 66 else 90.0
                        atoms.append({"resn": resn, "resi": resi, "x": x, "y": y, "z": z, "plddt": plddt})
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
                    "plddt": 94.0,
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


def draw_smart_title(draw: ImageDraw.ImageDraw, text: str, x: int, y: int, max_px: int, fill: Tuple[int, int, int]) -> None:
    for fnt, dy in ((FONT_H1, 0), (FONT_H1_MD, 3), (FONT_H2, 6), (FONT_H3, 8), (FONT_SM, 10)):
        if text_width(draw, text, fnt) <= max_px:
            draw.text((x, y + dy), text, fill=fill, font=fnt)
            return
    draw.text((x, y + 10), text, fill=fill, font=FONT_SM)


def draw_fitted_line(
    draw: ImageDraw.ImageDraw,
    text: str,
    x: int,
    y: int,
    max_px: int,
    fill: Tuple[int, int, int],
    preferred_fonts: Tuple[ImageFont.FreeTypeFont, ...] = (FONT_BODY, FONT_SM, FONT_XS, FONT_XXS),
) -> None:
    """Draws a single line of text by automatically stepping down font size instead of truncating with '...'."""
    for fnt in preferred_fonts:
        if text_width(draw, text, fnt) <= max_px:
            draw.text((x, y), text, fill=fill, font=fnt)
            return
    # Fallback: wrap onto 2 short lines in FONT_XXS so zero characters are lost
    draw_wrapped_box(draw, text, x, y - 4, max_px, FONT_XXS, fill, line_h=15, max_lines=2)


def fit_single_line(draw: ImageDraw.ImageDraw, text: str, max_px: int, font: ImageFont.FreeTypeFont) -> str:
    if text_width(draw, text, font) <= max_px:
        return text
    trimmed = text
    while len(trimmed) > 8 and text_width(draw, trimmed, font) > max_px:
        trimmed = trimmed[:-1]
    return trimmed


def draw_wrapped_box(
    draw: ImageDraw.ImageDraw,
    text: str,
    x: int,
    y: int,
    max_px: int,
    font: ImageFont.FreeTypeFont,
    fill: Tuple[int, int, int],
    line_h: int = 30,
    max_lines: int = 6,
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


def split_narration_chunks(narration: str, max_words_per_chunk: int = 15) -> List[str]:
    words = narration.split()
    if not words:
        return [""]
    chunks: List[str] = []
    cur: List[str] = []
    for w in words:
        cur.append(w)
        if len(cur) >= max_words_per_chunk or (len(cur) >= 8 and w.endswith((".", "!", "?", ";", ":"))):
            chunks.append(" ".join(cur))
            cur = []
    if cur:
        if chunks and len(cur) < 5:
            chunks[-1] = chunks[-1] + " " + " ".join(cur)
        else:
            chunks.append(" ".join(cur))
    return chunks


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
        chunks = split_narration_chunks(text, max_words_per_chunk=14)
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


def draw_header(
    draw: ImageDraw.ImageDraw,
    layer_str: str,
    act_idx: int,
    act_steps: List[str],
    progress: float,
    elapsed_sec: float,
    total_sec: float,
    subslide_idx: int,
    subslide_total: int,
) -> None:
    draw.rectangle([0, 0, W, 82], fill=(11, 13, 19))
    draw.line([(0, 82), (W, 82)], fill=(45, 55, 78), width=2)
    draw.rectangle([0, 79, int(W * progress), 82], fill=PRIMARY_BLUE)

    bw = min(460, text_width(draw, layer_str, FONT_XS) + 32)
    draw.rounded_rectangle([44, 18, 44 + bw, 60], radius=15, fill=PG_TOPIC_BLUE, outline=PRIMARY_BLUE, width=2)
    draw_fitted_line(draw, layer_str, 60, 29, bw - 28, (255, 255, 255), (FONT_XS, FONT_XXS))

    step_x = 530
    for idx, label in enumerate(act_steps[:5]):
        sw = text_width(draw, label, FONT_XS) + 22
        is_cur = idx == act_idx
        is_done = idx < act_idx
        fill_c = (37, 99, 235) if is_cur else PANEL_BG
        out_c = PRIMARY_BLUE if is_cur else PRIMARY_TEAL if is_done else (50, 60, 82)
        txt_c = (255, 255, 255) if is_cur else LIGHT_BLUE if is_done else TEXT_MUTED
        draw.rounded_rectangle([step_x, 19, step_x + sw, 59], radius=11, fill=fill_c, outline=out_c, width=2)
        draw.text((step_x + 11, 29), label, fill=txt_c, font=FONT_XS)
        step_x += sw + 8

    slide_badge = f"SLIDE {subslide_idx + 1}/{subslide_total}  |  {format_mm_ss(elapsed_sec)} / {format_mm_ss(total_sec)}"
    tw = text_width(draw, slide_badge, FONT_XS) + 28
    draw.rounded_rectangle([W - 44 - tw, 18, W - 44, 60], radius=15, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
    draw.text((W - 44 - tw + 14, 29), slide_badge, fill=TEXT_MAIN, font=FONT_XS)


def draw_caption_and_footer(draw: ImageDraw.ImageDraw, narration: str, act_local_progress: float, series_tag: str) -> Tuple[int, int, str]:
    chunks = split_narration_chunks(narration, max_words_per_chunk=15)
    c_idx = min(len(chunks) - 1, int(act_local_progress * len(chunks)))
    active_chunk = chunks[c_idx]

    draw.rounded_rectangle([96, 916, W - 96, 1012], radius=14, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
    dot_x = W - 125 - len(chunks) * 16
    for i in range(len(chunks)):
        c = PRIMARY_BLUE if i == c_idx else PRIMARY_TEAL if i < c_idx else (55, 65, 88)
        draw.ellipse([dot_x + i * 16, 926, dot_x + i * 16 + 8, 934], fill=c)

    wrapped = textwrap.wrap(active_chunk, width=88)[:2]
    cy = 934 if len(wrapped) > 1 else 950
    for line in wrapped:
        lw = text_width(draw, line, FONT_CAPTION)
        draw.text(((W - lw) // 2, cy), line, fill=TEXT_MAIN, font=FONT_CAPTION)
        cy += 30

    draw.rectangle([0, 1024, W, H], fill=(11, 13, 19))
    draw.line([(0, 1024), (W, 1024)], fill=BRAND_BLUE, width=2)

    clean_series = series_tag.replace("🧬 ", "")
    draw.text((52, 1040), f"Interactive 3D Studio: {BLOG_URL}/pg-bio-videos/", fill=PRIMARY_BLUE, font=FONT_MONO_SM)
    cw = text_width(draw, clean_series, FONT_MONO_SM)
    draw.text(((W - cw) // 2, 1040), clean_series, fill=TEXT_MUTED, font=FONT_MONO_SM)

    right_tag = f"GitHub: {REPO_URL}"
    rw = text_width(draw, right_tag, FONT_MONO_SM)
    draw.text((W - rw - 52, 1040), right_tag, fill=SECONDARY_PURPLE, font=FONT_MONO_SM)
    return c_idx, len(chunks), active_chunk


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
        px1, py1, z1, p1, a1, idx1 = proj[i]
        px2, py2, z2, p2, a2, idx2 = proj[i + 1]
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
    pocket2_pts: List[Tuple[int, int]] = []

    for px, py, rz2, persp, a, idx in proj:
        if idx in pocket_set:
            pocket_pts.append((px, py))
        if idx in pocket2_set:
            pocket2_pts.append((px, py))

    for px, py, rz2, persp, a, idx in sorted(proj, key=lambda item: item[2], reverse=True):
        is_pk1 = idx in pocket_set
        is_pk2 = idx in pocket2_set
        is_acid = a["resn"] in ("ASP", "GLU")

        if is_pk1:
            r = max(6, int(9 * persp))
            draw.ellipse([px - r - 4, py - r - 4, px + r + 4, py + r + 4], fill=(4, 78, 59))
            draw.ellipse([px - r, py - r, px + r, py + r], fill=EMERALD, outline=(255, 255, 255), width=2)
        elif is_pk2:
            r = max(6, int(9 * persp))
            draw.ellipse([px - r - 4, py - r - 4, px + r + 4, py + r + 4], fill=(120, 53, 15))
            draw.ellipse([px - r, py - r, px + r, py + r], fill=AMBER, outline=(255, 255, 255), width=2)
        elif highlight_mode == "brine_armor" and is_acid:
            r = max(4, int(6 * persp))
            draw.ellipse([px - r, py - r, px + r, py + r], fill=ORANGE)
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
        pulse = int(5 * math.sin(t * 3.8))
        pr = 44 + pulse
        draw.ellipse([pcx - pr, pcy - pr, pcx + pr, pcy + pr], outline=EMERALD, width=2)

    if pocket2_pts:
        pcx2 = sum(p[0] for p in pocket2_pts) // len(pocket2_pts)
        pcy2 = sum(p[1] for p in pocket2_pts) // len(pocket2_pts)
        pulse2 = int(5 * math.cos(t * 3.8))
        pr2 = 40 + pulse2
        draw.ellipse([pcx2 - pr2, pcy2 - pr2, pcx2 + pr2, pcy2 + pr2], outline=AMBER, width=2)

    return pocket_pts


def draw_dual_organism_comparison_stage(
    draw: ImageDraw.ImageDraw,
    box: Tuple[int, int, int, int],
    pair_cfg: Dict[str, Any],
    t: float,
    local_t: float,
    compact: bool = False,
) -> None:
    """
    Renders the focused 3D dual-organism comparison stage in Act IV ONLY.
    Notice the camera behavior:
      - Phase 1 (local_t < 0.55): Smooth 120-degree alignment sweep so the viewer sees the 3D fold.
      - Phase 2 (local_t >= 0.55): Camera LOCKS on the aligned catalytic pocket / SLiM groove
        (zero continuous dizzying rotation!) while pulsing the exact residue correspondence bridge.
    """
    x1, y1, x2, y2 = box
    w = x2 - x1
    mid_x = x1 + w // 2

    pdb1 = pair_cfg.get("pdb1", "AF-A0ABD5UZ94-F1-model_v4_aligned.pdb")
    pdb2 = pair_cfg.get("pdb2", "AF-D3SWG7-F1-model_v4.pdb")
    atoms1 = load_pdb_ca_trace(pdb1)
    atoms2 = load_pdb_ca_trace(pdb2)

    # Smooth eased rotation during first 55% of Act IV, then LOCKED camera focus!
    sweep_p = min(1.0, local_t / 0.55)
    eased_p = 1.0 - (1.0 - sweep_p) ** 3
    angle_y = -0.65 + eased_p * 1.45
    angle_x = 0.18 * (1.0 - eased_p) + 0.08
    is_locked = local_t >= 0.55

    vp_top = y1 + 60
    vp_bot = y2 - 88
    left_box = (x1 + 18, vp_top, mid_x - 10, vp_bot)
    right_box = (mid_x + 10, vp_top, x2 - 18, vp_bot)

    draw.rounded_rectangle(left_box, radius=14, fill=(9, 13, 22), outline=LIGHT_BLUE, width=2)
    draw.rounded_rectangle(right_box, radius=14, fill=(9, 13, 22), outline=ACCENT_PINK, width=2)

    header_title = pair_cfg.get("stage_header", "3D KABSCH-ALIGNED ORGANISM COMPARISON")
    draw_fitted_line(draw, header_title, x1 + 24, y1 + 18, w - 210, PRIMARY_BLUE, (FONT_H3, FONT_BODY, FONT_SM, FONT_XS))

    cam_status = "CAMERA LOCKED ON POCKET" if is_locked else "3D KABSCH ALIGNING..."
    cam_col = EMERALD if is_locked else AMBER
    cw = text_width(draw, cam_status, FONT_XXS) + 20
    draw.rounded_rectangle([x2 - cw - 18, y1 + 16, x2 - 18, y1 + 46], radius=8, fill=(10, 35, 30) if is_locked else (45, 30, 10), outline=cam_col, width=2)
    draw.text((x2 - cw - 8, y1 + 23), cam_status, fill=cam_col, font=FONT_XXS)

    # Organism Header Bars inside viewports — automatically sized so no '...' truncation occurs!
    draw.rectangle([left_box[0] + 2, left_box[1] + 2, left_box[2] - 2, left_box[1] + 44], fill=(22, 32, 52))
    draw_fitted_line(
        draw,
        pair_cfg.get("org1_label", "Bait Reference"),
        left_box[0] + 10,
        left_box[1] + 12,
        (left_box[2] - left_box[0]) - 20,
        LIGHT_BLUE,
        (FONT_XS, FONT_XXS),
    )

    draw.rectangle([right_box[0] + 2, right_box[1] + 2, right_box[2] - 2, right_box[1] + 44], fill=(42, 20, 46))
    draw_fitted_line(
        draw,
        pair_cfg.get("org2_label", "Extremophile Discovery"),
        right_box[0] + 10,
        right_box[1] + 12,
        (right_box[2] - right_box[0]) - 20,
        ACCENT_PINK,
        (FONT_XS, FONT_XXS),
    )

    pk1 = pair_cfg.get("pocket_indices", [22, 35, 48, 62, 78, 94])
    pk2 = pair_cfg.get("second_pocket_indices", [])
    hi_mode = "pocket" if is_locked else "brine_armor"

    c1 = ((left_box[0] + left_box[2]) // 2, (left_box[1] + left_box[3]) // 2 + 16)
    c2 = ((right_box[0] + right_box[2]) // 2, (right_box[1] + right_box[3]) // 2 + 16)
    rad = min((left_box[2] - left_box[0]), (left_box[3] - left_box[1])) * (0.40 if is_locked else 0.36)

    pts1 = render_3d_protein_trace(draw, atoms1, c1, rad, angle_y, angle_x, LIGHT_BLUE, pk1, hi_mode, t, pk2)
    pts2 = render_3d_protein_trace(draw, atoms2, c2, rad, angle_y, angle_x, ACCENT_PINK, pk1, hi_mode, t, pk2)

    if pts1 and pts2:
        active_bridge = 0 if is_locked else (int(t * 2.0) % min(len(pts1), len(pts2)))
        p_left = pts1[active_bridge]
        p_right = pts2[active_bridge]
        draw.line([p_left, p_right], fill=EMERALD, width=3 if is_locked else 2)
        mid_bx = (p_left[0] + p_right[0]) // 2
        mid_by = (p_left[1] + p_right[1]) // 2
        badge_txt = pair_cfg.get("bridge_badge", "0.44 Å Catalytic Match")
        bw = text_width(draw, badge_txt, FONT_XS) + 22
        draw.rounded_rectangle([mid_bx - bw // 2, mid_by - 17, mid_bx + bw // 2, mid_by + 17], radius=10, fill=(6, 78, 59), outline=EMERALD, width=2)
        draw.text((mid_bx - bw // 2 + 11, mid_by - 9), badge_txt, fill=(255, 255, 255), font=FONT_XS)

    strip_y = vp_bot + 10
    draw.rounded_rectangle([x1 + 18, strip_y, x2 - 18, y2 - 14], radius=12, fill=PANEL_BG, outline=EMERALD, width=2)
    m_line1 = pair_cfg.get("metric_line1", "3D Core RMSD: 0.44 Å (TM = 0.996)  |  9D Pocket Overlap: 95.51%")
    m_line2 = pair_cfg.get("metric_line2", "Orange = Acidic Brine Shield  |  Emerald = Locked 3D Catalytic Pocket")
    draw_fitted_line(draw, m_line1, x1 + 30, strip_y + 10, w - 60, EMERALD, (FONT_MONO_SM, FONT_XS, FONT_XXS))
    draw_fitted_line(draw, m_line2, x1 + 30, strip_y + 38, w - 60, TEXT_MAIN, (FONT_XS, FONT_XXS))


def draw_topic_2d_canvas(
    draw: ImageDraw.ImageDraw,
    visual_mode: str,
    box: Tuple[int, int, int, int],
    act_idx: int,
    step_idx: int,
    t: float,
    local_t: float,
    cfg: Dict[str, Any],
) -> None:
    """
    Renders a clean, low-density, animated 2D concept canvas on the right side of Acts I, II, III, and V
    that visually evolves with `step_idx` (the currently spoken bullet/step) instead of showing
    a spinning 3D protein all video long.
    """
    x1, y1, x2, y2 = box
    w = x2 - x1
    h = y2 - y1
    cx = (x1 + x2) // 2
    cy = (y1 + y2) // 2

    if visual_mode == "ipsae_minipae":
        # Visual Canvas for Layer 21: PAE Matrix (150 x 80), 10-aa SLiM vs 70-aa IDR Tail, 5x5 MiniPAE & 1,090 B PGPA CSR
        draw_fitted_line(
            draw,
            "VISUALIZING THE CO-FOLDING PAE MATRIX (150-aa Receptor x 80-aa Microprotein)",
            x1 + 28,
            y1 + 22,
            w - 56,
            PRIMARY_BLUE,
            (FONT_H3, FONT_BODY, FONT_SM),
        )

        # Draw 8x8 block schematic of the cross-chain PAE matrix (Receptor A rows vs Microprotein B cols)
        grid_x = x1 + 52
        grid_y = y1 + 86
        cell = 44
        draw.text((grid_x, grid_y - 24), "Chain B: 10-aa Docked SLiM [31..40]  +  70-aa Flexible IDR Linker", fill=TEXT_MUTED, font=FONT_XS)

        for r in range(8):
            for c in range(8):
                rx = grid_x + c * cell
                ry = grid_y + r * cell
                is_slim_diag = (2 <= r <= 4) and (1 <= c <= 3) and (r - 2 == c - 1)
                is_slim_box = (2 <= r <= 4) and (1 <= c <= 3)
                if is_slim_diag:
                    fill_c = (6, 95, 70)
                    out_c = EMERALD
                    txt = "0.45"
                elif is_slim_box:
                    fill_c = (15, 65, 55)
                    out_c = PRIMARY_TEAL
                    txt = "0.8"
                else:
                    fill_c = (48, 22, 36) if (act_idx in (0, 1) and step_idx == 0) else (24, 29, 42)
                    out_c = ACCENT_PINK if (act_idx in (0, 1) and step_idx == 0) else (45, 55, 75)
                    txt = "28.5"
                draw.rounded_rectangle([rx, ry, rx + cell - 4, ry + cell - 4], radius=6, fill=fill_c, outline=out_c, width=2)
                draw.text((rx + 6, ry + 13), txt, fill=(255, 255, 255) if is_slim_box else TEXT_MUTED, font=FONT_XXS)

        # Highlight the 5x5 Contiguous MiniPAE Window around cells (r=2..4, c=1..3)
        pulse = int(3 * math.sin(t * 4.0))
        win_x1 = grid_x + 1 * cell - 4 - pulse
        win_y1 = grid_y + 2 * cell - 4 - pulse
        win_x2 = grid_x + 4 * cell + pulse
        win_y2 = grid_y + 5 * cell + pulse
        draw.rounded_rectangle([win_x1, win_y1, win_x2, win_y2], radius=10, outline=AMBER, width=4)
        draw.text((grid_x, grid_y + 8 * cell + 10), "Amber Frame = 5x5 Contiguous MiniPAE Window (0.475 Å at Win 41x31)", fill=AMBER, font=FONT_XS)

        # Right column inside canvas: Equation & Score Meter synced to step_idx
        eq_x = grid_x + 8 * cell + 30
        eq_w = x2 - eq_x - 28

        # Score comparison card 1: Whole-Chain ipTM (Trap)
        c1_active = (act_idx in (0, 1) and step_idx in (0, 1))
        draw.rounded_rectangle([eq_x, grid_y, x2 - 28, grid_y + 160], radius=12, fill=PANEL_HIGHLIGHT if c1_active else PANEL_BG, outline=ACCENT_PINK, width=3 if c1_active else 2)
        draw.text((eq_x + 18, grid_y + 14), "[X] Global ipTM (Divides by L_B = 80 aa)", fill=ACCENT_PINK, font=FONT_SM)
        draw.text((eq_x + 18, grid_y + 48), "Denominator includes 70-aa IDR tail (28.5 Å)", fill=TEXT_MAIN, font=FONT_XS)
        draw.rounded_rectangle([eq_x + 18, grid_y + 88, eq_x + 18 + int((eq_w - 36) * 0.15), grid_y + 116], radius=6, fill=ACCENT_PINK)
        draw.text((eq_x + 18, grid_y + 124), "ipTM = 0.1250  (FALSE NEGATIVE!)", fill=ACCENT_PINK, font=FONT_MONO_SM)

        # Score comparison card 2: Contact-Restricted ipSAE (Solution)
        c2_active = (act_idx >= 2 or step_idx >= 2)
        cy2 = grid_y + 178
        draw.rounded_rectangle([eq_x, cy2, x2 - 28, cy2 + 174], radius=12, fill=PANEL_HIGHLIGHT if c2_active else PANEL_BG, outline=EMERALD, width=3 if c2_active else 2)
        draw.text((eq_x + 18, cy2 + 14), "[OK] Contact-Only ipSAE (PAE_ij < 10.0 Å)", fill=EMERALD, font=FONT_SM)
        draw.text((eq_x + 18, cy2 + 46), "d0(N_i) uses ONLY N_i = 10 contacting residues!", fill=TEXT_MAIN, font=FONT_XS)
        bar_fill = min(1.0, 0.8315)
        draw.rounded_rectangle([eq_x + 18, cy2 + 86, eq_x + 18 + int((eq_w - 36) * bar_fill), cy2 + 114], radius=6, fill=EMERALD)
        draw.text((eq_x + 18, cy2 + 126), "ipSAE_min = 0.8315 (+0.7065 IDR Rescue!)", fill=EMERALD, font=FONT_MONO_SM)

        # Bottom Bar: 5x5 MiniPAE vs Decoy & 1,090 B Zero-TOAST "PGPA" CSR Storage
        bot_y = grid_y + 8 * cell + 52
        draw.rounded_rectangle([x1 + 28, bot_y, x2 - 28, y2 - 24], radius=14, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
        draw.text((x1 + 48, bot_y + 16), "STORAGE & DECOY GATE: ZERO-TOAST 'PGPA' QUANTIZED CSR CONTAINER", fill=SECONDARY_PURPLE, font=FONT_SM)
        draw_fitted_line(
            draw,
            "• 5x5 MiniPAE Window: 0.475 Å (True 10-aa SLiM)  vs.  25.08 Å (1-Contact Random-Coil Decoy -> REJECTED)",
            x1 + 48,
            bot_y + 54,
            w - 96,
            EMERALD,
            (FONT_SM, FONT_XS),
        )
        draw_fitted_line(
            draw,
            "• PGPA Inline CSR: 211,600 B Dense float4[]  -->  1,090 B Inline Heap Tuple (194.13x Smaller, 0 B TOAST, 16B O(1) Header Peek)",
            x1 + 48,
            bot_y + 92,
            w - 96,
            LIGHT_BLUE,
            (FONT_SM, FONT_XS),
        )

    elif visual_mode == "enzref_pslim":
        # Visual Canvas for Layer 20: EnzRef Origin-Centered Vector Geometry & pSLiM 8-aa Window Linker Gate
        draw_fitted_line(
            draw,
            "ENZREF VECTOR GEOMETRY (WHY DEAD BINDERS COLLAPSE TO ZERO)",
            x1 + 28,
            y1 + 22,
            w - 56,
            PRIMARY_BLUE,
            (FONT_H3, FONT_BODY, FONT_SM),
        )
        top_box = (x1 + 28, y1 + 68, x2 - 28, y1 + 410)
        draw.rounded_rectangle(top_box, radius=14, fill=(11, 15, 25), outline=PRIMARY_BLUE, width=2)

        ox, oy = x1 + 135, y1 + 335
        rx_pt, ry_pt = x1 + 370, y1 + 170
        px_pt, py_pt = x1 + 690, y1 + 135
        mx_pt, my_pt = (rx_pt + px_pt) // 2, (ry_pt + py_pt) // 2

        draw.line([(ox, oy), (rx_pt, ry_pt)], fill=LIGHT_BLUE, width=4)
        draw.line([(ox, oy), (px_pt, py_pt)], fill=SECONDARY_PURPLE, width=4)
        draw.line([(ox, oy), (mx_pt, my_pt)], fill=AMBER, width=3)
        draw.line([(rx_pt, ry_pt), (px_pt, py_pt)], fill=EMERALD, width=6)

        draw.ellipse([ox - 12, oy - 12, ox + 12, oy + 12], fill=PRIMARY_BLUE, outline=(255, 255, 255), width=2)
        draw.text((ox - 30, oy + 18), "Enzyme Origin v_bg", fill=PRIMARY_BLUE, font=FONT_SM)

        draw.ellipse([rx_pt - 10, ry_pt - 10, rx_pt + 10, ry_pt + 10], fill=LIGHT_BLUE, outline=(255, 255, 255), width=2)
        draw.text((rx_pt - 140, ry_pt - 34), "d_react = r - v_bg", fill=LIGHT_BLUE, font=FONT_SM)

        draw.ellipse([px_pt - 10, py_pt - 10, px_pt + 10, py_pt + 10], fill=SECONDARY_PURPLE, outline=(255, 255, 255), width=2)
        draw.text((px_pt - 40, py_pt - 34), "d_prod = p - v_bg", fill=SECONDARY_PURPLE, font=FONT_SM)

        draw.text((mx_pt - 85, my_pt - 42), "c_fwd = p - r (Catalytic Vector)", fill=EMERALD, font=FONT_BODY)
        draw.text((x1 + 310, y1 + 295), "• Dead Binder (p == r)  -->  c_fwd = 0  -->  FluxRank = 0.0000", fill=ACCENT_PINK, font=FONT_SM)
        draw.text((x1 + 310, y1 + 335), "• Reverse Reaction (p <-> r)  -->  Flux(e,p,r) = -Flux(e,r,p)", fill=EMERALD, font=FONT_SM)

        bot_box = (x1 + 28, y1 + 430, x2 - 28, y2 - 24)
        draw.rounded_rectangle(bot_box, radius=14, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
        draw.text((x1 + 48, y1 + 448), "DISORDERED MICROPROTEIN pSLiM 8-AA WINDOW & LINKER SANITY GATE", fill=SECONDARY_PURPLE, font=FONT_SM)

        bar_y = y1 + 505
        draw.rounded_rectangle([x1 + 48, bar_y, x1 + 260, bar_y + 56], radius=10, fill=(6, 78, 59), outline=EMERALD, width=2)
        draw.text((x1 + 66, bar_y + 17), "SLiM A [1..8] (M1)", fill=EMERALD, font=FONT_SM)

        draw.rounded_rectangle([x1 + 268, bar_y + 8, x1 + 585, bar_y + 48], radius=8, fill=(30, 41, 59), outline=LIGHT_BLUE, width=2)
        draw.text((x1 + 284, bar_y + 17), "Linker: GGGSGGGGS (GRAVY -0.4)", fill=LIGHT_BLUE, font=FONT_XS)

        draw.rounded_rectangle([x1 + 593, bar_y, x2 - 48, bar_y + 56], radius=10, fill=(76, 29, 149), outline=SECONDARY_PURPLE, width=2)
        draw.text((x1 + 610, bar_y + 17), "SLiM B [28..35] (M2)", fill=SECONDARY_PURPLE, font=FONT_SM)

        draw_fitted_line(
            draw,
            "• 8-aa ZCA-Whitened Windows (||w||_2 = 1.0) prevent 10x whole-protein motif dilution (617,845 windows/sec)",
            x1 + 48,
            bar_y + 82,
            w - 96,
            TEXT_MAIN,
            (FONT_SM, FONT_XS),
        )
        draw_fitted_line(
            draw,
            "• Linker Gate: GRAVY <= 0.45 & HydroRun < 5 -> Single_Polypeptide_Fusion (else Bicistronic_Operon_Required)",
            x1 + 48,
            bar_y + 118,
            w - 96,
            EMERALD,
            (FONT_SM, FONT_XS),
        )

    elif visual_mode == "multipocket_posebusters":
        # Visual Canvas for Layer 19: Multi-Pocket SAS Voxel Hashing & 21-Check PoseBusters Dual-Ledger Gate
        draw_fitted_line(
            draw,
            "24-RAY FIBONACCI SAS CONCAVITY vs. SINGLE BLIND DOCKING BOX",
            x1 + 28,
            y1 + 22,
            w - 56,
            PRIMARY_BLUE,
            (FONT_H3, FONT_BODY, FONT_SM),
        )
        top_box = (x1 + 28, y1 + 68, x2 - 28, y1 + 365)
        draw.rounded_rectangle(top_box, radius=14, fill=(11, 15, 25), outline=PRIMARY_BLUE, width=2)

        lx, ly = x1 + 210, y1 + 215
        draw.text((x1 + 48, y1 + 82), "[X] Single Center-of-Mass Box", fill=ACCENT_PINK, font=FONT_SM)
        draw.ellipse([lx - 115, ly - 75, lx - 15, ly + 75], fill=(30, 58, 95), outline=LIGHT_BLUE, width=2)
        draw.ellipse([lx + 15, ly - 75, lx + 115, ly + 75], fill=(58, 28, 75), outline=SECONDARY_PURPLE, width=2)
        draw.rectangle([lx - 45, ly - 45, lx + 45, ly + 45], outline=ACCENT_PINK, width=3)
        draw.text((lx - 78, ly + 84), "Lands in empty hinge!", fill=ACCENT_PINK, font=FONT_XS)

        rx, ry = x1 + 610, y1 + 215
        draw.text((x1 + 430, y1 + 82), "[OK] 5.0 Å Voxel Grid + 24-Ray SAS (pk1 & pk2)", fill=EMERALD, font=FONT_SM)
        for gx in range(rx - 140, rx + 150, 35):
            draw.line([(gx, ry - 80), (gx, ry + 80)], fill=(25, 35, 52), width=1)
        for gy in range(ry - 80, ry + 90, 35):
            draw.line([(rx - 140, gy), (rx + 140, gy)], fill=(25, 35, 52), width=1)
        draw.ellipse([rx - 125, ry - 75, rx - 15, ry + 75], fill=(20, 48, 82), outline=LIGHT_BLUE, width=2)
        draw.ellipse([rx + 15, ry - 75, rx + 125, ry + 75], fill=(48, 22, 68), outline=SECONDARY_PURPLE, width=2)
        pulse = int(4 * math.sin(t * 4.0))
        draw.ellipse([rx - 85 - pulse, ry - 25 - pulse, rx - 45 + pulse, ry + 15 + pulse], fill=EMERALD, outline=(255, 255, 255), width=2)
        draw.text((rx - 115, ry + 28), "pk1 (MBL Zn2+)", fill=EMERALD, font=FONT_XS)
        draw.ellipse([rx + 45 - pulse, ry - 20 - pulse, rx + 85 + pulse, ry + 20 + pulse], fill=AMBER, outline=(255, 255, 255), width=2)
        draw.text((rx + 25, ry + 28), "pk2 (Rhodanese)", fill=AMBER, font=FONT_XS)

        bot_box = (x1 + 28, y1 + 382, x2 - 28, y2 - 24)
        draw.rounded_rectangle(bot_box, radius=14, fill=PANEL_BG, outline=EMERALD, width=2)
        draw.text((x1 + 48, y1 + 398), "21-CHECK POSEBUSTERS & NON-HIDING DUAL-LEDGER SALT RESCUE", fill=EMERALD, font=FONT_SM)

        checks = [
            ("1. Vendor Salt Stripping", "[Na+].Diclofenac -> Active Drug (>=3 heavy atoms; bare [Na+] rejected)", EMERALD),
            ("2. 3x3 Jacobi SVD Ring Planarity", "Best-fit aromatic plane RMSD <= 0.25 Å (catches warped rings)", LIGHT_BLUE),
            ("3. Pairwise VdW Steric Ratio", "d_ij / (r_vdW,i + r_vdW,j) >= 0.75 (blocks steric clash hallucinations)", AMBER),
            ("4. Non-Hiding Dual Ledger", "passes_all_21 = TRUE  |  raw_passes_all_21 = FALSE  |  'desalted_from_salt'", SECONDARY_PURPLE),
        ]
        cy_c = y1 + 438
        for idx_c, (c_title, c_desc, col) in enumerate(checks):
            is_spot = (idx_c == step_idx % len(checks))
            draw.rounded_rectangle([x1 + 48, cy_c, x2 - 48, cy_c + 56], radius=10, fill=PANEL_HIGHLIGHT if is_spot else PANEL_DARK, outline=col, width=3 if is_spot else 1)
            draw.text((x1 + 64, cy_c + 7), c_title, fill=col, font=FONT_SM)
            draw_fitted_line(draw, c_desc, x1 + 64, cy_c + 32, w - 128, TEXT_MAIN, (FONT_XS, FONT_XXS))
            cy_c += 64

    elif visual_mode == "morton_kabsch":
        # Visual Canvas for Ep 2: 64-Bit Morton Z-Order Curve, 21-Level Octree Pruning & 3x3 Jacobi Kabsch SVD
        draw_fitted_line(
            draw,
            "64-BIT MORTON Z-ORDER OCTREE PRUNING & 3x3 KABSCH SVD",
            x1 + 28,
            y1 + 22,
            w - 56,
            PRIMARY_BLUE,
            (FONT_H3, FONT_BODY, FONT_SM),
        )
        top_box = (x1 + 28, y1 + 68, x2 - 28, y1 + 420)
        draw.rounded_rectangle(top_box, radius=14, fill=(11, 15, 25), outline=PRIMARY_BLUE, width=2)

        # Draw 4x4 Octree / Quadtree grid with Z-order curve & 5.0 Å query sphere
        gx0, gy0 = x1 + 68, y1 + 100
        g_cell = 68
        z_pts = []
        for r in range(4):
            for c in range(4):
                rx, ry = gx0 + c * g_cell, gy0 + r * g_cell
                inside = (1 <= r <= 2) and (1 <= c <= 2)
                border = (0 <= r <= 3) and (0 <= c <= 3) and not inside and ((r in (1, 2)) or (c in (1, 2)))
                fc = (6, 78, 59) if inside else (45, 35, 15) if border else (18, 22, 33)
                oc = EMERALD if inside else AMBER if border else (40, 50, 70)
                draw.rounded_rectangle([rx, ry, rx + g_cell - 4, ry + g_cell - 4], radius=6, fill=fc, outline=oc, width=2)
                z_pts.append((rx + g_cell // 2, ry + g_cell // 2))
        for i in range(len(z_pts) - 1):
            draw.line([z_pts[i], z_pts[i + 1]], fill=LIGHT_BLUE, width=2)
        scx, scy = gx0 + 2 * g_cell, gy0 + 2 * g_cell
        draw.ellipse([scx - 108, scy - 108, scx + 108, scy + 108], outline=ACCENT_PINK, width=3)

        info_x = gx0 + 4 * g_cell + 36
        draw.text((info_x, y1 + 104), "Bit-Interleaving (X, Y, Z) -> 64-Bit BIGINT:", fill=LIGHT_BLUE, font=FONT_SM)
        draw.text((info_x, y1 + 140), "X = 101   Y = 011   Z = 110", fill=TEXT_MAIN, font=FONT_MONO)
        draw.text((info_x, y1 + 174), "Morton Key = [Z2 Y2 X2 | Z1 Y1 X1 | Z0 Y0 X0]", fill=EMERALD, font=FONT_MONO_SM)
        draw.text((info_x, y1 + 224), "• Emerald Cubes (fully_inside = TRUE):", fill=EMERALD, font=FONT_SM)
        draw.text((info_x + 18, y1 + 254), "Skip sqrt(dx^2 + dy^2 + dz^2) completely!", fill=TEXT_MAIN, font=FONT_XS)
        draw.text((info_x, y1 + 292), "• Dark Cubes (Disjoint Octants):", fill=TEXT_MUTED, font=FONT_SM)
        draw.text((info_x + 18, y1 + 322), "Pruned in O(log N) -> 74x fewer atoms scanned!", fill=AMBER, font=FONT_XS)

        bot_box = (x1 + 28, y1 + 440, x2 - 28, y2 - 24)
        draw.rounded_rectangle(bot_box, radius=14, fill=PANEL_BG, outline=EMERALD, width=2)
        draw.text((x1 + 48, y1 + 460), "3x3 JACOBI SVD KABSCH ALIGNMENT: H = P^T Q  -->  R = V * diag(1, 1, det(VU^T)) * U^T", fill=EMERALD, font=FONT_SM)
        draw_fitted_line(
            draw,
            "• Checking sign of det(V U^T) prevents mirror-image D-amino acid reflection (chirality bug)!",
            x1 + 48,
            y1 + 510,
            w - 96,
            TEXT_MAIN,
            (FONT_SM, FONT_XS),
        )
        draw_fitted_line(
            draw,
            "• Aligns catalytic triads across 2.34M proteins in microseconds inside PostgreSQL shared buffers.",
            x1 + 48,
            y1 + 555,
            w - 96,
            LIGHT_BLUE,
            (FONT_SM, FONT_XS),
        )

    elif visual_mode == "zca_whitening":
        # Visual Canvas for Ep 4: The 15-D Anisotropic SwiGLU Cone vs Full Symmetric ZCA Whitening Sphere
        draw_fitted_line(
            draw,
            "UNTANGLING THE ANISOTROPIC CONE WITH SYMMETRIC ZCA WHITENING",
            x1 + 28,
            y1 + 22,
            w - 56,
            PRIMARY_BLUE,
            (FONT_H3, FONT_BODY, FONT_SM),
        )
        top_box = (x1 + 28, y1 + 68, x2 - 28, y1 + 430)
        draw.rounded_rectangle(top_box, radius=14, fill=(11, 15, 25), outline=PRIMARY_BLUE, width=2)

        # Left: Narrow Anisotropic Cone (PR = 5.8 / 1280 dims, mean cos = 0.94)
        lx, ly = x1 + 215, y1 + 265
        draw.text((x1 + 50, y1 + 86), "[X] Raw 1,280-D ESM-2 Cone (PR = 5.8 dims)", fill=ACCENT_PINK, font=FONT_SM)
        draw.polygon([(lx - 130, ly + 95), (lx + 95, ly - 95), (lx + 135, ly - 45)], fill=(65, 20, 42), outline=ACCENT_PINK)
        for dx, dy in [(-20, 20), (15, -18), (45, -42), (70, -65), (90, -78), (30, -25)]:
            draw.ellipse([lx + dx - 5, ly + dy - 5, lx + dx + 5, ly + dy + 5], fill=AMBER)
        draw.text((x1 + 50, y1 + 382), "Unrelated proteins score 0.94 cosine!", fill=ACCENT_PINK, font=FONT_XS)

        # Right: Isotropic ZCA-Whitened Unit Sphere (W = V * Lambda^{-1/2} * V^T)
        rx, ry = x1 + 625, y1 + 255
        draw.text((x1 + 445, y1 + 86), "[OK] ZCA-Whitened Sphere (W = V Lambda^-1/2 V^T)", fill=EMERALD, font=FONT_SM)
        draw.ellipse([rx - 110, ry - 110, rx + 110, ry + 110], fill=(10, 48, 42), outline=EMERALD, width=3)
        for ang in range(0, 360, 45):
            rad_a = math.radians(ang)
            px = int(rx + 82 * math.cos(rad_a))
            py = int(ry + 82 * math.sin(rad_a))
            draw.line([(rx, ry), (px, py)], fill=PRIMARY_TEAL, width=2)
            draw.ellipse([px - 6, py - 6, px + 6, py + 6], fill=LIGHT_BLUE, outline=(255, 255, 255), width=1)
        draw.text((x1 + 445, y1 + 382), "SCOPe-40 Top-1 Superfamily: 58.55% -> 88.32%!", fill=EMERALD, font=FONT_XS)

        bot_box = (x1 + 28, y1 + 450, x2 - 28, y2 - 24)
        draw.rounded_rectangle(bot_box, radius=14, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
        draw.text((x1 + 48, y1 + 470), "1,304-BYTE ZERO-TOAST BioVecInt8G128 VARLENA & 1-BIT RABITQ PRISM", fill=SECONDARY_PURPLE, font=FONT_SM)
        draw_fitted_line(
            draw,
            "• RaBitQ Error Bound eps Tightens 2.78x (2.8108 -> 1.0118) once off-diagonal covariance is removed",
            x1 + 48,
            y1 + 516,
            w - 96,
            TEXT_MAIN,
            (FONT_SM, FONT_XS),
        )
        draw_fitted_line(
            draw,
            "• Group-wise int8 quantization (G=128) fits 1,280 dims in 1,304 B (< 2,000 B TOAST threshold, 0 B TOAST I/O!)",
            x1 + 48,
            y1 + 558,
            w - 96,
            EMERALD,
            (FONT_SM, FONT_XS),
        )

    elif visual_mode == "allostery_csr":
        # Visual Canvas for Ep 3: 500 ns 5-LOX Allosteric Triterpenoids & CSR Sparse Graph Diff
        draw_fitted_line(
            draw,
            "ALLOSTERIC WHISPERS: 20 Å DISTAL BINDING & CSR SPARSE GRAPH DIFF",
            x1 + 28,
            y1 + 22,
            w - 56,
            PRIMARY_BLUE,
            (FONT_H3, FONT_BODY, FONT_SM),
        )
        top_box = (x1 + 28, y1 + 68, x2 - 28, y1 + 420)
        draw.rounded_rectangle(top_box, radius=14, fill=(11, 15, 25), outline=PRIMARY_BLUE, width=2)

        # Draw 20 Å allosteric pathway from Triterpenoid (AKBA / Amyrin) -> Cryptic Hinges (Ile126, Val29, Leu66) -> Catalytic Fe2+ Core
        nodes = [
            (x1 + 145, y1 + 245, "Distal Allosteric Site\n(AKBA / alpha-Amyrin)", SECONDARY_PURPLE),
            (x1 + 430, y1 + 165, "Cryptic Hinge Ile126\n(& Val29 / Leu66)", AMBER),
            (x1 + 715, y1 + 245, "Catalytic Fe2+ Core\n(20 Å Away -> Frozen!)", EMERALD),
        ]
        draw.line([(nodes[0][0], nodes[0][1]), (nodes[1][0], nodes[1][1])], fill=AMBER, width=5)
        draw.line([(nodes[1][0], nodes[1][1]), (nodes[2][0], nodes[2][1])], fill=EMERALD, width=5)
        for nx, ny, lbl, col in nodes:
            draw.ellipse([nx - 36, ny - 36, nx + 36, ny + 36], fill=PANEL_DARK, outline=col, width=4)
            for idx_l, l_part in enumerate(lbl.split("\n")):
                draw.text((nx - 85, ny + 48 + idx_l * 24), l_part, fill=col if idx_l == 0 else TEXT_MAIN, font=FONT_XS)

        bot_box = (x1 + 28, y1 + 440, x2 - 28, y2 - 24)
        draw.rounded_rectangle(bot_box, radius=14, fill=PANEL_BG, outline=EMERALD, width=2)
        draw.text((x1 + 48, y1 + 460), "DOUBLE-ENTRY THERMODYNAMIC AUDIT & CSR SET-DIFFERENCING", fill=EMERALD, font=FONT_SM)
        draw_fitted_line(
            draw,
            "• Cycle Closure: Sum(Per-Residue Gas + Polar) vs. Global Delta-G isolates the +22 kcal/mol cavity shift",
            x1 + 48,
            y1 + 508,
            w - 96,
            TEXT_MAIN,
            (FONT_SM, FONT_XS),
        )
        draw_fitted_line(
            draw,
            "• Opposite Forces, Same Cryptic Locks: AKBA pulls ionically (+0.41) while Amyrin locks hydrophobically (-0.68)",
            x1 + 48,
            y1 + 552,
            w - 96,
            LIGHT_BLUE,
            (FONT_SM, FONT_XS),
        )

    else:
        # Default / Ep 1 (`central_dogma_gist`): Soda-Lake Charge Armor (pI 4.38) & 1D Genomic Safe-Harbor Intervals
        draw_fitted_line(
            draw,
            "EXTREMOPHILE BRINE ARMOR (pI 4.38) & 1D SAFE-HARBOR INTERVAL TREES",
            x1 + 28,
            y1 + 22,
            w - 56,
            PRIMARY_BLUE,
            (FONT_H3, FONT_BODY, FONT_SM),
        )
        top_box = (x1 + 28, y1 + 68, x2 - 28, y1 + 410)
        draw.rounded_rectangle(top_box, radius=14, fill=(11, 15, 25), outline=PRIMARY_BLUE, width=2)

        # Left: Mesophilic Enzyme in 3.5 M NaCl (Precipitates) vs Right: Lake Magadi D3SWG7 (Acidic Asp/Glu Hydration Shield)
        lx, ly = x1 + 215, y1 + 240
        draw.text((x1 + 52, y1 + 84), "[X] Standard Enzyme (pI 6.40 in 3.5 M Salt)", fill=ACCENT_PINK, font=FONT_SM)
        draw.ellipse([lx - 75, ly - 75, lx + 75, ly + 75], fill=(60, 25, 38), outline=ACCENT_PINK, width=3)
        draw.text((lx - 85, ly + 90), "Salt strips water -> Clumps!", fill=ACCENT_PINK, font=FONT_XS)

        rx, ry = x1 + 615, y1 + 240
        draw.text((x1 + 430, y1 + 84), "[OK] Lake Magadi D3SWG7 (pI 4.38 Brine Shield)", fill=EMERALD, font=FONT_SM)
        draw.ellipse([rx - 88, ry - 88, rx + 88, ry + 88], outline=ORANGE, width=3)
        draw.ellipse([rx - 65, ry - 65, rx + 65, ry + 65], fill=(10, 58, 48), outline=EMERALD, width=3)
        draw.text((rx - 48, ry - 10), "0.44 Å Triad", fill=EMERALD, font=FONT_XS)
        draw.text((rx - 115, ry + 98), "Asp-/Glu- Carboxylate Cage Holds H2O!", fill=ORANGE, font=FONT_XS)

        bot_box = (x1 + 28, y1 + 430, x2 - 28, y2 - 24)
        draw.rounded_rectangle(bot_box, radius=14, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
        draw.text((x1 + 48, y1 + 450), "1D GENOMIC SAFE-HARBOR INTERVAL TREE (int4range GiST && and <->)", fill=SECONDARY_PURPLE, font=FONT_SM)
        draw.line([(x1 + 58, y1 + 525), (x2 - 58, y1 + 525)], fill=TEXT_MUTED, width=3)
        draw.rounded_rectangle([x1 + 90, y1 + 503, x1 + 260, y1 + 547], radius=8, fill=(20, 48, 82), outline=LIGHT_BLUE, width=2)
        draw.text((x1 + 112, y1 + 515), "Promoter [8k..9.2k]", fill=LIGHT_BLUE, font=FONT_XS)
        draw.rounded_rectangle([x1 + 320, y1 + 500, x1 + 545, y1 + 550], radius=8, fill=(6, 78, 59), outline=EMERALD, width=3)
        draw.text((x1 + 340, y1 + 515), "Safe Harbor [10.5k..11.8k]", fill=EMERALD, font=FONT_XS)
        draw.rounded_rectangle([x1 + 610, y1 + 503, x2 - 70, y1 + 547], radius=8, fill=(75, 20, 38), outline=ACCENT_PINK, width=2)
        draw.text((x1 + 630, y1 + 515), "Oncogene (Zero Overlap!)", fill=ACCENT_PINK, font=FONT_XS)
        draw_fitted_line(
            draw,
            "• Henderson-Hasselbalch bisection + 1D GiST interval trees run directly on PostgreSQL 8 KB heap pages.",
            x1 + 48,
            y1 + 575,
            w - 96,
            TEXT_MAIN,
            (FONT_SM, FONT_XS),
        )


def draw_progressive_left_slide(
    draw: ImageDraw.ImageDraw,
    box: Tuple[int, int, int, int],
    section_badge: str,
    accent_color: Tuple[int, int, int],
    items: List[Tuple[str, str, Optional[str]]],
    active_idx: int,
    lead_summary: Optional[str] = None,
) -> None:
    """
    Renders a clean, low-density progressive slide on the left panel where ONLY `active_idx`
    is expanded as the main spotlight card (synchronized with the narrator's current sentence),
    while the other items appear as compact 1-line progress pills.
    """
    x1, y1, x2, y2 = box
    w = x2 - x1
    draw.rounded_rectangle([x1, y1, x2, y2], radius=18, fill=PANEL_BG, outline=accent_color, width=2)

    # Top header + Slide progress dots
    draw_fitted_line(draw, section_badge, x1 + 28, y1 + 22, w - 180, accent_color, (FONT_H2, FONT_H3, FONT_BODY))
    pill_txt = f"Focus Beat {active_idx + 1} of {len(items)}"
    pw = text_width(draw, pill_txt, FONT_XS) + 24
    draw.rounded_rectangle([x2 - pw - 24, y1 + 20, x2 - 24, y1 + 52], radius=12, fill=PANEL_DARK, outline=accent_color, width=2)
    draw.text((x2 - pw - 12, y1 + 28), pill_txt, fill=TEXT_MAIN, font=FONT_XS)

    cur_y = y1 + 68
    if lead_summary:
        cur_y = draw_wrapped_box(draw, lead_summary, x1 + 28, cur_y, w - 56, FONT_SM, TEXT_MUTED, line_h=26, max_lines=3)
        cur_y += 12

    for idx, (head, body, extra) in enumerate(items):
        is_active = idx == active_idx
        is_done = idx < active_idx
        if is_active:
            # Large, spacious spotlight card for the sentence currently being spoken!
            card_h = 286 if extra else 224
            card_bot = min(y2 - 24, cur_y + card_h)
            draw.rounded_rectangle([x1 + 24, cur_y, x2 - 24, card_bot], radius=16, fill=PANEL_HIGHLIGHT, outline=PRIMARY_BLUE, width=3)
            badge = f"NOW EXPLAINING ({idx + 1}/{len(items)})"
            draw.text((x1 + 48, cur_y + 16), badge, fill=AMBER, font=FONT_XS)
            draw_fitted_line(draw, head, x1 + 48, cur_y + 44, w - 96, (255, 255, 255), (FONT_H3, FONT_BODY, FONT_SM))
            by = draw_wrapped_box(draw, body, x1 + 48, cur_y + 84, w - 96, FONT_BODY, TEXT_MAIN, line_h=28, max_lines=3 if extra else 4)
            if extra:
                draw_wrapped_box(draw, extra, x1 + 48, by + 10, w - 96, FONT_SM, EMERALD, line_h=24, max_lines=2)
            cur_y = card_bot + 14
        else:
            # Compact 1-line progress bar pill (zero wall-of-text clutter!)
            pill_h = 54
            if cur_y + pill_h > y2 - 18:
                break
            status_tag = "[OK]" if is_done else f"[{idx + 1}]"
            st_col = EMERALD if is_done else TEXT_MUTED
            draw.rounded_rectangle([x1 + 24, cur_y, x2 - 24, cur_y + pill_h], radius=10, fill=PANEL_DARK, outline=EMERALD if is_done else (45, 55, 78), width=1)
            draw_fitted_line(draw, f"{status_tag}  {head}", x1 + 44, cur_y + 16, w - 88, st_col, (FONT_SM, FONT_XS))
            cur_y += pill_h + 12


def render_frame(
    frame_idx: int,
    total_frames: int,
    act_idx: int,
    local_t: float,
    t: float,
    total_sec: float,
    cfg: Dict[str, Any],
    narration: str,
) -> Image.Image:
    progress = frame_idx / max(1, total_frames)
    img = Image.new("RGB", (W, H), BG_COLOR)
    draw = ImageDraw.Draw(img)

    for gx in range(0, W, 80):
        draw.line([(gx, 82), (gx, 1024)], fill=(19, 23, 33), width=1)
    for gy in range(82, 1024, 80):
        draw.line([(0, gy), (W, gy)], fill=(19, 23, 33), width=1)

    display_title = cfg.get("display_title") or cfg["title"]
    layer = cfg["layer"]
    visual_mode = cfg.get("visual_mode", "multipocket_posebusters")
    act_steps = cfg.get("act_steps", DEFAULT_ACT_STEPS)
    pair_cfg = cfg.get("organism_pair", {})
    series_tag = cfg.get("series_tag", "pg_bio Science & Heuristics Series — Teaching Computational Biology")

    # Determine speech-synchronized chunk & sub-slide index
    c_idx, total_chunks, _ = draw_caption_and_footer(draw, narration, local_t, series_tag)

    if act_idx == 0:
        # ACT I: THE PARADOX — 4 Progressive Sub-Slides synced directly to the spoken sentences + 2D Concept Canvas (NO 3D spinning!)
        bullets = cfg.get("bullets", [])[:4]
        n_sub = max(1, len(bullets))
        sub_idx = min(n_sub - 1, int(local_t * n_sub))
        draw_header(draw, layer, act_idx, act_steps, progress, t, total_sec, sub_idx, n_sub)

        draw_smart_title(draw, display_title, 56, 96, W - 112, PRIMARY_BLUE)
        sub = cfg.get("subtitle", "A Veritasium-Style Journey into Molecular Heuristics")
        draw_fitted_line(draw, sub, 56, 144, W - 112, SECONDARY_PURPLE, (FONT_H3, FONT_BODY, FONT_SM))

        chunks = split_narration_chunks(narration, max_words_per_chunk=22)
        items = [
            (
                b.lstrip("• ").strip(),
                chunks[min(i, len(chunks) - 1)] if chunks else cfg["bio_problem"],
                f"Target Discovery Benchmark: {cfg['speedup']}" if i == sub_idx else None,
            )
            for i, b in enumerate(bullets)
        ]
        draw_progressive_left_slide(
            draw,
            (56, 188, 930, 898),
            "ACT I: THE BIOLOGICAL PARADOX",
            ACCENT_PINK,
            items,
            sub_idx,
        )

        draw.rounded_rectangle([960, 188, 1864, 898], radius=18, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw_topic_2d_canvas(draw, visual_mode, (960, 188, 1864, 898), act_idx, sub_idx, t, local_t, cfg)

    elif act_idx == 1:
        # ACT II: THE HIDDEN TRAP — 4 Progressive Sub-Slides (1 Trap Card spotlighted per slide in sync with speech!)
        trap_cards = cfg.get(
            "trap_cards",
            [
                ("The Naive Assumption", "Pooling or scoring the whole molecule as a single average number."),
                ("The Hidden Physical Trap", "High raw scores often come from steric clashes, salt counterions, or dead binders."),
                ("The Consequence in the Lab", "90% of top-ranked computer hits fail or precipitate the moment they enter a test tube."),
                ("The Heuristic Shift", "Decompose the geometry and center on the reference state before scoring."),
            ],
        )[:4]
        n_sub = max(1, len(trap_cards))
        sub_idx = min(n_sub - 1, int(local_t * n_sub))
        draw_header(draw, layer, act_idx, act_steps, progress, t, total_sec, sub_idx, n_sub)

        act2_title = cfg.get("act2_title", f"Act II: {trap_cards[sub_idx][0]}")
        draw_smart_title(draw, act2_title, 56, 98, W - 112, PRIMARY_BLUE)

        items = [(t_title, t_desc, None) for t_title, t_desc in trap_cards]
        draw_progressive_left_slide(
            draw,
            (56, 168, 930, 898),
            "WHERE NAIVE ALGORITHMS GO BLIND",
            ACCENT_PINK,
            items,
            sub_idx,
        )

        draw.rounded_rectangle([960, 168, 1864, 898], radius=18, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw_topic_2d_canvas(draw, visual_mode, (960, 168, 1864, 898), act_idx, sub_idx, t, local_t, cfg)

    elif act_idx == 2:
        # ACT III: THE HEURISTIC EUREKA — 3 Progressive Sub-Slides (1 Step spotlighted at a time!)
        mappings = cfg.get("rosetta_mappings", [])[:3]
        n_sub = max(1, len(mappings))
        sub_idx = min(n_sub - 1, int(local_t * n_sub))
        draw_header(draw, layer, act_idx, act_steps, progress, t, total_sec, sub_idx, n_sub)

        act3_title = cfg.get("act3_title", f"Act III: {mappings[sub_idx][0]}" if mappings else "Act III: Step-by-Step Heuristic")
        draw_smart_title(draw, act3_title, 56, 98, W - 112, PRIMARY_BLUE)

        items = [(m_title, f"Physical Intuition: {m_why}", f"Mathematical Invariant: {m_how}") for m_title, m_why, m_how in mappings]
        draw_progressive_left_slide(
            draw,
            (56, 168, 930, 898),
            "STEP-BY-STEP HEURISTIC BREAKDOWN",
            SECONDARY_PURPLE,
            items,
            sub_idx,
        )

        draw.rounded_rectangle([960, 168, 1864, 898], radius=18, fill=PANEL_DARK, outline=EMERALD, width=2)
        draw_topic_2d_canvas(draw, visual_mode, (960, 168, 1864, 898), act_idx, sub_idx, t, local_t, cfg)

    elif act_idx == 3:
        # ACT IV: THE FOCUSED 3D STRUCTURAL PROOF + DECLARATIVE SQL QUERY
        # (This is the ONLY act with 3D organism motion, and it locks onto the pocket at local_t >= 0.55!)
        sub_idx = 0 if local_t < 0.35 else (1 if local_t < 0.65 else 2)
        draw_header(draw, layer, act_idx, act_steps, progress, t, total_sec, sub_idx, 3)

        act4_title = cfg.get("act4_title", "Act IV: Focused 3D Structural Alignment & Declarative Query")
        draw_smart_title(draw, act4_title, 56, 98, W - 112, PRIMARY_BLUE)

        draw.rounded_rectangle([56, 168, 956, 898], radius=18, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw.ellipse([84, 192, 104, 212], fill=ACCENT_PINK)
        draw.ellipse([116, 192, 136, 212], fill=AMBER)
        draw.ellipse([148, 192, 168, 212], fill=PRIMARY_TEAL)
        draw.text((190, 190), "Expressing the Heuristic as a Declarative Query", fill=TEXT_MUTED, font=FONT_MONO_SM)
        draw.line([(56, 228), (956, 228)], fill=(45, 55, 78), width=2)

        sql_lines = cfg.get("sql_lines", [])
        visible_lines = max(1, min(len(sql_lines), int((local_t / 0.32) * len(sql_lines)) + 1))
        sy = 246
        for idx, line in enumerate(sql_lines[:visible_lines]):
            cursor_str = " █" if (idx == visible_lines - 1 and local_t < 0.32 and int(t * 6) % 2 == 0) else ""
            col = SECONDARY_PURPLE if any(k in line for k in ("SELECT", "FROM", "WHERE", "WITH", "JOIN", "LATERAL")) else TEXT_MAIN
            draw_fitted_line(draw, line + cursor_str, 78, sy, 856, col, (FONT_MONO_SM, FONT_XS, FONT_XXS))
            sy += 28

        if local_t >= 0.30:
            res_lines = cfg.get("sql_results", [])
            ry = 512
            for idx, rline in enumerate(res_lines):
                color = EMERALD if idx == len(res_lines) - 1 else AMBER if idx in (3, 4) else TEXT_MAIN
                draw_fitted_line(draw, rline, 78, ry, 856, color, (FONT_MONO_SM, FONT_XS, FONT_XXS))
                ry += 27

            draw.rounded_rectangle([78, 748, 934, 876], radius=12, fill=PANEL_HIGHLIGHT, outline=EMERALD, width=2)
            focus_note = (
                "3D Stage Locked on Conserved Pocket: Notice how the catalytic residues align within <0.6 Å!"
                if local_t >= 0.55
                else "3D Stage Aligning: Rotating Kabsch-superposed C-alpha backbones into shared pocket frame..."
            )
            draw.text((98, 762), focus_note, fill=EMERALD, font=FONT_SM)
            draw_wrapped_box(draw, cfg["cs_mechanism"], 98, 796, 816, FONT_XS, TEXT_MAIN, line_h=23, max_lines=3)

        # Focused 3D Dual-Organism Stage (rotates smoothly then locks onto pocket!)
        draw.rounded_rectangle([984, 168, 1864, 898], radius=18, fill=PANEL_DARK, outline=EMERALD, width=2)
        draw_dual_organism_comparison_stage(draw, (984, 168, 1864, 898), pair_cfg, t, local_t, compact=True)

    else:
        # ACT V: DISCOVERY IMPACT & RESEARCH TAKEAWAYS — 3 Progressive Sub-Slides + 2D Summary Canvas (NO 3D spinning!)
        takeaways = cfg.get(
            "takeaways",
            [
                ("1. Nature Conserves 3D Physics, Not 1D Spelling", "Across 3.5 billion years, extremophiles rewrite 70% of their surface sequence to survive salt, acid, or heat while locking the catalytic pocket within 0.5 Å."),
                ("2. Good Heuristics Expose Hidden Invariants", "Whether it is contact-restricted ipSAE, 24-ray SAS concavity, reference-centered EnzRef vectors, or ZCA whitening, the right math makes false positives collapse to zero."),
                ("3. Interactive Research in Your Browser", f"Rotate these aligned organisms in 3D at {BLOG_URL}/pg-bio-videos/ and explore the open-source research engine at {REPO_URL}."),
            ],
        )[:3]
        n_sub = max(1, len(takeaways))
        sub_idx = min(n_sub - 1, int(local_t * n_sub))
        draw_header(draw, layer, act_idx, act_steps, progress, t, total_sec, sub_idx, n_sub)

        act5_title = cfg.get("act5_title", "Act V: From Mathematical Insight to Extremophile Discovery")
        draw_smart_title(draw, act5_title, 56, 98, W - 112, PRIMARY_BLUE)

        items = [(t_head, t_body, f"Verified Result: {cfg['speedup']}") for t_head, t_body in takeaways]
        draw_progressive_left_slide(
            draw,
            (56, 168, 930, 898),
            "WHAT THIS HEURISTIC REVEALS IN NATURE",
            EMERALD,
            items,
            sub_idx,
        )

        draw.rounded_rectangle([960, 168, 1864, 898], radius=18, fill=PANEL_DARK, outline=SECONDARY_PURPLE, width=2)
        draw_topic_2d_canvas(draw, visual_mode, (960, 168, 1864, 898), act_idx, sub_idx, t, local_t, cfg)

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
            "Now let us bring in the 3D structures on the right. As they align and lock onto the active site, notice how their outer surfaces are completely different to survive extreme salt and pH, yet the glowing emerald catalytic pocket in the center locks in to less than half an Angstrom.",
            "When we understand the heuristic first, expressing it as a query becomes effortless—and suddenly we can explore millions of extremophile organisms in seconds. Dive into the interactive 3D models at ideia dot me.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Render Veritasium-style 1080p MP4 video + neural voiceover for a pg_bio episode.")
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
    print(f"Voiceover ready ({total_sec:.1f}s across 5 acts). Rendering {total_frames} synchronized 1080p multi-slide frames...")

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
    for act_idx, n_frames in enumerate(act_frame_counts):
        if broken_pipe:
            break
        for local_f in range(n_frames):
            local_t = local_f / max(1, n_frames - 1)
            t = global_frame / FPS
            img = render_frame(global_frame, total_frames, act_idx, local_t, t, total_sec, cfg, narrations[act_idx])
            try:
                proc.stdin.write(img.tobytes())
            except BrokenPipeError:
                broken_pipe = True
                break
            global_frame += 1

    try:
        proc.stdin.close()
    except BrokenPipeError:
        pass
    stderr_out = proc.stderr.read().decode("utf-8", errors="replace") if proc.stderr else ""
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"ffmpeg exited with {rc}: {stderr_out}")
    tmp_path.replace(out_path)
    print(f"Rendered 1080p MP4 with Multi-Slide Sync & Focused Act-IV 3D Stage ({total_frames} frames, {total_sec:.1f}s): {out_path}")


if __name__ == "__main__":
    main()
