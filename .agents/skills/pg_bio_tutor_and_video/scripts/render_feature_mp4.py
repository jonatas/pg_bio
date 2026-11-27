#!/usr/bin/env python3
"""
Veritasium-Style 1080p MP4 Video + Neural Voiceover Renderer for pg_bio Episodes.
Designed around curiosity-driven storytelling, clear teaching of the mathematical
and physical heuristics, and live 3D side-by-side organism comparison motion
(mirroring the website's 3Dmol.js synchronized dual-organism viewers).

Key Design Principles:
  1. Story & Heuristic First (Veritasium 5-Act Arc):
       - Act I:   The Paradox (A surprising biological mystery or failure mode)
       - Act II:  The Hidden Trap (Why naive intuition or raw scores fail — NO repetitive anti-Python sales pitch; only Episode 1 introduces the data-to-compute concept)
       - Act III: The Heuristic Eureka (Step-by-step visual derivation of WHY the heuristic works)
       - Act IV:  3D Organism Motion (Live Kabsch-aligned dual-organism 3D PDB comparison rotating in sync + live SQL probe)
       - Act V:   Real-World Discovery (What this unlocks in extreme biology & research)
  2. Real 3D Organism Comparison Motion:
       - Loads actual AlphaFold / experimental `.pdb` C-alpha backbones from `/Users/jonatas/code/ideia.me/assets/models/`,
         Kabsch-centers and scales them, and renders synchronized rotating 3D ribbons/tubes with glowing catalytic
         pocket spheres, acidic brine-armor halos, and live residue distance bridges!
  3. Zero Clipped Titles:
       - Header titles use clean, punchy display titles and dynamic font sizing so no title ever truncates with `...`.
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

# Exact `../ideia.me` Theme Palette
BG_COLOR = (15, 17, 23)            # --bg-color: #0f1117
PANEL_BG = (22, 27, 39)            # --card-bg-color: #161b27
PANEL_DARK = (26, 31, 46)          # --code-bg-color: #1a1f2e
PRIMARY_BLUE = (92, 158, 255)      # --primary-color: #5c9eff
BRAND_BLUE = (67, 97, 238)         # #4361ee
SECONDARY_PURPLE = (177, 88, 252)  # --secondary-color: #b158fc
ACCENT_PINK = (255, 90, 157)       # --accent-color: #ff5a9d
PRIMARY_TEAL = (13, 148, 136)      # --primary-teal: #0d9488
EMERALD = (16, 185, 129)           # #10b981
LIGHT_BLUE = (56, 189, 248)        # #38bdf8
PG_TOPIC_BLUE = (51, 103, 145)     # #336791
TEXT_MAIN = (204, 214, 246)        # #ccd6f6
TEXT_MUTED = (136, 146, 176)       # #8892b0
AMBER = (255, 193, 7)              # #FFC107
ORANGE = (249, 115, 22)            # #f97316

BLOG_URL = "https://ideia.me"
REPO_URL = "https://github.com/jonatas/pg_bio"
MODELS_DIR = Path("/Users/jonatas/code/ideia.me/assets/models")

DEFAULT_ACT_STEPS = [
    "I. The Paradox",
    "II. The Hidden Trap",
    "III. The Heuristic",
    "IV. 3D Organism Motion",
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
FONT_H1_MD = load_font(31)
FONT_H2 = load_font(25)
FONT_BODY = load_font(21)
FONT_SM = load_font(18)
FONT_XS = load_font(15)
FONT_CAPTION = load_font(22)
FONT_MONO = load_font(19, mono=True)
FONT_MONO_SM = load_font(16, mono=True)


# Cache of parsed & Kabsch-normalized 3D C-alpha traces from real PDB files
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
        # Fallback synthetic alpha/beta barrel trace if file missing
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

    # Subsample smoothly if very large so 3D rendering is crisp and fast
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
    if text_width(draw, text, FONT_H1) <= max_px:
        draw.text((x, y), text, fill=fill, font=FONT_H1)
    elif text_width(draw, text, FONT_H1_MD) <= max_px:
        draw.text((x, y + 3), text, fill=fill, font=FONT_H1_MD)
    else:
        draw.text((x, y + 6), text, fill=fill, font=FONT_H2)


def fit_single_line(draw: ImageDraw.ImageDraw, text: str, max_px: int, font: ImageFont.FreeTypeFont) -> str:
    if text_width(draw, text, font) <= max_px:
        return text
    trimmed = text
    while len(trimmed) > 8 and text_width(draw, trimmed + "…", font) > max_px:
        trimmed = trimmed[:-1]
    return trimmed + "…"


def draw_wrapped_box(
    draw: ImageDraw.ImageDraw,
    text: str,
    x: int,
    y: int,
    max_px: int,
    font: ImageFont.FreeTypeFont,
    fill: Tuple[int, int, int],
    line_h: int = 30,
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

    for idx, line in enumerate(lines[:max_lines]):
        if idx == max_lines - 1 and len(lines) > max_lines:
            line = fit_single_line(draw, line + "…", max_px, font)
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
) -> None:
    draw.rectangle([0, 0, W, 82], fill=(11, 13, 19))
    draw.line([(0, 82), (W, 82)], fill=(45, 55, 78), width=2)
    draw.rectangle([0, 79, int(W * progress), 82], fill=PRIMARY_BLUE)

    badge_text = fit_single_line(draw, layer_str, 460, FONT_XS)
    bw = text_width(draw, badge_text, FONT_XS) + 32
    draw.rounded_rectangle([44, 18, 44 + bw, 60], radius=15, fill=PG_TOPIC_BLUE, outline=PRIMARY_BLUE, width=2)
    draw.text((60, 29), badge_text, fill=(255, 255, 255), font=FONT_XS)

    step_x = 560
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

    timer_text = f"TIME {format_mm_ss(elapsed_sec)} / {format_mm_ss(total_sec)}"
    tw = text_width(draw, timer_text, FONT_XS) + 28
    draw.rounded_rectangle([W - 44 - tw, 18, W - 44, 60], radius=15, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
    draw.text((W - 44 - tw + 14, 29), timer_text, fill=TEXT_MAIN, font=FONT_XS)


def draw_caption_and_footer(draw: ImageDraw.ImageDraw, narration: str, act_local_progress: float, series_tag: str) -> None:
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
    """
    Renders a smooth 3D rotating C-alpha protein backbone with depth-shaded ribbon links,
    glowing catalytic pocket spheres, and extremophile surface charge highlights.
    Returns projected (px, py) coordinates for the primary pocket residues.
    """
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

    # 1. Draw smooth C-alpha backbone tube segments
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

    # 2. Draw residue spheres sorted back-to-front
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

    # 3. Draw glowing catalytic pocket bounding shell around pocket centroid
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
    Renders the website's signature side-by-side 3Dmol.js synchronized dual-organism
    comparison stage with mirrored camera rotation, Kabsch superposition bridge vectors,
    and organism habitat badges.
    """
    x1, y1, x2, y2 = box
    w = x2 - x1
    mid_x = x1 + w // 2

    pdb1 = pair_cfg.get("pdb1", "AF-A0ABD5UZ94-F1-model_v4_aligned.pdb")
    pdb2 = pair_cfg.get("pdb2", "AF-D3SWG7-F1-model_v4.pdb")
    atoms1 = load_pdb_ca_trace(pdb1)
    atoms2 = load_pdb_ca_trace(pdb2)

    # Synchronized mirrored 3D camera rotation (just like pg_bio_sync.js!)
    angle_y = t * 0.72
    angle_x = 0.28 * math.sin(t * 0.45)

    # Left Viewport (Bait / Reference Organism in Cyan) & Right Viewport (Extremophile Discovery in Magenta)
    vp_top = y1 + (58 if compact else 64)
    vp_bot = y2 - (84 if compact else 108)
    left_box = (x1 + 18, vp_top, mid_x - 10, vp_bot)
    right_box = (mid_x + 10, vp_top, x2 - 18, vp_bot)

    draw.rounded_rectangle(left_box, radius=14, fill=(9, 13, 22), outline=LIGHT_BLUE, width=2)
    draw.rounded_rectangle(right_box, radius=14, fill=(9, 13, 22), outline=ACCENT_PINK, width=2)

    # Top title banner
    header_title = pair_cfg.get("stage_header", "3D MIRRORED ORGANISM COMPARISON (KABSCH ALIGNED)")
    draw.text((x1 + 24, y1 + 18), fit_single_line(draw, header_title, w - 48, FONT_H2), fill=PRIMARY_BLUE, font=FONT_H2)

    # Organism Header Bars inside viewports
    draw.rectangle([left_box[0] + 2, left_box[1] + 2, left_box[2] - 2, left_box[1] + 40], fill=(22, 32, 52))
    draw.text(
        (left_box[0] + 12, left_box[1] + 10),
        fit_single_line(draw, pair_cfg.get("org1_label", "Bait Reference"), (left_box[2] - left_box[0]) - 24, FONT_XS),
        fill=LIGHT_BLUE,
        font=FONT_XS,
    )

    draw.rectangle([right_box[0] + 2, right_box[1] + 2, right_box[2] - 2, right_box[1] + 40], fill=(42, 20, 46))
    draw.text(
        (right_box[0] + 12, right_box[1] + 10),
        fit_single_line(draw, pair_cfg.get("org2_label", "Extremophile Discovery"), (right_box[2] - right_box[0]) - 24, FONT_XS),
        fill=ACCENT_PINK,
        font=FONT_XS,
    )

    pk1 = pair_cfg.get("pocket_indices", [22, 35, 48, 62, 78, 94])
    pk2 = pair_cfg.get("second_pocket_indices", [])
    hi_mode = "brine_armor" if (int(t / 3.5) % 2 == 1) else "pocket"

    c1 = ((left_box[0] + left_box[2]) // 2, (left_box[1] + left_box[3]) // 2 + 16)
    c2 = ((right_box[0] + right_box[2]) // 2, (right_box[1] + right_box[3]) // 2 + 16)
    rad = min((left_box[2] - left_box[0]), (left_box[3] - left_box[1])) * 0.36

    pts1 = render_3d_protein_trace(draw, atoms1, c1, rad, angle_y, angle_x, LIGHT_BLUE, pk1, hi_mode, t, pk2)
    pts2 = render_3d_protein_trace(draw, atoms2, c2, rad, angle_y, angle_x, ACCENT_PINK, pk1, hi_mode, t, pk2)

    # Draw animated catalytic equivalence bridge lines across the two organisms
    if pts1 and pts2:
        active_bridge = int(t * 3.0) % min(len(pts1), len(pts2))
        p_left = pts1[active_bridge]
        p_right = pts2[active_bridge]
        draw.line([p_left, p_right], fill=EMERALD, width=2)
        mid_bx = (p_left[0] + p_right[0]) // 2
        mid_by = (p_left[1] + p_right[1]) // 2
        badge_txt = pair_cfg.get("bridge_badge", "0.44 Å Catalytic Match")
        bw = text_width(draw, badge_txt, FONT_XS) + 20
        draw.rounded_rectangle([mid_bx - bw // 2, mid_by - 16, mid_bx + bw // 2, mid_by + 16], radius=10, fill=(6, 78, 59), outline=EMERALD, width=2)
        draw.text((mid_bx - bw // 2 + 10, mid_by - 8), badge_txt, fill=(255, 255, 255), font=FONT_XS)

    # Bottom Telemetry / Habitat Comparison Strip
    strip_y = vp_bot + 10
    draw.rounded_rectangle([x1 + 18, strip_y, x2 - 18, y2 - 14], radius=12, fill=PANEL_BG, outline=EMERALD, width=2)
    m_line1 = pair_cfg.get("metric_line1", "3D Core RMSD: 0.44 Å (TM = 0.996)  |  9D Pocket Overlap: 95.51%  |  Sequence Identity: 31.4%")
    m_line2 = pair_cfg.get("metric_line2", "Heuristic Proof: Surface mutates into acidic Asp/Glu brine shield while catalytic triad stays locked!")
    draw.text((x1 + 34, strip_y + 10), fit_single_line(draw, m_line1, w - 68, FONT_MONO_SM), fill=EMERALD, font=FONT_MONO_SM)
    if not compact:
        draw.text((x1 + 34, strip_y + 42), fit_single_line(draw, m_line2, w - 68, FONT_XS), fill=TEXT_MAIN, font=FONT_XS)


def draw_heuristic_diagram(
    draw: ImageDraw.ImageDraw,
    visual_mode: str,
    box: Tuple[int, int, int, int],
    t: float,
    local_t: float,
    pair_cfg: Dict[str, Any],
) -> None:
    """
    Renders the deep heuristic visualization in Act II (The Trap) and Act III (The Heuristic Eureka)
    customized to the exact scientific topic of each episode.
    """
    x1, y1, x2, y2 = box
    cx = (x1 + x2) // 2
    w = x2 - x1

    if visual_mode == "multipocket_posebusters":
        # Top: Single Blind Docking Box vs Multi-Pocket Connolly SAS Spheres (pk1..pk5)
        draw.text((x1 + 28, y1 + 20), "HEURISTIC 1: 24-RAY FIBONACCI SAS VOXEL HASHING", fill=PRIMARY_BLUE, font=FONT_H2)
        top_box = (x1 + 24, y1 + 64, x2 - 24, y1 + 356)
        draw.rounded_rectangle(top_box, radius=14, fill=(11, 15, 25), outline=PRIMARY_BLUE, width=2)

        # Left sub-diagram: Single Box Trap (misses pk1 & pk2)
        lx = x1 + 195
        ly = y1 + 210
        draw.text((x1 + 44, y1 + 78), "[X] Single Blind Box (Center of Mass)", fill=ACCENT_PINK, font=FONT_SM)
        draw.ellipse([lx - 115, ly - 75, lx - 15, ly + 75], fill=(30, 58, 95), outline=LIGHT_BLUE, width=2)
        draw.ellipse([lx + 15, ly - 75, lx + 115, ly + 75], fill=(58, 28, 75), outline=SECONDARY_PURPLE, width=2)
        draw.rectangle([lx - 45, ly - 45, lx + 45, ly + 45], outline=ACCENT_PINK, width=3)
        draw.text((lx - 72, ly + 82), "Lands in empty hinge!", fill=ACCENT_PINK, font=FONT_XS)

        # Right sub-diagram: 5.0 Å Voxel Grid + Fibonacci SAS Probes finding pk1 & pk2
        rx = x1 + 575
        ry = y1 + 210
        draw.text((x1 + 400, y1 + 78), "[OK] 5.0 Å Voxel Hash + SAS Concavity (pk1 & pk2)", fill=EMERALD, font=FONT_SM)
        for gx in range(rx - 140, rx + 150, 35):
            draw.line([(gx, ry - 80), (gx, ry + 80)], fill=(25, 35, 52), width=1)
        for gy in range(ry - 80, ry + 90, 35):
            draw.line([(rx - 140, gy), (rx + 140, gy)], fill=(25, 35, 52), width=1)
        draw.ellipse([rx - 125, ry - 75, rx - 15, ry + 75], fill=(20, 48, 82), outline=LIGHT_BLUE, width=2)
        draw.ellipse([rx + 15, ry - 75, rx + 125, ry + 75], fill=(48, 22, 68), outline=SECONDARY_PURPLE, width=2)
        pulse = int(4 * math.sin(t * 4.0))
        draw.ellipse([rx - 85 - pulse, ry - 25 - pulse, rx - 45 + pulse, ry + 15 + pulse], fill=EMERALD, outline=(255, 255, 255), width=2)
        draw.text((rx - 110, ry + 26), "★ pk1 (MBL Zn2+)", fill=EMERALD, font=FONT_XS)
        draw.ellipse([rx + 45 - pulse, ry - 20 - pulse, rx + 85 + pulse, ry + 20 + pulse], fill=AMBER, outline=(255, 255, 255), width=2)
        draw.text((rx + 30, ry + 26), "★ pk2 (Rhodanese)", fill=AMBER, font=FONT_XS)

        # Bottom: The 21-Check PoseBusters Physical Gate & Non-Hiding Dual-Ledger Salt Rescue
        bot_box = (x1 + 24, y1 + 374, x2 - 24, y2 - 20)
        draw.rounded_rectangle(bot_box, radius=14, fill=PANEL_BG, outline=EMERALD, width=2)
        draw.text((x1 + 44, y1 + 392), "HEURISTIC 2: 21-CHECK POSEBUSTERS & DUAL-LEDGER SALT RESCUE", fill=EMERALD, font=FONT_BODY)

        checks = [
            ("1. Vendor Salt Stripping", "[Na+].Diclofenac  ->  Active Diclofenac (>=3 heavy atoms; bare [Na+] rejected)", EMERALD),
            ("2. 3x3 Jacobi SVD Ring Planarity", "Best-fit plane RMSD <= 0.25 Å (catches warped aromatic rings instantly)", LIGHT_BLUE),
            ("3. Pairwise VdW Steric Ratio", "d_ij / (r_vdW,i + r_vdW,j) >= 0.75 (blocks atomic overlaps & clash hallucinations)", AMBER),
            ("4. Non-Hiding Dual-Ledger Contract", "passes_all_21 = TRUE  |  raw_passes_all_21 = FALSE  |  tag = 'desalted_from_salt'", SECONDARY_PURPLE),
        ]
        cy_c = y1 + 436
        for c_title, c_desc, col in checks:
            draw.rounded_rectangle([x1 + 44, cy_c, x2 - 44, cy_c + 56], radius=10, fill=PANEL_DARK, outline=col, width=2)
            draw.text((x1 + 60, cy_c + 8), c_title, fill=col, font=FONT_SM)
            draw.text((x1 + 60, cy_c + 32), fit_single_line(draw, c_desc, w - 125, FONT_XS), fill=TEXT_MAIN, font=FONT_XS)
            cy_c += 64

    elif visual_mode == "enzref_pslim":
        # Top: EnzRef Reference-Centered Vector Geometry (Why Dead Binders Collapse to Zero)
        draw.text((x1 + 28, y1 + 20), "HEURISTIC 1: ENZREF REFERENCE-CENTERED ROLE CONTRAST", fill=PRIMARY_BLUE, font=FONT_H2)
        top_box = (x1 + 24, y1 + 64, x2 - 24, y1 + 365)
        draw.rounded_rectangle(top_box, radius=14, fill=(11, 15, 25), outline=PRIMARY_BLUE, width=2)

        ox, oy = x1 + 125, y1 + 295
        rx_pt, ry_pt = x1 + 340, y1 + 165
        px_pt, py_pt = x1 + 610, y1 + 135
        mx_pt, my_pt = (rx_pt + px_pt) // 2, (ry_pt + py_pt) // 2

        # Vectors from Enzyme Background v_bg (Origin)
        draw.line([(ox, oy), (rx_pt, ry_pt)], fill=LIGHT_BLUE, width=4)
        draw.line([(ox, oy), (px_pt, py_pt)], fill=SECONDARY_PURPLE, width=4)
        draw.line([(ox, oy), (mx_pt, my_pt)], fill=AMBER, width=3)
        # Directed Forward Role Contrast c_fwd = d_prod - d_react
        draw.line([(rx_pt, ry_pt), (px_pt, py_pt)], fill=EMERALD, width=5)

        draw.ellipse([ox - 12, oy - 12, ox + 12, oy + 12], fill=PRIMARY_BLUE, outline=(255, 255, 255), width=2)
        draw.text((ox - 30, oy + 18), "Enzyme Origin v_bg", fill=PRIMARY_BLUE, font=FONT_XS)

        draw.ellipse([rx_pt - 10, ry_pt - 10, rx_pt + 10, ry_pt + 10], fill=LIGHT_BLUE, outline=(255, 255, 255), width=2)
        draw.text((rx_pt - 135, ry_pt - 28), "Reactant d_react = r - v_bg", fill=LIGHT_BLUE, font=FONT_XS)

        draw.ellipse([px_pt - 10, py_pt - 10, px_pt + 10, py_pt + 10], fill=SECONDARY_PURPLE, outline=(255, 255, 255), width=2)
        draw.text((px_pt + 16, py_pt - 10), "Product d_prod = p - v_bg", fill=SECONDARY_PURPLE, font=FONT_XS)

        draw.text((mx_pt - 75, my_pt - 36), "c_fwd = p - r (Turnover Vector!)", fill=EMERALD, font=FONT_SM)
        draw.text((x1 + 280, y1 + 265), "• Dead Binder (p == r)  -->  c_fwd = 0  -->  FluxRank = 0.000", fill=ACCENT_PINK, font=FONT_SM)
        draw.text((x1 + 280, y1 + 298), "• Swap(r <-> p)  -->  c_fwd flips sign  -->  Flux(e,p,r) = -Flux(e,r,p)", fill=EMERALD, font=FONT_SM)

        # Bottom: pSLiM Sliding Window ZCA Whitening & Hydrophobic-Core Linker Junction Gate
        bot_box = (x1 + 24, y1 + 382, x2 - 24, y2 - 20)
        draw.rounded_rectangle(bot_box, radius=14, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
        draw.text((x1 + 44, y1 + 400), "HEURISTIC 2: pSLiM 8-AA WINDOW ZCA + LINKER JUNCTION GATE", fill=SECONDARY_PURPLE, font=FONT_BODY)

        # Disordered microprotein bar with Motif A, Flexible Linker, Motif B
        bar_y = y1 + 460
        draw.rounded_rectangle([x1 + 54, bar_y, x1 + 254, bar_y + 52], radius=10, fill=(6, 78, 59), outline=EMERALD, width=2)
        draw.text((x1 + 72, bar_y + 16), "SLiM A [1..8] (M1)", fill=EMERALD, font=FONT_XS)

        draw.rounded_rectangle([x1 + 258, bar_y + 8, x1 + 538, bar_y + 44], radius=8, fill=(30, 41, 59), outline=LIGHT_BLUE, width=2)
        draw.text((x1 + 280, bar_y + 16), "Linker: GGGSGGGGS (GRAVY -0.4)", fill=LIGHT_BLUE, font=FONT_XS)

        draw.rounded_rectangle([x1 + 542, bar_y, x2 - 54, bar_y + 52], radius=10, fill=(76, 29, 149), outline=SECONDARY_PURPLE, width=2)
        draw.text((x1 + 558, bar_y + 16), "SLiM B [28..35] (M2)", fill=SECONDARY_PURPLE, font=FONT_XS)

        rules = (
            "Why 8-aa Windows? Full-protein averaging dilutes an 8-aa motif by 10x; Sine-bell + ZCA whitening restores unit vec16. "
            "Linker Gate: If linker GRAVY > 0.45 or has >=5 hydrophobic residues [IVLFWMC], fusion slices a buried core and precipitates—so "
            "we switch from Single_Polypeptide_Fusion to Bicistronic_Operon_Required."
        )
        draw_wrapped_box(draw, rules, x1 + 48, y1 + 530, w - 96, FONT_XS, TEXT_MAIN, line_h=26, max_lines=4)

    else:
        # For all other modes, render the dual-organism 3D comparison stage directly!
        draw_dual_organism_comparison_stage(draw, box, pair_cfg, t, local_t, compact=False)


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
    series_tag = cfg.get("series_tag", "🧬 pg_bio Science & Heuristics Series — Teaching Computational Biology")

    draw_header(draw, layer, act_idx, act_steps, progress, t, total_sec)

    if act_idx == 0:
        # ACT I: THE PARADOX / REAL-WORLD MYSTERY + LIVE 3D ORGANISM COMPARISON MOTION
        draw_smart_title(draw, display_title, 56, 98, W - 112, PRIMARY_BLUE)
        sub = cfg.get("subtitle", "A Veritasium-Style Journey into Molecular Heuristics")
        draw.text((56, 148), fit_single_line(draw, sub, W - 112, FONT_H2), fill=SECONDARY_PURPLE, font=FONT_H2)

        # Left Panel: The Biological Paradox & Intuition Hook
        draw.rounded_rectangle([56, 198, 930, 898], radius=18, fill=PANEL_BG, outline=ACCENT_PINK, width=2)
        paradox_heading = cfg.get("act1_heading", "ACT I: THE BIOLOGICAL PARADOX")
        draw.text((84, 222), paradox_heading, fill=ACCENT_PINK, font=FONT_H2)
        y = draw_wrapped_box(draw, cfg["bio_problem"], 84, 270, 816, FONT_BODY, TEXT_MAIN, line_h=33, max_lines=6)

        draw.rounded_rectangle([84, max(488, y + 14), 902, 874], radius=14, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        box_title = cfg.get("act1_box_title", "Key Questions We Will Solve Together:")
        draw.text((108, max(506, y + 30)), box_title, fill=PRIMARY_BLUE, font=FONT_BODY)
        by = max(548, y + 72)
        for b in cfg.get("bullets", [])[:4]:
            draw.text((108, by), fit_single_line(draw, b, 770, FONT_SM), fill=TEXT_MAIN, font=FONT_SM)
            by += 36

        # Right Panel: Live 3D Dual-Organism Comparison Stage (rotating in sync!)
        draw.rounded_rectangle([960, 198, 1864, 898], radius=18, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw_dual_organism_comparison_stage(draw, (960, 198, 1864, 898), pair_cfg, t, local_t, compact=False)

    elif act_idx == 1:
        # ACT II: THE HIDDEN TRAP (Why Naive Heuristics Fail — Story-Specific, NOT repetitive anti-Python!)
        act2_title = cfg.get("act2_title", "Act II: Why Naive Intuition & Raw Scores Fall Into a Trap")
        draw_smart_title(draw, act2_title, 56, 98, W - 112, PRIMARY_BLUE)

        # Left Panel: The Naive Trap vs The Physical Reality
        draw.rounded_rectangle([56, 168, 930, 898], radius=18, fill=PANEL_BG, outline=ACCENT_PINK, width=2)
        trap_cards = cfg.get(
            "trap_cards",
            [
                ("The Naive Assumption", "Pooling or scoring the whole molecule as a single average number."),
                ("The Hidden Physical Trap", "High raw scores often come from steric clashes, salt counterions, or dead binders."),
                ("The Consequence in the Lab", "90% of top-ranked computer hits fail or precipitate the moment they enter a test tube."),
                ("The Heuristic Shift", "Decompose the geometry and center on the reference state before scoring."),
            ],
        )
        draw.text((84, 194), "WHERE NAIVE ALGORITHMS GO BLIND", fill=ACCENT_PINK, font=FONT_H2)
        for idx, (t_title, t_desc) in enumerate(trap_cards[:4]):
            top = 248 + idx * 154
            col = EMERALD if idx == 3 else ACCENT_PINK if idx in (1, 2) else LIGHT_BLUE
            draw.rounded_rectangle([84, top, 902, top + 136], radius=12, fill=PANEL_DARK, outline=col, width=2)
            draw.text((108, top + 14), fit_single_line(draw, f"{idx+1}. {t_title}", 770, FONT_BODY), fill=col, font=FONT_BODY)
            draw_wrapped_box(draw, t_desc, 108, top + 52, 770, FONT_SM, TEXT_MAIN, line_h=26, max_lines=3)

        # Right Panel: Visual Heuristic Diagram
        draw.rounded_rectangle([960, 168, 1864, 898], radius=18, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw_heuristic_diagram(draw, visual_mode, (960, 168, 1864, 898), t, local_t, pair_cfg)

    elif act_idx == 2:
        # ACT III: THE HEURISTIC EUREKA (Step-by-Step Intuitive Breakdown of Why the Technique Works)
        act3_title = cfg.get("act3_title", "Act III: How the Heuristic Actually Works (Step by Step)")
        draw_smart_title(draw, act3_title, 56, 98, W - 112, PRIMARY_BLUE)

        draw.rounded_rectangle([56, 168, 980, 898], radius=18, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
        draw.text((84, 194), "STEP-BY-STEP HEURISTIC BREAKDOWN", fill=SECONDARY_PURPLE, font=FONT_H2)

        mappings = cfg.get("rosetta_mappings", [])
        for idx, (m_title, m_why, m_how) in enumerate(mappings[:3]):
            top = 248 + idx * 208
            draw.rounded_rectangle([84, top, 952, top + 190], radius=14, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
            draw.text((108, top + 14), fit_single_line(draw, m_title, 820, FONT_BODY), fill=PRIMARY_BLUE, font=FONT_BODY)
            draw_wrapped_box(draw, f"Intuition: {m_why}", 108, top + 52, 820, FONT_SM, TEXT_MAIN, line_h=25, max_lines=2)
            draw_wrapped_box(draw, f"Why It Works: {m_how}", 108, top + 114, 820, FONT_SM, EMERALD, line_h=25, max_lines=2)

        draw.rounded_rectangle([1008, 168, 1864, 898], radius=18, fill=PANEL_DARK, outline=EMERALD, width=2)
        draw_heuristic_diagram(draw, visual_mode, (1008, 168, 1864, 898), t, local_t, pair_cfg)

    elif act_idx == 3:
        # ACT IV: SIDE-BY-SIDE 3D ORGANISM COMPARISON MOTION + EXPRESSING THE HEURISTIC IN SQL
        act4_title = cfg.get("act4_title", "Act IV: Watching the Heuristic Align Two Extremophile Organisms in 3D")
        draw_smart_title(draw, act4_title, 56, 98, W - 112, PRIMARY_BLUE)

        # Left half: Live SQL Terminal expressing the heuristic cleanly
        draw.rounded_rectangle([56, 168, 980, 898], radius=18, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw.ellipse([84, 192, 104, 212], fill=ACCENT_PINK)
        draw.ellipse([116, 192, 136, 212], fill=AMBER)
        draw.ellipse([148, 192, 168, 212], fill=PRIMARY_TEAL)
        draw.text((190, 190), "Expressing the Heuristic as a Declarative Query", fill=TEXT_MUTED, font=FONT_MONO_SM)
        draw.line([(56, 228), (980, 228)], fill=(45, 55, 78), width=2)

        sql_lines = cfg.get("sql_lines", [])
        visible_lines = max(1, min(len(sql_lines), int((local_t / 0.28) * len(sql_lines)) + 1))
        sy = 246
        for idx, line in enumerate(sql_lines[:visible_lines]):
            cursor_str = " █" if (idx == visible_lines - 1 and local_t < 0.28 and int(t * 6) % 2 == 0) else ""
            line_fit = fit_single_line(draw, line + cursor_str, 880, FONT_MONO_SM)
            draw.text(
                (80, sy),
                line_fit,
                fill=SECONDARY_PURPLE if any(k in line for k in ("SELECT", "FROM", "WHERE", "WITH", "JOIN")) else TEXT_MAIN,
                font=FONT_MONO_SM,
            )
            sy += 28

        if local_t >= 0.25:
            res_lines = cfg.get("sql_results", [])
            ry = 520
            for idx, rline in enumerate(res_lines):
                color = EMERALD if idx == len(res_lines) - 1 else AMBER if idx in (3, 4) else TEXT_MAIN
                draw.text((80, ry), fit_single_line(draw, rline, 880, FONT_MONO_SM), fill=color, font=FONT_MONO_SM)
                ry += 27

            draw.rounded_rectangle([80, 765, 956, 874], radius=12, fill=PANEL_BG, outline=EMERALD, width=2)
            draw.text((100, 780), "Heuristic Realization:", fill=EMERALD, font=FONT_SM)
            draw_wrapped_box(draw, cfg["cs_mechanism"], 100, 810, 836, FONT_XS, TEXT_MAIN, line_h=22, max_lines=2)

        # Right half: Live Synchronized 3D Dual-Organism Motion
        draw.rounded_rectangle([1008, 168, 1864, 898], radius=18, fill=PANEL_DARK, outline=EMERALD, width=2)
        draw_dual_organism_comparison_stage(draw, (1008, 168, 1864, 898), pair_cfg, t, local_t, compact=True)

    else:
        # ACT V: REAL-WORLD SCIENTIFIC DISCOVERY & RESEARCH TAKEAWAYS
        act5_title = cfg.get("act5_title", "Act V: From Mathematical Insight to Extremophile Discovery")
        draw_smart_title(draw, act5_title, 56, 98, W - 112, PRIMARY_BLUE)

        # Left Panel: What the Heuristic Unlocked
        draw.rounded_rectangle([56, 168, 980, 898], radius=18, fill=PANEL_BG, outline=EMERALD, width=2)
        draw.text((84, 196), "WHAT THIS HEURISTIC REVEALS IN NATURE", fill=EMERALD, font=FONT_H2)

        takeaways = cfg.get(
            "takeaways",
            [
                ("1. Nature Conserves 3D Physics, Not 1D Spelling", "Across 3.5 billion years, extremophiles rewrite 70% of their surface sequence to survive salt, acid, or heat while locking the catalytic pocket within 0.5 Å."),
                ("2. Good Heuristics Expose Hidden Invariants", "Whether it is 24-ray SAS concavity, reference-centered reaction vectors, or ZCA whitening, the right math makes false positives collapse to zero."),
                ("3. Interactive Research in Your Browser", f"Rotate these aligned organisms in 3D at {BLOG_URL}/pg-bio-videos/ and explore the open-source research engine at {REPO_URL}."),
            ],
        )
        ty = 254
        for t_head, t_body in takeaways[:3]:
            draw.rounded_rectangle([84, ty, 952, ty + 190], radius=14, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
            draw.text((108, ty + 16), fit_single_line(draw, t_head, 820, FONT_BODY), fill=AMBER, font=FONT_BODY)
            draw_wrapped_box(draw, t_body, 108, ty + 56, 820, FONT_SM, TEXT_MAIN, line_h=28, max_lines=4)
            ty += 208

        # Right Panel: Dual Organism 3D Stage in Brine-Armor + Pocket Highlight Mode
        draw.rounded_rectangle([1008, 168, 1864, 898], radius=18, fill=PANEL_DARK, outline=SECONDARY_PURPLE, width=2)
        draw_dual_organism_comparison_stage(draw, (1008, 168, 1864, 898), pair_cfg, t, local_t, compact=False)

    draw_caption_and_footer(draw, narration, local_t, series_tag)
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
            "stage_header": "3D MIRRORED ORGANISM COMPARISON (0.44 Å KABSCH ALIGNED)",
            "org1_label": "Bait: A0ABD5UZ94 (Mesophilic Reference Hydrolase)",
            "org2_label": "Discovery: D3SWG7 (Natrialba magadii · Lake Magadi pH 10.5)",
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
            "Look at the two rotating organisms on the right. As they turn in lockstep, notice how their outer surfaces are completely different to survive extreme salt and pH, yet the glowing emerald catalytic pocket in the center aligns to less than half an Angstrom.",
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
    print(f"Voiceover ready ({total_sec:.1f}s across 5 acts). Rendering {total_frames} synchronized 1080p frames with 3D PDB organism motion...")

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
    print(f"Rendered 1080p MP4 with 3D Organism Motion ({total_frames} frames, {total_sec:.1f}s): {out_path}")


if __name__ == "__main__":
    main()
