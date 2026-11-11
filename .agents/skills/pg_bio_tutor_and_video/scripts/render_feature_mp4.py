#!/usr/bin/env python3
"""
Programmatic 1080p MP4 Video + Neural Voiceover Renderer for pg_bio Episodes.
Styled with the exact `ideia.me` dark theme palette and persistent footer bar
(`https://ideia.me` | `https://github.com/jonatas/pg_bio`).

Refined Flow Features:
  1. Dynamic `textbbox`-measured header pills + 5-Act Storyboard Stepper & `MM:SS` timer
     (zero text clipping or badge overflow).
  2. Full-coverage sentence-chunked captions synchronized to `act_local_progress`
     (never drops or truncates sentences).
  3. Episode-specific animated right-hand visual stages in Acts I & III:
       - `central_dogma_gist` : 1D Chromosome GiST Safe Harbor + Codon -> Peptide pI curve
       - `morton_3d_pocket`   : 3D Catalytic Zn2+ Pocket + Morton Z-Order bit interleaving
       - `csr_sparse_graph`   : Dense 90%-Zero Matrix -> CSR Arrays + Allosteric Network
       - `hnsw_twilight_zone` : 1280-D HNSW Graph Funnel + Twilight-Zone Homology Curve
  4. Animated packet flow in Act II and interactive psql cursor + EXPLAIN ANALYZE in Act IV.
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
from typing import Any, Dict, List, Tuple

from PIL import Image, ImageDraw, ImageFont

W, H = 1920, 1080
FPS = 12
DEFAULT_VOICE = "en-US-BrianMultilingualNeural"

# Exact `../ideia.me` Theme Palette
BG_COLOR = (15, 17, 23)           # --bg-color: #0f1117
PANEL_BG = (22, 27, 39)           # --card-bg-color: #161b27
PANEL_DARK = (26, 31, 46)         # --code-bg-color: #1a1f2e
PRIMARY_BLUE = (92, 158, 255)     # --primary-color (dark): #5c9eff
BRAND_BLUE = (67, 97, 238)        # --primary-color: #4361ee
SECONDARY_PURPLE = (177, 88, 252) # --secondary-color (dark): #b158fc
ACCENT_PINK = (255, 90, 157)      # --accent-color (dark): #ff5a9d
PRIMARY_TEAL = (13, 148, 136)     # --primary-teal: #0d9488
LIGHT_BLUE = (14, 165, 233)       # --primary-light-blue: #0ea5e9
PG_TOPIC_BLUE = (51, 103, 145)    # .topic-postgresql: #336791
TEXT_MAIN = (204, 214, 246)       # --text-color: #ccd6f6
TEXT_MUTED = (136, 146, 176)      # --text-muted / --gray-color: #8892b0
AMBER = (255, 193, 7)             # .topic-analytics: #FFC107

BLOG_URL = "https://ideia.me"
REPO_URL = "https://github.com/jonatas/pg_bio"

ACT_STEPS = [
    "I. Biotech Hook",
    "II. Python Wall",
    "III. Rosetta Bridge",
    "IV. Live SQL",
    "V. Benchmark",
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
            "/System/Library/Fonts/SFNSMono.ttf",
            "/System/Library/Fonts/SFNS.ttf",
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


FONT_H1 = load_font(38)
FONT_H2 = load_font(26)
FONT_BODY = load_font(21)
FONT_SM = load_font(17)
FONT_XS = load_font(15)
FONT_CAPTION = load_font(21)
FONT_MONO = load_font(20, mono=True)
FONT_MONO_SM = load_font(17, mono=True)


POCKET_ATOMS = [
    {"label": "ZN2+ (Cofactor)", "pos": (0.0, 0.0, 0.0), "color": AMBER, "r": 21},
    {"label": "HOH (2.05Å)", "pos": (-0.4, -1.6, 1.1), "color": ACCENT_PINK, "r": 15},
    {"label": "HIS94 ND1", "pos": (-1.5, 1.1, -0.8), "color": PRIMARY_BLUE, "r": 13},
    {"label": "HIS94 CA", "pos": (-2.9, 1.9, -1.4), "color": SECONDARY_PURPLE, "r": 10},
    {"label": "HIS96 ND1", "pos": (1.6, 0.9, -1.0), "color": PRIMARY_BLUE, "r": 13},
    {"label": "HIS96 CA", "pos": (2.8, 1.8, -1.6), "color": SECONDARY_PURPLE, "r": 10},
    {"label": "GLU106 OE1", "pos": (0.2, -2.3, -1.2), "color": PRIMARY_TEAL, "r": 14},
    {"label": "GLU106 CA", "pos": (0.5, -3.6, -2.0), "color": SECONDARY_PURPLE, "r": 10},
    {"label": "Substrate C=O", "pos": (-1.4, -2.4, 2.2), "color": LIGHT_BLUE, "r": 14},
]

POCKET_BONDS = [
    (0, 1, True),
    (0, 2, True),
    (0, 4, True),
    (1, 6, True),
    (1, 8, True),
    (2, 3, False),
    (4, 5, False),
    (6, 7, False),
]


def text_width(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.FreeTypeFont) -> int:
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0]


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


def default_payload(args: argparse.Namespace) -> Dict[str, Any]:
    return {
        "title": args.title,
        "layer": args.layer,
        "visual_mode": "morton_3d_pocket",
        "bio_problem": args.bio_problem,
        "cs_mechanism": args.cs_mechanism,
        "speedup": args.speedup,
        "subtitle": "Why Software & Database Engineers Already Know the Math of Biotech",
        "bullets": [
            "• Zn2+ cofactor must coordinate HOH within <= 2.40 Å",
            "• Bürgi-Dunitz nucleophilic attack angle: 100° to 115°",
            "• Zero Buried Unsatisfied Polar Atoms (BUNS) in core",
            "• Kekulé, Aromatic & Salt forms must hash identically",
        ],
        "rosetta_mappings": [
            [
                "64-Bit Morton Z-Order Curve (z_order_encode)",
                "Interleaves 21-bit (X, Y, Z) floats into a single BIGINT B-Tree key",
                "Extracts 3.5 Å 3D protein binding pockets in O(log N) without KD-trees",
            ],
            [
                "Weisfeiler-Lehman 64-Bit Hash + VF2 (@>)",
                "Hückel 4n+2 π-electron ring perception + canonical neighbor hashing",
                "Unifies Kekulé, aromatic & salted SMILES with 100% recall in 0.67 ms",
            ],
            [
                "3x3 Jacobi SVD + Fibonacci Sphere Grid",
                "Zero-allocation Kabsch RMSD + 5 Å voxel hash SASA/BUNS audit",
                "Filters RFdiffusion3 / Boltz-2 designs for catalytic water & zero BUNS",
            ],
        ],
        "sql_lines": [
            "bio_demo=# SELECT",
            "    design_id,",
            "    molecule_wl_hash('C1=CC=CC=C1') = molecule_wl_hash('c1ccccc1') AS kekule_eq_aromatic,",
            "    theozyme_rmsd(active_site_coords, '4qhp_zn_pentad')            AS zn_pentad_rmsd,",
            "    bio_angle(water_o, substrate_c, substrate_o)                   AS burgi_dunitz_deg,",
            "    has_coordinated_water_shell(active_site_coords, 'ZN', 2.4)     AS water_shell_ok",
            "FROM metalloprotease_designs",
            "WHERE theozyme_rmsd(active_site_coords, '4qhp_zn_pentad') <= 0.90;",
        ],
        "sql_results": [
            "┌──────────────────┬────────────────────┬────────────────┬──────────────────┬────────────────┐",
            "│    design_id     │ kekule_eq_aromatic │ zn_pentad_rmsd │ burgi_dunitz_deg │ water_shell_ok │",
            "├──────────────────┼────────────────────┼────────────────┼──────────────────┼────────────────┤",
            "│ Design_ZnO_042   │ t                  │          0.412 │            108.4 │ t              │",
            "│ Design_ZnO_119   │ t                  │          0.589 │            111.2 │ t              │",
            "└──────────────────┴────────────────────┴────────────────┴──────────────────┴────────────────┘",
            f"Time: 0.673 ms  ({args.speedup})",
        ],
        "narrations": [
            (
                f"Welcome to pg_bio: {args.title}. "
                f"In wet-lab synthetic biology, {args.bio_problem} "
                "If you already understand B-trees, spatial indexes, and graph hashes, you already know the math required to fix this."
            ),
            (
                "Why do traditional Python pipelines hit a wall? Because pulling hundreds of thousands of biological structures "
                "and sequences over a TCP socket into Python spends 99 percent of your CPU time on serialization and memory allocation, "
                "whereas pg_bio executes compiled Rust directly on PostgreSQL 8-kilobyte shared buffer pages."
            ),
            (
                f"Here is the systems-to-biotech Rosetta bridge: {args.cs_mechanism} "
                "Every biological invariant maps directly to a high-performance PostgreSQL index or zero-allocation Rust loop."
            ),
            (
                "Watch what happens inside psql. In a single declarative SQL query, we push multi-constraint biological "
                "filtering and physical verification directly to the database engine in sub-millisecond time."
            ),
            (
                f"The result is verified across the benchmark: {args.speedup}. "
                "Explore all discoveries at ideia dot me and star the pg_bio repository on GitHub."
            ),
        ],
    }


def sanitize_for_tts(text: str) -> str:
    """
    Normalizes arrows, dashes, and compound tokens into natural spoken punctuation
    so TTS engines never pronounce 'dash', 'hyphen', or 'greater than' (e.g.,
    'Dock -> Cut -> Transform' or 'Dock → Cut → Transform' -> 'Dock, Cut, Transform').
    """
    spoken = text
    # 1. Convert ASCII/Unicode arrows and em-dashes into comma pauses
    spoken = re.sub(r"\s*(?:<->|<=>|-->|->|<-|→|←|↔|⇒|⟶|—|--)\s*", ", ", spoken)
    # 2. Convert spaced hyphens / en-dashes (' - ' or ' – ') into comma pauses
    spoken = re.sub(r"\s+[-–]\s+", ", ", spoken)
    # 3. Convert numeric ranges ('100-115' or '3.5–5.0') into '100 to 115'
    spoken = re.sub(r"(\d)\s*[-–]\s*(\d)", r"\1 to \2", spoken)
    # 4. Replace remaining intra-word hyphens ('wet-lab', 'B-trees', 'multi-stage') with spaces
    spoken = re.sub(r"(?<=\w)[-–](?=\w)", " ", spoken)
    # 5. Clean up any stray standalone hyphens/dashes
    spoken = re.sub(r"[-–—]+", ", ", spoken)
    # 6. Normalize duplicate commas or colon-comma sequences and whitespace
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

        communicate = edge_tts.Communicate(spoken_text, voice, rate="+6%")
        await communicate.save(str(out_mp3))
    except Exception:
        aiff_path = out_mp3.with_suffix(".aiff")
        subprocess.run(["say", "-v", "Samantha", "-r", "195", "-o", str(aiff_path), spoken_text], check=True)
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
    progress: float,
    elapsed_sec: float,
    total_sec: float,
) -> None:
    draw.rectangle([0, 0, W, 84], fill=(11, 13, 19))
    draw.line([(0, 84), (W, 84)], fill=(45, 55, 78), width=2)
    draw.rectangle([0, 81, int(W * progress), 84], fill=PRIMARY_BLUE)

    # Left Badge: dynamically measured so text never overflows
    badge_text = fit_single_line(draw, f"pg_bio • {layer_str}", 560, FONT_XS)
    bw = text_width(draw, badge_text, FONT_XS) + 36
    draw.rounded_rectangle([48, 18, 48 + bw, 62], radius=16, fill=PG_TOPIC_BLUE, outline=PRIMARY_BLUE, width=2)
    draw.text((66, 29), badge_text, fill=(255, 255, 255), font=FONT_XS)

    # Center 5-Act Storyboard Stepper
    step_x = 670
    for idx, label in enumerate(ACT_STEPS):
        sw = text_width(draw, label, FONT_XS) + 24
        is_cur = idx == act_idx
        is_done = idx < act_idx
        fill_c = PANEL_BG if not is_cur else (37, 99, 235)
        out_c = PRIMARY_BLUE if is_cur else PRIMARY_TEAL if is_done else (50, 60, 82)
        txt_c = (255, 255, 255) if is_cur else LIGHT_BLUE if is_done else TEXT_MUTED
        draw.rounded_rectangle([step_x, 20, step_x + sw, 60], radius=12, fill=fill_c, outline=out_c, width=2)
        draw.text((step_x + 12, 30), label, fill=txt_c, font=FONT_XS)
        step_x += sw + 10

    # Right Timer Pill
    timer_text = f"🔊 {format_mm_ss(elapsed_sec)} / {format_mm_ss(total_sec)}"
    tw = text_width(draw, timer_text, FONT_XS) + 32
    draw.rounded_rectangle([W - 48 - tw, 18, W - 48, 62], radius=16, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
    draw.text((W - 48 - tw + 16, 29), timer_text, fill=TEXT_MAIN, font=FONT_XS)


def draw_caption_and_footer(draw: ImageDraw.ImageDraw, narration: str, act_local_progress: float) -> None:
    chunks = split_narration_chunks(narration, max_words_per_chunk=15)
    c_idx = min(len(chunks) - 1, int(act_local_progress * len(chunks)))
    active_chunk = chunks[c_idx]

    # Caption container (y: 916 .. 1012)
    draw.rounded_rectangle([110, 916, W - 110, 1012], radius=14, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
    # Chunk progress dots inside caption box top-right
    dot_x = W - 135 - len(chunks) * 16
    for i in range(len(chunks)):
        c = PRIMARY_BLUE if i == c_idx else PRIMARY_TEAL if i < c_idx else (55, 65, 88)
        draw.ellipse([dot_x + i * 16, 926, dot_x + i * 16 + 8, 934], fill=c)

    wrapped = textwrap.wrap(active_chunk, width=88)[:2]
    cy = 934 if len(wrapped) > 1 else 950
    for line in wrapped:
        lw = text_width(draw, line, FONT_CAPTION)
        draw.text(((W - lw) // 2, cy), line, fill=TEXT_MAIN, font=FONT_CAPTION)
        cy += 30

    # Persistent ideia.me + pg_bio Footer Bar (y: 1024 .. 1080)
    draw.rectangle([0, 1024, W, H], fill=(11, 13, 19))
    draw.line([(0, 1024), (W, 1024)], fill=BRAND_BLUE, width=2)

    draw.text((56, 1040), f"🌐 Blog & 3D Video Library: {BLOG_URL}/pg-bio-videos/", fill=PRIMARY_BLUE, font=FONT_MONO_SM)
    center_tag = "🧬 Hello World, pg_bio! — Multiomics in PostgreSQL"
    cw = text_width(draw, center_tag, FONT_MONO_SM)
    draw.text(((W - cw) // 2, 1040), center_tag, fill=TEXT_MUTED, font=FONT_MONO_SM)

    right_tag = f"🐙 Repo: {REPO_URL}"
    rw = text_width(draw, right_tag, FONT_MONO_SM)
    draw.text((W - rw - 56, 1040), right_tag, fill=SECONDARY_PURPLE, font=FONT_MONO_SM)


def draw_visual_stage(
    draw: ImageDraw.ImageDraw,
    visual_mode: str,
    box: Tuple[int, int, int, int],
    t: float,
    local_t: float,
    show_morton_sphere: bool,
) -> None:
    x1, y1, x2, y2 = box
    cx = (x1 + x2) // 2
    cy = (y1 + y2) // 2 + 20

    if visual_mode == "central_dogma_gist":
        # 1. Top half: 1D Chromosome int4range GiST Safe Harbor Diagram
        draw.text((x1 + 32, y1 + 26), "1D GENOMIC SAFE HARBOR (int4range GiST)", fill=PRIMARY_BLUE, font=FONT_H2)
        track_y = y1 + 125
        draw.line([(x1 + 48, track_y), (x2 - 48, track_y)], fill=TEXT_MUTED, width=4)

        # Promoter block
        p_label = "Promoter [8k..9.2k]"
        draw.rounded_rectangle([x1 + 52, track_y - 26, x1 + 266, track_y + 26], radius=10, fill=PRIMARY_TEAL, outline=(255, 255, 255), width=2)
        draw.text((x1 + 52 + (214 - text_width(draw, p_label, FONT_XS)) // 2, track_y - 10), p_label, fill=(255, 255, 255), font=FONT_XS)

        # Animated Safe Harbor insertion block
        sh_label = "SAFE HARBOR (< 5kb)"
        pulse = int(4 * math.sin(t * 4.0))
        draw.rounded_rectangle(
            [x1 + 286 - pulse, track_y - 30 - pulse, x1 + 522 + pulse, track_y + 30 + pulse],
            radius=10,
            fill=BRAND_BLUE,
            outline=AMBER,
            width=3,
        )
        draw.text((x1 + 286 + (236 - text_width(draw, sh_label, FONT_XS)) // 2, track_y - 10), sh_label, fill=AMBER, font=FONT_XS)

        # Oncogene / Tumor Suppressor block (avoided!)
        on_label = "Oncogene (NOT &&)"
        draw.rounded_rectangle([x1 + 542, track_y - 26, x2 - 48, track_y + 26], radius=10, fill=(90, 24, 48), outline=ACCENT_PINK, width=2)
        draw.text((x1 + 542 + ((x2 - 48 - (x1 + 542)) - text_width(draw, on_label, FONT_XS)) // 2, track_y - 10), on_label, fill=ACCENT_PINK, font=FONT_XS)

        # 2. Middle: DNA -> Codon -> Amino Acid Compiler Pipeline
        comp_y = y1 + 230
        draw.rounded_rectangle([x1 + 36, comp_y, x2 - 36, comp_y + 185], radius=14, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
        draw.text((x1 + 60, comp_y + 16), "CENTRAL DOGMA COMPILER (bio_translate_dna)", fill=SECONDARY_PURPLE, font=FONT_BODY)

        codons = [("ATG", "Met (M)"), ("TTC", "Phe (F)"), ("GAA", "Glu (E)"), ("GGT", "Gly (G)"), ("TTC", "Phe (F)")]
        active_c = int(t * 2.5) % len(codons)
        for idx, (dna_c, aa_c) in enumerate(codons):
            bx = x1 + 60 + idx * 145
            is_hi = idx == active_c
            draw.rounded_rectangle([bx, comp_y + 62, bx + 130, comp_y + 106], radius=8, fill=PANEL_DARK, outline=AMBER if is_hi else PRIMARY_BLUE, width=2)
            draw.text((bx + 36, comp_y + 74), dna_c, fill=AMBER if is_hi else TEXT_MAIN, font=FONT_MONO)
            draw.rounded_rectangle([bx, comp_y + 122, bx + 130, comp_y + 166], radius=8, fill=PG_TOPIC_BLUE if is_hi else PANEL_DARK, outline=PRIMARY_TEAL, width=2)
            draw.text((bx + 22, comp_y + 134), aa_c, fill=(255, 255, 255), font=FONT_XS)

        # 3. Bottom: Henderson-Hasselbalch Isoelectric Point (pI = 4.68) Titration Curve
        pi_y = y1 + 440
        draw.rounded_rectangle([x1 + 36, pi_y, x2 - 36, y2 - 28], radius=14, fill=PANEL_BG, outline=PRIMARY_BLUE, width=2)
        draw.text((x1 + 60, pi_y + 16), "ISOELECTRIC POINT BISECTION (pI = 4.68 in 200 ns)", fill=PRIMARY_BLUE, font=FONT_BODY)
        zero_y = pi_y + 130
        draw.line([(x1 + 80, zero_y), (x2 - 80, zero_y)], fill=TEXT_MUTED, width=2)
        draw.text((x1 + 60, zero_y - 24), "Net Charge = 0", fill=TEXT_MUTED, font=FONT_XS)

        pts = []
        for i in range(40):
            ph = 2.0 + (10.0 * i / 39.0)
            px = int(x1 + 80 + (x2 - x1 - 160) * (i / 39.0))
            charge = -math.tanh((ph - 4.68) * 0.65)
            py = int(zero_y - charge * 65)
            pts.append((px, py))
        draw.line(pts, fill=ACCENT_PINK, width=4)
        pi_px = int(x1 + 80 + (x2 - x1 - 160) * ((4.68 - 2.0) / 10.0))
        draw.ellipse([pi_px - 9, zero_y - 9, pi_px + 9, zero_y + 9], fill=AMBER, outline=(255, 255, 255), width=2)
        draw.text((pi_px + 14, zero_y - 28), "pI = 4.68 (Bisection Root)", fill=AMBER, font=FONT_SM)

    elif visual_mode == "csr_sparse_graph":
        draw.text((x1 + 32, y1 + 26), "DENSE 90%-ZERO MATRIX vs. CSR GRAPH", fill=PRIMARY_BLUE, font=FONT_H2)
        # Left mini-grid: 8x8 matrix mostly 0s
        gx0, gy0 = x1 + 48, y1 + 88
        draw.text((gx0, gy0 - 26), "Dense scRNA / Contact Matrix (90% 0s):", fill=ACCENT_PINK, font=FONT_XS)
        non_zero_coords = {(0, 2), (1, 5), (2, 0), (2, 6), (3, 4), (4, 3), (5, 1), (6, 2), (7, 7)}
        cell_s = 34
        for r in range(8):
            for c in range(8):
                cx0 = gx0 + c * cell_s
                cy0 = gy0 + r * cell_s
                nz = (r, c) in non_zero_coords
                draw.rectangle(
                    [cx0, cy0, cx0 + cell_s - 3, cy0 + cell_s - 3],
                    fill=PRIMARY_TEAL if nz else (20, 25, 36),
                    outline=PRIMARY_BLUE if nz else (40, 48, 68),
                )
                draw.text((cx0 + 10, cy0 + 7), "1" if nz else "0", fill=(255, 255, 255) if nz else (70, 80, 105), font=FONT_XS)

        # Right side of top: CSR Compact Arrays
        csr_x = gx0 + 8 * cell_s + 32
        draw.rounded_rectangle([csr_x, gy0, x2 - 36, gy0 + 8 * cell_s], radius=12, fill=PANEL_BG, outline=PRIMARY_TEAL, width=2)
        draw.text((csr_x + 20, gy0 + 16), "SparseAttentionMap (CSR):", fill=PRIMARY_TEAL, font=FONT_BODY)
        csr_lines = [
            "row_ptr: [0, 1, 2, 4, 5, 6, 7, 8, 9]",
            "col_idx: [2, 5, 0, 6, 4, 3, 1, 2, 7]",
            "values:  [-2.4, -1.8, -3.1, -2.9...]",
            "",
            "• Memory Slashed by >90%",
            "• O(nnz) Slice Scan in < 1 µs",
            "• Uncovers Cryptic Pocket Anchors:",
            "  Ile126, Val29, Leu66 (5-LOX)",
        ]
        cy_line = gy0 + 58
        for cl in csr_lines:
            draw.text((csr_x + 20, cy_line), cl, fill=AMBER if "Ile126" in cl else TEXT_MAIN, font=FONT_MONO_SM)
            cy_line += 26

        # Bottom: Allosteric Hotspot Network Graph
        net_y = y1 + 390
        draw.rounded_rectangle([x1 + 36, net_y, x2 - 36, y2 - 28], radius=14, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
        draw.text((x1 + 60, net_y + 16), "ALLOSTERIC HOTSPOT GRAPH (discover_cryptic_hotspots)", fill=SECONDARY_PURPLE, font=FONT_BODY)
        nodes = [
            ("AKBA Ligand", (cx, net_y + 145), AMBER, 24),
            ("Arg68 (Pocket)", (cx - 210, net_y + 95), PRIMARY_BLUE, 18),
            ("Glu134 (Lock)", (cx - 190, net_y + 215), PRIMARY_BLUE, 18),
            ("Ile126 (Cryptic!)", (cx + 195, net_y + 90), ACCENT_PINK, 20),
            ("Val29 (Cryptic!)", (cx + 225, net_y + 165), ACCENT_PINK, 20),
            ("Leu66 (Cryptic!)", (cx + 175, net_y + 230), ACCENT_PINK, 20),
        ]
        for idx in range(1, len(nodes)):
            nx, ny = nodes[idx][1]
            lx, ly = nodes[0][1]
            is_cryptic = idx >= 3
            draw.line([(lx, ly), (nx, ny)], fill=ACCENT_PINK if is_cryptic else PRIMARY_BLUE, width=3)
        for label, (nx, ny), col, rad in nodes:
            draw.ellipse([nx - rad, ny - rad, nx + rad, ny + rad], fill=col, outline=(255, 255, 255), width=2)
            draw.text((nx + rad + 10, ny - 10), label, fill=TEXT_MAIN, font=FONT_XS)

    elif visual_mode == "hnsw_twilight_zone":
        draw.text((x1 + 32, y1 + 26), "1280-D HNSW GRAPH & TWILIGHT-ZONE FUNNEL", fill=PRIMARY_BLUE, font=FONT_H2)
        # Top: 3-stage Query Planner Funnel
        stages = [
            ("1. 1280-D halfvec HNSW Index (<=>)", "2,340,573 Dark Proteins -> Top 15 in 0.4 ms", PRIMARY_BLUE, 720),
            ("2. Rolling-Buffer Smith-Waterman DP", "Scores only 15 survivors in O(N) L1 cache", SECONDARY_PURPLE, 540),
            ("3. LATERAL bio_fetch_uniprot() REST", "Live organism & taxonomy enrichment in SQL", PRIMARY_TEAL, 390),
        ]
        fy = y1 + 88
        for s_title, s_sub, col, w_box in stages:
            bx0 = cx - w_box // 2
            bx1 = cx + w_box // 2
            draw.rounded_rectangle([bx0, fy, bx1, fy + 82], radius=12, fill=PANEL_BG, outline=col, width=2)
            draw.text((bx0 + 20, fy + 12), s_title, fill=col, font=FONT_BODY)
            draw.text((bx0 + 20, fy + 44), s_sub, fill=TEXT_MAIN, font=FONT_XS)
            fy += 98

        # Bottom: Twilight-Zone Bell Curve (Sequence Identity vs 3D Structural Homology)
        bell_y = y1 + 395
        draw.rounded_rectangle([x1 + 36, bell_y, x2 - 36, y2 - 28], radius=14, fill=PANEL_BG, outline=AMBER, width=2)
        draw.text((x1 + 60, bell_y + 16), "TWILIGHT-ZONE BELL CURVE (15%–75% Sequence Identity)", fill=AMBER, font=FONT_BODY)
        base_y = y2 - 62
        draw.line([(x1 + 80, base_y), (x2 - 80, base_y)], fill=TEXT_MUTED, width=2)
        draw.text((x1 + 80, base_y + 10), "0% Identity (Noise)", fill=TEXT_MUTED, font=FONT_XS)
        draw.text((cx - 95, base_y + 10), "25%–50% TWILIGHT ZONE", fill=AMBER, font=FONT_XS)
        draw.text((x2 - 240, base_y + 10), "100% Clones", fill=TEXT_MUTED, font=FONT_XS)

        curve_pts = []
        for i in range(50):
            s = i / 49.0
            bell = 4.0 * s * (1.0 - s)
            px = int(x1 + 80 + (x2 - x1 - 160) * s)
            py = int(base_y - bell * 165)
            curve_pts.append((px, py))
        draw.line(curve_pts, fill=ACCENT_PINK, width=4)
        peak_x = cx - 40
        peak_y = int(base_y - 0.96 * 165)
        draw.ellipse([peak_x - 10, peak_y - 10, peak_x + 10, peak_y + 10], fill=AMBER, outline=(255, 255, 255), width=2)
        draw.text((peak_x + 16, peak_y - 12), "A0A5B9DCV6 (Cosine Dist: 0.068)", fill=AMBER, font=FONT_SM)

    else:
        # Default / Episode 2: 3D Rotating Catalytic Pocket + Morton Z-Order Bit Interleaving
        draw.text((x1 + 32, y1 + 26), "3D CATALYTIC POCKET & 64-BIT MORTON CURVE", fill=PRIMARY_BLUE, font=FONT_H2)
        if show_morton_sphere:
            r_bound = 195
            draw.ellipse([cx - r_bound, cy - 30 - r_bound, cx + r_bound, cy - 30 + r_bound], outline=PRIMARY_BLUE, width=3)
            draw.text((cx - r_bound + 14, cy - 30 - r_bound - 28), "Morton Z-Order 3.5 Å B-Tree Bounding Sphere", fill=PRIMARY_BLUE, font=FONT_XS)

        cos_a, sin_a = math.cos(t * 0.85), math.sin(t * 0.85)
        scale = 66.0
        projected = []
        for atom in POCKET_ATOMS:
            x, y, z = atom["pos"]
            rx = x * cos_a - z * sin_a
            rz = x * sin_a + z * cos_a
            ry = y * math.cos(0.35) - rz * math.sin(0.35)
            rz2 = y * math.sin(0.35) + rz * math.cos(0.35)
            persp = 1.0 / (1.0 + rz2 * 0.08)
            px = int(cx + rx * scale * persp)
            py = int((cy - 30) + ry * scale * persp)
            projected.append((px, py, rz2, persp, atom))

        for i, j, is_coord in POCKET_BONDS:
            px1, py1, _, _, _ = projected[i]
            px2, py2, _, _, _ = projected[j]
            draw.line([(px1, py1), (px2, py2)], fill=AMBER if is_coord else TEXT_MUTED, width=3 if is_coord else 5)

        for px, py, _, persp, atom in sorted(projected, key=lambda item: item[2], reverse=True):
            r = max(6, int(atom["r"] * persp))
            c = atom["color"]
            draw.ellipse([px - r - 4, py - r - 4, px + r + 4, py + r + 4], fill=(c[0] // 4, c[1] // 4, c[2] // 4))
            draw.ellipse([px - r, py - r, px + r, py + r], fill=c, outline=(255, 255, 255), width=2)
            draw.text((px + r + 8, py - 10), atom["label"], fill=TEXT_MAIN, font=FONT_XS)

        # Bottom Bit-Interleaving Strip
        draw.rounded_rectangle([x1 + 36, y2 - 86, x2 - 36, y2 - 22], radius=10, fill=PANEL_BG, outline=AMBER, width=2)
        draw.text(
            (x1 + 54, y2 - 66),
            "z_order_encode(3.95, -1.65, -2.15) -> Interleaved [Z|Y|X] BIGINT Key",
            fill=AMBER,
            font=FONT_MONO_SM,
        )


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
        draw.line([(gx, 84), (gx, 1024)], fill=(20, 24, 34), width=1)
    for gy in range(84, 1024, 80):
        draw.line([(0, gy), (W, gy)], fill=(20, 24, 34), width=1)

    title = fit_single_line(draw, cfg["title"], W - 128, FONT_H1)
    layer = cfg["layer"]
    bio_problem = cfg["bio_problem"]
    speedup = cfg["speedup"]
    visual_mode = cfg.get("visual_mode", "morton_3d_pocket")

    draw_header(draw, layer, act_idx, progress, t, total_sec)

    if act_idx == 0:
        draw.text((64, 102), title, fill=PRIMARY_BLUE, font=FONT_H1)
        sub = fit_single_line(draw, cfg.get("subtitle", "Hello World, pg_bio!"), W - 128, FONT_H2)
        draw.text((64, 154), sub, fill=SECONDARY_PURPLE, font=FONT_H2)

        # Left Card
        draw.rounded_rectangle([64, 206, 956, 898], radius=20, fill=PANEL_BG, outline=ACCENT_PINK, width=2)
        draw.text((96, 232), "🧬 THE BIOLOGICAL PROBLEM (WHY BIOTECH CARES)", fill=ACCENT_PINK, font=FONT_H2)
        y = draw_wrapped_box(draw, bio_problem, 96, 282, 825, FONT_BODY, TEXT_MAIN, line_h=34, max_lines=6)

        draw.rounded_rectangle([96, max(495, y + 16), 924, 872], radius=14, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw.text((120, max(512, y + 32)), "Core Multiomics Invariants in SQL:", fill=PRIMARY_BLUE, font=FONT_BODY)
        by = max(554, y + 74)
        for b in cfg.get("bullets", [])[:4]:
            b_fit = fit_single_line(draw, b, 780, FONT_SM)
            draw.text((120, by), b_fit, fill=TEXT_MAIN, font=FONT_SM)
            by += 34

        # Right Card: Episode-Specific Visual Stage
        draw.rounded_rectangle([996, 206, 1856, 898], radius=20, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw_visual_stage(draw, visual_mode, (996, 206, 1856, 898), t, local_t, show_morton_sphere=False)

    elif act_idx == 1:
        draw.text((64, 102), "Why Duct-Tape Bioinformatics Pipelines Hit a Wall", fill=PRIMARY_BLUE, font=FONT_H1)

        # Left Panel: Python Bottleneck
        draw.rounded_rectangle([64, 176, 920, 898], radius=20, fill=PANEL_BG, outline=ACCENT_PINK, width=3)
        draw.text((96, 204), "🐍 Traditional Python Stack (Socket + RAM Wall)", fill=ACCENT_PINK, font=FONT_H2)
        py_boxes = [
            ("1. PostgreSQL / Disk Files", "Millions of FASTA, VCF, PDB & CSR rows"),
            ("2. TCP Socket Serialization", "psycopg wire copy + UTF-8 string decoding"),
            ("3. Python Object Boxing", "PyObject + Pandas / NumPy heap allocations"),
            ("4. O(N²) Interpreter Loops", "GIL contention & CPU cache thrashing"),
        ]
        for idx, (b_title, b_sub) in enumerate(py_boxes):
            top = 272 + idx * 148
            draw.rounded_rectangle([104, top, 880, top + 114], radius=14, fill=PANEL_DARK, outline=ACCENT_PINK, width=2)
            draw.text((132, top + 18), b_title, fill=TEXT_MAIN, font=FONT_BODY)
            draw.text((132, top + 58), fit_single_line(draw, b_sub, 720, FONT_SM), fill=TEXT_MUTED, font=FONT_SM)

        # Right Panel: pg_bio Engine
        draw.rounded_rectangle([984, 176, 1856, 898], radius=20, fill=PANEL_BG, outline=PRIMARY_TEAL, width=3)
        draw.text((1016, 204), "🦀 pg_bio Native Rust Engine (Compute at Data)", fill=PRIMARY_TEAL, font=FONT_H2)
        pg_boxes = [
            ("1. Zero-IPC Shared Buffer Pointers", "pgrx reads 8KB Postgres pages directly in Rust"),
            ("2. Query-Planner Index Sieve", "Prunes 99.9% via GiST, B-Tree, CSR or HNSW"),
            ("3. Stack-Allocated LLVM SIMD Math", fit_single_line(draw, cfg["cs_mechanism"], 720, FONT_SM)),
            ("4. Only Winning Hits Cross the Wire", fit_single_line(draw, f"Verified: {speedup}", 720, FONT_SM)),
        ]
        for idx, (b_title, b_sub) in enumerate(pg_boxes):
            top = 272 + idx * 148
            active = local_t * 4 >= idx
            outline_c = PRIMARY_TEAL if active else (55, 68, 95)
            draw.rounded_rectangle([1024, top, 1816, top + 114], radius=14, fill=PANEL_DARK, outline=outline_c, width=3)
            draw.text((1052, top + 18), b_title, fill=LIGHT_BLUE if active else TEXT_MAIN, font=FONT_BODY)
            draw.text((1052, top + 58), b_sub, fill=TEXT_MAIN, font=FONT_SM)

    elif act_idx == 2:
        draw.text((64, 102), "Translating Database Internals into Molecular Biology", fill=PRIMARY_BLUE, font=FONT_H1)

        # Left side: Rosetta Mapping Cards with wrapped 2-line descriptions (no truncation!)
        draw.rounded_rectangle([64, 176, 1044, 898], radius=20, fill=PANEL_BG, outline=SECONDARY_PURPLE, width=2)
        draw.text((96, 204), "⚙️ SYSTEMS ENGINEERING   -->   🧬 BIOTECH REALITY", fill=SECONDARY_PURPLE, font=FONT_H2)

        mappings = cfg.get("rosetta_mappings", [])
        for idx, (m_title, m_cs, m_bio) in enumerate(mappings[:3]):
            top = 264 + idx * 202
            draw.rounded_rectangle([96, top, 1012, top + 186], radius=14, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
            draw.text((120, top + 14), fit_single_line(draw, m_title, 860, FONT_BODY), fill=PRIMARY_BLUE, font=FONT_BODY)
            draw_wrapped_box(draw, f"CS / Rust: {m_cs}", 120, top + 52, 865, FONT_SM, TEXT_MAIN, line_h=24, max_lines=2)
            draw_wrapped_box(draw, f"Biotech:   {m_bio}", 120, top + 112, 865, FONT_SM, LIGHT_BLUE, line_h=24, max_lines=2)

        # Right side: Episode-Specific Visual Stage with Bounding Sphere active
        draw.rounded_rectangle([1080, 176, 1856, 898], radius=20, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw_visual_stage(draw, visual_mode, (1080, 176, 1856, 898), t, local_t, show_morton_sphere=True)

    elif act_idx == 3:
        draw.text((64, 102), "Declarative Multiomics inside psql (PostgreSQL 18)", fill=PRIMARY_BLUE, font=FONT_H1)

        draw.rounded_rectangle([64, 176, 1856, 898], radius=20, fill=PANEL_DARK, outline=PRIMARY_BLUE, width=2)
        draw.ellipse([96, 202, 118, 224], fill=ACCENT_PINK)
        draw.ellipse([132, 202, 154, 224], fill=AMBER)
        draw.ellipse([168, 202, 190, 224], fill=PRIMARY_TEAL)
        draw.text((220, 200), "psql — pg_bio (PostgreSQL 18.0 + pgrx Rust Engine) — \\timing on", fill=TEXT_MUTED, font=FONT_MONO_SM)
        draw.line([(64, 242), (1856, 242)], fill=(45, 55, 78), width=2)

        sql_lines = cfg.get("sql_lines", [])
        # Type all SQL lines smoothly during the first 28% of Act IV
        visible_lines = max(1, min(len(sql_lines), int((local_t / 0.28) * len(sql_lines)) + 1))
        sy = 264
        for idx, line in enumerate(sql_lines[:visible_lines]):
            cursor_str = " █" if (idx == visible_lines - 1 and local_t < 0.28 and int(t * 6) % 2 == 0) else ""
            draw.text(
                (104, sy),
                line + cursor_str,
                fill=SECONDARY_PURPLE if any(k in line for k in ("SELECT", "FROM", "WHERE", "WITH", "JOIN")) else TEXT_MAIN,
                font=FONT_MONO,
            )
            sy += 32

        if local_t >= 0.26:
            res_lines = cfg.get("sql_results", [])
            ry = 535
            for idx, rline in enumerate(res_lines):
                color = PRIMARY_BLUE if idx == len(res_lines) - 1 else AMBER if idx in (3, 4) else TEXT_MAIN
                draw.text((104, ry), rline, fill=color, font=FONT_MONO)
                ry += 31

            # Query Planner / Shared Buffer Telemetry Box at the bottom of the terminal
            draw.rounded_rectangle([104, 780, 1816, 874], radius=12, fill=PANEL_BG, outline=PRIMARY_TEAL, width=2)
            draw.text(
                (130, 796),
                "⚡ EXPLAIN (ANALYZE, BUFFERS): Shared Hit Blocks = 12  |  Zero TCP Socket Copy  |  SIMD Rust Execution",
                fill=PRIMARY_TEAL,
                font=FONT_MONO_SM,
            )
            draw.text(
                (130, 832),
                fit_single_line(draw, f"Algorithmic Engine: {cfg['cs_mechanism']}", 1650, FONT_MONO_SM),
                fill=TEXT_MAIN,
                font=FONT_MONO_SM,
            )

    else:
        h1_bench = fit_single_line(draw, f"Verified Progress: {speedup}", W - 128, FONT_H1)
        draw.text((64, 102), h1_bench, fill=PRIMARY_BLUE, font=FONT_H1)

        draw.rounded_rectangle([64, 176, 1856, 540], radius=20, fill=PANEL_BG, outline=PRIMARY_TEAL, width=2)
        draw.text((108, 206), "WALL-CLOCK EXECUTION LATENCY (LOWER IS BETTER)", fill=TEXT_MAIN, font=FONT_H2)

        draw.text((108, 272), "🐍 Traditional Python Pipeline (Socket Pull + User-Space Math)", fill=ACCENT_PINK, font=FONT_BODY)
        py_w = int(1420 * min(1.0, local_t * 1.8))
        draw.rounded_rectangle([108, 312, 108 + max(20, py_w), 360], radius=12, fill=ACCENT_PINK)
        draw.text((1560, 318), "100% Baseline", fill=ACCENT_PINK, font=FONT_H2)

        draw.text((108, 394), "🦀 pg_bio Native Rust inside PostgreSQL 18", fill=PRIMARY_BLUE, font=FONT_BODY)
        draw.rounded_rectangle([108, 434, 136, 482], radius=12, fill=PRIMARY_TEAL)
        draw.text((160, 440), fit_single_line(draw, speedup, 1620, FONT_H2), fill=LIGHT_BLUE, font=FONT_H2)

        draw.rounded_rectangle([64, 568, 1856, 898], radius=20, fill=PANEL_DARK, outline=SECONDARY_PURPLE, width=2)
        draw.text((108, 600), "🚀 HELLO WORLD, pg_bio! — SERIES TAKEAWAYS", fill=SECONDARY_PURPLE, font=FONT_H2)
        laws = [
            "1. Move Compute to 8KB Shared-Buffer Pages — Never serialize millions of molecules over sockets.",
            "2. Use the Query Planner as a Biological Sieve — Prune with B-Tree/GiST/HNSW before O(N) physics.",
            f"3. Explore 3D views at {BLOG_URL}/pg-bio-videos/ & contribute at {REPO_URL}",
        ]
        ly = 668
        for law in laws:
            draw.text((108, ly), law, fill=TEXT_MAIN, font=FONT_BODY)
            ly += 48

    draw_caption_and_footer(draw, narration, local_t)
    return img


def main() -> None:
    parser = argparse.ArgumentParser(description="Render 1080p MP4 video + neural voiceover for a pg_bio episode.")
    parser.add_argument("--out", required=True, help="Output MP4 file path")
    parser.add_argument("--title", default="Hello World, pg_bio!")
    parser.add_argument("--layer", default="Multiomics Stack")
    parser.add_argument("--bio-problem", default="")
    parser.add_argument("--cs-mechanism", default="")
    parser.add_argument("--speedup", default="100x faster than Python")
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
    print(f"Checking / synthesizing 5-Act neural voiceover ({args.voice}) for {cfg['title']}...")
    full_vo_mp3, act_durations = generate_voiceover_and_srt(out_dir, narrations, args.voice, ffmpeg_bin, ffprobe_bin)

    act_frame_counts = [max(FPS * 3, int(math.ceil(d * FPS))) for d in act_durations]
    total_frames = sum(act_frame_counts)
    total_sec = total_frames / FPS
    print(f"Voiceover ready ({total_sec:.1f}s across 5 acts). Rendering {total_frames} synchronized 1080p frames...")

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
    print(f"Rendered 1080p MP4 with Refined Flow ({total_frames} frames, {total_sec:.1f}s): {out_path}")


if __name__ == "__main__":
    main()
