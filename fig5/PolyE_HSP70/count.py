#!/usr/bin/env python3
"""Interactive cell phenotype counting and per-class export tool.

This script aggregates single-cell crops from microscopy images that share the
same sample name, filters edge-touching cells, randomly samples a fixed number
per sample, and presents them in a blinded order for manual phenotype
classification. Labeled cells can be exported to an Excel summary and optional
colorized PNGs grouped by sample and class.

Input files must follow the convention:
    <sample>_<image_id>_<channel>.tif
where <channel> is chXX (e.g., ch02) or masks.
"""

from __future__ import annotations

import argparse
import os
import random
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import pandas as pd
except Exception as exc:  # pragma: no cover - dependency guard
    raise SystemExit("Please install pandas (pip install pandas)") from exc

try:
    import tifffile
except Exception as exc:  # pragma: no cover - dependency guard
    raise SystemExit("Please install tifffile (pip install tifffile)") from exc

try:
    from skimage.measure import regionprops
    from skimage.morphology import dilation, square
except Exception as exc:  # pragma: no cover - dependency guard
    raise SystemExit("Please install scikit-image (pip install scikit-image)") from exc

try:
    from PIL import Image
except Exception as exc:  # pragma: no cover - dependency guard
    raise SystemExit("Please install Pillow (pip install pillow)") from exc

import matplotlib.pyplot as plt


# ------------------------------ Data classes ------------------------------


@dataclass
class ChannelDisplay:
    name: str
    color: str


@dataclass
class ImageSet:
    sample: str
    image_id: str
    mask_path: str
    channel_paths: Dict[str, str]


@dataclass
class CellInfo:
    sample: str
    image_id: str
    label_id: int
    bbox: Tuple[int, int, int, int]
    mask_path: str
    channel_paths: Dict[str, str]


@dataclass
class LabelRecord:
    sample: str
    image_id: str
    label_id: int
    phenotype: str
    cell_index: int


@dataclass
class SaveSpec:
    base_dir: Path
    fmt: str = "png"
    panel_gap: int = 4
    mode: str = "panel"
    include_composite: bool = True


# ------------------------------ Constants ---------------------------------

_COLOR_TABLE: Dict[str, Tuple[float, float, float]] = {
    "gray": (1.0, 1.0, 1.0),
    "white": (1.0, 1.0, 1.0),
    "red": (1.0, 0.0, 0.0),
    "green": (0.0, 1.0, 0.0),
    "blue": (0.0, 0.0, 1.0),
    "cyan": (0.0, 1.0, 1.0),
    "magenta": (1.0, 0.0, 1.0),
    "yellow": (1.0, 1.0, 0.0),
    "orange": (1.0, 0.5, 0.0),
}


_EDGE_MARGIN = 2  # pixels from the image border required to keep a cell


# ------------------------------ Utilities ---------------------------------


class TiffCache:
    """Cache TIFF images in memory to avoid repeated disk reads."""

    def __init__(self) -> None:
        self._cache: Dict[str, np.ndarray] = {}

    def get(self, path: str) -> np.ndarray:
        if path not in self._cache:
            arr = tifffile.imread(path)
            self._cache[path] = arr
        return self._cache[path]


def parse_channel_specs(specs: Optional[Sequence[str]]) -> List[ChannelDisplay]:
    if not specs:
        return [ChannelDisplay(name="ch00", color="gray")]

    parsed: List[ChannelDisplay] = []
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"Channel spec '{spec}' must be formatted as <channel>=<color>")
        name, color = spec.split("=", 1)
        name = name.strip().lower()
        color = color.strip().lower()
        if not name:
            raise ValueError(f"Channel spec '{spec}' is missing a channel name")
        if color not in _COLOR_TABLE:
            raise ValueError(
                f"Unknown color '{color}' in spec '{spec}'. Supported colors: {', '.join(sorted(_COLOR_TABLE))}"
            )
        parsed.append(ChannelDisplay(name=name, color=color))
    return parsed


def discover_image_sets(
    root: Path,
    samples: Optional[Sequence[str]],
    required_channels: Iterable[ChannelDisplay],
) -> List[ImageSet]:
    sample_filter = set(s.lower() for s in samples) if samples else None
    required = {ch.name.lower() for ch in required_channels}
    grouped: Dict[Tuple[str, str], Dict[str, Dict[str, str]]] = {}

    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if not fname.lower().endswith(".tif"):
                continue
            stem = Path(fname).stem
            parts = stem.split("_")
            if len(parts) < 3:
                continue
            channel = parts[-1].lower()
            image_id = parts[-2]
            sample = "_".join(parts[:-2])
            if sample_filter and sample.lower() not in sample_filter:
                continue

            key = (sample, image_id)
            entry = grouped.setdefault(key, {"channels": {}, "mask": None})
            path = os.path.join(dirpath, fname)
            if channel == "masks":
                entry["mask"] = path
            elif channel.startswith("ch"):
                entry["channels"][channel] = path

    image_sets: List[ImageSet] = []
    for (sample, image_id), data in grouped.items():
        mask_path = data.get("mask")
        if not mask_path:
            continue
        channel_paths = data.get("channels", {})
        if not required.issubset(channel_paths.keys()):
            missing = sorted(required - channel_paths.keys())
            print(
                f"[WARN] Skipping {sample}_{image_id}: missing channels {missing}",
                file=sys.stderr,
            )
            continue
        image_sets.append(
            ImageSet(
                sample=sample,
                image_id=image_id,
                mask_path=mask_path,
                channel_paths=channel_paths,
            )
        )
    return image_sets


def ensure_2d(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 2:
        return arr
    if arr.ndim == 3:
        if arr.shape[0] == 1:
            return arr[0]
        if arr.shape[-1] == 1:
            return arr[..., 0]
        return arr.max(axis=0)
    raise ValueError(f"Unsupported array ndim={arr.ndim}")


def robust_percentile_norm(arr: np.ndarray, low: float, high: float) -> np.ndarray:
    lo = np.percentile(arr, low)
    hi = np.percentile(arr, high)
    if hi <= lo:
        hi = lo + 1e-6
    norm = (arr.astype(np.float32) - lo) / (hi - lo)
    return np.clip(norm, 0.0, 1.0)


def colorize(img01: np.ndarray, color: str) -> np.ndarray:
    rgb = np.zeros(img01.shape + (3,), dtype=np.float32)
    weights = _COLOR_TABLE[color]
    for idx in range(3):
        rgb[..., idx] = img01 * weights[idx]
    return rgb


def overlay_outline(base_rgb: np.ndarray, mask: np.ndarray, width: int = 1) -> np.ndarray:
    edges = mask.astype(bool) ^ dilation(mask.astype(bool), square(3))
    if width > 1:
        for _ in range(width - 1):
            edges = dilation(edges, square(3))
    outlined = base_rgb.copy()
    outlined[edges, 1] = 1.0
    outlined[edges, 2] = 1.0
    return outlined


def touches_border(coords: np.ndarray, shape: Tuple[int, int], margin: int = _EDGE_MARGIN) -> bool:
    if coords.size == 0:
        return False
    margin = max(0, margin)
    rows = coords[:, 0]
    cols = coords[:, 1]
    h, w = shape
    if margin == 0:
        return (rows == 0).any() or (rows == h - 1).any() or (cols == 0).any() or (cols == w - 1).any()
    return (rows.min() < margin or rows.max() >= h - margin or
            cols.min() < margin or cols.max() >= w - margin)


def compute_crop_bounds(
    bbox: Tuple[int, int, int, int],
    shape: Tuple[int, int],
    pad: int,
) -> Tuple[int, int, int, int]:
    minr, minc, maxr, maxc = bbox
    h, w = shape
    r0 = max(0, minr - pad)
    c0 = max(0, minc - pad)
    r1 = min(h, maxr + pad)
    c1 = min(w, maxc + pad)
    return r0, r1, c0, c1


def extract_cells(
    image_sets: Sequence[ImageSet],
    cache: TiffCache,
    channel_configs: Sequence[ChannelDisplay],
    pad: int,
) -> Tuple[Dict[str, List[CellInfo]], Dict[str, Dict[str, int]]]:
    per_sample_cells: Dict[str, List[CellInfo]] = {}
    stats: Dict[str, Dict[str, int]] = {}

    for image_set in image_sets:
        seg = ensure_2d(cache.get(image_set.mask_path))
        sample_stats = stats.setdefault(image_set.sample, {"total": 0, "kept": 0, "edge": 0})
        for region in regionprops(seg):
            if region.label == 0:
                continue
            sample_stats["total"] += 1
            if touches_border(region.coords, seg.shape):
                sample_stats["edge"] += 1
                continue
            bbox = region.bbox
            cell = CellInfo(
                sample=image_set.sample,
                image_id=image_set.image_id,
                label_id=region.label,
                bbox=bbox,
                mask_path=image_set.mask_path,
                channel_paths={cfg.name: image_set.channel_paths[cfg.name] for cfg in channel_configs},
            )
            per_sample_cells.setdefault(image_set.sample, []).append(cell)
            sample_stats["kept"] += 1

    for cells in per_sample_cells.values():
        cells.sort(key=lambda c: (c.image_id, c.label_id))
    return per_sample_cells, stats


def render_cell(
    cell: CellInfo,
    cache: TiffCache,
    channel_configs: Sequence[ChannelDisplay],
    pad: int,
    perc_low: float,
    perc_high: float,
    outline: bool,
) -> Tuple[np.ndarray, List[np.ndarray], np.ndarray]:
    seg_full = ensure_2d(cache.get(cell.mask_path))
    r0, r1, c0, c1 = compute_crop_bounds(cell.bbox, seg_full.shape, pad)
    seg_crop = seg_full[r0:r1, c0:c1]
    cell_mask = seg_crop == cell.label_id

    composite = np.zeros((seg_crop.shape[0], seg_crop.shape[1], 3), dtype=np.float32)
    per_channel_norms: List[np.ndarray] = []

    for cfg in channel_configs:
        channel_path = cell.channel_paths[cfg.name]
        channel_img = ensure_2d(cache.get(channel_path))
        channel_crop = channel_img[r0:r1, c0:c1]
        norm = robust_percentile_norm(channel_crop, perc_low, perc_high)
        per_channel_norms.append(norm)
        composite = np.maximum(composite, colorize(norm, cfg.color))

    composite = np.clip(composite, 0.0, 1.0)
    display = overlay_outline(composite, cell_mask) if outline else composite
    return display, per_channel_norms, composite


def select_cells_per_sample(
    per_sample_cells: Dict[str, List[CellInfo]],
    per_sample_count: int,
    seed: Optional[int],
) -> List[CellInfo]:
    rng = random.Random(seed)
    selected: List[CellInfo] = []

    for sample, cells in sorted(per_sample_cells.items()):
        if not cells:
            continue
        if len(cells) <= per_sample_count:
            selected.extend(cells)
        else:
            selected.extend(rng.sample(cells, per_sample_count))
    return selected


def sanitize_token(token: str) -> str:
    token = token.strip()
    if not token:
        return "unknown"
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", token)
    return safe or "unknown"


def build_panels(
    per_channel_norms: Sequence[np.ndarray],
    composite: np.ndarray,
    channel_configs: Sequence[ChannelDisplay],
    gap: int,
    include_composite: bool = True,
) -> np.ndarray:
    panels: List[np.ndarray] = [colorize(norm, cfg.color) for cfg, norm in zip(channel_configs, per_channel_norms)]
    if include_composite:
        panels.append(composite)

    heights = [p.shape[0] for p in panels]
    if len(set(heights)) != 1:
        raise ValueError("Panels must have equal heights")
    h = heights[0]
    sep = np.zeros((h, gap, 3), dtype=np.float32)

    strip_parts: List[np.ndarray] = []
    for idx, panel in enumerate(panels):
        strip_parts.append(panel)
        if gap > 0 and idx < len(panels) - 1:
            strip_parts.append(sep)
    strip = np.concatenate(strip_parts, axis=1)
    return strip


def save_cell_view(
    cell: CellInfo,
    phenotype: str,
    per_channel_norms: Sequence[np.ndarray],
    composite: np.ndarray,
    channel_configs: Sequence[ChannelDisplay],
    save_spec: SaveSpec,
    cell_index: int,
) -> Optional[Path]:
    sample_dir = sanitize_token(cell.sample)
    class_dir = sanitize_token(phenotype)
    out_dir = save_spec.base_dir / sample_dir / class_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    base_name = f"{sample_dir}_{cell.image_id}_label{cell.label_id:04d}_idx{cell_index:04d}"
    saved_paths: List[Path] = []

    if save_spec.mode in {"panel", "both"}:
        panels = build_panels(
            per_channel_norms,
            composite,
            channel_configs,
            save_spec.panel_gap,
            include_composite=save_spec.include_composite,
        )
        img_uint8 = (np.clip(panels, 0.0, 1.0) * 255).astype(np.uint8)
        image = Image.fromarray(img_uint8)
        suffix = "" if save_spec.mode == "panel" else "_panel"
        panel_path = out_dir / f"{base_name}{suffix}.{save_spec.fmt}"
        image.save(panel_path)
        saved_paths.append(panel_path)

    if save_spec.mode in {"separate", "both"}:
        for cfg, norm in zip(channel_configs, per_channel_norms):
            colored = colorize(norm, cfg.color)
            chan_uint8 = (np.clip(colored, 0.0, 1.0) * 255).astype(np.uint8)
            chan_img = Image.fromarray(chan_uint8)
            chan_name = sanitize_token(cfg.name)
            chan_path = out_dir / f"{base_name}_{chan_name}.{save_spec.fmt}"
            chan_img.save(chan_path)
            saved_paths.append(chan_path)
        if save_spec.include_composite:
            comp_uint8 = (np.clip(composite, 0.0, 1.0) * 255).astype(np.uint8)
            comp_img = Image.fromarray(comp_uint8)
            comp_path = out_dir / f"{base_name}_composite.{save_spec.fmt}"
            comp_img.save(comp_path)
            saved_paths.append(comp_path)

    return saved_paths[-1] if saved_paths else None


def backend_supports_interaction() -> bool:
    backend = plt.get_backend().lower()
    return not any(backend.startswith(prefix) for prefix in ("agg", "pdf", "svg", "ps", "template"))


def capture_keypress(fig: plt.Figure) -> Optional[str]:
    """Wait for a key press on the matplotlib figure and return a normalized value."""
    selection: Dict[str, Optional[str]] = {"key": None}

    def normalize(key: Optional[str]) -> Optional[str]:
        if not key:
            return None
        key = key.lower()
        if "+" in key:
            key = key.split("+")[-1]
        mapping = {
            " " : "s",
            "space": "s",
            "escape": "q",
            "esc": "q",
            "left": "b",
            "backspace": "b",
            "delete": "b",
        }
        key = mapping.get(key, key)
        if key in {"enter", "return"}:
            return None
        return key

    def on_key(event) -> None:
        key = normalize(getattr(event, "key", None))
        if key is not None:
            selection["key"] = key

    cid = fig.canvas.mpl_connect("key_press_event", on_key)
    try:
        while selection["key"] is None and plt.fignum_exists(fig.number):
            plt.pause(0.05)
    finally:
        fig.canvas.mpl_disconnect(cid)
    return selection["key"]



def run_labeling_session(
    cells: Sequence[CellInfo],
    cache: TiffCache,
    channel_configs: Sequence[ChannelDisplay],
    classes: Sequence[str],
    pad: int,
    perc_low: float,
    perc_high: float,
    outline: bool,
    save_spec: Optional[SaveSpec],
) -> List[LabelRecord]:
    if not classes:
        raise ValueError("At least one phenotype class must be provided")
    if len(classes) > 9:
        raise ValueError("Currently limited to 9 classes")

    key_to_class = {str(idx + 1): cls for idx, cls in enumerate(classes)}
    valid_inputs = set(key_to_class.keys()) | {"s", "b", "q"}
    records: List[LabelRecord] = []
    total = len(cells)
    idx = 0
    interactive_supported = backend_supports_interaction()

    ui_help = " / ".join(f"{digit}:{cls}" for digit, cls in key_to_class.items())
    ui_help += " | commands: s=skip, b=back, q=quit"
    console_help = ui_help
    if interactive_supported:
        console_help += " | 図ウィンドウ上でキーを押すと即時入力されます"

    while 0 <= idx < total:
        cell = cells[idx]
        _display_img, per_channel, composite = render_cell(
            cell, cache, channel_configs, pad, perc_low, perc_high, outline
        )

        tinted_channels = [colorize(norm, cfg.color) for cfg, norm in zip(channel_configs, per_channel)]

        choice: Optional[str] = None
        while True:
            ncols = max(1, len(tinted_channels))
            fig, axes = plt.subplots(1, ncols, figsize=(4 * ncols, 4))
            axes = np.atleast_1d(axes)
            for ax, cfg, img in zip(axes, channel_configs, tinted_channels):
                ax.imshow(img, vmin=0.0, vmax=1.0)
                ax.axis("off")
                ax.set_title(f"{cfg.name}\n({cfg.color})")
            fig.suptitle(f"Cell {idx + 1}/{total}\nChoose: {ui_help}")
            plt.tight_layout()
            plt.show(block=False)
            plt.pause(0.5)

            if interactive_supported:
                choice = capture_keypress(fig)
            else:
                choice = None

            if choice is None:
                try:
                    choice = input(f"Cell {idx + 1}/{total} > ").strip().lower()
                except (EOFError, KeyboardInterrupt):
                    print("\nSession interrupted. Saving progress...")
                    choice = "q"

            if plt.fignum_exists(fig.number):
                plt.close(fig)

            if not choice:
                print("入力が認識されませんでした。もう一度キーを押してください。")
                continue

            if choice in valid_inputs:
                break

            print(f"Unrecognized input '{choice}'. Type one of: {console_help}")
            choice = None

        if choice in key_to_class:
            phenotype = key_to_class[choice]
            record = LabelRecord(
                sample=cell.sample,
                image_id=cell.image_id,
                label_id=cell.label_id,
                phenotype=phenotype,
                cell_index=idx + 1,
            )
            records.append(record)
            if save_spec:
                save_cell_view(
                    cell,
                    phenotype,
                    per_channel,
                    composite,
                    channel_configs,
                    save_spec,
                    cell_index=idx + 1,
                )
            idx += 1
        elif choice == "s":
            idx += 1
        elif choice == "b":
            if records:
                records.pop()
                idx = max(0, idx - 1)
            else:
                print("Nothing to go back to.")
        elif choice == "q":
            print("Session ended early by user.")
            break

    plt.close("all")
    return records

def summarize_and_export(
    records: Sequence[LabelRecord],
    output_path: Path,
    include_raw: bool = True,
) -> None:
    if not records:
        print("No labeled cells to export.")
        return

    df = pd.DataFrame([
        {
            "sample": rec.sample,
            "image_id": rec.image_id,
            "label_id": rec.label_id,
            "phenotype": rec.phenotype,
            "cell_index": rec.cell_index,
        }
        for rec in records
    ])

    counts = (
        df.pivot_table(index="sample", columns="phenotype", values="label_id", aggfunc="count", fill_value=0)
        .sort_index()
        .sort_index(axis=1)
    )
    counts["total"] = counts.sum(axis=1)

    with pd.ExcelWriter(output_path) as writer:
        counts.to_excel(writer, sheet_name="summary")
        if include_raw:
            df.to_excel(writer, sheet_name="raw", index=False)

    print(f"Exported summary to {output_path}")
    print(counts)


# ------------------------------ CLI ---------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Blinded manual phenotype counting tool with per-class image export",
    )
    parser.add_argument("--data-root", required=True, help="Directory containing the TIFF files")
    parser.add_argument(
        "--samples",
        nargs="*",
        help="Sample names to include (defaults to all detected samples)",
    )
    parser.add_argument(
        "--channels",
        nargs="*",
        help="Channel display specification, e.g. ch02=green ch06=magenta",
    )
    parser.add_argument(
        "--classes",
        nargs="+",
        help="Phenotype class names, order determines hotkeys 1..N",
    )
    parser.add_argument(
        "--per-sample",
        type=int,
        default=200,
        help="Number of cells to randomly sample per sample",
    )
    parser.add_argument(
        "--pad",
        type=int,
        default=10,
        help="Padding (pixels) around each bounding box",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed for reproducible sampling",
    )
    parser.add_argument(
        "--percentile-low",
        type=float,
        default=1.0,
        help="Lower percentile for intensity normalization",
    )
    parser.add_argument(
        "--percentile-high",
        type=float,
        default=99.0,
        help="Upper percentile for intensity normalization",
    )
    parser.add_argument(
        "--outline",
        action="store_true",
        help="Overlay segmentation outline on the displayed composite view",
    )
    parser.add_argument(
        "--output",
        default="cell_counts.xlsx",
        help="Path to the Excel file where results will be saved",
    )
    parser.add_argument(
        "--save-dir",
        help="Directory where labeled cell images will be stored by sample/class",
    )
    parser.add_argument(
        "--save-format",
        default="png",
        help="Image format for exports (default: png)",
    )
    parser.add_argument(
        "--save-mode",
        choices=["panel", "separate", "both"],
        default="panel",
        help="Export mode: panel strip, per-channel images, or both",
    )
    parser.add_argument(
        "--no-save-composite",
        action="store_true",
        help="Skip saving composite RGB image when using separate/both modes",
    )
    parser.add_argument(
        "--save-gap",
        type=int,
        default=4,
        help="Gap (pixels) between channel panels (panel/both modes)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    channel_configs = parse_channel_specs(args.channels)

    if args.per_sample <= 0:
        parser.error("--per-sample must be positive")
    if args.pad < 0:
        parser.error("--pad must be non-negative")
    if args.percentile_low < 0 or args.percentile_high > 100 or args.percentile_low >= args.percentile_high:
        parser.error("Invalid percentile bounds")
    if args.save_gap < 0:
        parser.error("--save-gap must be non-negative")

    data_root = Path(args.data_root)
    if not data_root.exists():
        parser.error(f"Data root '{data_root}' does not exist")

    image_sets = discover_image_sets(data_root, args.samples, channel_configs)
    if not image_sets:
        print("No matching image sets found. Check data root, sample names, and channel specs.")
        return 1

    print(f"Discovered {len(image_sets)} image sets.")

    cache = TiffCache()
    per_sample_cells, stats = extract_cells(image_sets, cache, channel_configs, args.pad)

    for sample, sample_stats in sorted(stats.items()):
        print(
            f"Sample {sample}: total={sample_stats['total']}, edge_removed={sample_stats['edge']}, usable={sample_stats['kept']}"
        )

    sampled_cells = select_cells_per_sample(per_sample_cells, args.per_sample, args.seed)
    if not sampled_cells:
        print("No usable cells after filtering.")
        return 1

    rng = random.Random(args.seed)
    rng.shuffle(sampled_cells)

    print(f"Prepared {len(sampled_cells)} cells for labeling.")
    print("Cells will be shown in random order without sample names.")
    print("Press 1..N to assign phenotypes, or use commands s/b/q as listed.")

    save_spec = None
    if args.save_dir:
        base_dir = Path(args.save_dir)
        base_dir.mkdir(parents=True, exist_ok=True)
        save_spec = SaveSpec(
            base_dir=base_dir,
            fmt=args.save_format,
            panel_gap=args.save_gap,
            mode=args.save_mode,
            include_composite=not args.no_save_composite,
        )
        mode_desc = {
            "panel": "panel strips",
            "separate": "per-channel images",
            "both": "panel strips and per-channel images",
        }
        print(f"Labeled cell images will be saved under {base_dir} ({mode_desc[args.save_mode]}).")
        if args.save_mode in {"separate", "both"}:
            if save_spec.include_composite:
                print("Composite RGB images will also be saved.")
            else:
                print("Composite RGB images will be skipped (per --no-save-composite).")

    records = run_labeling_session(
        sampled_cells,
        cache,
        channel_configs,
        classes=args.classes or ["phenotype"],
        pad=args.pad,
        perc_low=args.percentile_low,
        perc_high=args.percentile_high,
        outline=args.outline,
        save_spec=save_spec,
    )

    summarize_and_export(records, Path(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
