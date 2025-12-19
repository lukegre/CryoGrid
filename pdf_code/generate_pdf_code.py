#!/usr/bin/env python3
"""Generate searchable PDFs containing CryoGrid MATLAB source code.

Spec (from pdf_code/README.md):
- Each subfolder in /source gets its own PDF containing the source code as text.
- Special case: for /source/IO, each of its subfolders gets its own PDF.

This script walks the MATLAB source tree and produces PDFs with embedded text
(using reportlab), so the output is searchable and copy/paste friendly.

Typical usage:
  python pdf_code/generate_pdf_code.py

"""

from __future__ import annotations

import argparse
import datetime as _dt
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


def _iter_matlab_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        # MATLAB + common text docs in this repo
        if p.suffix.lower() in {".m", ".md", ".txt"}:
            files.append(p)
    # Stable ordering for reproducible PDFs
    files.sort(key=lambda x: str(x).lower())
    return files


def _safe_pdf_name(name: str) -> str:
    # Keep it readable and filesystem-safe.
    out = []
    for ch in name:
        if ch.isalnum() or ch in {"-", "_", "."}:
            out.append(ch)
        else:
            out.append("_")
    return "".join(out)


@dataclass(frozen=True)
class Unit:
    label: str
    source_dir: Path
    out_pdf: Path


def _build_units(source_root: Path, out_dir: Path) -> list[Unit]:
    """Determine which PDFs to generate."""
    units: list[Unit] = []

    for child in sorted([p for p in source_root.iterdir() if p.is_dir()], key=lambda p: p.name.lower()):
        label = child.name
        units.append(
            Unit(
                label=label,
                source_dir=child,
                out_pdf=out_dir / f"cryogrid-{_safe_pdf_name(label)}.pdf",
            )
        )

    return units


def _render_pdf(
    *,
    unit: Unit,
    inputs: list[Path],
    source_root: Path,
    title_prefix: str,
    overwrite: bool,
) -> None:
    """Create a single PDF with embedded text."""

    if unit.out_pdf.exists() and not overwrite:
        return

    # Import lazily so `--dry-run` works without reportlab installed.
    try:
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.units import inch
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont
        from reportlab.pdfgen.canvas import Canvas
    except ImportError as e:
        raise SystemExit(
            "Missing dependency: reportlab. Install it with `pip install reportlab` "
            "(ideally in a virtualenv)."
        ) from e

    unit.out_pdf.parent.mkdir(parents=True, exist_ok=True)

    # Use a mono font for code (DejaVuSansMono is often available; fall back to Courier).
    # Registering is optional; reportlab always has built-in 'Courier'.
    font_name = "Courier"
    try:
        # If system font exists, use it for nicer rendering.
        # Common on many systems, but not guaranteed.
        for probe in [
            "/Library/Fonts/DejaVuSansMono.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
            str(Path.home() / ".fonts" / "DejaVuSansMono.ttf"),
        ]:
            if Path(probe).exists():
                pdfmetrics.registerFont(TTFont("DejaVuSansMono", probe))
                font_name = "DejaVuSansMono"
                break
    except Exception:
        font_name = "Courier"

    pagesize = letter
    width, height = pagesize

    left = 0.65 * inch
    right = 0.55 * inch
    top = 0.6 * inch
    bottom = 0.6 * inch

    font_size = 8.5
    line_height = font_size * 1.25

    canvas = Canvas(str(unit.out_pdf), pagesize=pagesize)
    canvas.setTitle(f"{title_prefix} - {unit.label}")
    canvas.setAuthor("CryoGrid pdf_code generator")

    y = height - top

    def new_page() -> None:
        nonlocal y
        canvas.showPage()
        canvas.setFont(font_name, font_size)
        y = height - top

    canvas.setFont(font_name, font_size)

    header_lines = [
        f"{title_prefix}",
        f"Unit: {unit.label}",
        f"Generated: {_dt.datetime.now().isoformat(timespec='seconds')}",
        f"Source root: {source_root}",
        "",
    ]

    for hl in header_lines:
        if y < bottom + line_height:
            new_page()
        canvas.drawString(left, y, hl)
        y -= line_height

    for file_path in inputs:
        rel = file_path.relative_to(source_root)
        banner = f"===== {rel} ====="

        for chunk_line in ["", banner, ""]:
            if y < bottom + line_height:
                new_page()
            canvas.drawString(left, y, chunk_line)
            y -= line_height

        try:
            text = file_path.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            text = f"<ERROR reading file: {e}>\n"

        # Draw line-by-line, wrapping long lines to remain inside margins.
        max_chars = max(20, int((width - left - right) / (font_size * 0.60)))
        for raw_line in text.splitlines():
            # Keep tabs readable.
            line = raw_line.replace("\t", "    ")
            if line == "":
                wrapped = [""]
            else:
                wrapped = [line[i : i + max_chars] for i in range(0, len(line), max_chars)]

            for wline in wrapped:
                if y < bottom + line_height:
                    new_page()
                canvas.drawString(left, y, wline)
                y -= line_height

    canvas.save()


def _format_unit_summary(units: list[Unit]) -> str:
    return "\n".join(f"- {u.label}: {u.out_pdf.name}" for u in units)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate searchable PDFs containing MATLAB source code from /source.",
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "source",
        help="Path to the CryoGrid source directory (default: <repo>/source).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Output directory for PDFs (default: pdf_code/).",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="CryoGrid source code",
        help="PDF title prefix.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing PDFs.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be generated, but don't write PDFs.",
    )
    parser.add_argument(
        "--include",
        action="append",
        default=[],
        help="Optional: only generate units whose label contains this substring (repeatable).",
    )

    args = parser.parse_args(argv)

    source_root: Path = args.source.resolve()
    out_dir: Path = args.out.resolve()

    if not source_root.exists() or not source_root.is_dir():
        raise SystemExit(f"Source directory not found: {source_root}")

    units = _build_units(source_root, out_dir)

    if args.include:
        inc = [s.lower() for s in args.include]
        units = [u for u in units if any(s in u.label.lower() for s in inc)]

    if args.dry_run:
        print("Units to generate:")
        print(_format_unit_summary(units))
        return 0

    for unit in units:
        inputs = _iter_matlab_files(unit.source_dir)
        _render_pdf(
            unit=unit,
            inputs=inputs,
            source_root=source_root,
            title_prefix=args.title,
            overwrite=args.overwrite,
        )

    print(f"Generated {len(units)} PDF(s) into: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
