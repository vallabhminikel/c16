#!/usr/bin/env python3
"""Export every slide of a PowerPoint deck to high-resolution PNG images.

LibreOffice's own PNG filter rasterizes at a fixed, low resolution, so this
goes through PDF (vector, resolution-independent) and rasterizes that instead:

    deck.pptx  --[LibreOffice headless]-->  deck.pdf  --[renderer]-->  slide-01.png

Two renderers are supported and auto-detected:
  * PyMuPDF  (pip install pymupdf)   -- no external binary needed
  * pdftoppm (Poppler)               -- `poppler-utils` / `brew install poppler`

Requires LibreOffice (`soffice`) on PATH either way.

Usage:
    python pptx_to_png.py deck.pptx
    python pptx_to_png.py deck.pptx -o renders/ --dpi 600
    python pptx_to_png.py *.pptx --dpi 300 --renderer pdftoppm
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# A 13.333 x 7.5 in widescreen slide at 600 dpi is 8000 x 4500 px (~15 MB
# per PNG). Above this the intermediate bitmaps get unwieldy and some
# renderers start failing on allocation rather than telling you why.
MAX_DPI = 1200


def find_soffice() -> str:
    """Locate the LibreOffice binary, including the usual macOS/Windows spots."""
    candidates = [
        "soffice",
        "libreoffice",
        "/Applications/LibreOffice.app/Contents/MacOS/soffice",
        r"C:\Program Files\LibreOffice\program\soffice.exe",
    ]
    for c in candidates:
        found = shutil.which(c) or (c if Path(c).exists() else None)
        if found:
            return found
    sys.exit(
        "LibreOffice not found. Install it and make sure `soffice` is on PATH.\n"
        "  macOS:  brew install --cask libreoffice\n"
        "  Debian: sudo apt install libreoffice-impress\n"
        "  Fedora: sudo dnf install libreoffice-impress"
    )


def pptx_to_pdf(pptx: Path, outdir: Path, timeout: int = 300) -> Path:
    """Convert a deck to PDF with headless LibreOffice."""
    soffice = find_soffice()
    # -env:UserInstallation isolates the profile so this works even when a
    # desktop LibreOffice is already running (otherwise the call silently
    # hands off to that instance and exits without converting).
    profile = (outdir / "lo-profile").resolve().as_uri()
    cmd = [
        soffice,
        f"-env:UserInstallation={profile}",
        "--headless",
        "--norestore",
        "--invisible",
        "--convert-to", "pdf:impress_pdf_Export",
        "--outdir", str(outdir),
        str(pptx),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    pdf = outdir / (pptx.stem + ".pdf")
    if not pdf.exists():
        sys.exit(
            f"PDF conversion failed for {pptx.name}\n"
            f"exit {proc.returncode}\nstdout: {proc.stdout}\nstderr: {proc.stderr}"
        )
    return pdf


def render_pymupdf(pdf: Path, outdir: Path, stem: str, dpi: int) -> list[Path]:
    import fitz  # PyMuPDF

    written = []
    with fitz.open(pdf) as doc:
        width = len(str(doc.page_count))
        for i, page in enumerate(doc, start=1):
            pix = page.get_pixmap(dpi=dpi, alpha=False)
            path = outdir / f"{stem}-{i:0{width}d}.png"
            pix.save(path)
            written.append(path)
    return written


def render_pdftoppm(pdf: Path, outdir: Path, stem: str, dpi: int) -> list[Path]:
    prefix = outdir / stem
    # -r sets both axes; pdftoppm zero-pads the index to the page count.
    subprocess.run(
        ["pdftoppm", "-png", "-r", str(dpi), str(pdf), str(prefix)],
        check=True,
        capture_output=True,
    )
    return sorted(outdir.glob(f"{stem}-*.png"))


def pick_renderer(requested: str) -> str:
    have_pymupdf = False
    try:
        import fitz  # noqa: F401
        have_pymupdf = True
    except ImportError:
        pass
    have_pdftoppm = shutil.which("pdftoppm") is not None

    if requested == "pymupdf" and not have_pymupdf:
        sys.exit("PyMuPDF not installed. Run: pip install pymupdf")
    if requested == "pdftoppm" and not have_pdftoppm:
        sys.exit("pdftoppm not found. Install Poppler (apt install poppler-utils).")
    if requested != "auto":
        return requested

    if have_pymupdf:
        return "pymupdf"
    if have_pdftoppm:
        return "pdftoppm"
    sys.exit(
        "No PDF renderer available. Install one:\n"
        "  pip install pymupdf\n"
        "  or: apt install poppler-utils / brew install poppler"
    )


def export(pptx: Path, outdir: Path, dpi: int, renderer: str,
           keep_pdf: bool = False) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        pdf = pptx_to_pdf(pptx, work)
        if keep_pdf:
            shutil.copy2(pdf, outdir / pdf.name)
        if renderer == "pymupdf":
            return render_pymupdf(pdf, outdir, pptx.stem, dpi)
        return render_pdftoppm(pdf, outdir, pptx.stem, dpi)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pptx", nargs="+", type=Path, help="one or more .pptx/.ppt files")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("slides"),
                    help="output directory (default: ./slides)")
    ap.add_argument("--dpi", type=int, default=600, help="resolution (default: 600)")
    ap.add_argument("--renderer", choices=["auto", "pymupdf", "pdftoppm"],
                    default="auto", help="PDF rasterizer (default: auto)")
    ap.add_argument("--keep-pdf", action="store_true",
                    help="also save the intermediate PDF")
    args = ap.parse_args()

    if not 1 <= args.dpi <= MAX_DPI:
        sys.exit(f"--dpi must be between 1 and {MAX_DPI}")

    renderer = pick_renderer(args.renderer)

    for deck in args.pptx:
        if not deck.exists():
            print(f"skipping {deck}: not found", file=sys.stderr)
            continue
        pages = export(deck, args.outdir, args.dpi, renderer, args.keep_pdf)
        print(f"{deck.name}: {len(pages)} slide(s) at {args.dpi} dpi via {renderer}")
        for p in pages:
            print(f"  {p}")


if __name__ == "__main__":
    main()
