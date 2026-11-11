"""
docs_hooks.py - MkDocs build hook: bring the notebooks and the explorer into the website.

What this file is FOR
---------------------
The notebooks live in notebooks/ and the explorer in explorer/, outside docs/
(the folder MkDocs builds from). Keeping them where they are is better for
users and tests, so at build time this hook COPIES them into docs/:

    notebooks/*.ipynb   -> docs/notebooks/      (rendered by mkdocs-jupyter, not re-executed)
    explorer/index.html -> docs/explorer/index.html   (served as-is: the live explorer)

Both copies are git-ignored and rebuilt every time, so the originals stay the
only source. Image links inside the notebooks that point at ../docs/figures/
are rewritten so they still work from the notebook's page on the site.

Registered in mkdocs.yml under `hooks:` (MkDocs >= 1.4).
Reference: https://www.mkdocs.org/user-guide/configuration/#hooks
"""
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent          # repository root (next to mkdocs.yml)
DOCS = ROOT / "docs"


def _copy_notebooks():
    """Copy every notebook into docs/notebooks/, fixing links to docs/figures/."""
    out = DOCS / "notebooks"
    out.mkdir(exist_ok=True)
    for src in sorted((ROOT / "notebooks").glob("*.ipynb")):
        nb = json.loads(src.read_text(encoding="utf-8"))
        for cell in nb["cells"]:
            if cell["cell_type"] != "markdown":
                continue
            text = "".join(cell["source"])
            # a notebook page is served at /notebooks/<name>/, so docs/figures/ is two levels up
            text = text.replace("](../docs/figures/", "](../../figures/")
            cell["source"] = text
        dst = out / src.name
        new = json.dumps(nb, ensure_ascii=False, indent=1)
        # only write when changed, so `mkdocs serve` does not rebuild in a loop
        if not dst.exists() or dst.read_text(encoding="utf-8") != new:
            dst.write_text(new, encoding="utf-8")


def _copy_explorer():
    """Copy the built explorer page so the site serves it at /explorer/."""
    src = ROOT / "explorer" / "index.html"
    if src.exists():
        dst = DOCS / "explorer" / "index.html"
        dst.parent.mkdir(exist_ok=True)
        if not dst.exists() or dst.read_bytes() != src.read_bytes():
            shutil.copyfile(src, dst)


def on_pre_build(config, **kwargs):
    """MkDocs calls this once before every build (and before every reload of `mkdocs serve`)."""
    _copy_notebooks()
    _copy_explorer()
