"""
Explorer Unit 2: the browser decoders (explorer/decoders.js) must give the same
answers as the Python decoders in qdec, so the pictures in the explorer are the
library's real behaviour.

How: we export the codes exactly as the explorer sees them (code_to_dict),
run decoders.js in Node on a batch of random errors, and compare with
  * MWPMDecoder on MatchingGraph.from_check_matrix  - same matched pairs, same weight
  * MinSumBP                                         - same beliefs at every iteration
  * BPOSDDecoder(osd_order=0)                        - same correction
Skipped when Node is not installed.
"""
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from qdec.families import catalogue
from qdec.export import code_to_dict, build_atlas_data
from qdec.dem import MatchingGraph, DecodingProblem
from qdec.decoders import MWPMDecoder, BPOSDDecoder
from qdec.bp import MinSumBP

ROOT = Path(__file__).resolve().parents[1]
JS = ROOT / "explorer" / "decoders.js"
NODE = shutil.which("node") or ("/opt/node22/bin/node" if Path("/opt/node22/bin/node").exists() else None)
pytestmark = pytest.mark.skipif(NODE is None, reason="Node.js not installed")

P, ITERS, ALPHA = 0.05, 30, 0.625

# Node driver: reads {code, cases:[{err:[...], kind}], what} on stdin, prints results as JSON
DRIVER = r"""
const Q = require(process.argv[process.argv.length - 1]);   // with -e, the path is the last argument
let buf = ""; process.stdin.on("data", (d) => (buf += d)).on("end", () => {
  const job = JSON.parse(buf), out = [];
  for (const c of job.cases) {
    const e = new Set(c.err);
    if (job.what === "mwpm") {
      const r = Q.mwpm(job.code, e, c.kind);
      out.push({ cost: r.cost, pairs: r.pairs.map((p) => [p.a, p.b]),
                 correction: [...r.correction].sort((a, b) => a - b) });
    } else if (job.what === "bp") {
      const r = Q.bp(job.code, e, c.kind, job.p, job.iters, job.alpha);
      out.push({ llr: r.llr, converged: r.converged, iterations: r.iterations });
    } else {
      const r = Q.bposd(job.code, e, c.kind, job.p, job.iters, job.alpha);
      out.push({ correction: [...r.correction].sort((a, b) => a - b), osd: !!r.osd });
    }
  }
  process.stdout.write(JSON.stringify(out));
});
"""


def run_js(code_dict, cases, what):
    """Run decoders.js in Node on a batch of cases; returns the list of results."""
    job = {"code": code_dict, "cases": cases, "what": what, "p": P, "iters": ITERS, "alpha": ALPHA}
    r = subprocess.run([NODE, "-e", DRIVER, str(JS)], input=json.dumps(job),
                       capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


CODES = {c.name: c for c in catalogue()}


def H_of(code, kind):
    """Check matrix of the generators of one kind, in generator order (as the JS sees them)."""
    return code.HX if kind == "X" else code.HZ


def random_cases(code, kind, n_cases, max_w, seed):
    rng = np.random.default_rng(seed)
    cases = []
    for _ in range(n_cases):
        w = int(rng.integers(1, max_w + 1))
        cases.append({"err": sorted(int(q) for q in rng.choice(code.n, w, replace=False)), "kind": kind})
    return cases


# ---------------------------------------------------------------- MWPM
@pytest.mark.parametrize("name", ["surface (rotated) d=5", "surface (unrotated) d=5", "toric 4x4", "repetition d=5"])
@pytest.mark.parametrize("kind", ["X", "Z"])
def test_mwpm_matches_python(name, kind):
    """Same matched pairs and same total weight as qdec.MWPMDecoder, on 40 random errors."""
    code = CODES[name]
    H = H_of(code, kind)
    if H.shape[0] == 0 or H.sum(0).max() > 2:
        pytest.skip("no graph for this check type")
    rows = list(code.index_of(kind))                      # JS reports generator indices
    cases = random_cases(code, kind, 40, 5, seed=1)
    js = run_js(code_to_dict(code, dist={"d": None, "how": ""}), cases, "mwpm")
    g = MatchingGraph.from_check_matrix(H, np.zeros((1, code.n), np.uint8), P)
    dec = MWPMDecoder(g)
    for case, out in zip(cases, js):
        e = np.zeros(code.n, np.uint8); e[case["err"]] = 1
        s = H.astype(int) @ e % 2
        dec.decode(s)
        py_pairs = sorted((rows[a], None if b == g.boundary else rows[b]) for a, b in dec.last_matching)
        js_pairs = sorted((a, b) for a, b in out["pairs"])
        assert js_pairs == py_pairs, case
        # total log-weight of the matched paths, from the Python graph's edge weights
        py_cost = sum(g.edges[(min(a, b), max(a, b))][2]
                      for path in dec.last_paths.values() for a, b in zip(path[:-1], path[1:]))
        assert out["cost"] == pytest.approx(py_cost, abs=1e-9), case
        # and the JS correction really clears the syndrome
        c = np.zeros(code.n, np.uint8); c[out["correction"]] = 1
        assert not (H.astype(int) @ ((e + c) % 2) % 2).any()


def test_mwpm_refuses_hypergraph_codes():
    """Steane's checks share qubits three ways: matching must say 'not possible' (null)."""
    code = CODES["Steane [[7,1,3]]"]
    r = subprocess.run([NODE, "-e", f"const Q=require({json.dumps(str(JS))});"
                        "let b='';process.stdin.on('data',d=>b+=d).on('end',()=>{const c=JSON.parse(b);"
                        "console.log(JSON.stringify([Q.isGraphLike(c,'Z'),Q.mwpm(c,new Set([0]),'Z')]))})"],
                       input=json.dumps(code_to_dict(code, dist={"d": None, "how": ""})),
                       capture_output=True, text=True, check=True)
    assert json.loads(r.stdout) == [False, None]


# ---------------------------------------------------------------- BP
@pytest.mark.parametrize("name", ["surface (rotated) d=5", "toric 4x4", "BB [[72,12,6]]", "gross [[144,12,12]]"])
def test_bp_beliefs_match_python(name):
    """Same LLR at every iteration as qdec.MinSumBP (X errors seen by Z checks)."""
    code = CODES[name]
    cases = random_cases(code, "Z", 15, 6, seed=2)
    js = run_js(code_to_dict(code, dist={"d": None, "how": ""}), cases, "bp")
    bp = MinSumBP(code.HZ, P, max_iter=ITERS, alpha=ALPHA)
    for case, out in zip(cases, js):
        e = np.zeros(code.n, np.uint8); e[case["err"]] = 1
        res = bp.decode(code.HZ.astype(int) @ e % 2)
        assert out["converged"] == bool(res.converged)
        assert out["iterations"] == res.iterations
        np.testing.assert_allclose(np.array(out["llr"]), res.llr_history, atol=1e-9)


# ---------------------------------------------------------------- BP + OSD-0
@pytest.mark.parametrize("name", ["toric 4x4", "BB [[72,12,6]]", "gross [[144,12,12]]", "HGP(Hamming, Hamming)"])
def test_bposd_matches_python(name):
    """Same correction as qdec.BPOSDDecoder(osd_order=0), including the shots where OSD runs."""
    code = CODES[name]
    cases = random_cases(code, "Z", 30, 8, seed=3)
    if name == "gross [[144,12,12]]":
        cases.append({"err": [18, 79, 94, 117], "kind": "Z"})   # the explorer's oscillation preset: OSD must run
    js = run_js(code_to_dict(code, dist={"d": None, "how": ""}), cases, "bposd")
    dec = BPOSDDecoder(DecodingProblem.code_capacity(code.HZ, np.zeros((1, code.n), np.uint8), P),
                       max_iter=ITERS, osd_order=0, alpha=ALPHA)
    osd_seen = 0
    for case, out in zip(cases, js):
        e = np.zeros(code.n, np.uint8); e[case["err"]] = 1
        corr = dec.decode_error(code.HZ.astype(int) @ e % 2)
        assert out["correction"] == [int(q) for q in np.nonzero(corr)[0]], case
        assert out["osd"] == dec.last["osd"] or out["osd"] == (not dec.last["bp"].converged)
        osd_seen += out["osd"]
    if name == "gross [[144,12,12]]":
        assert osd_seen > 0, "the gross test should exercise OSD at least once"


# ---------------------------------------------------------------- the explorer build
def test_template_inlines_decoders(tmp_path):
    """The built page carries decoders.js inline (single self-contained file)."""
    from qdec.export import build_explorer_html
    out = tmp_path / "x.html"
    build_explorer_html(ROOT / "explorer" / "template.html", out, use_stim=False)
    html = out.read_text(encoding="utf-8")
    assert "function mwpm(code, errSet, kind, p = 0.05)" in html
    assert "/*__DECODERS_JS__*/" not in html


# ---------------------------------------------------------------- the Decode tab in a real browser
EXPECTED_PRESETS = {                       # label prefix -> (decoder used, final verdict)
    "Surface d=5 · three X errors": ("mwpm", "exact"),
    "Surface d=5 · a different but equivalent": ("mwpm", "stabilizer"),      # degeneracy
    "Surface d=5 · three on a line": ("mwpm", "logical"),
    "Toric 4×4 · a tie that goes the wrong way": ("mwpm", "logical"),
    "Toric 4×4 · BP freezes": ("bposd", "exact"),
    "Steane · one Y error": ("lookup", "exact"),
    "Five-qubit · non-CSS": ("lookup", "exact"),
    "Gross · BP oscillates": ("bposd", "stabilizer"),
}


def test_decode_tab_in_browser(tmp_path):
    """Every example gives the expected verdict; every code steps through without a JS error."""
    pw = pytest.importorskip("playwright.sync_api")
    from qdec.export import build_explorer_html
    out = tmp_path / "x.html"
    build_explorer_html(ROOT / "explorer" / "template.html", out, use_stim=False)
    with pw.sync_playwright() as p:
        try:
            browser = p.chromium.launch()
        except Exception as e:                                   # no browser installed
            pytest.skip(f"no Chromium: {e}")
        page = browser.new_page()
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.goto(out.as_uri())
        page.click("#tab-decode")
        labels = page.evaluate("PRESETS.map(p => p.label)")
        assert len(labels) == len(EXPECTED_PRESETS)
        for i, label in enumerate(labels):
            page.click(f'button[data-preset="{i}"]')
            dec, verdict, osd = page.evaluate("""() => { const r = state.dec.run;
                return [r.dec, r.steps[r.steps.length - 1].verdict, r.steps.some(s => s.title.includes('OSD'))]; }""")
            want = next(v for k, v in EXPECTED_PRESETS.items() if label.startswith(k))
            assert (dec, verdict) == want, label
            if "OSD" in label:
                assert osd, "the OSD examples must show the OSD step"
        # every code, a Y and an X error, every step drawn (with the BP slider moved)
        for ci in range(page.evaluate("ATLAS.codes.length")):
            page.evaluate(f"selectCode({ci}); state.err = {{1: 'Y', 3: 'X'}}; render(); syncPlayer();")
            for s in range(page.evaluate("state.dec.run.steps.length")):
                page.evaluate(f"state.dec.step = {s}; state.iter = 3; render(); syncPlayer();")
            assert page.evaluate("state.dec.run.steps[state.dec.run.steps.length - 1].title") == "Result"
        browser.close()
    assert not errors, errors
