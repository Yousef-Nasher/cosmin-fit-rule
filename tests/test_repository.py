"""Checks on the repository as a whole. Run from the repository root:  python -m pytest tests

They cover the integrity of the stored data, the consistency of the Study 3 tables and draw,
the agreement between the R results and the independent Python results, and that every Python
script reproduces its stored output."""
import csv
import hashlib
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def read_csv(path):
    with open(ROOT / path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def text(path):
    return (ROOT / path).read_text(encoding="utf-8")


def run_script(folder, script, *args):
    done = subprocess.run([sys.executable, script, *args], cwd=ROOT / folder, capture_output=True,
                          text=True, encoding="utf-8", check=True)
    return done.stdout


# ---- stored data ---------------------------------------------------------------------
def test_checksums():
    lines = text("SHA256SUMS").strip().splitlines()
    assert len(lines) >= 5
    for line in lines:
        digest, name = line.split(maxsplit=1)
        assert hashlib.sha256((ROOT / name.strip()).read_bytes()).hexdigest() == digest, name


def test_simulation_archives():
    with zipfile.ZipFile(ROOT / "study2/data/study2_final_chunks.zip") as z:
        names = z.namelist()
        assert sum(n.startswith("final/") and n.endswith(".csv") for n in names) == 152 * 10
        assert "population.csv" in names
    with zipfile.ZipFile(ROOT / "study2b/data/study2b_chunks.zip") as z:
        names = z.namelist()
        assert sum(n.startswith("chunks/") and n.endswith(".csv") for n in names) == 288 * 4
        assert "population2b.csv" in names


# ---- Study 3 tables --------------------------------------------------------------------
def test_study3_tables_are_consistent():
    frame = read_csv("study3/data/study3_frame.csv")
    outside = read_csv("study3/data/study3_frame_outside_years.csv")
    search = {r["stage"] + "|" + r["source"]: int(r["n"]) for r in read_csv("study3/data/study3_search_counts.csv")}
    screen = read_csv("study3/data/study3_fulltext_screening.csv")
    reviews = read_csv("study3/data/study3_reviews.csv")
    rows = read_csv("study3/data/study3_cfa_results.csv")
    instruments = read_csv("study3/data/study3_instruments.csv")
    ids = {r["frame_id"] for r in frame}
    assert len(frame) == len(ids) == 1120
    assert not ids & {r["frame_id"] for r in outside} and len(outside) == 57
    found = search["records identified|PubMed search"] + search["records identified|COSMIN database export"]
    screened = search["records screened on title and abstract|both sources"]
    assert found - search["duplicates removed|both sources"] == screened
    excluded = sum(n for k, n in search.items() if k.startswith("excluded:") and k.endswith("|title and abstract"))
    assert screened - excluded == search["records passed to the sampling frame|title and abstract"]
    assert search["records passed to the sampling frame|title and abstract"] - len(outside) == len(frame)
    assert {r["frame_id"] for r in screen} == ids and len(screen) == len(frame)
    assert sum(r["draw_list"] == "drawn" for r in screen) == 100
    included = {r["frame_id"] for r in screen if r["fulltext_status"] == "included"}
    assert included == {r["review_id"] for r in reviews} and len(reviews) == 97
    assert all(r["reason"] for r in screen if r["fulltext_status"] == "excluded")
    assert {r["review_id"] for r in rows} <= included and {r["review_id"] for r in instruments} <= included
    assert len({r["row_id"] for r in rows}) == len(rows)
    assert all(r["wording_code"] in "12345" and r["cosmin_version"] in "12" for r in reviews)
    assert all(r["level"] in ("study", "pooled", "summary") for r in rows)
    assert all(r["rating"] in ("", "+", "-", "+/-", "?") for r in rows)
    assert [r["n_cfa_rows"] for r in reviews] == [
        str(sum(x["review_id"] == r["review_id"] for x in rows)) for r in reviews]
    for r in rows:   # every single printed value is stored as the same number
        for k in ["cfi", "tli", "rmsea", "srmr"]:
            if re.fullmatch(r"-?\d*\.?\d+", r[k + "_printed"]):
                assert float(r[k]) == float(r[k + "_printed"]), (r["row_id"], k)


def test_study3_screening_follows_the_draw():
    screen = read_csv("study3/data/study3_fulltext_screening.csv")
    drawn = run_script("study3", "code/draw_sample.py").splitlines()
    assert drawn[0] == "frame_id,year_stratum,draw_list,rank"
    stored = {r["frame_id"]: (r["year_stratum"], r["draw_list"], r["rank"]) for r in screen}
    assert {line.split(",")[0]: tuple(line.split(",")[1:]) for line in drawn[1:]} == stored
    # within a stratum, records are reached in queue order: nothing after the first record not reached
    for stratum in {r["year_stratum"] for r in screen}:
        queue = sorted((int(r["rank"]), r["fulltext_status"]) for r in screen
                       if r["year_stratum"] == stratum and r["draw_list"] == "replacement queue")
        status = [s for _, s in queue]
        if "not reached" in status:
            assert all(s == "not reached" for s in status[status.index("not reached"):])


# ---- R results against the independent Python results ---------------------------------
def test_study2_r_and_python_agree():
    r_out = text("study2/results/study2_results.txt")
    py = text("study2/results/verify_independent.txt")
    primary = r_out.split("SENSITIVITY")[0]
    for arm in ["ml_cont", "ml_ord_as_cont", "wlsmv_ord"]:
        block = primary.split(f"--- arm: {arm} ---")[1].split("--- arm:")[0]
        for n in [100, 200, 500, 1000]:
            m = re.search(rf"^\s*{n} \|\s*([\d.]+)\([\d.]+\) /\s*([\d.]+)\(", block, re.M)
            ch = re.search(rf"N=\s*{n}\s.*dropping RMSEA changes (\d+) of", block)
            p = re.search(rf"^{arm}\s+N=\s*{n}\s+false-fail\s+([\d.]+)\s+false-pass\s+([\d.]+)\s+RMSEA-changes (\d+)", py, re.M)
            assert (m.group(1), m.group(2), ch.group(1)) == p.groups(), (arm, n)


def test_study2b_r_and_python_agree():
    r_out = text("study2b/results/study2b_results.txt")
    py = text("study2b/results/verify_independent_2b.txt")
    e1 = re.search(r"E1 SRMR false-fail > 20%: holds in (\d+) of (\d+)", r_out).groups()
    e2 = re.search(r"E2 RMSEA rescue > 5%: holds in (\d+) of (\d+)", r_out).groups()
    e3 = re.search(r"E3 RMSEA rescue < 1% for 5-cat symmetric: holds in (\d+) of (\d+)", r_out).groups()
    assert re.search(r"E1 (\d+) of (\d+)", py).groups() == e1
    assert re.search(r"E2 (\d+) of (\d+)", py).groups() == e2
    assert re.search(r"E3 (\d+) of (\d+)", py).groups() == e3


def test_study3_r_and_python_agree():
    r_out = text("study3/results/study3_results.txt")
    py = text("study3/results/verify_independent.txt")
    q1 = re.search(r"Q1 wording: (\{.*\})", py).group(1)
    for label in ["disjunctive", "conjunctive", "mixed", "availability", "not stated"]:
        n = re.search(rf"^  {label}\s+(\d+)", r_out, re.M).group(1)
        assert f"'{label}': {n}" in q1, label
    for rule in ["cosmin", "no_rmsea", "srmr_df", "two_of_three"]:
        r = re.search(rf"^  {rule}\s+sufficient (\d+)\s+insufficient (\d+)\s+indeterminate\s+(\d+)", r_out, re.M).groups()
        p = re.search(rf"^Q4 {rule}\s+\{{'sufficient': (\d+), 'insufficient': (\d+), 'indeterminate': (\d+)\}}", py, re.M).groups()
        assert r == p, rule
    r = re.search(r"At least one reported index fails its own threshold: (\d+)/(\d+)", r_out).groups()
    assert re.search(r"^Q3 (\d+)/(\d+)", py, re.M).groups() == r
    r = (re.search(r"Rated \+ although no reported index meets its threshold: (\d+)", r_out).group(1),
         re.search(r"Rated - although a reported index meets its threshold: (\d+)", r_out).group(1))
    assert re.search(r"no passing index: (\d+) \| printed - with a passing index: (\d+)", py).groups() == r
    for rule in ["no_rmsea", "srmr_df", "two_of_three"]:
        r = re.search(rf"^  {rule}\s+A -> B: (\d+); B -> C: (\d+); changed (\d+)/(\d+)", r_out, re.M).groups()
        p = re.search(rf"^{rule}\s+.*A -> B (\d+); B -> C (\d+); recommendations changed (\d+)/(\d+)", py, re.M).groups()
        assert r == p, rule


# ---- every Python script reproduces its stored output -----------------------------------
@pytest.mark.parametrize("folder, script, stored", [
    ("study2", "code/verify_independent.py", "results/verify_independent.txt"),
    ("study2b", "code/verify_independent_2b.py", "results/verify_independent_2b.txt"),
    ("study3", "code/verify_independent.py", "results/verify_independent.txt"),
])
def test_python_scripts_reproduce(folder, script, stored):
    assert run_script(folder, script).splitlines() == text(f"{folder}/{stored}").splitlines()
