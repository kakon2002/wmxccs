"""Attempt machine resolution of the 87 steroid names to InChIKeys. READ-ONLY.

Run: python tools/resolve_steroid_names.py

WHAT THIS IS AND IS NOT
-----------------------
It is a REPORT. It writes one file, a report, outside anything the model reads, and it
never touches data/seed. No assignment it produces enters the corpus, so running it cannot
move `corpus_sha256` and cannot change any answer the API gives. That is deliberate: a
structure assignment is a curation act with a provenance trail, and a name that resolved
today through a public service is not a provenance trail.

It exists to size the manual job, not to do it.

WHAT IT TRIES
-------------
Two public resolvers, in order, for each name AND for the commercial name beside it:

  PubChem PUG REST   https://pubchem.ncbi.nlm.nih.gov/rest/pug/...
  NIH CIR            https://cactus.nci.nih.gov/chemical/structure/...

A name is reported RESOLVED only where a resolver returns exactly one InChIKey. Anything
else - no answer, several answers, or two resolvers disagreeing - is reported as AMBIGUOUS
or FAILED with which.

THE AMBIGUITY FLAGS, and why they matter more than the successes
----------------------------------------------------------------
  stereo      the InChIKey's second block differs between resolvers or between the
              systematic and commercial name, meaning they agree on the skeleton and not
              on the stereochemistry. For steroids this is the difference between two real
              compounds with different cross sections.
  collision   two DIFFERENT names in this corpus resolved to the same InChIKey. One of
              them is wrong, and nothing here can say which.
  name-only   only the commercial name resolved, not the systematic one. The commercial
              name is the weaker identifier: "alpha-trenbolone/epitrenbolone" names two
              things, and several rows here carry a slash.

A resolution that is merely "found something" is not a resolution. The corpus holds epimer
pairs - 4,9,11-estratiene-17beta-ol-3-one and its 17alpha partner, same formula, same m/z,
one stereo descriptor apart - and a resolver that returns a stereo-free key for both has
merged two compounds while reporting success.
"""

from __future__ import annotations

import csv
import json
import pathlib
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

SEED = REPO / "data" / "seed" / "steroid_jasms2022.csv"
# Written OUTSIDE data/, so nothing here is covered by corpus_sha256 and nothing can be
# mistaken for an assignment the corpus holds.
REPORT = REPO / "reports" / "steroid_name_resolution.md"

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/InChIKey,IsomericSMILES/JSON"
CIR = "https://cactus.nci.nih.gov/chemical/structure/{}/stdinchikey"

TIMEOUT = 20
PAUSE = 0.25  # courtesy to two public services; neither is being hammered


def _get(url: str) -> str | None:
    context = ssl.create_default_context()
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    request = urllib.request.Request(url, headers={"User-Agent": "wmxccs-name-resolution/0.5.0"})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT, context=context) as response:
            return response.read().decode("utf-8", errors="replace")
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError):
        return None


def pubchem(name: str) -> list[str]:
    body = _get(PUBCHEM.format(urllib.parse.quote(name, safe="")))
    if not body:
        return []
    try:
        payload = json.loads(body)
    except json.JSONDecodeError:
        return []
    keys = [
        entry.get("InChIKey")
        for entry in payload.get("PropertyTable", {}).get("Properties", [])
        if entry.get("InChIKey")
    ]
    return sorted(set(keys))


def cir(name: str) -> list[str]:
    body = _get(CIR.format(urllib.parse.quote(name, safe="")))
    if not body:
        return []
    keys = {
        line.strip().replace("InChIKey=", "")
        for line in body.splitlines()
        if line.strip().startswith("InChIKey=")
    }
    return sorted(keys)


def skeleton(key: str) -> str:
    """The first block of an InChIKey: the connectivity, without the stereochemistry."""
    return key.split("-")[0]


def compounds() -> list[tuple[str, str]]:
    """(systematic name, commercial name) for every distinct compound, read from the corpus.

    Read-only. The file is opened and never written.
    """
    seen: dict[str, str] = {}
    with SEED.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = row["analyte_dataset_compound_id"].removeprefix("steroid_jasms2022:")
            seen.setdefault(name, row["analyte_display_name"])
    return sorted(seen.items())


def main() -> int:
    names = compounds()
    print(f"{len(names)} distinct compounds read from {SEED.relative_to(REPO)} (read-only)")

    rows = []
    for index, (systematic, commercial) in enumerate(names, 1):
        result = {
            "systematic": systematic,
            "commercial": commercial,
            "pubchem_systematic": pubchem(systematic),
            "cir_systematic": cir(systematic),
            "pubchem_commercial": pubchem(commercial) if commercial else [],
            "cir_commercial": cir(commercial) if commercial else [],
        }
        rows.append(result)
        time.sleep(PAUSE)
        if index % 10 == 0:
            print(f"   {index}/{len(names)} queried")

    # --- classify -------------------------------------------------------------------------
    for row in rows:
        systematic_keys = sorted(set(row["pubchem_systematic"]) | set(row["cir_systematic"]))
        commercial_keys = sorted(set(row["pubchem_commercial"]) | set(row["cir_commercial"]))
        flags = []
        if len(systematic_keys) == 1:
            row["key"] = systematic_keys[0]
            row["source"] = "PubChem" if row["pubchem_systematic"] else "CIR"
            if row["pubchem_systematic"] and row["cir_systematic"]:
                row["source"] = "PubChem+CIR"
            row["status"] = "RESOLVED"
        elif len(systematic_keys) > 1:
            row["key"] = ", ".join(systematic_keys)
            row["source"] = "PubChem/CIR disagree"
            row["status"] = "AMBIGUOUS"
            if len({skeleton(k) for k in systematic_keys}) == 1:
                flags.append("stereo")
        elif len(commercial_keys) == 1:
            row["key"] = commercial_keys[0]
            row["source"] = "commercial name only"
            row["status"] = "AMBIGUOUS"
            flags.append("name-only")
        elif len(commercial_keys) > 1:
            row["key"] = ", ".join(commercial_keys)
            row["source"] = "commercial name, resolvers disagree"
            row["status"] = "AMBIGUOUS"
            flags.append("name-only")
            if len({skeleton(k) for k in commercial_keys}) == 1:
                flags.append("stereo")
        else:
            row["key"] = "FAILED"
            row["source"] = "-"
            row["status"] = "FAILED"
        if commercial and "/" in (row["commercial"] or ""):
            flags.append("commercial name names two things")
        row["flags"] = flags

    # collisions: two different corpus names landing on one key
    by_key = defaultdict(list)
    for row in rows:
        if row["status"] == "RESOLVED":
            by_key[row["key"]].append(row["systematic"])
    for key, holders in by_key.items():
        if len(holders) > 1:
            for row in rows:
                if row["status"] == "RESOLVED" and row["key"] == key:
                    row["flags"].append(f"collision with {len(holders) - 1} other name(s)")
                    row["status"] = "AMBIGUOUS"

    # skeleton collisions: same connectivity, different stereochemistry - the epimer case
    by_skeleton = defaultdict(set)
    for row in rows:
        if row["status"] in {"RESOLVED", "AMBIGUOUS"} and row["key"] != "FAILED":
            for key in row["key"].split(", "):
                by_skeleton[skeleton(key)].add(row["systematic"])
    for core, holders in by_skeleton.items():
        if len(holders) > 1:
            for row in rows:
                if row["systematic"] in holders and "stereo-pair" not in row["flags"]:
                    row["flags"].append(f"stereo-pair: {len(holders)} names share connectivity {core}")

    # A clean resolution that shares connectivity with another name here is still a clean
    # resolution - and it is the one place a wrong answer would look right, because the only
    # thing separating the two compounds is the stereo block nobody downstream inspects.
    # Counted as resolved, and marked so a reader cannot miss it.
    for row in rows:
        row["needs_confirmation"] = row["status"] == "RESOLVED" and any(
            flag.startswith("stereo-pair") for flag in row["flags"]
        )

    resolved = [r for r in rows if r["status"] == "RESOLVED"]
    ambiguous = [r for r in rows if r["status"] == "AMBIGUOUS"]
    failed = [r for r in rows if r["status"] == "FAILED"]
    confirm = [r for r in resolved if r["needs_confirmation"]]

    # --- report ---------------------------------------------------------------------------
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with REPORT.open("w", encoding="utf-8") as out:
        out.write("# Steroid name resolution: report only, nothing assigned\n\n")
        out.write(
            "Produced by `tools/resolve_steroid_names.py` against the 87 distinct compound names in\n"
            "`data/seed/steroid_jasms2022.csv`, read-only. **No assignment here has entered the corpus.**\n"
            "Nothing in this file is covered by `corpus_sha256`, and running the script cannot change\n"
            "any answer the API gives.\n\n"
            "A structure assignment is a curation act needing a provenance trail. A name that resolved\n"
            "through a public service on one day is not one.\n\n"
        )
        out.write(
            f"**{len(rows)} names, {len(resolved)} resolved clean, {len(ambiguous)} ambiguous,"
            f" {len(failed)} failed - manual chemist review needed for"
            f" {len(ambiguous) + len(failed)}.**\n\n"
        )
        out.write(
            f"{len(confirm)} of the {len(resolved)} clean resolutions are additionally marked"
            " **NEEDS CONFIRMATION**: they resolved from their systematic name AND share connectivity"
            " with another name in this corpus. See below.\n\n"
        )
        out.write("| name | InChIKey or FAILED | source | ambiguity |\n")
        out.write("| --- | --- | --- | --- |\n")
        for row in rows:
            flags = "; ".join(row["flags"]) if row["flags"] else ""
            if row["needs_confirmation"]:
                flags = f"**NEEDS CONFIRMATION** - {flags}"
            out.write(
                f"| `{row['systematic']}` | {row['key']} | {row['source']} | {flags} |\n"
            )
        out.write("\n## What \"resolved clean\" means, and what NEEDS CONFIRMATION means\n\n")
        out.write(
            "A name is RESOLVED only where a resolver returned exactly one InChIKey FOR THE SYSTEMATIC\n"
            "NAME. Resolving from the commercial name instead is not a resolution and is reported as\n"
            "ambiguous, because a commercial name is the weaker identifier and several here name two\n"
            "compounds at once.\n\n"
            "**NEEDS CONFIRMATION** marks a clean resolution that shares connectivity with another name\n"
            "in this corpus, differing only in stereochemistry. In these the resolvers DID distinguish\n"
            "the pair, returning different stereo blocks, so the answer looks right - and that is the\n"
            "point. It is the one case where a WRONG answer would also look right, because the only thing\n"
            "separating the two compounds is a stereo block nothing downstream inspects, and the two have\n"
            "different cross sections. These do not need re-resolving. They need a chemist to confirm\n"
            "which stereoisomer is which, once, and record it.\n\n"
            "REPRODUCIBILITY. This report is a SNAPSHOT of two live public services. Unlike the model's\n"
            "own digests, re-running it on another day may give different answers if PubChem or CIR change\n"
            "what they return. That is a further reason no assignment here belongs in the corpus without a\n"
            "person putting it there.\n\n"
        )
        out.write("## Flags\n\n")
        counts = Counter(flag.split(":")[0].split(" with ")[0] for row in rows for flag in row["flags"])
        for flag, count in counts.most_common():
            out.write(f"- **{flag}** - {count}\n")

    print()
    print(f"report written to {REPORT.relative_to(REPO)}")
    print(
        f"{len(rows)} names, {len(resolved)} resolved clean, {len(ambiguous)} ambiguous,"
        f" {len(failed)} failed - manual chemist review needed for {len(ambiguous) + len(failed)}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
