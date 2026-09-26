"""The Bush Lab conversion, and the five structural traps it exists to survive.

Split the way test_ingest_steroid.py is split, and for the same reason: data/raw is
gitignored, so the committed conversion in data/seed/as_delivered IS the reference these
tests run against. Any clone can run them.

- THE PURE GUARDS, on synthetic input. Each of the five traps has a refusal behind it, and
  each refusal is exercised in both directions: once where it must fire and once where it
  must not. A converter that refused everything would pass half of these.
- WHAT THE COMMITTED CONVERSION HOLDS. Counts, per-sheet gases, the unit arithmetic, and
  the blockers. If the workbook is ever re-read differently, one of these goes red with a
  number a reviewer can compare against the file.
- THAT NO SEED FILE EXISTS. This is the assertion that matters most here. The workbook
  names no platform, a seed row needs one, and the whole point of this adapter is that it
  refuses rather than inventing one. A seed file appearing without a platform map would mean
  somebody had guessed.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import pytest

from wmxccs.identity import adduct_carrier_is_unstated, parse_adduct

import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools"))
import ingest_bushlab as adapter  # noqa: E402

CONVERSION = REPO / "data" / "seed" / "as_delivered" / "bushlab_ccs_database_converted.csv"
SEED = REPO / "data" / "seed" / "bushlab_microsource.csv"
REFERENCE_PAPER = "10.1021/acs.analchem.7b01709"

EXPECTED_DATA_ROWS = 1804
EXPECTED_CCS_VALUES = 2045


@pytest.fixture(scope="module")
def rows() -> list[dict]:
    assert CONVERSION.is_file(), (
        f"{CONVERSION} is missing. It is committed precisely because data/raw is not, so"
        " these tests do not need the workbook."
    )
    with CONVERSION.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


# --- the pure guards ----------------------------------------------------------------------


class FakeSheet:
    """Enough of an openpyxl sheet for `rows_of`. Title included: refusals name it."""

    def __init__(self, rows, title="fake"):
        self._rows = rows
        self.title = title

    def iter_rows(self, values_only=True):
        return iter(self._rows)


def test_a_sheet_whose_header_depth_changed_is_refused_not_reinterpreted():
    """Trap 4. Three sheets put a banner above their column names.

    If the file gains or loses a header row, every column index in LAYOUT shifts and the
    converter would read a charge column as a cross section. Refusing is the only safe
    answer, and the message says not to adjust the number to make it pass.
    """
    two_row_header = [("", "z = 1", ""), ("n", "A", "B"), (3.0, 89.0, 151.0)]
    assert adapter.rows_of(FakeSheet(two_row_header), 2)[1] == [(3.0, 89.0, 151.0)]
    with pytest.raises(adapter.SheetRefused, match="declared 1 header"):
        adapter.rows_of(FakeSheet(two_row_header), 1)


def test_a_unit_is_read_from_the_header_and_a_wrong_declaration_is_refused():
    """Trap 1, and the direction that matters: the unit is never inferred from magnitude."""
    header = ["Protein", "n", "m / kDa", "z", "Ω(He) / nm^2", "Ω(N2) / nm^2"]
    assert adapter.unit_from_header(header, 4, "nm^2") == "nm^2"
    with pytest.raises(adapter.SheetRefused, match="does not state the declared unit"):
        adapter.unit_from_header(header, 4, "A^2")


def test_an_angstrom_header_is_read_as_angstroms_however_it_is_spelled():
    """The other direction. A converter that refused every unit would pass the test above."""
    header = ["n", "Ω(He) / Å^2", "Ω(N2) / Å^2"]
    assert adapter.unit_from_header(header, 1, "A^2") == "A^2"
    assert adapter.unit_from_header(header, 2, "A^2") == "A^2"


def test_the_two_identically_labelled_helium_columns_are_refused_without_their_banner():
    """Trap 2, the worst of the five.

    Columns 4 and 5 are both headed as helium and are told apart only by a merged banner
    above them. Read as helium-against-nitrogen, an anion measurement silently becomes a
    nitrogen measurement of a cation - a wrong gas AND a wrong charge sign, on a number that
    still looks like a plausible cross section.
    """
    adapter.check_banner([[None, None, None, None, "Cations", "Anions"], ["Sample"]])
    with pytest.raises(adapter.SheetRefused, match="cannot be told apart"):
        adapter.check_banner([[None, None, None, None, None, None], ["Sample"]])
    with pytest.raises(adapter.SheetRefused, match="cannot be told apart"):
        adapter.check_banner([[None, None, None, None, "Cations", "Cations"], ["Sample"]])


def test_the_unstated_carrier_adduct_is_well_formed_and_keys_per_charge():
    """Trap 3. A charge with no adduct is an ion the source did not define."""
    for charge in (1, 2, 3, 24, -1, -2):
        adduct = adapter.unstated_carrier_adduct(charge)
        assert parse_adduct(adduct), adduct
        assert adduct_carrier_is_unstated(adduct), adduct
    # distinct charges give distinct adducts, so they cannot collapse into one key
    made = {adapter.unstated_carrier_adduct(z) for z in (1, 2, 3, 24)}
    assert len(made) == 4
    assert adapter.unstated_carrier_adduct(3) == "[M+3?]3+"
    assert adapter.unstated_carrier_adduct(-2) == "[M-2?]2-"


def test_a_named_adduct_is_not_reported_as_an_unstated_carrier():
    """The other direction: MicroSource states real adducts and they must read as stated."""
    assert not adduct_carrier_is_unstated("[M+H]+")
    assert not adduct_carrier_is_unstated("[M+Na]+")


def test_an_adduct_the_schema_refuses_is_flagged_and_never_repaired():
    """The CEO's instruction: handle ambiguous adducts explicitly, do not normalise them.

    Three MicroSource values carry strings this schema refuses. `[M+3H]+3` is almost
    certainly `[M+3H]3+` - and "almost certainly" is exactly why the adapter does not
    rewrite it. A repaired adduct is a fact invented in the column that decides which ions
    match each other.
    """
    assert adapter.adduct_problem("[M+H]+") == ""
    assert "bracket notation" in adapter.adduct_problem("M+")
    assert "bracket notation" in adapter.adduct_problem("[M+3H]+3")


# --- what the committed conversion holds ---------------------------------------------------


def test_the_conversion_holds_every_value_the_workbook_states(rows):
    assert len(rows) == EXPECTED_CCS_VALUES
    located = {r["source_locator"] for r in rows}
    assert len(located) == EXPECTED_CCS_VALUES, "two values share a locator, so one is untraceable"


def test_values_outnumber_rows_because_a_row_can_hold_several_cross_sections(rows):
    """1,804 data rows, 2,045 values. A row carries a helium AND a nitrogen column, and the
    polymer sheets carry one column pair per charge state."""
    assert EXPECTED_CCS_VALUES > EXPECTED_DATA_ROWS
    per_sheet = Counter(r["sheet"] for r in rows)
    assert per_sheet["MicroSource Collection"] == 1440, "one value per row on the only 1:1 sheet"
    assert per_sheet["Polyalaine Cations"] == 86, "31 rows, up to six values each"


def test_every_converted_value_is_its_native_value_times_its_stated_factor(rows):
    """The conversion is auditable rather than assumed: the arithmetic is in the file."""
    for r in rows:
        native = float(r["native_value"])
        factor = float(r["conversion_factor"])
        assert float(r["ccs_a2"]) == pytest.approx(native * factor, rel=0, abs=1e-9)


def test_the_factor_is_a_hundred_for_nanometres_and_one_for_angstroms(rows):
    """Trap 1's arithmetic, asserted per row rather than trusted once."""
    factors = {r["native_unit"]: {float(x["conversion_factor"])
                                  for x in rows if x["native_unit"] == r["native_unit"]}
               for r in rows}
    assert factors["nm^2"] == {100.0}
    assert factors["A^2"] == {1.0}
    assert adapter.A2_PER_NM2 == 100.0


def test_the_nanometre_sheets_are_the_protein_and_peptide_ones(rows):
    """A unit belongs to a sheet, and mixing them up is a two-order-of-magnitude error."""
    units = {}
    for r in rows:
        units.setdefault(r["sheet"], set()).add(r["native_unit"])
    assert units["Native-Like Protein Cations"] == {"nm^2"}
    assert units["Denatured Protein Cations"] == {"nm^2"}
    assert units["Other Peptides"] == {"nm^2"}
    assert units["MicroSource Collection"] == {"A^2"}
    assert units["Polyalaine Cations"] == {"A^2"}
    assert units["Small Molecular Ions"] == {"A^2"}


def test_the_banner_sheet_is_helium_only_and_the_banner_became_a_charge_sign(rows):
    """Trap 2, on the real data. Both its columns are helium; the banner set the polarity."""
    banner = [r for r in rows if r["sheet"] == "Native-Like Protein Cations and"]
    assert banner, "the banner sheet contributed nothing, so this test proves nothing"
    assert {r["drift_gas"] for r in banner} == {"He"}, "a column was read as nitrogen"
    assert {r["polarity"] for r in banner} == {"positive", "negative"}


def test_no_other_sheet_invented_a_negative_charge(rows):
    """The other direction: only the sheets whose banners say so carry anions."""
    negative = {r["sheet"] for r in rows if r["polarity"] == "negative"}
    assert negative == {"Native-Like Protein Cations and", "Anionic Homopolymers"}


def test_every_protein_peptide_and_polymer_value_carries_an_unstated_carrier(rows):
    """Trap 3 on the real data, and the count is the finding: 559 values."""
    # the three unusable adducts cannot be parsed at all, so they are excluded here and
    # asserted separately in test_the_other_two_blockers_are_counted_separately
    usable = [r for r in rows if r["adduct"] and not adapter.adduct_problem(r["adduct"])]
    unstated = [r for r in usable if adduct_carrier_is_unstated(r["adduct"])]
    assert len(unstated) == 559
    assert {r["sheet"] for r in unstated} <= set(adapter.NO_ADDUCT_SHEETS)
    stated = [r for r in usable if not adduct_carrier_is_unstated(r["adduct"])]
    assert {r["sheet"] for r in stated} == {"MicroSource Collection"}, (
        "MicroSource is the only sheet with an adduct column"
    )


def test_only_hines_2017_is_resolved_and_the_other_eight_papers_are_still_blocked(rows):
    """The outcome of this ingest, asserted rather than described.

    The workbook names no platform anywhere. One paper of the nine has been read - Hines 2017,
    which governs MicroSource - so its 1,440 values carry a platform and the other 605 do not.
    The eight are unread deliberately: they govern 559 unstated-carrier values that can never
    match whatever platform they turn out to be on. If this ever passes with a different
    split, somebody read another paper or guessed one, and either is a deliberate act.
    """
    resolved = [r for r in rows if r["cited_paper"] == "Hines 2017"]
    unresolved = [r for r in rows if r["cited_paper"] != "Hines 2017"]
    assert len(resolved) == 1440
    assert len(unresolved) == 605

    assert all(r["ims_type"] == "TWIMS" for r in resolved)
    assert all(r["drift_gas"] == "N2" for r in resolved)
    assert all(REFERENCE_PAPER in r["platform_evidence"] for r in resolved)

    assert all("no platform stated" in r["blocked_because"] for r in unresolved)
    assert all(r["ims_type"] == "" for r in unresolved)


def test_the_three_unusable_adducts_stay_blocked_even_with_a_platform(rows):
    """A resolved platform does not resolve an adduct the schema refuses."""
    still_blocked = [
        r for r in rows if r["cited_paper"] == "Hines 2017" and r["blocked_because"]
    ]
    assert len(still_blocked) == 3
    assert {r["adduct"] for r in still_blocked} == {"M+", "[M+3H]+3"}


def test_the_other_two_blockers_are_counted_separately(rows):
    """Small Molecular Ions states no charge; three MicroSource adducts are unusable."""
    no_charge = [r for r in rows if "no charge stated" in r["blocked_because"]]
    assert len(no_charge) == 46
    assert {r["sheet"] for r in no_charge} == {"Small Molecular Ions"}
    bad_adduct = [r for r in rows if "adduct not in a form" in r["blocked_because"]]
    assert len(bad_adduct) == 3
    assert {r["adduct"] for r in bad_adduct} == {"M+", "[M+3H]+3"}


def test_the_nine_cited_papers_are_the_bounded_task(rows):
    """Every value names the paper that would resolve its platform."""
    papers = {r["cited_paper"] for r in rows}
    assert "" not in papers, "a value with no cited paper has nothing to resolve it"
    assert len(papers) == 9


def seeded_rows() -> list[dict]:
    with SEED.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_the_seed_holds_only_the_paper_that_was_read():
    """A seed row needs an ims_type, and this one has evidence behind it.

    The previous version of this test asserted that NO seed file existed, because no platform
    was known. It went red when Hines 2017 was read, which is what it was for: a seed file
    appearing without a resolution is an invented platform, and it would enter the fit
    silently because build_default_model globs every CSV in data/seed.
    """
    seeded = seeded_rows()
    assert len(seeded) == 1437, "1,440 MicroSource values less the 3 unusable adducts"
    assert {r["ims_type"] for r in seeded} == {"TWIMS"}
    assert {r["drift_gas"] for r in seeded} == {"N2"}
    assert {r["source"] for r in seeded} == {"Bush Lab CCS database"}


def test_every_seeded_record_carries_the_calibration_lineage():
    """CONDITION 1, and the reason for the other three.

    These TWIMS values were calibrated against their own laboratory nitrogen DTIMS
    measurements. Leaving the lineage blank would lose the one fact that makes them dangerous
    to compare, and the comparison guard in test_circularity.py has nothing to read without
    it. A blank calref column here would disarm that guard silently.
    """
    seeded = seeded_rows()
    assert all(r["calref_reference_set"] for r in seeded)
    assert {r["calref_doi"] for r in seeded} == {REFERENCE_PAPER}
    assert {r["calref_platform"] for r in seeded} == {"DTIMS"}
    assert {r["calref_method"] for r in seeded} == {"stepped_field"}
    assert all("RF-confining drift cell" in r["calref_reference_set"] for r in seeded)


def test_the_publication_is_in_the_locator_and_not_in_the_doi_column():
    """Permission comes from the database page, not from the ACS paper.

    The paper licence has not been read - Europe PMC flags it open access and names no
    licence, and an open-access flag is not a licence - so a reuse claim attached to its DOI
    would be unbacked and the gate refuses it. This is the decision already taken for
    CCSbase. The publication still has to be findable, so it lives in the locator, and the
    circularity guard reads it from there.
    """
    seeded = seeded_rows()
    assert all(r["doi"] == "" for r in seeded)
    assert all(f"doi:{REFERENCE_PAPER}" in r["source_locator"] for r in seeded)


def test_the_twelve_protomer_pairs_are_conformers_and_none_is_dropped_or_averaged():
    """CONDITION 3. Two peaks of one ion are a legitimate pair, not a duplicate."""
    seeded = seeded_rows()
    conformers = [r for r in seeded if r["conformer"]]
    assert len(conformers) == 24, "twelve pairs, both members kept"
    assert {r["conformers_total"] for r in conformers} == {"2"}
    assert sorted(r["conformer"] for r in conformers) == ["1"] * 12 + ["2"] * 12
    by_compound: dict[str, list[dict]] = {}
    for r in conformers:
        by_compound.setdefault(r["analyte_display_name"], []).append(r)
    assert len(by_compound) == 12
    for name, pair in by_compound.items():
        assert len({r["adduct"] for r in pair}) == 1, f"{name}: a conformer pair shares its adduct"
        assert len({r["ccs"] for r in pair}) == 2, f"{name}: two conformers with one CCS"


def test_the_six_pairs_the_paper_does_not_name_carry_a_curation_flag():
    """The honest half of condition 3.

    Hines 2017 says the sixteen "display two peaks, had two major adducts, or were mixtures".
    Four differ by adduct and are ordinary separate ions. Of the twelve that share an adduct,
    the paper NAMES only the fluoroquinolone protomers and its cephalosporin finding. The
    other six are recorded as conformers and flagged, because a mixture two components are
    two analytes and calling them conformers of one ion would merge two compounds.
    """
    seeded = seeded_rows()
    flagged = [r for r in seeded if r["curation_flag"]]
    assert len(flagged) == 12, "six pairs, both members flagged"
    named = {
        r["analyte_display_name"].casefold()
        for r in seeded
        if r["conformer"] and not r["curation_flag"]
    }
    assert "ciprofloxacin" in named and "norfloxacin" in named
    assert "cefpodoxime proxetil" in named
    assert all("mixture" in r["curation_flag"] for r in flagged)


def test_the_four_pairs_differing_by_adduct_are_not_conformers():
    """The other direction: a different adduct is a different ion, not a second peak."""
    seeded = seeded_rows()
    by_name: dict[str, list[dict]] = {}
    for r in seeded:
        by_name.setdefault(r["analyte_display_name"].casefold(), []).append(r)
    two_adducts = {n: rs for n, rs in by_name.items() if len({r["adduct"] for r in rs}) > 1}
    assert len(two_adducts) == 4
    assert set(two_adducts) == {
        "glycocholic acid", "methyldopa", "methoxamine hydrochloride", "podofilox",
    }
    for name, rs in two_adducts.items():
        assert all(not r["conformer"] for r in rs), f"{name} was labelled a conformer"


def test_the_conversion_is_not_picked_up_by_the_model_loader():
    """It lives in as_delivered/ on purpose: build_default_model globs data/seed/*.csv only."""
    seeds = sorted(p.name for p in (REPO / "data" / "seed").glob("*.csv"))
    assert "bushlab_ccs_database_converted.csv" not in seeds
    assert "bushlab_microsource.csv" in seeds, "the seed itself IS meant to be loaded"
    assert CONVERSION.parent.name == "as_delivered"
