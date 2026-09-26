"""The glycoenzyme table: which residue-and-linkage combinations a human enzyme makes."""

import pytest

from wmxglycan.enzymes import EnzymeCatalogue

HEADER = "glycoenzyme,glycoclass,monolink,acceptor,glycan_class,species\n"
ROWS = (
    "MGAT1,MGAT,GlcNAc(b1-2),,N,Homo_sapiens\n"
    "MGAT2,MGAT,GlcNAc(b1-2),,N,Homo_sapiens\n"
    "FUT8,FUT,Fuc(a1-6),,N,Homo_sapiens\n"
    "B4GALT1,B4GALT,Gal(b1-4),,,Homo_sapiens\n"
)


def test_a_monolink_carries_every_enzyme_that_makes_it():
    catalogue = EnzymeCatalogue.from_csv(HEADER + ROWS)
    assert catalogue.makes("GlcNAc", "b1-2")
    assert catalogue.enzymes_for("GlcNAc", "b1-2") == ("MGAT1", "MGAT2")
    assert catalogue.enzymes_for("Fuc", "a1-6") == ("FUT8",)
    assert len(catalogue) == 3  # three distinct monolinks from four rows
    assert catalogue.enzyme_count == 4


def test_a_linkage_no_enzyme_makes_is_not_claimed():
    catalogue = EnzymeCatalogue.from_csv(HEADER + ROWS)
    assert not catalogue.makes("Gal", "a1-3")  # humans have no alpha1-3 galactosyltransferase
    assert catalogue.enzymes_for("Gal", "a1-3") == ()


def test_monolinks_for_a_residue():
    catalogue = EnzymeCatalogue.from_csv(HEADER + ROWS)
    assert [m.text for m in catalogue.monolinks_for("GlcNAc")] == ["GlcNAc(b1-2)"]
    assert catalogue.residues == ("Fuc", "Gal", "GlcNAc")


def test_a_leading_byte_order_mark_is_tolerated():
    catalogue = EnzymeCatalogue.from_csv("﻿" + HEADER + ROWS)
    assert catalogue.makes("Fuc", "a1-6")


# --- rows that state no attachment are counted, not dropped quietly ----------


def test_every_row_lands_in_exactly_one_bucket():
    catalogue = EnzymeCatalogue.from_csv(
        HEADER + ROWS + "AGA,AGA,,,,Homo_sapiens\nCHST1,CHST,6S,,,Homo_sapiens\n"
    )
    assert catalogue.rows_with_monolink == 4
    assert catalogue.rows_without_monolink == 1
    assert catalogue.rows_with_modification == 1
    assert catalogue.rows_read == 6


def test_an_enzyme_with_no_monolink_is_counted():
    # Glycosidases and chaperones are in the table with nothing to attach.
    catalogue = EnzymeCatalogue.from_csv(HEADER + ROWS + "AGA,AGA,,,,Homo_sapiens\n")
    assert catalogue.rows_without_monolink == 1
    assert len(catalogue) == 3  # the vocabulary is unchanged
    assert "no monolink" in catalogue.summary()


def test_modification_rows_are_counted_apart_from_the_marks_they_spell():
    # "6S" is a sulfation mark, not a residue arriving on a linkage. Three rows, two marks.
    catalogue = EnzymeCatalogue.from_csv(
        HEADER + ROWS + "CHST1,CHST,6S,,,Homo_sapiens\nCHST2,CHST,6S,,,Homo_sapiens\nCHST3,CHST,4S,,,Homo_sapiens\n"
    )
    assert catalogue.rows_with_modification == 3
    assert catalogue.modifications == ("4S", "6S")
    assert len(catalogue) == 3
    assert "modification" in catalogue.summary()


def test_a_row_naming_no_enzyme_stops_the_load():
    with pytest.raises(ValueError, match="names no enzyme"):
        EnzymeCatalogue.from_csv(HEADER + ",MGAT,GlcNAc(b1-2),,N,Homo_sapiens\n")


def test_missing_columns_are_an_error():
    with pytest.raises(ValueError, match="missing the columns"):
        EnzymeCatalogue.from_csv("glycoenzyme,monolink\nMGAT1,GlcNAc(b1-2)\n")


def test_a_header_with_no_rows_is_empty_not_malformed():
    # The columns are all present; blaming them would send the reader after the wrong fault.
    with pytest.raises(ValueError, match="no glycoenzymes"):
        EnzymeCatalogue.from_csv(HEADER)


def test_a_table_with_no_usable_monolink_is_an_error():
    with pytest.raises(ValueError, match="no glycoenzymes"):
        EnzymeCatalogue.from_csv(HEADER + "AGA,AGA,,,,Homo_sapiens\n")


# --- the table that actually ships -------------------------------------------


def test_the_shipped_table():
    # Tied to the pinned glycowork release.
    catalogue = EnzymeCatalogue.from_glycowork()
    assert catalogue.rows_read == 355
    assert catalogue.rows_with_monolink == 161
    assert catalogue.rows_without_monolink == 161
    assert catalogue.rows_with_modification == 33
    assert catalogue.modifications == ("3S", "4S", "6S", "GlcA2S", "GlcNAc3S", "OS")
    assert len(catalogue) == 35
    assert catalogue.enzyme_count == 150
    for residue, linkage, enzyme in [
        ("Fuc", "a1-6", "FUT8"),
        ("GlcNAc", "b1-2", "MGAT1"),
        ("GlcNAc", "b1-4", "MGAT3"),
        ("GlcNAc", "b1-6", "MGAT5"),
        ("Gal", "b1-4", "B4GALT1"),
        ("Neu5Ac", "a2-6", "ST6GAL1"),
    ]:
        assert catalogue.makes(residue, linkage)
        assert enzyme in catalogue.enzymes_for(residue, linkage)


def test_the_shipped_table_has_no_human_alpha1_3_galactosyltransferase():
    # Humans carry no GGTA1, which is why anti-Gal antibodies exist; the table reflects that.
    assert not EnzymeCatalogue.from_glycowork().makes("Gal", "a1-3")


def test_a_monolink_list_mixes_contexts_and_must_not_be_read_as_bond_attribution():
    # GlcNAc(b1-4) is made by MGAT3 (bisecting), MGAT4A-D (the alpha1-3 arm branch),
    # the hyaluronan synthases and a chitinase that only cleaves. The table records
    # no acceptor, so the list belongs to the monolink, never to one bond.
    enzymes = EnzymeCatalogue.from_glycowork().enzymes_for("GlcNAc", "b1-4")
    assert {"MGAT3", "MGAT4A", "CHIA"} <= set(enzymes)
