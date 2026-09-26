"""The store: frozen means frozen, appending never merges, and a restart changes nothing.

Every test here drives the real SQLite code path. An in-memory database is still SQLite, and the
one property a file gives that memory does not - surviving a restart - is tested against a real
file in a temporary directory, because that is the whole reason the owner asked for persistence.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

from wmxglycan.store import (
    ALREADY_FROZEN,
    NO_SUCH_PREDICTION,
    PAYLOAD_MOVED,
    AlreadyFrozen,
    Store,
    UnknownPrediction,
    digest_of,
    new_id,
    now,
)

BODY = {"candidates_total": 6, "classes_total": 6, "is_a_ranking": False}
OTHER = {"candidates_total": 999, "classes_total": 1, "is_a_ranking": True}


def freeze(store: Store, **overrides):
    fields = dict(
        prediction_id="pred_one",
        created_at=now(),
        composition="Hex5HexNAc2",
        adduct="[M-H]-",
        charge=-1,
        pipeline_version="0.1.0",
        fingerprint="params/corpus",
        snapshot_digest="snapshot",
        payload=BODY,
    )
    fields.update(overrides)
    return store.freeze(**fields)


@pytest.fixture
def store():
    made = Store()
    yield made
    made.close()


# --- a frozen prediction is immutable ------------------------------------------------------------


def test_a_second_write_to_a_frozen_prediction_is_refused(store):
    first = freeze(store)
    with pytest.raises(AlreadyFrozen) as caught:
        freeze(store, payload=OTHER)
    assert "already frozen and was not modified" in str(caught.value)
    assert "REFUSED rather than merged" in str(caught.value)
    # AND NOTHING WAS MERGED. The distinction the owner asked for: a refusal that left the row
    # half-updated would satisfy "refused" and break the thing refusing protects.
    assert store.get("pred_one").payload == first.payload
    assert store.get("pred_one").payload_digest == first.payload_digest
    assert json.loads(store.get("pred_one").payload)["candidates_total"] == 6


def test_a_second_write_under_the_same_client_reference_is_refused_and_names_what_stands(store):
    first = freeze(store, client_reference="run-A")
    with pytest.raises(AlreadyFrozen) as caught:
        freeze(store, prediction_id="pred_two", payload=OTHER, client_reference="run-A")
    # The message names the prediction that STANDS, not the one that was refused, so a caller
    # can fetch it instead of guessing what happened.
    assert first.prediction_id in str(caught.value)
    assert store.get("pred_two") is None
    assert store.by_client_reference("run-A").payload == first.payload


def test_the_refusal_wording_is_importable_rather_than_matched_as_prose():
    assert "{prediction_id!r}" in ALREADY_FROZEN
    assert "REFUSED rather than merged" in ALREADY_FROZEN


def test_the_store_has_no_update_statement_anywhere():
    # The design, asserted. A frozen payload is never rewritten, and the way that is guaranteed
    # is that no code path exists to rewrite it - not that callers are careful.
    source = Path(__import__("wmxglycan.store", fromlist=["x"]).__file__).read_text(encoding="utf-8")
    statements = source.upper()
    assert "UPDATE PREDICTIONS" not in statements
    assert "UPDATE VALIDATIONS" not in statements
    assert "DELETE FROM" not in statements
    # INSERT is the only write verb used against predictions.
    assert "INSERT INTO PREDICTIONS" in statements


def test_two_predictions_of_the_same_composition_are_both_kept(store):
    # Creating twice is not a mutation: two ids, two rows, both readable. What is refused is a
    # second write to ONE prediction, not a second prediction.
    one = freeze(store, prediction_id="pred_a")
    two = freeze(store, prediction_id="pred_b")
    assert one.payload_digest == two.payload_digest
    assert {row.prediction_id for row in store.predictions()} == {"pred_a", "pred_b"}


# --- attaching appends and never touches the prediction --------------------------------------------


def test_attaching_a_validation_leaves_the_prediction_byte_identical(store):
    frozen = freeze(store)
    store.attach(
        prediction_id="pred_one",
        validation_id="val_one",
        created_at=now(),
        fingerprint="params/corpus",
        snapshot_digest="snapshot",
        payload={"measurements": [{"ccs": 500.0}]},
    )
    after = store.get("pred_one")
    assert after.payload == frozen.payload
    assert after.payload_digest == frozen.payload_digest


def test_a_second_attach_appends_rather_than_replacing(store):
    freeze(store)
    for n in (1, 2, 3):
        store.attach(
            prediction_id="pred_one",
            validation_id=f"val_{n}",
            created_at=f"2026-09-27T00:00:0{n}+00:00",
            fingerprint="params/corpus",
            snapshot_digest="snapshot",
            payload={"measurements": [{"ccs": 500.0 + n}]},
        )
    records = store.validations("pred_one")
    assert len(records) == 3, "attaching is append-only; a second attach must not replace the first"
    assert [one.validation_id for one in records] == ["val_1", "val_2", "val_3"]
    assert [one.body["measurements"][0]["ccs"] for one in records] == [501.0, 502.0, 503.0]


def test_each_validation_records_the_prediction_digest_it_was_attached_to(store):
    frozen = freeze(store)
    record = store.attach(
        prediction_id="pred_one",
        validation_id="val_one",
        created_at=now(),
        fingerprint="params/corpus",
        snapshot_digest="snapshot",
        payload={"measurements": []},
    )
    # So that if the prediction were ever modified, this would no longer match and the mismatch
    # is discoverable rather than lost.
    assert record.predicted_digest == frozen.payload_digest


def test_attaching_to_a_prediction_that_does_not_exist_is_refused(store):
    with pytest.raises(UnknownPrediction, match="no prediction"):
        store.attach(
            prediction_id="nope",
            validation_id="val_one",
            created_at=now(),
            fingerprint="f",
            snapshot_digest="s",
            payload={},
        )
    assert "{prediction_id!r}" in NO_SUCH_PREDICTION


def test_the_payload_moved_guard_exists_and_names_both_digests():
    # It can only fire if a future edit starts touching the prediction inside `attach`, which is
    # the change this module exists to make impossible - so the wording is pinned even though no
    # input can reach it today. An unreachable guard with no test is the shape this project
    # keeps meeting; an unreachable guard whose CONTRACT is tested is the honest version.
    assert "{before}" in PAYLOAD_MOVED and "{after}" in PAYLOAD_MOVED
    assert "must never touch the prediction" in PAYLOAD_MOVED


def test_a_validation_cannot_reference_a_missing_prediction_at_the_schema_level(store):
    # Foreign keys are OFF by default in SQLite. The schema declares the reference; the PRAGMA
    # makes it real, and this proves the PRAGMA took rather than trusting it.
    with pytest.raises(sqlite3.IntegrityError):
        store._db.execute(
            "INSERT INTO validations (validation_id, prediction_id, created_at, fingerprint,"
            " snapshot_digest, payload, payload_digest, predicted_digest)"
            " VALUES ('v', 'ghost', 'now', 'f', 's', '{}', 'd', 'p')"
        )


# --- it survives a restart -------------------------------------------------------------------------


def test_a_frozen_prediction_survives_a_restart(tmp_path: Path):
    database = tmp_path / "nested" / "glycan.sqlite3"
    first = Store(database)
    frozen = freeze(first, client_reference="run-A")
    first.attach(
        prediction_id="pred_one",
        validation_id="val_one",
        created_at=now(),
        fingerprint="params/corpus",
        snapshot_digest="snapshot",
        payload={"measurements": [{"ccs": 500.0}]},
    )
    first.record_run(
        run_id="run_one",
        created_at=now(),
        kind="prediction.created",
        fingerprint="params/corpus",
        snapshot_digest="snapshot",
        detail="six candidates",
        prediction_id="pred_one",
        composition="Hex5HexNAc2",
    )
    first.close()

    # A new process would open the file exactly like this.
    second = Store(database)
    assert second.get("pred_one").payload == frozen.payload
    assert second.get("pred_one").payload_digest == frozen.payload_digest
    assert len(second.validations("pred_one")) == 1
    assert len(second.runs()) == 1
    # And it is STILL frozen after the restart, which is the property that matters: immutability
    # that only held in memory would be no immutability at all.
    with pytest.raises(AlreadyFrozen):
        freeze(second, payload=OTHER)
    assert second.get("pred_one").payload == frozen.payload
    second.close()


def test_the_database_directory_is_created_if_it_is_missing(tmp_path: Path):
    database = tmp_path / "a" / "b" / "c" / "glycan.sqlite3"
    store = Store(database)
    freeze(store)
    store.close()
    assert database.is_file()


# --- the run log -----------------------------------------------------------------------------------


def _log(store: Store, n: int, kind: str, composition: str, prediction_id: str | None = None):
    return store.record_run(
        run_id=f"run_{n}",
        created_at=f"2026-09-27T00:00:{n:02d}+00:00",
        kind=kind,
        fingerprint="params/corpus",
        snapshot_digest="snapshot",
        detail=f"detail {n}",
        prediction_id=prediction_id,
        composition=composition,
    )


def test_runs_are_returned_newest_first(store):
    for n in range(1, 4):
        _log(store, n, "prediction.created", "Hex5HexNAc2")
    assert [one.run_id for one in store.runs()] == ["run_3", "run_2", "run_1"]


def test_every_run_filter_composes_with_and(store):
    _log(store, 1, "prediction.created", "Hex5HexNAc2")
    _log(store, 2, "validation.attached", "Hex5HexNAc2")
    _log(store, 3, "prediction.created", "Hex5HexNAc4Fuc1")
    assert len(store.runs(kind="prediction.created")) == 2
    assert len(store.runs(composition="Hex5HexNAc2")) == 2
    assert len(store.runs(kind="prediction.created", composition="Hex5HexNAc2")) == 1
    assert len(store.runs(since="2026-09-27T00:00:02+00:00")) == 2
    assert len(store.runs(kind="nothing.like.this")) == 0


def test_a_filter_value_cannot_reach_the_query_as_syntax(store):
    # Parameterised throughout. If a filter were interpolated, this would drop the table and the
    # next assertion would raise instead of returning nothing.
    _log(store, 1, "prediction.created", "Hex5HexNAc2")
    assert store.runs(composition="'; DROP TABLE runs; --") == ()
    assert len(store.runs()) == 1


def test_paging_reports_the_total_before_the_page(store):
    for n in range(1, 8):
        _log(store, n, "prediction.created", "Hex5HexNAc2")
    page = store.runs(limit=3, offset=0)
    assert len(page) == 3
    assert store.count_runs() == 7
    assert store.count_runs(kind="prediction.created") == 7
    assert store.count_runs(kind="validation.attached") == 0
    # The last page must be distinguishable from a full one, which is what the total is for.
    assert len(store.runs(limit=3, offset=6)) == 1


# --- digests ----------------------------------------------------------------------------------------


def test_the_digest_does_not_depend_on_key_order(store):
    assert digest_of({"a": 1, "b": 2}) == digest_of({"b": 2, "a": 1})


def test_a_string_payload_is_digested_as_given(store):
    # Because a string is already the bytes that were served: re-rendering it could change them.
    assert digest_of('{"a":1}') == digest_of('{"a":1}')
    assert digest_of('{"a":1}') != digest_of('{"a": 1}')


def test_an_id_carries_no_ordering(store):
    # uuid4, so a caller cannot read sequence or volume out of an id.
    ids = {new_id("pred") for _ in range(50)}
    assert len(ids) == 50
    assert all(one.startswith("pred_") for one in ids)
