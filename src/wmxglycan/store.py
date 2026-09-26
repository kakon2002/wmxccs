"""Persistence for predictions, validations and runs. SQLite, and immutable after freeze.

WHY SQLITE AND NOT A DATABASE SERVICE. The owner's instruction: whatever is simplest that
survives a restart, and do not add a service three days out. SQLite is in the standard
library, it is a single file, and it gives the one property that matters here - a frozen
prediction that is still frozen after the process dies.

IMMUTABILITY IS ENFORCED BY THE SCHEMA, NOT BY THE CALLER.

A frozen prediction is a row whose payload is never updated. That is not a convention here:

  - `predictions.payload` is written once, by an INSERT, and there is no UPDATE statement
    anywhere in this module. A second freeze of the same id raises `AlreadyFrozen` and
    writes nothing.
  - `validations` is APPEND ONLY. Attaching experimental measurements to a prediction adds
    a row; it does not touch the prediction. `attach` re-reads the payload digest afterwards
    and raises if it moved, so a future edit that did touch it cannot pass silently.
  - every row carries the pipeline fingerprint and the data snapshot digest it was made
    under, so a prediction frozen against one corpus is distinguishable from one frozen
    against another even if the composition and the answer are identical.

A CORRECTION IS A NEW ROW WITH ITS OWN PROVENANCE. There is no field anywhere that a
correction overwrites, which is the CCS core's rule three carried across a package boundary:
the original and the corrected value are both returned, and which is which is never in doubt
because they are separate rows with separate timestamps and separate sources.

WHY THE PAYLOAD IS STORED AS TEXT AND DIGESTED

The frozen payload is the JSON the API served, byte for byte, with its digest beside it. Two
reasons, and the second is the one that matters: a caller can be handed exactly what was
served rather than a re-rendering that a schema change might have altered; and the digest
makes "this prediction has not been modified" a checkable claim rather than a promise about
code nobody re-reads.

TIME IS PASSED IN, NEVER READ HERE. `created_at` is a parameter with no default, so a test
can write a row at a fixed instant and a caller cannot accidentally get two different clocks.
That also keeps this module free of the one thing that makes a store untestable.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

# Refusal wording is an interface: written once, imported by tests rather than matched as prose.
ALREADY_FROZEN = (
    "prediction {prediction_id!r} is already frozen and was not modified. A frozen prediction is"
    " immutable: the second write is REFUSED rather than merged, because merging would leave no"
    " way to tell which parts of the answer were produced by which pipeline. If the inputs have"
    " changed, create a new prediction; the two will differ in their fingerprint and their"
    " snapshot and both remain readable"
)
PAYLOAD_MOVED = (
    "the frozen payload of prediction {prediction_id!r} changed while a validation was being"
    " attached: its digest was {before} and is now {after}. Attaching experimental measurements"
    " must never touch the prediction, so this is refused and the validation is not recorded"
)
NO_SUCH_PREDICTION = "no prediction {prediction_id!r} exists"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS predictions (
    prediction_id     TEXT PRIMARY KEY NOT NULL,
    created_at        TEXT NOT NULL,
    composition       TEXT NOT NULL,
    adduct            TEXT NOT NULL,
    charge            INTEGER NOT NULL,
    client_reference  TEXT,
    pipeline_version  TEXT NOT NULL,
    fingerprint       TEXT NOT NULL,
    snapshot_digest   TEXT NOT NULL,
    payload           TEXT NOT NULL,
    payload_digest    TEXT NOT NULL
);
-- A client reference is the caller's own idempotency handle. UNIQUE so that a repeated
-- create is REFUSED at the storage layer rather than silently producing a second prediction
-- the caller believes is the same one.
CREATE UNIQUE INDEX IF NOT EXISTS predictions_client_reference
    ON predictions (client_reference) WHERE client_reference IS NOT NULL;
CREATE INDEX IF NOT EXISTS predictions_composition ON predictions (composition);
CREATE INDEX IF NOT EXISTS predictions_created_at ON predictions (created_at);

CREATE TABLE IF NOT EXISTS validations (
    validation_id     TEXT PRIMARY KEY NOT NULL,
    prediction_id     TEXT NOT NULL REFERENCES predictions (prediction_id),
    created_at        TEXT NOT NULL,
    fingerprint       TEXT NOT NULL,
    snapshot_digest   TEXT NOT NULL,
    payload           TEXT NOT NULL,
    payload_digest    TEXT NOT NULL,
    -- The prediction payload digest AS IT WAS when this validation was attached. If the
    -- prediction were ever modified, this would no longer match and the mismatch is
    -- discoverable rather than lost.
    predicted_digest  TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS validations_prediction ON validations (prediction_id, created_at);

CREATE TABLE IF NOT EXISTS runs (
    run_id            TEXT PRIMARY KEY NOT NULL,
    created_at        TEXT NOT NULL,
    kind              TEXT NOT NULL,
    prediction_id     TEXT,
    composition       TEXT,
    fingerprint       TEXT NOT NULL,
    snapshot_digest   TEXT NOT NULL,
    detail            TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS runs_created_at ON runs (created_at);
CREATE INDEX IF NOT EXISTS runs_kind ON runs (kind);
CREATE INDEX IF NOT EXISTS runs_composition ON runs (composition);
"""


class StoreError(RuntimeError):
    """Something the store refuses to do."""


class AlreadyFrozen(StoreError):
    """A second write to a frozen prediction. Refused, never merged."""


class UnknownPrediction(StoreError):
    """A prediction id nothing was ever frozen under."""


def digest_of(payload: Mapping | str) -> str:
    """sha256 over the payload as it will be stored.

    A mapping is rendered with sorted keys so the digest does not depend on insertion order;
    a string is digested as given, because a string is already the bytes that were served.
    """
    text = (
        payload
        if isinstance(payload, str)
        else json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def new_id(prefix: str) -> str:
    """A fresh identifier. uuid4, so ids carry no ordering a caller could read meaning into."""
    return f"{prefix}_{uuid.uuid4().hex}"


def now() -> str:
    """An ISO-8601 instant in UTC. The only place this module looks at a clock, and callers
    pass the result in, so nothing inside the store reads one."""
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


@dataclass(frozen=True)
class StoredPrediction:
    prediction_id: str
    created_at: str
    composition: str
    adduct: str
    charge: int
    client_reference: str | None
    pipeline_version: str
    fingerprint: str
    snapshot_digest: str
    payload: str
    payload_digest: str

    @property
    def body(self) -> dict:
        """The payload as it was served, parsed. Never re-rendered from live objects."""
        return json.loads(self.payload)


@dataclass(frozen=True)
class StoredValidation:
    validation_id: str
    prediction_id: str
    created_at: str
    fingerprint: str
    snapshot_digest: str
    payload: str
    payload_digest: str
    predicted_digest: str

    @property
    def body(self) -> dict:
        return json.loads(self.payload)


@dataclass(frozen=True)
class StoredRun:
    run_id: str
    created_at: str
    kind: str
    prediction_id: str | None
    composition: str | None
    fingerprint: str
    snapshot_digest: str
    detail: str


class Store:
    """The SQLite-backed store. One file, or ":memory:" for a test."""

    def __init__(self, path: Path | str = ":memory:") -> None:
        self.path = str(path)
        if self.path != ":memory:":
            Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        self._db = sqlite3.connect(self.path, check_same_thread=False)
        self._db.row_factory = sqlite3.Row
        # Foreign keys are OFF by default in SQLite, which would let a validation reference a
        # prediction that does not exist. The schema declares the reference; this makes it real.
        self._db.execute("PRAGMA foreign_keys = ON")
        self._db.executescript(_SCHEMA)
        self._db.commit()

    def close(self) -> None:
        self._db.close()

    # --- freezing, which happens exactly once per prediction ----------------------------------

    def freeze(
        self,
        *,
        prediction_id: str,
        created_at: str,
        composition: str,
        adduct: str,
        charge: int,
        pipeline_version: str,
        fingerprint: str,
        snapshot_digest: str,
        payload: Mapping,
        client_reference: str | None = None,
    ) -> StoredPrediction:
        """Write a prediction once. A second write with the same id or client reference raises.

        There is no `update` counterpart anywhere in this module, and that is the design rather
        than an omission.
        """
        text = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
        row = StoredPrediction(
            prediction_id=prediction_id,
            created_at=created_at,
            composition=composition,
            adduct=adduct,
            charge=charge,
            client_reference=client_reference,
            pipeline_version=pipeline_version,
            fingerprint=fingerprint,
            snapshot_digest=snapshot_digest,
            payload=text,
            payload_digest=digest_of(text),
        )
        try:
            self._db.execute(
                "INSERT INTO predictions (prediction_id, created_at, composition, adduct, charge,"
                " client_reference, pipeline_version, fingerprint, snapshot_digest, payload,"
                " payload_digest) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    row.prediction_id,
                    row.created_at,
                    row.composition,
                    row.adduct,
                    row.charge,
                    row.client_reference,
                    row.pipeline_version,
                    row.fingerprint,
                    row.snapshot_digest,
                    row.payload,
                    row.payload_digest,
                ),
            )
        except sqlite3.IntegrityError as clash:
            # NOTHING WAS WRITTEN. The message names which prediction stands, so a caller can
            # fetch it rather than guessing what happened.
            self._db.rollback()
            existing = self.by_client_reference(client_reference) if client_reference else None
            standing = existing.prediction_id if existing else prediction_id
            raise AlreadyFrozen(ALREADY_FROZEN.format(prediction_id=standing)) from clash
        self._db.commit()
        return row

    # --- attaching, which never touches the prediction ------------------------------------------

    def attach(
        self,
        *,
        prediction_id: str,
        validation_id: str,
        created_at: str,
        fingerprint: str,
        snapshot_digest: str,
        payload: Mapping,
    ) -> StoredValidation:
        """Append experimental measurements to a prediction without modifying it.

        The prediction's payload digest is read before and after and compared. That check can
        only fail if some future edit starts touching the prediction here, which is exactly the
        change this module exists to make impossible - so it is asserted rather than trusted.
        """
        before = self.get(prediction_id)
        if before is None:
            raise UnknownPrediction(NO_SUCH_PREDICTION.format(prediction_id=prediction_id))
        text = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
        row = StoredValidation(
            validation_id=validation_id,
            prediction_id=prediction_id,
            created_at=created_at,
            fingerprint=fingerprint,
            snapshot_digest=snapshot_digest,
            payload=text,
            payload_digest=digest_of(text),
            predicted_digest=before.payload_digest,
        )
        self._db.execute(
            "INSERT INTO validations (validation_id, prediction_id, created_at, fingerprint,"
            " snapshot_digest, payload, payload_digest, predicted_digest)"
            " VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (
                row.validation_id,
                row.prediction_id,
                row.created_at,
                row.fingerprint,
                row.snapshot_digest,
                row.payload,
                row.payload_digest,
                row.predicted_digest,
            ),
        )
        after = self.get(prediction_id)
        assert after is not None
        if after.payload_digest != before.payload_digest:
            self._db.rollback()
            raise StoreError(
                PAYLOAD_MOVED.format(
                    prediction_id=prediction_id,
                    before=before.payload_digest[:12],
                    after=after.payload_digest[:12],
                )
            )
        self._db.commit()
        return row

    # --- reading -------------------------------------------------------------------------------

    def get(self, prediction_id: str) -> StoredPrediction | None:
        found = self._db.execute(
            "SELECT * FROM predictions WHERE prediction_id = ?", (prediction_id,)
        ).fetchone()
        return None if found is None else StoredPrediction(**dict(found))

    def by_client_reference(self, client_reference: str) -> StoredPrediction | None:
        found = self._db.execute(
            "SELECT * FROM predictions WHERE client_reference = ?", (client_reference,)
        ).fetchone()
        return None if found is None else StoredPrediction(**dict(found))

    def validations(self, prediction_id: str) -> tuple[StoredValidation, ...]:
        """Every validation attached to a prediction, oldest first. Append-only, so this grows."""
        rows = self._db.execute(
            "SELECT * FROM validations WHERE prediction_id = ? ORDER BY created_at, validation_id",
            (prediction_id,),
        ).fetchall()
        return tuple(StoredValidation(**dict(row)) for row in rows)

    def predictions(self, limit: int = 50) -> tuple[StoredPrediction, ...]:
        rows = self._db.execute(
            "SELECT * FROM predictions ORDER BY created_at DESC, prediction_id DESC LIMIT ?",
            (limit,),
        ).fetchall()
        return tuple(StoredPrediction(**dict(row)) for row in rows)

    # --- the run log, append-only -------------------------------------------------------------

    def record_run(
        self,
        *,
        run_id: str,
        created_at: str,
        kind: str,
        fingerprint: str,
        snapshot_digest: str,
        detail: str,
        prediction_id: str | None = None,
        composition: str | None = None,
    ) -> StoredRun:
        row = StoredRun(
            run_id=run_id,
            created_at=created_at,
            kind=kind,
            prediction_id=prediction_id,
            composition=composition,
            fingerprint=fingerprint,
            snapshot_digest=snapshot_digest,
            detail=detail,
        )
        self._db.execute(
            "INSERT INTO runs (run_id, created_at, kind, prediction_id, composition, fingerprint,"
            " snapshot_digest, detail) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (
                row.run_id,
                row.created_at,
                row.kind,
                row.prediction_id,
                row.composition,
                row.fingerprint,
                row.snapshot_digest,
                row.detail,
            ),
        )
        self._db.commit()
        return row

    def runs(
        self,
        *,
        kind: str | None = None,
        composition: str | None = None,
        prediction_id: str | None = None,
        since: str | None = None,
        limit: int = 50,
        offset: int = 0,
    ) -> tuple[StoredRun, ...]:
        """Search the history. Every filter is optional and they compose with AND.

        Parameterised throughout: the filters are values, never interpolated SQL, so a
        composition string cannot reach the query as syntax.
        """
        where: list[str] = []
        values: list[object] = []
        for column, value in (
            ("kind", kind),
            ("composition", composition),
            ("prediction_id", prediction_id),
        ):
            if value is not None:
                where.append(f"{column} = ?")
                values.append(value)
        if since is not None:
            where.append("created_at >= ?")
            values.append(since)
        clause = f" WHERE {' AND '.join(where)}" if where else ""
        rows = self._db.execute(
            f"SELECT * FROM runs{clause} ORDER BY created_at DESC, run_id DESC LIMIT ? OFFSET ?",
            (*values, limit, offset),
        ).fetchall()
        return tuple(StoredRun(**dict(row)) for row in rows)

    def count_runs(self, **filters) -> int:
        """How many runs match, so a paged response can say what it is a page of."""
        return len(self.runs(limit=1_000_000, offset=0, **filters))

    def __iter__(self) -> Iterator[StoredPrediction]:
        return iter(self.predictions(limit=1_000_000))


def frozen_payloads(rows: Iterable[StoredPrediction]) -> Sequence[str]:
    """The stored payloads, for a test that wants to assert none of them moved."""
    return [row.payload for row in rows]
