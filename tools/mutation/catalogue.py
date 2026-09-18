"""The mutations. One entry per behaviour that is supposed to have a test behind it.

Hand-written, never generated. Each entry was authored against a specific guard,
and the label says what would go WRONG rather than which operator was flipped.

Anchors are literal strings and must identify exactly ONE site in their file.
When an ordinary edit moves the line an anchor names, the mutation must be
RE-ANCHORED, not deleted: tests/test_mutation_catalogue.py fails on every pytest
run until it is. Deleting one is a deliberate act, and the comment left behind
must say what now covers its intent.

THIS CATALOGUE IS NEW, AND SMALLER THAN THE ONE IT REPLACES.
The glycan platform's catalogue held 154 mutations. Forty-six of those anchored
into modules that do not come across at all - the isomer enumerator, the glycan
graph, the featuriser, the enzyme table - and another twenty-six into splits.py,
evaluation.py and contracts.py, which are not in this milestone. Counting a
mutation that cannot exist here as a loss would be as misleading as counting one
that never applies as a pass. The floor in test_mutation_catalogue.py is set to
what this catalogue actually holds, and it is a floor: it goes up as modules
arrive, never quietly down.

Sections:
  [I] identity.py   - the analyte union and the matched-ion key
  [M] models.py     - the measurement record and its validators
  [L] licensing.py  - the gate
  [R] reuse.py      - the status tiers
  [S] sources.py    - the registry and the claim checks
  [D] loader.py     - reading a file without losing or inventing a row
  [K] readiness.py  - what refuses and what merely warns
  [C] models.py     - cyclic ion mobility (M1)
  [G] identity.py   - the glycopeptide (M1)
  [B] identity.py   - antibody, ADC and protein subunit identity (M1)
  [X] matching.py   - matched-ion construction (M2)
"""

from __future__ import annotations

from . import Expect, Mutation

MUTATIONS: tuple[Mutation, ...] = (
    # --- [I] identity: the matched-ion key -----------------------------------
    Mutation(
        label="[I] adduct components are not sorted, so two spellings of one ion give two keys",
        file="identity.py",
        find='    parts.sort(key=lambda part: (part[0] != "+", part[2]))',
        replace="    pass",
    ),
    Mutation(
        label="[I] an adduct naming one species twice is merged instead of refused",
        file="identity.py",
        find="    if repeated:\n        raise ValueError(\n            f\"adduct {adduct!r} names {', '.join(repeated)} more than once; write each species once, with its count\"\n        )",
        replace="    if False:\n        raise ValueError(\n            f\"adduct {adduct!r} names {', '.join(repeated)} more than once; write each species once, with its count\"\n        )",
    ),
    Mutation(
        label="[I] a charge carrier the source never named is reported as named",
        file="identity.py",
        find="    return any(species == UNSTATED_CARRIER for _sign, _count, species in adduct_components(adduct))",
        replace="    return False",
    ),
    Mutation(
        label="[I] the placeholder list stops rejecting fillers",
        file="identity.py",
        find='    return value.strip().casefold() in _PLACEHOLDERS or not any(ch.isalnum() for ch in value)',
        replace="    return False",
    ),
    Mutation(
        label="[I] a record built without validation is trusted anyway",
        file="identity.py",
        find="    if not unchanged:\n        return f\"it differs from its validated form, so it was built or changed without validation; {advice}\"\n    return None",
        replace="    return None",
    ),
    Mutation(
        label="[I] a glycan identified by nothing at all is accepted",
        file="identity.py",
        find="        if self.composition is None and self.wurcs is None and self.glytoucan_ac is None and self.iupac_condensed is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[I] a glycan keys on its composition even when its structure is resolved",
        file="identity.py",
        find='        if self.wurcs is not None:\n            return (AnalyteKind.GLYCAN.value, "wurcs", self.wurcs)',
        replace='        if False:\n            return (AnalyteKind.GLYCAN.value, "wurcs", self.wurcs)',
        also=(
            (
                '        if self.iupac_condensed is not None and _IUPAC_LINKAGE.search(self.iupac_condensed) is not None:\n            return (AnalyteKind.GLYCAN.value, "iupac", self.iupac_condensed)',
                '        if False:\n            return (AnalyteKind.GLYCAN.value, "iupac", self.iupac_condensed)',
            ),
        ),
    ),
    Mutation(
        label="[I] a glycan's structural state drops the derivatisation",
        file="identity.py",
        find="        return (str(self.reducing_end_label), str(self.derivatisation))",
        replace="        return (str(self.reducing_end_label),)",
    ),
    Mutation(
        label="[I] an open reducing-end label no longer blocks training",
        file="identity.py",
        find="        if is_one_of(label, OPEN_LABELS):",
        replace="        if False:",
    ),
    Mutation(
        label="[I] a sialic-acid derivatisation is accepted on a composition with no sialic acid",
        file="identity.py",
        find="            and self.composition.neuac + self.composition.neugc == 0\n        ):",
        replace="            and False\n        ):",
    ),
    Mutation(
        label="[I] a sialic-acid derivatisation is accepted where no composition can check it",
        file="identity.py",
        find="        if self.derivatisation in SIALIC_ACID_DERIVATISATIONS and self.composition is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[I] an unstated folding state stops blocking training",
        file="identity.py",
        find='        if is_one_of(getattr(self, "folding_state", None), UNDEFINED_FOLDING_STATES):',
        replace="        if False:",
    ),
    Mutation(
        label="[I] an ambiguity code is accepted in a peptide sequence",
        file="identity.py",
        find='_PEPTIDE_SEQUENCE = re.compile(r"[ACDEFGHIKLMNPQRSTVWYUO]+")',
        replace='_PEPTIDE_SEQUENCE = re.compile(r"[A-Z]+")',
    ),
    Mutation(
        label="[I] anything is accepted as an InChIKey",
        file="identity.py",
        find="    if not _INCHIKEY.fullmatch(value):",
        replace="    if False:",
    ),
    Mutation(
        label="[I] modifications are left in the order written, so one peptide gets two keys",
        file="identity.py",
        find="    return tuple(sorted(values))",
        replace="    return tuple(values)",
    ),
    Mutation(
        label="[I] a composition is left in the order written, so one glycan gets two keys",
        file="identity.py",
        find='        return "".join(f"{name}{n}" for name, n in self.counts().items())',
        replace='        return "".join(f"{name}{n}" for name, n in sorted(self.counts().items(), reverse=True))',
    ),
    Mutation(
        label="[I] an antibody INN is not normalised, so two spellings give two keys",
        file="identity.py",
        find="        if self.inn is not None and self.inn != self.inn.casefold():",
        replace="        if False:",
    ),
    Mutation(
        label="[I] an antibody identified by nothing at all is accepted",
        file="identity.py",
        find="        if self.inn is None and self.accession is None and self.sequence is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[I] a protein identified only by a display name is accepted",
        file="identity.py",
        find="        if self.accession is None and self.sequence is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[I] an ADC whose loading is unrecorded is accepted",
        file="identity.py",
        find="        if self.dar is None and self.conjugation_state is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[I] an ADC's key drops the payload class, pooling two conjugates of one antibody",
        file="identity.py",
        find="            *self.antibody.key(),\n            self.linker_payload_class,\n            *loading,",
        replace="            *self.antibody.key(),\n            *loading,",
    ),
    Mutation(
        label="[I] an antibody's key drops the glycoform",
        file="identity.py",
        find="        return (AnalyteKind.INTACT_ANTIBODY.value, *self.antibody.key(), self.glycoform)",
        replace="        return (AnalyteKind.INTACT_ANTIBODY.value, *self.antibody.key())",
    ),
    Mutation(
        label="[I] an antibody's structural state drops the CIU state",
        file="identity.py",
        # Re-anchored once already: component_records was removed from the two
        # antibody variants when it was found to refuse every antibody record.
        find="    def structural_state(self) -> tuple:\n        return (str(self.folding_state), self.ciu_state)\n\n\nclass ADCAnalyte",
        replace="    def structural_state(self) -> tuple:\n        return (str(self.folding_state),)\n\n\nclass ADCAnalyte",
    ),
    Mutation(
        label="[I] the analyte union loses its discriminator, so an arm can silently change",
        file="identity.py",
        find='    Field(discriminator="kind_tag"),',
        replace="    Field(),",
    ),
    # --- [M] models: the measurement and its key ------------------------------
    Mutation(
        label="[M] the matched-ion key carries the platform, so no pair can ever be found",
        file="models.py",
        find="        return MatchedIonKey(\n            analyte=self.analyte.identity_key(),",
        replace="        return MatchedIonKey(\n            analyte=(self.ims_type, *self.analyte.identity_key()),",
    ),
    Mutation(
        label="[M] the matched-ion key uses the gas in the cell, not the gas the value refers to",
        file="models.py",
        find="            drift_gas=self.drift_gas,\n            state=self.analyte.structural_state(),",
        replace="            drift_gas=self.cell_gas,\n            state=self.analyte.structural_state(),",
    ),
    Mutation(
        label="[M] the matched-ion key drops the charge",
        file="models.py",
        find="            adduct=self.adduct,\n            charge=self.charge,",
        replace="            adduct=self.adduct,\n            charge=0,",
    ),
    Mutation(
        label="[M] the matched-ion key drops the structural state, pooling native with denatured",
        file="models.py",
        # Re-anchored once already: the key gained its `unmatchable` component.
        find="            state=self.analyte.structural_state(),\n            unmatchable=unmatchable,",
        replace="            state=(),\n            unmatchable=unmatchable,",
    ),
    Mutation(
        label="[M] the calibration group carries the analyte, so the shared-peak check can never fire",
        file="models.py",
        # Rewritten after the first sweep, where this survived for a reason that
        # had nothing to do with test coverage: it read
        # `self.analyte.identity_key() and self.ims_type`, and since an identity
        # key is always a non-empty tuple that expression just returns ims_type.
        # The text changed and the behaviour did not, so the harness called it
        # applied. Anchor.INERT catches a mutation that leaves the FILE
        # unchanged; nothing catches one that leaves the BEHAVIOUR unchanged,
        # which is worth knowing when writing a replacement that looks clever.
        find="        return CalibrationGroup(\n            self.ims_type,",
        replace="        return CalibrationGroup(\n            self.analyte.identity_key(),",
    ),
    Mutation(
        label="[M] the calibration group drops the calibrant",
        file="models.py",
        find="            self.calibrant,\n            self.adduct,\n            self.analyte.structural_state(),",
        replace="            None,\n            self.adduct,\n            self.analyte.structural_state(),",
    ),
    Mutation(
        label="[M] a spread is stored with no statement of what it is",
        file="models.py",
        find="        if self.ccs_uncertainty is not None and self.uncertainty_type is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] a type that describes a number is accepted with no number to describe",
        file="models.py",
        find="            and self.uncertainty_type is not UncertaintyType.UNKNOWN\n            and self.ccs_uncertainty is None",
        replace="            and self.uncertainty_type is UncertaintyType.UNKNOWN\n            and self.ccs_uncertainty is None",
    ),
    Mutation(
        label="[M] an uncertainty of unknown kind stops blocking training",
        file="models.py",
        find="        if kind is not None and not is_one_of(kind, TRAINABLE_UNCERTAINTY_TYPES):",
        replace="        if False:",
    ),
    Mutation(
        label="[M] an undefined drift gas stops blocking training",
        file="models.py",
        find='        if is_one_of(getattr(self, "drift_gas", None), UNDEFINED_GASES):',
        replace="        if False:",
    ),
    Mutation(
        label="[M] an adduct that does not name its charge carrier stops blocking training",
        file="models.py",
        find="            if unstated_carrier:",
        replace="            if False:",
    ),
    Mutation(
        label="[M] an adduct may carry a charge the record disagrees with",
        file="models.py",
        find="        if adduct_charge != self.charge:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] a charge may disagree with the stated polarity",
        file="models.py",
        find="        if (self.charge > 0) != (self.polarity is Polarity.POSITIVE):",
        replace="        if False:",
    ),
    Mutation(
        label="[M] a DTIMS record need not say whether it is stepped-field or single-field",
        file="models.py",
        find="        if self.ims_type is IMSType.DTIMS and self.dtims_method is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] a calibrated platform may record no calibrant",
        file="models.py",
        find="        if self.ccs_is_calibrated and self.calibrant is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] a primary stepped-field value may name a calibrant",
        file="models.py",
        find="        if not self.ccs_is_calibrated and self.calibrant is not None:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] single-field DTIMS is treated as primary",
        file="models.py",
        find="            return self.dtims_method is not DTIMSMethod.STEPPED_FIELD",
        replace="            return False",
    ),
    Mutation(
        label="[M] a primary value may refer to a gas it was not measured in",
        file="models.py",
        find="        if self.cell_gas is not None and not self.ccs_is_calibrated and self.cell_gas is not self.drift_gas:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] a conformer index may exceed the total the source resolved",
        file="models.py",
        find="        if self.conformer is not None and self.conformers_total is not None and self.conformer > self.conformers_total:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] a conformer index is accepted with no total beside it",
        file="models.py",
        find="        if self.conformer is not None and self.conformers_total is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] a row claiming several conformers need not say which one it is",
        file="models.py",
        find="        if self.conformers_total is not None and self.conformers_total > 1 and self.conformer is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] a primary value may carry a calibration reference it cannot have",
        file="models.py",
        find="        if self.calibration_reference is not None and not self.ccs_is_calibrated:",
        replace="        if False:",
    ),
    Mutation(
        label="[M] an unrecorded calibration chain is reported as reaching first principles",
        file="models.py",
        find="        if self.calibration_reference is None:\n            return None",
        replace="        if self.calibration_reference is None:\n            return True",
    ),
    Mutation(
        label="[M] a reference set of unstated platform is reported as primary",
        file="models.py",
        find="        if self.platform is None:\n            return None",
        replace="        if self.platform is None:\n            return True",
    ),
    Mutation(
        label="[M] single-field DTIMS reference values are reported as primary",
        file="models.py",
        find="        return self.method is DTIMSMethod.STEPPED_FIELD",
        replace="        return True",
    ),
    Mutation(
        label="[M] the analyte is not offered to the licence gate as a component",
        file="models.py",
        find="        parts: list = [self.analyte]",
        replace="        parts: list = []",
    ),
    # --- [L] licensing: the gate ----------------------------------------------
    Mutation(
        label="[L] the gate reports a refusal instead of raising one",
        file="licensing.py",
        find="    if licence or unbacked or other:",
        replace="    if False:",
    ),
    Mutation(
        label="[L] the records a measurement is built from are never checked",
        file="licensing.py",
        find='    components = getattr(subject, "component_records", None)',
        replace="    components = None",
    ),
    Mutation(
        label="[L] a record's own training blockers are never asked for",
        file="licensing.py",
        find='    blockers = getattr(subject, "training_blockers", None)',
        replace="    blockers = None",
    ),
    Mutation(
        label="[L] a record whose components cannot be listed is cleared instead of refused",
        file="licensing.py",
        find="        except Exception as exc:  # fail closed: licences that cannot be read are not cleared\n            licence.append(",
        replace="        except Exception as exc:  # fail closed: licences that cannot be read are not cleared\n            [].append(",
    ),
    Mutation(
        label="[L] an unbacked claim is reported as a licence refusal, pointing at the wrong remedy",
        file="licensing.py",
        find="        error = LicenceGateError if licence else UnbackedClaimError if unbacked else TrainingGateError",
        replace="        error = LicenceGateError",
    ),
    Mutation(
        label="[L] a synthetic declaration on a component is not counted",
        file="licensing.py",
        find='        components = getattr(record, "component_records", None)\n        if components is None:\n            return False',
        replace="        components = None\n        if components is None:\n            return False",
    ),
    Mutation(
        label="[L] a claim that cannot be checked is waved through instead of refused",
        file="licensing.py",
        find='        return f"its reuse-status claim could not be checked ({type(exc).__name__}: {exc})"',
        replace="        return None",
    ),
    # --- [R] reuse: the tiers --------------------------------------------------
    Mutation(
        label="[R] an unverified record becomes trainable",
        file="reuse.py",
        find="_TRAINABLE = frozenset(\n    {ReuseStatus.OPEN_ATTRIBUTION, ReuseStatus.INTERNAL_PROPRIETARY, ReuseStatus.SYNTHETIC_FIXTURE}\n)",
        replace="_TRAINABLE = frozenset(\n    {\n        ReuseStatus.OPEN_ATTRIBUTION,\n        ReuseStatus.INTERNAL_PROPRIETARY,\n        ReuseStatus.SYNTHETIC_FIXTURE,\n        ReuseStatus.UNVERIFIED,\n    }\n)",
    ),
    Mutation(
        label="[R] a no-derivatives record becomes usable at inference",
        file="reuse.py",
        find="_INFERENCE_ONLY = frozenset(\n    {ReuseStatus.OPEN_SHARE_ALIKE, ReuseStatus.NON_COMMERCIAL, ReuseStatus.ACADEMIC_ONLY}\n)",
        replace="_INFERENCE_ONLY = frozenset(\n    {\n        ReuseStatus.OPEN_SHARE_ALIKE,\n        ReuseStatus.NON_COMMERCIAL,\n        ReuseStatus.ACADEMIC_ONLY,\n        ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES,\n    }\n)",
    ),
    Mutation(
        label="[R] an unrecognised status string is read as trainable rather than refused",
        file="reuse.py",
        find='            raise ValueError(f"unknown reuse status {status!r}; expected one of: {known}") from None',
        replace="            return ReuseStatus.OPEN_ATTRIBUTION",
    ),
    # --- [S] sources: the registry and the claim checks ------------------------
    Mutation(
        label="[S] a claim that exceeds its licence record is accepted",
        file="sources.py",
        find="    if entry.reuse_status is not claimed:",
        replace="    if False:",
    ),
    Mutation(
        label="[S] a trainable claim on a DOI with no licence record is accepted",
        file="sources.py",
        find="    if entry is None:\n        return NO_RECORD_FOR_DOI.format(claimed=claimed, doi=doi)",
        replace="    if entry is None:\n        return None",
    ),
    Mutation(
        label="[S] a published claim with no DOI at all is accepted",
        file="sources.py",
        find="    if not doi:\n        return NO_DOI_FOR_PUBLISHED_CLAIM.format(claimed=claimed)",
        replace="    if not doi:\n        return None",
    ),
    Mutation(
        label="[S] a row's DOI backs an analyte identity that names a different source",
        file="sources.py",
        find="    names_the_paper = (\n        bool(doi)\n        and component_source is not None\n        and measurement_source is not None\n        and component_source.strip() == measurement_source.strip()\n    )",
        replace="    names_the_paper = bool(doi)",
    ),
    Mutation(
        label="[S] an analyte claim backed by nothing at all is accepted",
        file="sources.py",
        find="    return COMPONENT_CLAIM_UNBACKED.format(claimed=claimed, source=component_source)",
        replace="    return None",
    ),
    Mutation(
        label="[S] a dataset provenance is matched loosely, so any spelling backs a claim",
        file="sources.py",
        find="    return DATASETS.get(provenance.strip())",
        replace="    return next((entry for key, entry in DATASETS.items() if provenance.strip() in key), None)",
    ),
    Mutation(
        label="[S] an analyte claim exceeding its dataset's licence is accepted",
        file="sources.py",
        find="        if dataset.reuse_status is claimed:\n            return None",
        replace="        if True:\n            return None",
    ),
    Mutation(
        label="[S] an unverified source is recorded without saying what would settle it",
        file="sources.py",
        find="        if self.reuse_status is ReuseStatus.UNVERIFIED and not self.what_to_check.strip():",
        replace="        if False:",
    ),
    Mutation(
        label="[S] a licence that permits training is recorded with no attribution to discharge",
        file="sources.py",
        find="        if can_train_commercial(self.reuse_status) and not self.attribution.strip():",
        replace="        if False:",
    ),
    Mutation(
        label="[S] an entry is accepted with no evidence behind it",
        file="sources.py",
        find="        if not self.evidence or not all(item.strip() for item in self.evidence):",
        replace="        if False:",
    ),
    Mutation(
        label="[S] an entry is accepted with nobody named as having read it",
        file="sources.py",
        find="        if not self.reported_by.strip():",
        replace="        if False:",
    ),
    Mutation(
        label="[S] an entry is accepted that nobody can navigate to",
        file="sources.py",
        find="        if not (self.doi or self.url):",
        replace="        if False:",
    ),
    Mutation(
        label="[S] evidence is left as whatever was passed, so a frozen record holds a mutable field",
        file="sources.py",
        find='        object.__setattr__(self, "evidence", tuple(self.evidence))',
        replace="        pass",
    ),
    Mutation(
        label="[S] the registry becomes a plain dict a caller can append to",
        file="sources.py",
        find="REGISTRY: Mapping[str, SourceLicence] = MappingProxyType({entry.key: entry for entry in _ENTRIES})",
        replace="REGISTRY: Mapping[str, SourceLicence] = {entry.key: entry for entry in _ENTRIES}",
    ),
    Mutation(
        label="[S] a DOI is matched case-sensitively, so one spelling of a DOI loses its record",
        file="sources.py",
        find="    return REGISTRY.get(doi.strip().casefold())",
        replace="    return REGISTRY.get(doi.strip())",
    ),
    # --- [D] loader: reading a file without losing or inventing a row ----------
    Mutation(
        label="[D] a row with the wrong number of cells is skipped without being counted",
        file="loader.py",
        find="            if len(cells) != len(fields):\n                fault = f\"{len(cells)} cells against {len(fields)} header columns\"\n                rows.append(Row(line=start, cells={}, fault=fault))\n                continue",
        replace="            if len(cells) != len(fields):\n                continue",
    ),
    Mutation(
        label="[D] a row folded across two lines by a stray quote is read as one row",
        file="loader.py",
        find="            if end != start:\n                raise ValueError(",
        replace="            if False:\n                raise ValueError(",
    ),
    Mutation(
        label="[D] a header naming one column twice is accepted, so one column's cells are lost",
        file="loader.py",
        find="        if repeated:\n            raise ValueError(f\"{label}: the header names a column more than once: {repeated}\")",
        replace="        if False:\n            raise ValueError(f\"{label}: the header names a column more than once: {repeated}\")",
    ),
    Mutation(
        label="[D] a file missing a required column is read as a clean run over nothing",
        file="loader.py",
        find="    if missing:\n        raise ValueError(f\"{label}: the header is missing the required columns {missing}\")",
        replace="    if missing:\n        pass",
    ),
    Mutation(
        label="[D] a synthetic-fixture status is accepted from a file",
        file="loader.py",
        find="        if column in _STATUS_COLUMNS and text.casefold() == ReuseStatus.SYNTHETIC_FIXTURE.value:",
        replace="        if False:",
    ),
    Mutation(
        label="[D] a placeholder cell is mapped to null instead of reaching validation",
        file="loader.py",
        find="        if not text:\n            continue\n        value: object = text",
        replace="        if not text or text.casefold() in {'n/a', 'n.d.', 'unknown', '-'}:\n            continue\n        value: object = text",
    ),
    Mutation(
        label="[D] a suspected shared peak is cleared for training",
        file="loader.py",
        find="        if where in shared:",
        replace="        if False:",
    ),
    Mutation(
        label="[D] a conformer whose siblings are missing from the file is cleared",
        file="loader.py",
        find="        if where in lone:",
        replace="        if False:",
    ),
    Mutation(
        label="[D] a row a person held for curation review is cleared anyway",
        file="loader.py",
        find="        if flag:\n            _note(gate, gate_examples, GATE_CURATION, f\"{where}: {flag}\")",
        replace="        if False:\n            _note(gate, gate_examples, GATE_CURATION, f\"{where}: {flag}\")",
    ),
    Mutation(
        label="[D] a row refused by the gate is dropped from the file rather than kept",
        file="loader.py",
        find="        records.append(record)\n        built.append((where, record, flag))",
        replace="        built.append((where, record, flag))",
    ),
    Mutation(
        label="[D] a shared peak against a row that did not validate is no longer held",
        file="loader.py",
        find="    for where, record, _ in built:\n        if where in flagged:\n            continue",
        replace="    for where, record, _ in built:\n        if True:\n            continue",
    ),
    Mutation(
        label="[D] conformer siblings pool across sample origins, so a lost sibling is not noticed",
        file="loader.py",
        find="        key = (record.calibration_group, record.analyte.identity_key(), total, record.source_locator)",
        replace="        key = (record.calibration_group, record.analyte.identity_key(), total)",
    ),
    Mutation(
        label="[D] a number written in a second spelling is accepted",
        file="loader.py",
        find="            if not _FLOAT.fullmatch(text) or not math.isfinite(float(text)):",
        replace="            if False:",
    ),
    Mutation(
        label="[D] a date written as a run of digits is read as a Unix timestamp",
        file="loader.py",
        find="            if not _ISO_DATE.fullmatch(text):\n                raise _Coercion(NOT_A_DATE, f\"{column}={text!r}\")",
        replace="            if False:\n                raise _Coercion(NOT_A_DATE, f\"{column}={text!r}\")",
    ),
    Mutation(
        label="[D] the analyte kind is guessed from the columns rather than read from the row",
        file="loader.py",
        find='    raw_kind = fields.pop("kind", None)',
        replace='    raw_kind = fields.pop("kind", None) or ("glycan" if "composition" in fields else "small_molecule")',
    ),
    Mutation(
        label="[D] an unidentified failed row is joined to every other unidentified one",
        file="loader.py",
        find='        for atom in atoms or (f"unidentified:{where}",):',
        replace='        for atom in atoms or ("unidentified",):',
    ),
    Mutation(
        label="[G] a glycopeptide's structural state drops the derivatisation",
        file="identity.py",
        find="    def structural_state(self) -> tuple:\n        return (str(self.derivatisation),)",
        replace="    def structural_state(self) -> tuple:\n        return ()",
    ),
    Mutation(
        label="[G] a glycopeptide's sialic-acid rule stops asking for a composition to check against",
        file="identity.py",
        find="        if self.derivatisation in SIALIC_ACID_DERIVATISATIONS and self.glycan_composition is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[G] a glycopeptide accepts a sialic-acid derivatisation on a glycan with no sialic acid",
        file="identity.py",
        find="            and self.glycan_composition.neuac + self.glycan_composition.neugc == 0",
        replace="            and False",
    ),
    Mutation(
        label="[G] a glycopeptide's identity atoms drop the attachment site",
        file="identity.py",
        find="        stem = f\"glycopeptide:{self.sequence}|{'+'.join(self.modifications)}|{self.attachment_site or ''}\"",
        replace="        stem = f\"glycopeptide:{self.sequence}|{'+'.join(self.modifications)}\"",
    ),
    Mutation(
        label="[B] an ADC's key drops the DAR, pooling two drug loads of one conjugate",
        file="identity.py",
        find='        loading = ("dar", self.dar) if self.dar is not None else ("conjugation", self.conjugation_state)',
        replace='        loading = ("dar",) if self.dar is not None else ("conjugation",)',
    ),
    Mutation(
        label="[B] an ADC's identity atoms drop the DAR",
        file="identity.py",
        find='        loading = f"dar={self.dar}" if self.dar is not None else f"conjugation={self.conjugation_state}"',
        replace='        loading = "dar" if self.dar is not None else "conjugation"',
    ),
    Mutation(
        label="[B] an ADC's structural state drops the CIU state",
        file="identity.py",
        find="    def structural_state(self) -> tuple:\n        return (str(self.folding_state), self.ciu_state)\n\n\n# The tagged union.",
        replace="    def structural_state(self) -> tuple:\n        return (str(self.folding_state),)\n\n\n# The tagged union.",
    ),
    Mutation(
        label="[B] a protein's key drops the subunit, pooling a light chain with a heavy chain",
        file="identity.py",
        find="        return (AnalyteKind.PROTEIN.value, *which, self.subunit)",
        replace="        return (AnalyteKind.PROTEIN.value, *which)",
    ),
    Mutation(
        label="[B] a protein's sequence atom drops the subunit, so atoms and the key disagree",
        file="identity.py",
        find="            atoms.add(f\"sequence:{self.sequence}|{self.subunit or ''}\")",
        replace='            atoms.add(f"sequence:{self.sequence}")',
    ),
    # --- [C] models.py: cyclic ion mobility (M1) -------------------------------
    Mutation(
        label="[C] a cyclic record may omit its cyclic settings, so it cannot be told from single-pass TWIMS",
        file="models.py",
        find="        if self.ims_type is IMSType.CYCLIC and self.cyclic is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[C] cyclic settings are accepted on a platform that has no passes",
        file="models.py",
        find="        if self.ims_type is not IMSType.CYCLIC and self.cyclic is not None:",
        replace="        if False:",
    ),
    Mutation(
        label="[C] the pass count leaves the calibration group, so six passes pool with one",
        file="models.py",
        find="            None if self.cyclic is None else (self.cyclic.passes, str(self.cyclic.pass_mode)),",
        replace="            None,",
    ),
    Mutation(
        label="[C] a pass count and a pass mode that contradict each other are accepted",
        file="models.py",
        find="        if self.pass_mode is not PassMode.UNSTATED and self.pass_mode is not expected:",
        replace="        if False:",
    ),
    Mutation(
        label="[C] a multipass value that never says whether ions lapped stops blocking training",
        file="models.py",
        find="        if self.is_multipass and self.wrap_around is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[C] wrap-around with no correction recorded stops blocking training",
        file="models.py",
        find="        if self.wrap_around is True and self.arrival_time_correction is None:",
        replace="        if False:",
    ),
    Mutation(
        label="[C] a cyclic record stops reporting its cyclic settings' blockers",
        file="models.py",
        find='        cyclic = getattr(self, "cyclic", None)\n        if cyclic is not None:',
        replace='        cyclic = getattr(self, "cyclic", None)\n        if False:',
    ),
    Mutation(
        label="[C] a single pass is read as multipass",
        file="models.py",
        find="            return self.passes > 1",
        replace="            return self.passes >= 1",
    ),
    # --- [G] identity.py: the glycopeptide (M1) --------------------------------
    Mutation(
        label="[G] a glycopeptide's key drops the attachment site, pooling two sites of one backbone",
        file="identity.py",
        find="            self.modifications,\n            self.attachment_site,\n            *self._glycan_key(),",
        replace="            self.modifications,\n            *self._glycan_key(),",
    ),
    Mutation(
        label="[G] a glycopeptide's key drops the glycan, pooling every glycoform of one backbone",
        file="identity.py",
        find="            self.attachment_site,\n            *self._glycan_key(),\n        )",
        replace="            self.attachment_site,\n        )",
    ),
    Mutation(
        label="[G] a glycopeptide emits a bare peptide atom, so it merges with its own backbone",
        file="identity.py",
        find='        stem = f"glycopeptide:{self.sequence}|{\'+\'.join(self.modifications)}|{self.attachment_site or \'\'}"',
        replace='        stem = f"peptide:{self.sequence}|{\'+\'.join(self.modifications)}"',
    ),
    Mutation(
        label="[G] a glycopeptide with no glycan at all is accepted",
        file="identity.py",
        find="            self.glycan_composition is None\n            and self.glycan_wurcs is None",
        replace="            False\n            and self.glycan_wurcs is None",
    ),
    # --- [X] matching.py: matched-ion construction (M2) ------------------------
    Mutation(
        label="[X] stepped-field and single-field DTIMS are read as one platform",
        file="matching.py",
        # A primary value and a calibrated one are not the same measurement in the
        # way that matters here, and folding them together loses a real comparison
        # AND hides the primary-versus-derived distinction from any later fit.
        find='    return str(ims) if method is None else f"{ims}/{method}"',
        replace="    return str(ims)",
    ),
    Mutation(
        label="[X] a matched set counts measurements rather than platforms, so replicates become matches",
        file="matching.py",
        find="        return len(self.platforms) >= 2",
        replace="        return len(self.measurements) >= 2",
    ),
    Mutation(
        label="[X] deduplication keys on the DOI, so one measurement republished becomes two",
        file="matching.py",
        find="    return (\n        getattr(record, \"matched_ion_key\", None),\n        _platform_of(record),",
        replace="    return (\n        getattr(record, \"matched_ion_key\", None),\n        _doi_of(record),\n        _platform_of(record),",
    ),
    Mutation(
        label="[X] deduplication drops the platform, so two laboratories agreeing are read as one measurement",
        file="matching.py",
        find="        _platform_of(record),\n        getattr(record, \"ccs\", None),",
        replace="        getattr(record, \"ccs\", None),",
    ),
    Mutation(
        label="[X] deduplication drops the conformer index, so two conformers collapse into one",
        file="matching.py",
        find='        getattr(record, "ccs", None),\n        getattr(record, "conformer", None),\n    )',
        replace='        getattr(record, "ccs", None),\n    )',
    ),
    Mutation(
        label="[X] grouping drops the conformer index, so two conformers become one matched set",
        file="matching.py",
        find="        group_key = (key, conformer)",
        replace="        group_key = (key,)",
    ),
    Mutation(
        label="[X] an ion whose charge carrier is unstated is allowed to pair",
        file="matching.py",
        find='        if key is not None and not getattr(key, "matchable", True):',
        replace="        if False:",
    ),
    Mutation(
        label="[X] a synthetic matched set may be quoted as a result",
        file="matching.py",
        find="        return not self.synthetic_groups",
        replace="        return True",
    ),
    Mutation(
        label="[X] assert_quotable stops refusing a synthetic report",
        file="matching.py",
        find="        raise NotQuotableError(refusal)",
        replace="        pass",
    ),
    Mutation(
        label="[X] a set counts only its outermost synthetic declarations",
        file="matching.py",
        find="        return sum(1 for m in self.measurements if declares_synthetic(m.record))",
        replace="        return sum(1 for m in self.measurements if str(getattr(m.record, 'reuse_status', '')) == 'synthetic_fixture' and False)",
    ),
    Mutation(
        label="[X] a member nobody may use stops blocking its matched set",
        file="matching.py",
        find="            elif not can_train_commercial(status):",
        replace="            elif False:",
    ),
    Mutation(
        label="[X] a member with an unreadable status is waved through instead of refused",
        file="matching.py",
        find="            if status is None:",
        replace="            if False:",
    ),
    Mutation(
        label="[X] conformer indices are paired across sources as though they corresponded",
        file="matching.py",
        find="        if self.conformer is not None and self.is_cross_platform:",
        replace="        if False:",
    ),
    Mutation(
        label="[X] the report calls itself quotable when only its matched sets are real",
        file="matching.py",
        find="        return tuple(group for group in self.all_groups if group.declares_synthetic)",
        replace="        return tuple(group for group in self.matched if group.declares_synthetic)",
    ),
    Mutation(
        label="[X] a blocked matched set is reported as usable",
        file="matching.py",
        find="        return tuple(group for group in self.matched if group.usable)",
        replace="        return self.matched",
    ),
    Mutation(
        label="[X] the republished citations are dropped rather than kept beside the measurement",
        file="matching.py",
        find="                also_published_as=tuple((_source_of(other), _doi_of(other)) for other in rest),",
        replace="                also_published_as=(),",
    ),
    # --- [K] readiness: what refuses and what merely warns ---------------------
    Mutation(
        label="[K] a single-platform corpus is reported as ready to compare",
        file="readiness.py",
        find="        if len(self.platforms) < MIN_PLATFORMS:",
        replace="        if False:",
    ),
    Mutation(
        label="[K] one platform is enough for a comparison",
        file="readiness.py",
        find="MIN_PLATFORMS = 2",
        replace="MIN_PLATFORMS = 1",
    ),
    Mutation(
        label="[K] a corpus with no ion on two platforms is reported as ready",
        file="readiness.py",
        find="        if self.matched_ions_multi_platform == 0:",
        replace="        if False:",
    ),
    Mutation(
        label="[K] a platform pair too small to cross-validate is reported as scorable",
        file="readiness.py",
        find="        if not any(count.scorable for count in self.by_pair.values()):",
        replace="        if False:",
    ),
    Mutation(
        label="[K] an empty corpus is reported as ready rather than as having nothing in it",
        file="readiness.py",
        find="        if self.records_cleared == 0:\n            return (NOTHING_HELD,)",
        replace="        if False:\n            return (NOTHING_HELD,)",
    ),
    Mutation(
        label="[K] the conformal calibration floor is off by one, so an unformable interval is allowed",
        file="readiness.py",
        find="    return math.ceil(1 / alpha) - 1",
        replace="    return math.ceil(1 / alpha) - 2",
    ),
    Mutation(
        label="[K] a calibration set too small for any finite interval is not refused",
        file="readiness.py",
        find="    if held < needed:",
        replace="    if False:",
    ),
    Mutation(
        label="[K] a matched-ion set accepts records that never went through the gate",
        file="readiness.py",
        find="            try:\n                assert_trainable(record)\n            except TrainingGateError as exc:",
        replace="            try:\n                pass\n            except TrainingGateError as exc:",
    ),
    Mutation(
        label="[K] a bare list is accepted where a gated set is required",
        file="readiness.py",
        find="    if not isinstance(matched_ion_set, MatchedIonSet):",
        replace="    if False:",
    ),
    Mutation(
        label="[K] synthetic fixtures are counted as real data",
        file="readiness.py",
        find="            synthetic_records=sum(1 for record in self.records if declares_synthetic(record)),",
        replace="            synthetic_records=0,",
    ),
    Mutation(
        label="[K] the matched-ion count is taken from keys held rather than keys on two platforms",
        file="readiness.py",
        # Re-anchored once already: the count runs over matchable keys only.
        find="            matched_ions_multi_platform=sum(1 for measured_on in matchable.values() if len(measured_on) > 1),",
        replace="            matched_ions_multi_platform=len(matchable),",
    ),
    # --- added after the first sweep, for behaviour that did not exist then ----
    Mutation(
        label="[I] two ions whose charge carrier is unstated silently match each other",
        file="identity.py",
        find="        return self.unmatchable is None",
        replace="        return True",
    ),
    Mutation(
        label="[M] an unstated charge carrier no longer makes its key unique to its record",
        file="models.py",
        find="        if adduct_carrier_is_unstated(self.adduct):",
        replace="        if False:",
    ),
    Mutation(
        label="[I] a glycan with no keyable identifier crashes instead of keying on what it states",
        file="identity.py",
        find='        if self.glytoucan_ac is not None:\n            return (AnalyteKind.GLYCAN.value, "glytoucan_unverified_form", self.glytoucan_ac)',
        replace='        if False:\n            return (AnalyteKind.GLYCAN.value, "glytoucan_unverified_form", self.glytoucan_ac)',
    ),
    Mutation(
        label="[S] a bare string is coerced into one piece of evidence per character",
        file="sources.py",
        find="            if isinstance(getattr(self, name), str):",
        replace="            if False:",
    ),
    Mutation(
        label="[D] an analyte kind the union does not hold is counted as bad data",
        file="loader.py",
        find="        except _Coercion as exc:\n            # Its own bucket.",
        replace="        except () as exc:\n            # Its own bucket.",
    ),
    Mutation(
        label="[D] a fully held file stops saying what its uncertainty types are",
        file="loader.py",
        find="        return _uncertainty_tally(self.records)",
        replace="        return {}",
    ),
    Mutation(
        label="[D] unmatchable records are folded back into the held matched-ion count",
        file="loader.py",
        find="        return len({record.matched_ion_key for record in self.records if record.matched_ion_key.matchable})",
        replace="        return len({record.matched_ion_key for record in self.records})",
    ),
    Mutation(
        label="[K] the maturity stamp reports a validated result",
        file="readiness.py",
        find="            data_maturity=DataMaturity.PROVISIONAL, matched_ion_count=self.matched_ions_multi_platform",
        replace="            data_maturity=DataMaturity.VALIDATED, matched_ion_count=self.matched_ions_multi_platform",
    ),
)
