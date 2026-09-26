"""Request and response bodies for the HTTP API."""

from enum import StrEnum
from typing import Literal, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .models import CompositionField, MeasurementConditions, check_derivatisation_fits


class DataMaturity(StrEnum):
    """How far the basis of a prediction has been taken.

    V1 ships "provisional": no model has been trained on experimental CCS, let
    alone checked against held-out measurements. The stamp travels with every
    prediction output so the caveat cannot be separated from the number.
    """

    PROVISIONAL = "provisional"
    VALIDATED = "validated"


class MaturityStamp(BaseModel):
    """The maturity of the data behind a prediction, carried in the prediction itself.

    Frozen, like every other record here. A stamp is routinely shared between
    reports, so a mutable one could be flipped to validated after the check that
    refused it, and every report holding it would re-render the new claim.
    """

    model_config = ConfigDict(frozen=True)

    data_maturity: DataMaturity = Field(
        description="'provisional' until a model trained on experimental CCS has been checked against held-out"
        " measurements; 'validated' only after that."
    )
    training_record_count: int = Field(
        ge=0,
        description="Experimental CCS measurements the loaded model was trained on. 0 while no model exists,"
        " which is why V1 is provisional.",
    )


class HealthResponse(BaseModel):
    status: Literal["ok"] = "ok"
    version: str
    model_loaded: bool = Field(
        description="Whether a trained CCS model is loaded. False until one has been trained on real data."
    )


class CompositionCheckRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    composition: CompositionField


class CompositionCheckResponse(BaseModel):
    canonical: str
    counts: dict[str, int]
    monoisotopic_mass: float = Field(
        description="Neutral monoisotopic mass of the free glycan (residue masses plus one water), Da."
    )
    has_full_man3glcnac2_core: bool = Field(
        description="The counts cover the full Man3GlcNAc2 core: at least 3 Hex and 2 HexNAc. False is a warning,"
        " not a rejection: truncated N-glycans exist."
    )
    warnings: list[str] = Field(description="Worth a look, but not grounds to rule the composition out.")
    plausible_n_glycan: bool = Field(
        description="No N-glycan rule is broken. Every N-glycan keeps at least one GlcNAc; the other rules assume"
        " mammalian biosynthesis and a glycan released with both core GlcNAc."
    )
    n_glycan_reasons: list[str] = Field(description="The N-glycan rules the composition breaks. Empty when plausible.")


class PredictRequest(MeasurementConditions):
    """A composition, and the conditions (calibration group) to predict its CCS under."""

    composition: CompositionField

    @model_validator(mode="after")
    def derivatisation_fits_composition(self) -> Self:
        check_derivatisation_fits(self.composition, self.derivatisation)
        return self


class PredictionUnavailable(BaseModel):
    """The 501 body. It carries the maturity stamp too, so the caveat ships with every prediction output."""

    error: Literal["not_implemented"] = "not_implemented"
    detail: str
    maturity: MaturityStamp
