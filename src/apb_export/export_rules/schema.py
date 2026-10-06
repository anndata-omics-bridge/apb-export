"""Export rule documents: apb2's rule schema, read from APB results instead of vendor tables.

Names follow apb2's ``rules.json``. Where a key keeps its name but its direction flips, its
description says so: an entry's ``name`` is the target field and its ``source`` the APB field.
"""

from __future__ import annotations

from typing import Annotated, Literal, Self

from pydantic import BaseModel, ConfigDict, Field, JsonValue, model_validator

type SchemaVersion = Literal["0.1"]
type ColumnType = Literal["string", "integer", "number", "boolean"]
type FloatWidth = Literal[32, 64]
type Scalar = str | int | float | bool | None
type ProformaColumn = Literal[
    "ProForma_peptide", "ProForma_peptidoform", "ProForma_ion", "ProForma_fragment"
]
"""The column names apb2's rule schema fixes for every vendor."""


class ModelBase(BaseModel):
    """Strict base: unknown keys are errors, as in apb2's rules."""

    model_config = ConfigDict(extra="forbid", frozen=True)


class RoleSource(ModelBase):
    """The var column carrying an apb2 role."""

    role: Literal["fasta_accessions", "protein_assignment"]


class AxisSource(ModelBase):
    """The level's observation or variable key column; a composite key is refused."""

    axis: Literal["obs_keys", "var_keys"]


class MeasurementSource(ModelBase):
    """The level's primary layer, or the abundance layer the caller picks."""

    measurements: Literal["primary_layer"]


class RuleSource(ModelBase):
    """A top-level field of the apb2 rule stored with the level, such as ``software_name``."""

    rule: str = Field(min_length=1)


class CatalogueSource(ModelBase):
    """A field found by meaning through apb-catalog."""

    catalogue: str = Field(min_length=1)
    concept: str = Field(min_length=1)
    kind: str = Field(min_length=1)
    location: Literal["layers", "var"] = "layers"
    level: Literal["protein"] | None = Field(
        default=None,
        description="Read the field from this level through each variable's protein_assignment, "
        "as a layer: a cell takes its protein group's value in its run; a var column's value "
        "repeats in every run. The assignment names its group whole, or else the group it "
        "leads, as MaxQuant's razor protein does. A group the level lacks reads as missing.",
    )
    accept: Literal["present", "nonzero"] = Field(
        default="present",
        description="nonzero rejects a layer holding no finite nonzero value, such as "
        "DIA-NN's library q-value without match-between-runs.",
    )


type Reference = (
    ProformaColumn | RoleSource | AxisSource | MeasurementSource | RuleSource | CatalogueSource
)
type Source = Reference | list[Reference]
"""One APB reference, or candidates in order: the first one that resolves is used."""


class EntryBase(ModelBase):
    """Facts shared by every column and ``uns`` entry."""

    name: str = Field(min_length=1, description="The target field; apb2's rules name APB fields.")
    type: ColumnType = "string"
    width: FloatWidth | None = Field(
        default=None, description="Float width of a number entry; 64 when omitted."
    )
    required: bool = Field(
        default=True,
        description="An absent source is an error; when false, the entry is skipped or filled.",
    )
    fill_value: Scalar = Field(
        default=None,
        description='Written when an optional source is absent; "NaN" for numbers.',
    )
    index_only: bool = Field(
        default=False,
        description="An axis entry that becomes the index without being written as a column.",
    )
    named_index: bool = Field(
        default=False,
        description="The index takes this axis entry's name, as AlphaPeptTools' levels need.",
    )

    @model_validator(mode="after")
    def _width_is_for_numbers(self) -> Self:
        if self.width is not None and self.type != "number":
            raise ValueError(f"entry {self.name!r}: width applies only to type number")
        return self


class SourcedEntry(EntryBase):
    """A target field read from one APB source."""

    source: Source = Field(description="The APB field; apb2's rules name a vendor column here.")


class Constant(EntryBase):
    """The same value everywhere."""

    how: Literal["constant"]
    value: JsonValue = Field(description="Any JSON in uns; a scalar in columns.")


class JoinNonempty(EntryBase):
    """Non-empty input values joined with a separator, as apb2's ``join_nonempty``."""

    how: Literal["join_nonempty"]
    inputs: list[Reference] = Field(min_length=2, description="APB references, not entries.")
    separator: str = Field(min_length=1)


class ProformaCharge(EntryBase):
    """The charge of a ProForma ion, ``SEQ/2`` giving 2."""

    how: Literal["proforma_charge"]
    inputs: list[Reference] = Field(min_length=1, max_length=1)


class SequenceLength(EntryBase):
    """The number of characters of a sequence."""

    how: Literal["sequence_length"]
    inputs: list[Reference] = Field(min_length=1, max_length=1)


class Accessions(EntryBase):
    """Protein-group members rewritten by a named ``accession_syntax``."""

    how: Literal["accessions"]
    inputs: list[Reference] = Field(min_length=1, max_length=1)
    syntax: str = Field(min_length=1)


class Contaminant(EntryBase):
    """Whether any protein-group member carries a contaminant marker of the syntax."""

    how: Literal["contaminant"]
    inputs: list[Reference] = Field(min_length=1, max_length=1)
    syntax: str = Field(min_length=1)


class Lowercase(EntryBase):
    """The input text in lower case."""

    how: Literal["lowercase"]
    inputs: list[Reference] = Field(min_length=1, max_length=1)


class AnnotationMap(EntryBase):
    """Every obs annotation column mapped to its own name; ``uns`` only."""

    how: Literal["annotation_map"]


type ComputedEntry = Annotated[
    Constant
    | JoinNonempty
    | ProformaCharge
    | SequenceLength
    | Accessions
    | Contaminant
    | Lowercase
    | AnnotationMap,
    Field(discriminator="how"),
]
type ColumnEntry = SourcedEntry | ComputedEntry


class AllValues(ModelBase):
    """Everything the APB level carries that the rule does not name.

    In ``varm``: every var column and layer value at the written rows, the level's key first.
    In ``columns.obs``: every obs annotation column, under its own name.
    """

    name: str = Field(min_length=1)
    how: Literal["all_values"]


class LayerEntry(ModelBase):
    """A target matrix read from one APB layer."""

    name: str = Field(min_length=1)
    source: Source = Field(description="An APB layer: the primary layer or a catalogued one.")
    width: FloatWidth = 64
    required: bool = False
    missing_values: list[float] = Field(
        default_factory=list,
        description="APB values treated as absent; apb2 lists vendor values here.",
    )


class Measurements(ModelBase):
    """The target's matrices; the primary layer becomes X."""

    primary_layer: str
    layers: list[LayerEntry] = Field(min_length=1)

    @model_validator(mode="after")
    def _primary_is_declared(self) -> Self:
        names = [layer.name for layer in self.layers]
        if len(names) != len(set(names)):
            raise ValueError("measurement layer names must be unique")
        if self.primary_layer not in names:
            raise ValueError(
                f"measurements.primary_layer={self.primary_layer!r} matches no layer; "
                f"available: {sorted(names)}"
            )
        return self


class Axis(ModelBase):
    """The entries that become the target's observation and variable index."""

    obs_keys: list[str] = Field(min_length=1, max_length=1)
    var_keys: list[str] = Field(min_length=1, max_length=1)


class Columns(ModelBase):
    """The target's observation and variable columns, and var-aligned tables."""

    obs: list[ColumnEntry | AllValues] = Field(default_factory=list)
    var: list[ColumnEntry] = Field(min_length=1)
    varm: list[AllValues] = Field(default_factory=list)


class AccessionSyntax(ModelBase):
    """How a target writes protein-group members, like apb2's ``sequence_syntax``."""

    separator: str = Field(min_length=1)
    contaminant_markers: list[str] = Field(default_factory=list)
    contaminant_prefix: str = ""


type Entries = list[ColumnEntry] | list[ColumnEntry | AllValues]


def _unique(group: str, entries: Entries) -> set[str]:
    names = [entry.name for entry in entries]
    if len(names) != len(set(names)):
        raise ValueError(f"{group} entry names must be unique")
    return set(names)


def _check_axis_flags(group: str, entries: Entries, keys: list[str]) -> None:
    for entry in entries:
        for flag in ("index_only", "named_index"):
            if getattr(entry, flag, False) and entry.name not in keys:
                raise ValueError(f"{group} entry {entry.name!r} is {flag} but no axis key")


def _check_uns_only(group: str, entries: Entries) -> None:
    for entry in entries:
        if isinstance(entry, AnnotationMap) or (
            isinstance(entry, Constant) and isinstance(entry.value, dict | list)
        ):
            raise ValueError(f"{group} entry {entry.name!r} can only be written to uns")


class LevelRule(ModelBase):
    """One level's export, base and level merged: apb2's effective rule, reversed."""

    modality: str | None = Field(default=None, description="The MuData modality it fills.")
    required: bool = Field(default=True, description="When false, a result without it is fine.")
    software: list[str] | None = Field(
        default=None,
        min_length=1,
        description="The software, as apb2's rule names it, whose results this level exports; "
        "every software when omitted. Another software's result skips the level.",
    )
    axis: Axis
    measurements: Measurements
    columns: Columns
    uns: list[ColumnEntry] = Field(default_factory=list)
    accession_syntax: dict[str, AccessionSyntax] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _references_are_declared(self) -> Self:
        obs = _unique("columns.obs", self.columns.obs)
        var = _unique("columns.var", self.columns.var)
        _unique("uns", self.uns)
        if self.axis.obs_keys[0] not in obs:
            raise ValueError(f"axis.obs_keys={self.axis.obs_keys} matches no columns.obs entry")
        if self.axis.var_keys[0] not in var:
            raise ValueError(f"axis.var_keys={self.axis.var_keys} matches no columns.var entry")
        _check_axis_flags("columns.obs", self.columns.obs, self.axis.obs_keys)
        _check_axis_flags("columns.var", self.columns.var, self.axis.var_keys)
        _check_axis_flags("uns", self.uns, [])
        _check_uns_only("columns.obs", self.columns.obs)
        _check_uns_only("columns.var", self.columns.var)
        for entry in (*self.columns.obs, *self.columns.var, *self.uns):
            syntax = getattr(entry, "syntax", None)
            if syntax is not None and syntax not in self.accession_syntax:
                raise ValueError(f"entry {entry.name!r} names no accession_syntax {syntax!r}")
        return self


class Output(ModelBase):
    """What one table writes, like apb2's ``input``.

    ``wide`` keeps APB2's samples-by-features matrices; ``long`` writes one row per measured
    cell. ``.h5ad`` holds one level: the first listed that the result's software and levels
    allow. ``.h5mu`` holds one modality per level.
    """

    shape: Literal["long", "wide"]
    extensions: list[Literal[".h5mu", ".h5ad"]] = Field(min_length=1, max_length=1)


class Table(ModelBase):
    """One output and the levels it reads; each level's block overrides ``base``."""

    output: Output
    base: dict[str, JsonValue] = Field(default_factory=dict)
    levels: dict[str, dict[str, JsonValue]] = Field(min_length=1)


class ExportRuleDocument(ModelBase):
    """One target's export rule file."""

    schema_version: SchemaVersion
    file_version: str = Field(min_length=1)
    target_name: str = Field(min_length=1)
    target_version_pattern: str = Field(
        min_length=1, description="Target releases the rule was checked against; recorded only."
    )
    tables: list[Table] = Field(min_length=1)
