"""``apb-export``: vendor output in, the AnnData or MuData a downstream tool opens.

One subcommand per target, each reading the APB2 levels its export rule names: APB2 converts
the vendor files, apb2's sample annotation attaches ``--annotation``, apb-fasta checks the
peptides and protein groups against an optional ``--fasta``, and
:class:`~apb_export.api.Exporter` builds the file. Existing files are never overwritten.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Annotated, Any, cast

import anndata as ad
import mudata as md
from apb2.api import (
    AnnotationCompiler,
    ParsedLevels,
    ParseRuleCompiler,
    QuantificationLevel,
    UnsJsonCodec,
    write_container_representation,
)
from apb_catalog.api import UnresolvedField
from apb_fasta.api import FastaAnnotator
from cyclopts import App, Parameter
from loguru import logger

from apb_export.api import Exporter

app = App(
    name="apb-export",
    help="Convert vendor output into the AnnData or MuData a downstream tool opens.",
)

Data = Annotated[Path, Parameter(help="Vendor table or vendor-result directory")]
Params = Annotated[Path | None, Parameter(help="Vendor parameter file")]
Software = Annotated[
    str | None, Parameter(help="Restrict recognition to one vendor, as in apb2 convert")
]
Abundance = Annotated[
    str | None,
    Parameter(help="Layer of the target's level to become X; APB2's primary layer by default"),
]
SampleAnnotation = Annotated[
    Path | None,
    Parameter(help="SDRF or prolfquapp-style CSV/TSV sample annotation, as apb2 annotate reads"),
]
Strict = Annotated[
    bool, Parameter(negative=False, help="Promote APB2 layer-contract warnings to errors")
]
Fasta = Annotated[
    tuple[Path, ...],
    Parameter(
        negative=False,
        help="Optional FASTA files, or their protein-fasta database: check the peptides and "
        "annotate the protein groups first, and log how many the FASTA holds; msmu also flags "
        "the FASTA's contaminants",
    ),
]


def _parse(
    data: Path,
    params: Path | None,
    levels: tuple[QuantificationLevel, ...],
    *,
    software: str | None,
    strict: bool,
) -> ParsedLevels:
    """Read the vendor files into the levels the target reads, those the vendor has.

    Without a parameter file the software's columns alone choose the rule, as ``apb2
    convert --software`` does.
    """
    checks = "strict" if strict else "standard"
    compiler = (
        ParseRuleCompiler.from_software(data, software, levels, checks)
        if params is None and software is not None
        else ParseRuleCompiler(data, params, levels, checks, software)
    )
    return compiler.compile().parse()


def _exported(
    target: str,
    data: Path,
    output: Path,
    *,
    params: Path | None,
    software: str | None,
    abundance: str | None,
    annotation: Path | None,
    strict: bool,
    fasta: tuple[Path, ...] = (),
) -> ad.AnnData | md.MuData | None:
    """Convert, annotate, FASTA-check and export; ``None`` after logging why not."""
    exporter = Exporter(target, abundance)
    if output.suffix != exporter.extension:
        logger.error(f"{target} opens {exporter.extension} files, not {output}")
        return None
    if output.exists():
        logger.error(f"refusing to overwrite {output}")
        return None
    try:
        parsed = _parse(data, params, exporter.levels, software=software, strict=strict)
        if annotation is not None:
            parsed = AnnotationCompiler().compile(annotation).parse(parsed).annotate().parsed
        if fasta:
            parsed = _fasta_checked(parsed, fasta)
        return exporter.export(parsed)
    except (OSError, ValueError) as error:
        logger.error(str(error))
    except UnresolvedField as error:
        # Unknown is not missing: the catalogue has not reviewed this exact vendor rule, so it
        # cannot say whether a q-value or PEP exists. Guessing would export a silent gap.
        logger.error(
            f"apb-catalog cannot answer for this vendor rule ({error}); re-review the "
            "catalogues against this APB2 revision"
        )
    return None


def _fasta_checked(parsed: ParsedLevels, fasta: tuple[Path, ...]) -> ParsedLevels:
    """Verify the peptides and annotate the protein groups the result holds against the FASTA,
    and log the coverage; protein groups need apb2's ``fasta_accessions`` role."""
    annotator = FastaAnnotator.read(fasta)
    if set(parsed.levels) - {"protein"}:
        verified = annotator.verify_peptides(parsed)
        for level, coverage in verified.reports.peptide_levels.items():
            logger.info(
                f"fasta level={level} peptides_in_fasta={coverage.matched_feature_count}/"
                f"{coverage.feature_count - coverage.decoy_feature_count} "
                f"unmatched={coverage.unmatched_feature_count} "
                f"decoys={coverage.decoy_feature_count}"
            )
        parsed = verified.parsed
    protein = parsed.levels.get("protein")
    if protein is not None and "fasta_accessions" in protein.var.roles:
        annotated = annotator.merge_annotations(parsed)
        groups = annotated.reports.protein_groups
        if groups is not None:
            logger.info(
                f"fasta level=protein members_in_fasta={groups.matched_member_count}/"
                f"{groups.member_count} unmatched={groups.unmatched_member_count} "
                f"ambiguous={groups.ambiguous_member_count}"
            )
        parsed = annotated.parsed
    return parsed


def _summary(adata: ad.AnnData) -> str:
    """Its shape and the APB layer behind each written matrix."""
    sources = _level_part(adata)["export"]["provenance"]["sources"]
    layers = [
        f"{key.split('.', 1)[1]}={value.get('name', value.get('location'))}"
        for key, value in sources.items()
        if key.startswith("layers.")
    ]
    return f"{adata.n_obs} x {adata.n_vars} ({', '.join(layers)})"


def _level_part(adata: ad.AnnData) -> dict[str, Any]:
    """One exported AnnData's own APB part: ``uns[<level>]["apb"]`` standalone, else ``uns["apb"]``."""
    parts = [
        value["apb"] for value in adata.uns.values() if isinstance(value, dict) and "apb" in value
    ]
    namespace = parts[0] if parts else adata.uns["apb"]
    codec = UnsJsonCodec()
    return cast(dict[str, Any], codec.decode(namespace, codec.storage(namespace)))


def _write(result: ad.AnnData | md.MuData | None, output: Path) -> int:
    if result is None:
        return 1
    output.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(result, md.MuData):
        result.write_h5mu(output)
        shapes = "; ".join(
            f"{name} {_summary(modality)}"
            for name, modality in result.mod.items()
            if isinstance(modality, ad.AnnData)
        )
    else:
        result.write_h5ad(output)
        shapes = _summary(result)
    write_container_representation(result, output)
    logger.info(f"wrote {output} and its sidecar: {shapes}")
    return 0


@app.command
def msmu(
    data: Data,
    output: Annotated[Path, Parameter(help="New .h5mu file for msmu.read_h5mu")],
    /,
    *,
    params: Params = None,
    software: Software = None,
    abundance: Abundance = None,
    strict: Strict = False,
    fasta: Fasta = (),
) -> int:
    """APB2's ion level as msmu's psm modality: one feature per run and precursor."""
    result = _exported(
        "msmu",
        data,
        output,
        params=params,
        software=software,
        abundance=abundance,
        annotation=None,
        strict=strict,
        fasta=fasta,
    )
    if isinstance(result, md.MuData):
        # An APB2 result does not record its source files; this command read them.
        result["psm"].uns["identification_file"] = str(data)
    return _write(result, output)


@app.command
def prolfqua(
    data: Data,
    output: Annotated[Path, Parameter(help="New .h5ad file for LFQData_from_anndata")],
    /,
    *,
    params: Params = None,
    software: Software = None,
    abundance: Abundance = None,
    annotation: SampleAnnotation = None,
    strict: Strict = False,
    fasta: Fasta = (),
) -> int:
    """APB2's ion level for prolfquapp; the annotation's columns become its factors."""
    return _write(
        _exported(
            "prolfqua",
            data,
            output,
            params=params,
            software=software,
            abundance=abundance,
            annotation=annotation,
            strict=strict,
            fasta=fasta,
        ),
        output,
    )


@app.command
def proteopy(
    data: Data,
    output: Annotated[Path, Parameter(help="New .h5ad file for ProteoPy")],
    /,
    *,
    params: Params = None,
    software: Software = None,
    abundance: Abundance = None,
    annotation: SampleAnnotation = None,
    strict: Strict = False,
    fasta: Fasta = (),
) -> int:
    """APB2's protein level for ProteoPy, with sample_id and protein_id."""
    return _write(
        _exported(
            "proteopy",
            data,
            output,
            params=params,
            software=software,
            abundance=abundance,
            annotation=annotation,
            strict=strict,
            fasta=fasta,
        ),
        output,
    )


@app.command
def alphapepttools(
    data: Data,
    output: Annotated[Path, Parameter(help="New .h5mu file for mulink_from_anndatas")],
    /,
    *,
    params: Params = None,
    software: Software = None,
    annotation: SampleAnnotation = None,
    strict: Strict = False,
    fasta: Fasta = (),
) -> int:
    """APB2's ion, peptide and protein levels as AlphaPeptTools modalities, each primary layer."""
    return _write(
        _exported(
            "alphapepttools",
            data,
            output,
            params=params,
            software=software,
            abundance=None,
            annotation=annotation,
            strict=strict,
            fasta=fasta,
        ),
        output,
    )


def main() -> None:
    """Run the command line and exit with its status."""
    sys.exit(app(sys.argv[1:]))
