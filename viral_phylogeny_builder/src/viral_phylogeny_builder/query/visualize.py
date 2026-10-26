# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Build query tree (IQ-TREE) and render with query + reference metadata."""

from __future__ import annotations

import importlib.resources
import logging
from pathlib import Path

from viral_phylogeny_builder._commands import QUERY_COMMAND, QUERY_VISUALIZE_SUBCOMMAND
from viral_phylogeny_builder._pipeline import (
    run_subprocess,
    should_skip,
    staged_output,
    step_workdir,
    validate_input_file,
    validate_step_output,
    validate_virus_family,
)
from viral_phylogeny_builder._tool_config import resolve_tools
from viral_phylogeny_builder.query._concat_metadata import (
    concat_metadata,
    concat_metadata_filename,
)
from viral_phylogeny_builder.query._query_metadata import (
    prepare_query_metadata,
    query_metadata_filename,
)

logger = logging.getLogger(__name__)

_IQTREE_REQUIRED_TOOLS = ("iqtree",)
_TAXONKIT_REQUIRED_TOOLS = ("taxonkit",)
_RSCRIPT_REQUIRED_TOOLS = ("Rscript",)

_IQTREE_SUBSTITUTION_MODEL = "LG+G"
_IQTREE_PREFIX_SUFFIX = ".iqtree"
_VIROME_TREE_R_SCRIPT = "virome_ggtree_add.R"

QUERY_VISUALIZE_TREEFILE_FILENAME_SUFFIX = ".iqtree.treefile"
QUERY_VISUALIZE_TREE_PDF_FILENAME_SUFFIX = ".tree.pdf"
QUERY_VISUALIZE_TREE_LEGEND_PDF_FILENAME_SUFFIX = ".tree_legend.pdf"


def query_visualize_treefile_filename(virus_family: str) -> str:
    return f"{virus_family}{QUERY_VISUALIZE_TREEFILE_FILENAME_SUFFIX}"


def query_visualize_tree_pdf_filename(virus_family: str) -> str:
    return f"{virus_family}{QUERY_VISUALIZE_TREE_PDF_FILENAME_SUFFIX}"


def query_visualize_tree_legend_pdf_filename(virus_family: str) -> str:
    return f"{virus_family}{QUERY_VISUALIZE_TREE_LEGEND_PDF_FILENAME_SUFFIX}"


def _resolve_r_script() -> Path:
    return Path(
        str(importlib.resources.files("viral_phylogeny_builder.scripts") / _VIROME_TREE_R_SCRIPT),
    )


def _iqtree_prefix(workdir: Path, virus_family: str) -> Path:
    return workdir / f"{virus_family}{_IQTREE_PREFIX_SUFFIX}"


def _run_rscript(
    treefile: Path,
    concat_metadata_csv: Path,
    output_pdf: Path,
    reference_taxon_csv: Path,
    *,
    midpoint_rooting: bool = False,
    rscript_cmd: str = "Rscript",
) -> None:
    r_script = _resolve_r_script()
    cmd = [
        rscript_cmd,
        str(r_script),
        str(treefile),
        str(concat_metadata_csv),
        str(output_pdf),
        str(reference_taxon_csv),
        str(midpoint_rooting).upper(),
    ]
    run_subprocess(cmd, "Rscript virome_ggtree_add.R")


def _run_iqtree(
    alignment_fasta: Path,
    prefix: Path,
    *,
    threads: int,
    iqtree_cmd: str = "iqtree",
) -> None:
    run_subprocess(
        [
            iqtree_cmd,
            "-s",
            str(alignment_fasta),
            "-m",
            _IQTREE_SUBSTITUTION_MODEL,
            "-T",
            str(threads),
            "--prefix",
            str(prefix),
            "--redo",
        ],
        "iqtree",
    )


def run_query_visualize(
    virus_family: str,
    alignment_fasta: Path,
    name_map_tsv: Path,
    host_csv: Path,
    reference_taxon_csv: Path,
    output_dir: Path,
    *,
    threads: int = 4,
    force: bool = False,
    midpoint_rooting: bool = False,
    tool_config: dict[str, str] | None = None,
) -> Path:
    """Run query visualize pipeline and return path to treefile."""
    resolved = resolve_tools(
        _IQTREE_REQUIRED_TOOLS + _TAXONKIT_REQUIRED_TOOLS + _RSCRIPT_REQUIRED_TOOLS,
        tool_config,
    )
    validate_virus_family(virus_family)
    validate_input_file(alignment_fasta, "Alignment FASTA")
    validate_input_file(name_map_tsv, "Name map TSV")
    validate_input_file(host_csv, "Host CSV")
    validate_input_file(reference_taxon_csv, "Reference taxon CSV")

    workdir = step_workdir(output_dir, QUERY_COMMAND, QUERY_VISUALIZE_SUBCOMMAND, virus_family)
    workdir.mkdir(parents=True, exist_ok=True)

    treefile = workdir / query_visualize_treefile_filename(virus_family)

    if not should_skip(treefile, force=force, step_name="iqtree"):
        prefix = _iqtree_prefix(workdir, virus_family)
        _run_iqtree(alignment_fasta, prefix, threads=threads, iqtree_cmd=resolved["iqtree"])
        validate_step_output(treefile, "IQ-TREE treefile")

    query_metadata_csv = workdir / query_metadata_filename(virus_family)

    if not should_skip(query_metadata_csv, force=force, step_name="query metadata"):
        with staged_output(query_metadata_csv) as staging_path:
            prepare_query_metadata(
                name_map_tsv=name_map_tsv,
                alignment_fasta=alignment_fasta,
                host_csv=host_csv,
                virus_family=virus_family,
                output_csv=staging_path,
                output_dir=output_dir,
                taxonkit_cmd=resolved["taxonkit"],
            )
            validate_step_output(staging_path, "query metadata CSV")

    concat_csv = workdir / concat_metadata_filename(virus_family)

    if not should_skip(concat_csv, force=force, step_name="concat metadata"):
        with staged_output(concat_csv) as staging_path:
            concat_metadata(
                ref_csv=reference_taxon_csv,
                add_csv=query_metadata_csv,
                alignment_fasta=alignment_fasta,
                output_csv=staging_path,
            )
            validate_step_output(staging_path, "concat metadata CSV")

    tree_pdf = workdir / query_visualize_tree_pdf_filename(virus_family)
    tree_legend_pdf = workdir / query_visualize_tree_legend_pdf_filename(virus_family)

    if not should_skip(tree_pdf, force=force, step_name="Rscript virome_ggtree_add.R"):
        _run_rscript(
            treefile,
            concat_csv,
            tree_pdf,
            reference_taxon_csv,
            midpoint_rooting=midpoint_rooting,
            rscript_cmd=resolved["Rscript"],
        )
        validate_step_output(tree_pdf, "tree PDF")
        validate_step_output(tree_legend_pdf, "tree legend PDF")

    logger.info("Query visualize complete: %s", tree_pdf)
    return treefile
