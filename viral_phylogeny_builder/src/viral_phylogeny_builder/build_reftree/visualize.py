# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

"""Build phylogenetic tree with IQ-TREE and render visualization."""

from __future__ import annotations

import importlib.resources
import logging
from pathlib import Path

from viral_phylogeny_builder._commands import (
    BUILD_REFTREE_COMMAND,
    BUILD_REFTREE_VISUALIZE_SUBCOMMAND,
)
from viral_phylogeny_builder._pipeline import (
    mark_step_done,
    run_subprocess,
    should_skip,
    staged_output,
    step_done_marker,
    step_is_done,
    validate_input_file,
    validate_step_output,
    validate_virus_family,
)
from viral_phylogeny_builder._tool_config import resolve_tools
from viral_phylogeny_builder.build_reftree._tree_metadata import prepare_tree_metadata

logger = logging.getLogger(__name__)

_IQTREE_REQUIRED_TOOLS = ("iqtree",)
_RSCRIPT_REQUIRED_TOOLS = ("Rscript",)
_VIROME_TREE_R_SCRIPT = "virome_ggtree_ref.R"

_IQTREE_SUBSTITUTION_MODEL = "LG+G"
_IQTREE_PREFIX_SUFFIX = ".iqtree"

BUILD_REFTREE_VISUALIZE_TREEFILE_FILENAME_SUFFIX = ".iqtree.treefile"
BUILD_REFTREE_VISUALIZE_TREE_METADATA_FILENAME_SUFFIX = ".tree_metadata.csv"
BUILD_REFTREE_VISUALIZE_TREE_PDF_FILENAME_SUFFIX = ".tree.pdf"
BUILD_REFTREE_VISUALIZE_TREE_LEGEND_PDF_FILENAME_SUFFIX = ".tree_legend.pdf"


def build_reftree_visualize_treefile_filename(virus_family: str) -> str:
    return f"{virus_family}{BUILD_REFTREE_VISUALIZE_TREEFILE_FILENAME_SUFFIX}"


def build_reftree_visualize_tree_metadata_filename(virus_family: str) -> str:
    return f"{virus_family}{BUILD_REFTREE_VISUALIZE_TREE_METADATA_FILENAME_SUFFIX}"


def build_reftree_visualize_tree_pdf_filename(virus_family: str) -> str:
    return f"{virus_family}{BUILD_REFTREE_VISUALIZE_TREE_PDF_FILENAME_SUFFIX}"


def build_reftree_visualize_tree_legend_pdf_filename(virus_family: str) -> str:
    return f"{virus_family}{BUILD_REFTREE_VISUALIZE_TREE_LEGEND_PDF_FILENAME_SUFFIX}"


def _iqtree_prefix(workdir: Path, virus_family: str) -> Path:
    return workdir / f"{virus_family}{_IQTREE_PREFIX_SUFFIX}"


def _resolve_r_script() -> Path:
    return Path(
        str(importlib.resources.files("viral_phylogeny_builder.scripts") / _VIROME_TREE_R_SCRIPT),
    )


def _run_rscript(
    treefile: Path,
    tree_metadata_csv: Path,
    output_pdf: Path,
    *,
    midpoint_rooting: bool = False,
    rscript_cmd: str = "Rscript",
) -> None:
    r_script = _resolve_r_script()
    run_subprocess(
        [
            rscript_cmd,
            str(r_script),
            str(treefile),
            str(tree_metadata_csv),
            str(output_pdf),
            "",
            str(midpoint_rooting).upper(),
        ],
        "Rscript virome_ggtree_ref.R",
    )


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
            "-fast",
            "-T",
            str(threads),
            "--prefix",
            str(prefix),
            "--redo",
        ],
        "iqtree",
    )


def run_build_reftree_visualize(
    virus_family: str,
    alignment_fasta: Path,
    rdrp_metadata_csv: Path,
    output_dir: Path,
    *,
    threads: int = 4,
    force: bool = False,
    midpoint_rooting: bool = False,
    tool_config: dict[str, str] | None = None,
) -> Path:
    """Run build-reftree visualize pipeline and return path to treefile."""
    resolved = resolve_tools(_IQTREE_REQUIRED_TOOLS + _RSCRIPT_REQUIRED_TOOLS, tool_config)
    validate_virus_family(virus_family)
    validate_input_file(alignment_fasta, "Trimmed alignment FASTA")
    validate_input_file(rdrp_metadata_csv, "RdRp metadata CSV")

    workdir = output_dir / BUILD_REFTREE_COMMAND / BUILD_REFTREE_VISUALIZE_SUBCOMMAND / virus_family
    workdir.mkdir(parents=True, exist_ok=True)

    treefile = workdir / build_reftree_visualize_treefile_filename(virus_family)
    tree_metadata_csv = workdir / build_reftree_visualize_tree_metadata_filename(virus_family)
    tree_pdf = workdir / build_reftree_visualize_tree_pdf_filename(virus_family)
    tree_legend_pdf = workdir / build_reftree_visualize_tree_legend_pdf_filename(virus_family)
    marker = step_done_marker(
        workdir,
        BUILD_REFTREE_COMMAND,
        BUILD_REFTREE_VISUALIZE_SUBCOMMAND,
    )
    step_outputs = [treefile, tree_metadata_csv, tree_pdf, tree_legend_pdf]

    if not force and step_is_done(marker, step_outputs):
        logger.info("Skipping %s (cached): %s", BUILD_REFTREE_VISUALIZE_SUBCOMMAND, workdir)
        return treefile

    marker.unlink(missing_ok=True)
    if force:
        logger.info("Force flag set; removing cached outputs in %s", workdir)
        for output in step_outputs:
            output.unlink(missing_ok=True)

    if not should_skip(treefile, force=force, step_name="iqtree"):
        prefix = _iqtree_prefix(workdir, virus_family)
        _run_iqtree(alignment_fasta, prefix, threads=threads, iqtree_cmd=resolved["iqtree"])
        validate_step_output(treefile, "IQ-TREE treefile")

    if not should_skip(tree_metadata_csv, force=force, step_name="tree metadata"):
        with staged_output(tree_metadata_csv) as staging_path:
            prepare_tree_metadata(alignment_fasta, rdrp_metadata_csv, staging_path)
            validate_step_output(staging_path, "tree metadata CSV")

    if not should_skip(tree_pdf, force=force, step_name="Rscript virome_ggtree_ref.R"):
        _run_rscript(
            treefile,
            tree_metadata_csv,
            tree_pdf,
            midpoint_rooting=midpoint_rooting,
            rscript_cmd=resolved["Rscript"],
        )
        validate_step_output(tree_pdf, "tree PDF")

    mark_step_done(marker, step_outputs)
    logger.info("Build-reftree visualize complete: %s", treefile)
    return treefile
