# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript virome_ggtree_ref.R <tree_file> <node_file> <output_file> [ref_taxon_csv]")
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grep("^--file=", cmd_args)][1]
script_path <- normalizePath(sub("^--file=", "", file_arg))
script_dir <- dirname(script_path)
options(
  virome_script_dir = script_dir,
  virome_project_dir = normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
)
source(file.path(script_dir, "virome_tree_plot_common.R"))

ref_taxon_csv <- if (length(args) >= 4) args[4] else NULL
midpoint_rooting <- length(args) >= 5 && toupper(args[5]) == "TRUE"
render_virome_tree(args[1], args[2], args[3], label_mode = "ref", ref_taxon_csv = ref_taxon_csv, midpoint_rooting = midpoint_rooting)
