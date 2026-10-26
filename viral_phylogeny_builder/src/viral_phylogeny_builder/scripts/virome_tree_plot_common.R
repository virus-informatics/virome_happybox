# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright 2026 Junna Kawasaki & Jumpei Ito

load_virome_plot_packages <- function(midpoint_rooting = FALSE) {
  packages <- c(
    "ape",
    "ggtree",
    "dplyr",
    "ggplot2",
    "ggnewscale",
    "RColorBrewer"
  )
  if (midpoint_rooting) {
    packages <- c(packages, "phangorn")
  }

  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(paste("Missing R packages:", paste(missing_packages, collapse = ", ")))
  }

  suppressPackageStartupMessages({
    invisible(lapply(packages, library, character.only = TRUE))
  })
}

normalise_label <- function(x) {
  trimws(gsub("\\.[0-9]+$", "", as.character(x)))
}

normalize_plot_na <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "NaN", "Nan", "nan")] <- NA_character_
  x
}

make_palette <- function(levels, palette_name, max_base) {
  n <- length(levels)
  if (n == 0) {
    return(character(0))
  }
  if (n <= max_base) {
    colors <- RColorBrewer::brewer.pal(max(3, n), palette_name)[seq_len(n)]
  } else {
    colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(max_base, palette_name))(n)
  }
  stats::setNames(colors, levels)
}

legend_output_file <- function(output_file) {
  if (grepl("_final\\.pdf$", output_file)) {
    return(sub("_final\\.pdf$", "_legend.pdf", output_file))
  }
  sub("\\.pdf$", "_legend.pdf", output_file)
}

write_pdf_plot <- function(plot, output_file, width_cm, height_cm) {
  grDevices::pdf(
    file = output_file,
    width = width_cm / 2.54,
    height = height_cm / 2.54,
    family = "Helvetica",
    version = "1.4",
    useDingbats = FALSE
  )
  print(plot)
  invisible(grDevices::dev.off())
}

write_pdf_legend <- function(plot, output_file, width_cm = 18, height_cm = 14) {
  legend_file <- legend_output_file(output_file)
  grDevices::pdf(
    file = legend_file,
    width = width_cm / 2.54,
    height = height_cm / 2.54,
    family = "Helvetica",
    version = "1.4",
    useDingbats = FALSE
  )
  guide_box <- ggplot2::ggplotGrob(plot)
  guide_index <- which(vapply(guide_box$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(guide_index) == 0) {
    return(invisible(grDevices::dev.off()))
  }

  grid::grid.draw(guide_box$grobs[[guide_index[1]]])
  invisible(grDevices::dev.off())
}

read_node_data <- function(node_file, label_mode) {
  node_data <- read.csv(node_file, header = TRUE, stringsAsFactors = FALSE, sep = ",")
  required_columns <- c("Accession", "Genus", "manual_host", "Organism_Name", "GenBank_RefSeq")
  for (column in required_columns) {
    if (!column %in% colnames(node_data)) {
      node_data[[column]] <- ""
    }
  }

  if (label_mode == "add" && "node" %in% colnames(node_data)) {
    node_data$label <- normalise_label(node_data$node)
  } else {
    node_data$label <- normalise_label(node_data$Accession)
  }
  node_data$Accession <- normalise_label(node_data$Accession)

  node_data %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(
      Accession = dplyr::first(Accession),
      Organism_Name = dplyr::first(Organism_Name),
      Genus = paste(unique(Genus), collapse = ", "),
      manual_host = paste(unique(manual_host), collapse = ", "),
      GenBank_RefSeq = paste(unique(GenBank_RefSeq), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::distinct(label, .keep_all = TRUE)
}

infer_family <- function(tree_file, node_data) {
  if ("VirusFamily" %in% colnames(node_data)) {
    families <- unique(node_data$VirusFamily[!is.na(node_data$VirusFamily) & node_data$VirusFamily != ""])
    if (length(families) > 0) {
      return(families[[1]])
    }
  }
  sub("\\..*$", "", basename(tree_file))
}

project_dir <- function() {
  configured <- Sys.getenv("VIROME_PROJECT_DIR", unset = "")
  if (configured != "") {
    return(normalizePath(configured, mustWork = FALSE))
  }
  script_dir <- getOption("virome_script_dir", default = "")
  if (script_dir != "") {
    return(normalizePath(file.path(script_dir, ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

read_reference_csv <- function(family, filename_suffix) {
  root <- project_dir()
  candidates <- c(
    file.path(root, "reference", family, paste0(family, filename_suffix))
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    return(NULL)
  }
  read.csv(existing[[1]], header = TRUE, stringsAsFactors = FALSE)
}

read_reference_taxon <- function(ref_taxon_csv = NULL) {
  if (is.null(ref_taxon_csv) || !file.exists(ref_taxon_csv)) {
    return(NULL)
  }
  ref <- read.csv(ref_taxon_csv, header = TRUE, stringsAsFactors = FALSE)
  for (column in c("Accession", "Genus", "manual_host", "Organism_Name", "GenBank_RefSeq")) {
    if (!column %in% colnames(ref)) {
      ref[[column]] <- NA_character_
    }
  }
  ref %>%
    dplyr::mutate(label = normalise_label(Accession)) %>%
    dplyr::select(
      label,
      ref_Genus = Genus,
      ref_manual_host = manual_host,
      ref_Organism_Name = Organism_Name,
      ref_GenBank_RefSeq = GenBank_RefSeq
    ) %>%
    dplyr::distinct(label, .keep_all = TRUE) %>%
    dplyr::mutate(ref_label_exists = TRUE)
}

read_species_labels <- function(family) {
  ictv <- read_reference_csv(family, ".ictv.csv")
  if (is.null(ictv) || !"Accession" %in% colnames(ictv) || !"Species" %in% colnames(ictv)) {
    return(NULL)
  }
  ictv %>%
    dplyr::mutate(label = normalise_label(Accession)) %>%
    dplyr::select(label, Species_Label = Species) %>%
    dplyr::filter(!is.na(Species_Label) & Species_Label != "") %>%
    dplyr::distinct(label, .keep_all = TRUE)
}

fill_blank_with_ref <- function(plot_data, ref_taxon_csv = NULL) {
  ref_taxon <- read_reference_taxon(ref_taxon_csv)
  if (is.null(ref_taxon)) {
    return(plot_data)
  }

  plot_data %>%
    dplyr::left_join(ref_taxon, by = "label") %>%
    dplyr::mutate(
      has_ref_metadata = !is.na(ref_label_exists),
      Genus = dplyr::if_else(has_ref_metadata & !is.na(ref_Genus) & ref_Genus != "", ref_Genus, Genus),
      manual_host = dplyr::if_else(has_ref_metadata & !is.na(ref_manual_host) & ref_manual_host != "", ref_manual_host, manual_host),
      Organism_Name = dplyr::if_else(has_ref_metadata & !is.na(ref_Organism_Name) & ref_Organism_Name != "", ref_Organism_Name, Organism_Name),
      GenBank_RefSeq = dplyr::if_else(has_ref_metadata & !is.na(ref_GenBank_RefSeq) & ref_GenBank_RefSeq != "", ref_GenBank_RefSeq, GenBank_RefSeq)
    ) %>%
    dplyr::select(-ref_Genus, -ref_manual_host, -ref_Organism_Name, -ref_GenBank_RefSeq, -ref_label_exists)
}

add_species_labels <- function(plot_data, family) {
  species_labels <- read_species_labels(family)
  if (is.null(species_labels)) {
    plot_data$Species_Label <- NA_character_
    return(plot_data)
  }
  plot_data %>% dplyr::left_join(species_labels, by = "label")
}

add_display_labels <- function(plot_data) {
  if (!"has_ref_metadata" %in% colnames(plot_data)) {
    plot_data$has_ref_metadata <- FALSE
  }

  plot_data %>%
    dplyr::mutate(
      is_added = !is.na(GenBank_RefSeq) & GenBank_RefSeq == "Input",
      legacy_add_reference = !has_ref_metadata & !is_added & GenBank_RefSeq == "add",
      uploaded_label = substr(dplyr::if_else(!is.na(Accession) & Accession != "", as.character(Accession), as.character(label)), 1, 10),
      display_label = dplyr::case_when(
        is_added ~ uploaded_label,
        legacy_add_reference ~ as.character(label),
        !is.na(Species_Label) & Species_Label != "" ~ Species_Label,
        !is.na(Organism_Name) & Organism_Name != "" ~ as.character(Organism_Name),
        TRUE ~ as.character(label)
      ),
      point_type = dplyr::if_else(is_added, "added", "reference")
    )
}

render_virome_tree <- function(tree_file, node_file, output_file, label_mode = "ref", ref_taxon_csv = NULL, midpoint_rooting = FALSE) {
  load_virome_plot_packages(midpoint_rooting)

  tree <- ape::read.tree(tree_file)
  tree$tip.label <- normalise_label(tree$tip.label)

  if (midpoint_rooting) {
    tree <- phangorn::midpoint(tree)
  }

  node_data <- read_node_data(node_file, label_mode)
  family <- infer_family(tree_file, node_data)

  tree_data <- suppressWarnings(ggtree::fortify(tree))
  tree_data$label <- normalise_label(tree_data$label)

  plot_data <- tree_data %>%
    dplyr::left_join(node_data, by = "label") %>%
    fill_blank_with_ref(ref_taxon_csv) %>%
    add_species_labels(family) %>%
    add_display_labels()

  tip_data <- plot_data %>%
    dplyr::filter(isTip) %>%
    dplyr::mutate(
      Genus = normalize_plot_na(Genus),
      manual_host = normalize_plot_na(manual_host)
    )

  tip_data$Genus <- as.factor(tip_data$Genus)
  tip_data$manual_host <- as.factor(tip_data$manual_host)

  palette_colors_genus <- make_palette(levels(tip_data$Genus), "Set3", 12)
  palette_colors_host <- make_palette(levels(tip_data$manual_host), "Spectral", 11)

  legend_levels <- c(levels(tip_data$Genus), levels(tip_data$manual_host), "reference", "added")
  legend_levels <- legend_levels[!is.na(legend_levels) & legend_levels != ""]

  legend_width_cm <- 18 #max(18, min(36, 7 + max(nchar(legend_levels), na.rm = TRUE) * 0.18))
  legend_height_cm <- 20 #max(10, min(60, 4 + length(legend_levels) * 0.35))

  max_x <- max(plot_data$x, na.rm = TRUE)
  max_y <- max(plot_data$y, na.rm = TRUE)

  node_height_cm <- as.numeric(Sys.getenv("VIROME_NODE_HEIGHT_CM", unset = "0.2"))
  if (is.na(node_height_cm) || node_height_cm <= 0) {
    node_height_cm <- 0.2
  }

  label_size_pt <- max(8, as.numeric(Sys.getenv("VIROME_TIP_LABEL_PT", unset = "8")))
  if (is.na(label_size_pt)) {
    label_size_pt <- 8
  }

  label_size <- label_size_pt / ggplot2::.pt
  point_size <- max(1.8, label_size * 0.75)

  n_tips <- nrow(tip_data)
  output_height_cm <- max(12, n_tips * node_height_cm + 3)

  output_width_cm <- as.numeric(Sys.getenv("VIROME_TREE_WIDTH_CM", unset = "30"))
  if (is.na(output_width_cm) || output_width_cm <= 0) {
    output_width_cm <- 30
  }

  max_label_chars <- max(nchar(as.character(tip_data$display_label)), na.rm = TRUE)
  if (is.na(max_label_chars)) {
    max_label_chars <- 0
  }

  auto_width_cm <- 30 + max(0, max_label_chars - 30) * 0.18
  output_width_cm <- max(output_width_cm, min(52, auto_width_cm))

  tree_linewidth <- as.numeric(Sys.getenv("VIROME_TREE_LINEWIDTH", unset = "0.9"))
  if (is.na(tree_linewidth) || tree_linewidth <= 0) {
    tree_linewidth <- 0.9
  }

  tip_guide_linewidth <- as.numeric(Sys.getenv("VIROME_TIP_GUIDE_LINEWIDTH", unset = "0.08"))
  if (is.na(tip_guide_linewidth) || tip_guide_linewidth < 0) {
    tip_guide_linewidth <- 0.08
  }

  label_x <- max_x + 0.12
  label_space <- max(3.2, min(8.5, max_label_chars * 0.08))

  tile_width <- 0.45
  tile_height <- 0.82

  heatmap_start <- label_x + label_space + 0.35
  heatmap_x <- c(
    Genus = heatmap_start,
    manual_host = heatmap_start + 0.55
  )

  column_labels <- data.frame(
    x = unname(heatmap_x),
    y = max_y + 2.5,
    label = names(heatmap_x)
  )

  p <- ggtree::ggtree(
    tree,
    layout = "rectangular",
    size = tree_linewidth,
    color = "black"
  ) +
    ggplot2::geom_segment(
      data = tip_data,
      ggplot2::aes(
        x = x,
        xend = label_x,
        y = y,
        yend = y
      ),
      linewidth = tip_guide_linewidth,
      color = "gray88",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = tip_data,
      ggplot2::aes(
        x = label_x,
        y = y,
        label = display_label
      ),
      hjust = 0,
      size = label_size,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = tip_data,
      ggplot2::aes(
        x = x,
        y = y,
        color = point_type
      ),
      size = point_size,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        reference = "gray60",
        added = "orange"
      ),
      breaks = c(
        "reference",
        "added"
      ),
      name = "tip"
    ) +
    ggplot2::geom_tile(
      data = tip_data,
      ggplot2::aes(
        x = heatmap_x[["Genus"]],
        y = y,
        fill = Genus
      ),
      width = tile_width,
      height = tile_height,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = palette_colors_genus,
      na.value = "gray80"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(title = "Genus", ncol = 4)
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      data = tip_data,
      ggplot2::aes(
        x = heatmap_x[["manual_host"]],
        y = y,
        fill = manual_host
      ),
      width = tile_width,
      height = tile_height,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = palette_colors_host,
      na.value = "gray80"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(title = "manual_host")
    ) +
    ggplot2::geom_text(
      data = column_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      angle = 90,
      hjust = 0,
      size = 3,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = 0.5,
        y = -1,
        yend = -1
      ),
      linewidth = tree_linewidth,
      lineend = "round",
      color = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = 0.25,
      y = -2.2,
      label = "0.5",
      size = 3
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, max(heatmap_x) + 0.75),
      ylim = c(-3, max_y + 5),
      clip = "off"
    ) +
    ggplot2::labs(
      x = "Branch Length",
      y = NULL
    ) +
    ggtree::theme_tree2(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(5.5, 16, 5.5, 5.5),
      text = ggplot2::element_text(family = "Helvetica")
    )

  write_pdf_legend(p, output_file, legend_width_cm, legend_height_cm)

  write_pdf_plot(
    p + ggplot2::theme(legend.position = "none"),
    output_file,
    output_width_cm,
    output_height_cm
  )
}
