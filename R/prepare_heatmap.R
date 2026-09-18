prepare_heatmap <- function(long_data, regulated_sites, time_point) {

  heatmap_data <- long_data %>%
    filter(
      PTM.CollapseKey %in% regulated_sites,
      tp == time_point
    ) %>%
    mutate(
      id = paste0(
        PG.Genes, "_",
        PTM.SiteAA, PTM.SiteLocation, "_",
        PTM.Multiplicity
      )
    ) %>%
    group_by(PTM.CollapseKey, id) %>%
    mutate(
      z = as.numeric(scale(quantity))
    ) %>%
    ungroup() %>%
    group_by(
      PTM.CollapseKey,
      id,
      genotype,
      treatment,
      rep
    ) %>%
    summarise(
      z = mean(z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      condition = paste0(genotype, "_", treatment, "_", rep)
    )

  wide <- data.table::dcast(
    as.data.table(heatmap_data),
    id ~ condition,
    value.var = "z"
  )

  # Order columns by genotype
  column_order <- c(
    grep("^KO_", colnames(wide), value = TRUE),
    grep("^WT_", colnames(wide), value = TRUE),
    grep("^ND_", colnames(wide), value = TRUE)
  )

  # Convert to data.frame before assigning row names
  wide <- as.data.frame(wide)

  rownames(wide) <- wide$id

  wide <- wide[, column_order, drop = FALSE]

  # Column annotations
  annotation <- data.frame(
    Treatment = ifelse(
      grepl("_TNF_", colnames(wide)),
      "TNF",
      "Ctrl"
    ),
    Genotype = case_when(
      grepl("^KO_", colnames(wide)) ~ "KO",
      grepl("^WT_", colnames(wide)) ~ "WT",
      grepl("^ND_", colnames(wide)) ~ "ND"
    ),
    row.names = colnames(wide)
  )

  annotation$Treatment <- factor(
    annotation$Treatment,
    levels = c("Ctrl", "TNF")
  )

  annotation$Genotype <- factor(
    annotation$Genotype,
    levels = c("KO", "WT", "ND")
  )

  list(
    matrix = as.matrix(wide),
    annotation = annotation
  )
}
