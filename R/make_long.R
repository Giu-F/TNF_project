make_long <- function(data) {

  long <- reshape2::melt(
    data,
    id.vars = c(
      "PG.Genes",
      "PG.ProteinDescriptions",
      "PG.ProteinNames",
      "PTM.ProteinId",
      "PTM.CollapseKey",
      "PTM.Multiplicity",
      "PTM.SiteAA",
      "PTM.SiteLocation",
      "PTM.FlankingRegion"
    ),
    variable.name = "Run",
    value.name = "quantity"
  )

  long <- long[!is.na(long$quantity), ]
  long$Run <- as.character(long$Run)

  long$genotype <- ifelse(
    grepl("WT", long$Run), "WT",
    ifelse(grepl("KO", long$Run), "KO", "ND")
  )

  long$treatment <- ifelse(grepl("Ctrl", long$Run), "Ctrl", "TNF")
  long$tp <- ifelse(grepl("15min", long$Run), "15min", "5min")

  long$rep <- ifelse(
    grepl("d1", long$Run), "d1",
    ifelse(grepl("d2", long$Run), "d2", "d3")
  )

  long$genotype <- factor(long$genotype, levels = c("KO", "WT", "ND"))
  long$treatment <- factor(long$treatment, levels = c("Ctrl", "TNF"))
  long$tp <- factor(long$tp, levels = c("5min", "15min"))
  long$rep <- factor(long$rep, levels = c("d1", "d2", "d3"))

  long <- long[
    order(long$genotype, long$treatment, long$tp, long$rep),
  ]

  long$Run <- factor(long$Run, levels = unique(long$Run))

  return(long)
}
