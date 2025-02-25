
# Load data ---------------------------------------------------------------
setwd("N:/SUN-CPR-proteomics_jvo/GIULIA/N. COLLAB_TNF/20241128_60spd_TNF_phos/data_analysis_no_norm_no_imp")

library(data.table)
library(tidyverse)
data <- fread("20241126_Report_PTM_pivot_GF (Pivot).tsv")
data <- as.data.frame(data)

# Filter data --------------------------------------------------------------

#Remove contaminants
unique(data$PG.FastaFiles)
data <- data[data$PG.FastaFiles == "human_uniprotkb_reviewed_2024_05_24",]
data <- data %>% select(-PG.FastaFiles)

#Remove rows wo Gene name
data <- data[data$PG.Genes != "",]

#Keep phos
unique(data$PTM.ModificationTitle)
data <- data[data$PTM.ModificationTitle == "Phospho (STY)",]
data <- data %>% select(-PTM.ModificationTitle)

# Rename samples ----------------------------------------------------------

#remove ".raw.PTM.Quantity"
colnames(data) <- gsub(".raw.PTM.Quantity", "", colnames(data))

#remove text before "GF"
colnames(data) <- gsub(".*HCT116_", "", colnames(data))

# Change quantity cols to numerical ---------------------------------------
data[,10:81] <- sapply(data[,10:81], as.numeric)

# Check if there are any values <1 in the entire data frame ----------------
any(data[,10:81] < 1)

# Replace all values <1 with NA -------------------------------------------
data[,10:81][data[,10:81] < 1] <- NA

#Log2 transformation ------------------------------------------------------
data[,10:81] <- log2(data[,10:81])

# Remove duplicates -----------------------------------------------------
#keep rows with unique values in all columns except PG.Genes, PG.ProteinDescriptions, PG.ProteinNames, PTM.ProteinId, PTM.CollapseKey, PTM.FlankingRegion
data <- data %>% distinct_at(vars(-PG.Genes, -PG.ProteinDescriptions, -PG.ProteinNames, -PTM.ProteinId, -PTM.CollapseKey, -PTM.FlankingRegion), .keep_all = TRUE)

# Remove outlier ----------------------------------------------------------
data <- data %>% select(-KO_Ctrl_05min_d1_1)

#export data cleaned
write.table(data, "data_cleaned.txt", sep = "\t", quote = F, row.names = F)

data <- fread("data_cleaned.txt")
data <- as.data.frame(data)

# Filtering --------------------------------------------------------------
#keep rows with more than 50% values present
data_f <- data[rowSums(!is.na(data[,10:80])) > 0.5 * ncol(data[,10:80]),]

col_KO_Ctrl <- grep("KO_Ctrl", colnames(data))
col_KO_TNF <- grep("KO_TNF", colnames(data))
col_WT_Ctrl <- grep("WT_Ctrl", colnames(data))
col_WT_TNF <- grep("WT_TNF", colnames(data))
col_ND_Ctrl <- grep("ND_Ctrl", colnames(data))
col_ND_TNF <- grep("ND_TNF", colnames(data))

#keep rows with at least 4 values in one condition
data_f <- data_f[rowSums(!is.na(data_f[, col_KO_Ctrl])) >= 4 |
                   rowSums(!is.na(data_f[, col_KO_TNF])) >= 4 |
                   rowSums(!is.na(data_f[, col_WT_Ctrl])) >= 4 |
                   rowSums(!is.na(data_f[, col_WT_TNF])) >= 4 |
                   rowSums(!is.na(data_f[, col_ND_Ctrl])) >= 4 |
                   rowSums(!is.na(data_f[, col_ND_TNF])) >= 4, ]

#export data filtered
write.table(data_f, "data_filtered.txt", sep = "\t", quote = F, row.names = F)

data_f <- fread("data_filtered.txt")
data_f <- as.data.frame(data_f)

# Long format with reshape2 after filtering -------------------------------
long <- reshape2::melt(
  data_f,
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
long <- long[!is.na(long$quantity),]
long$Run <- as.character(long$Run)
length(unique(long$Run))
#If Run contains "WT", write "WT", otherwise write "KO", otherwise write "ND"
long$genotype <- ifelse(grepl("WT", long$Run), "WT", ifelse(grepl("KO", long$Run), "KO", "ND"))
unique(long$genotype)
long$genotype <- factor(long$genotype, levels = c("WT", "KO", "ND"))
#if Run contains "Ctrl", write "Ctrl", otherwise write "TNF"
long$treatment <- ifelse(grepl("Ctrl", long$Run), "Ctrl", "TNF")
unique(long$treatment)
long$treatment <- factor(long$treatment, levels = c("Ctrl", "TNF"))
#if Run contains "15min", write "15min", otherwise write "5min"
long$tp <- ifelse(grepl("15min", long$Run), "15min", "5min")
unique(long$tp)
long$tp <- factor(long$tp, levels = c("5min", "15min"))
#if Run contains "d1", write "d1", otherwise write "d2", otherwise write "d3"
long$rep <- ifelse(grepl("d1", long$Run), "d1", ifelse(grepl("d2", long$Run), "d2", "d3"))
unique(long$rep)
long$rep <- factor(long$rep, levels = c("d1", "d2", "d3"))

#order by condition and time point and rep
long <- long[order(long$genotype, long$treatment, long$tp, long$rep),]
long$Run <- factor(long$Run, levels = unique(long$Run))

#boxplot
ggplot(long, aes(x = Run, y = quantity, fill = treatment, color = genotype)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9), alpha = 0.5, lwd = 0.5) +
  scale_color_manual(values = c("black", "darkgrey", "blue")) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Log2 transformed quantity", fill = "Time point", x = NULL)

# Count psites -----------------------------------------------------------

#total
totcount <- long %>%
  group_by(Run, genotype, treatment, tp) %>%
  summarise(Count = n_distinct(PTM.CollapseKey))
write.table(totcount, "id_total_count.txt", sep = "\t", quote = F, row.names = F)

#mean
meancount <- long %>%
  group_by(Run, genotype, treatment, tp) %>%
  summarise(C = n_distinct(PTM.CollapseKey)) %>%
  group_by(genotype, treatment, tp) %>%
  summarise(Count = mean(C, na.rm = T), sd = sd(C, na.rm = T))
write.table(meancount, "id_mean_count.txt", sep = "\t", quote = F, row.names = F)

#plot
ggplot(meancount, aes(x = treatment, y = Count, fill = tp)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  geom_errorbar(aes(ymin = Count - sd, ymax = Count + sd), position = position_dodge(width = 0.9), width = 0.25) +
  geom_text(aes(label = round(Count, 1)), vjust = -3, position = position_dodge(width = 0.9), size=2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Number of localized phosphosites", fill = "Time point") +
  facet_grid(~ genotype) +
  scale_y_continuous(limits = c(0, 14000), breaks = seq(0, 15000, 2000))

# Normalize data after filtering ------------------------------------------
library(DAPAR)
norm <- QuantileCentering(
  data_f[,10:80],
  conds = NULL,
  type = "overall",
  quantile = 0.5
)
norm <- as.data.frame(norm)
norm <- cbind(data_f[,1:9], norm)

#export normalized data
write.table(norm, "data_filtered_norm.txt", sep = "\t", quote = F, row.names = F)
write.table(norm, "data_filtered_norm_2.txt", sep = "\t", quote = F, row.names = F)

norm <- fread("data_filtered_norm.txt")
norm <- as.data.frame(norm)

norm <- fread("data_filtered_norm_2.txt")
norm <- as.data.frame(norm)

# Long format with reshape2 after normalization ----------------------------
long <- reshape2::melt(
  norm,
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
long <- long[!is.na(long$quantity),]
long$Run <- as.character(long$Run)
#If Run contains "WT", write "WT", otherwise write "KO", otherwise write "ND"
long$genotype <- ifelse(grepl("WT", long$Run), "WT", ifelse(grepl("KO", long$Run), "KO", "ND"))
unique(long$genotype)
long$genotype <- factor(long$genotype, levels = c("WT", "KO", "ND"))
#if Run contains "Ctrl", write "Ctrl", otherwise write "TNF"
long$treatment <- ifelse(grepl("Ctrl", long$Run), "Ctrl", "TNF")
unique(long$treatment)
long$treatment <- factor(long$treatment, levels = c("Ctrl", "TNF"))
#if Run contains "15min", write "15min", otherwise write "5min"
long$tp <- ifelse(grepl("15min", long$Run), "15min", "5min")
unique(long$tp)
long$tp <- factor(long$tp, levels = c("5min", "15min"))
#if Run contains "d1", write "d1", otherwise write "d2", otherwise write "d3"
long$rep <- ifelse(grepl("d1", long$Run), "d1", ifelse(grepl("d2", long$Run), "d2", "d3"))
unique(long$rep)
long$rep <- factor(long$rep, levels = c("d1", "d2", "d3"))

#order by condition and time point and rep
long <- long[order(long$genotype, long$treatment, long$tp, long$rep),]
long$Run <- factor(long$Run, levels = unique(long$Run))

#boxplot
ggplot(long, aes(x = Run, y = quantity, fill = treatment, color = genotype)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9), alpha = 0.5, lwd = 0.5) +
  scale_color_manual(values = c("black", "darkgrey", "blue")) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Log2 transformed normalized quantity", fill = "Time point", x = NULL)

# Limma -------------------------------------------------------------------
library(limma)

#remove the last 5 characters from column names
group <- substr(colnames(norm[,10:80]), 1, nchar(colnames(norm[,10:80]))-5)
unique(group)

#design matrix
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))

fit <- lmFit(norm[,10:80], design)

cont.matrix <- makeContrasts(
  WT_5min = WT_TNF_05min - WT_Ctrl_05min,
  WT_15min = WT_TNF_15min - WT_Ctrl_15min,
  KO_5min = KO_TNF_05min - KO_Ctrl_05min,
  KO_15min = KO_TNF_15min - KO_Ctrl_15min,
  ND_5min = ND_TNF_05min - ND_Ctrl_05min,
  ND_15min = ND_TNF_15min - ND_Ctrl_15min,
  KO_WT_5min = KO_Ctrl_05min - WT_Ctrl_05min,
  KO_WT_15min = KO_Ctrl_15min - WT_Ctrl_15min,
  ND_WT_5min = ND_Ctrl_05min - WT_Ctrl_05min,
  ND_WT_15min = ND_Ctrl_15min - WT_Ctrl_15min,
  KO_WT_TNF_5min = KO_TNF_05min - WT_TNF_05min,
  KO_WT_TNF_15min = KO_TNF_15min - WT_TNF_15min,
  ND_WT_TNF_5min = ND_TNF_05min - WT_TNF_05min,
  ND_WT_TNF_15min = ND_TNF_15min - WT_TNF_15min,
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2, trend = T)

coefficients <- c(
  "WT_5min",
  "WT_15min",
  "KO_5min",
  "KO_15min",
  "ND_5min",
  "ND_15min",
  "KO_WT_5min",
  "KO_WT_15min",
  "ND_WT_5min",
  "ND_WT_15min",
  "KO_WT_TNF_5min",
  "KO_WT_TNF_15min",
  "ND_WT_TNF_5min",
  "ND_WT_TNF_15min"
)

tabs <- list()

for (i in coefficients) {
  tab <- topTable(fit3, coef = i, number = Inf, sort.by = "none")
  tab$contrast <- i
  tab <- cbind(norm[,1:9], tab)
  tab <- tab %>%
    mutate(
      Regulation = case_when(
        logFC > 2*sd(tab$logFC[tab$logFC > 0], na.rm = T) & adj.P.Val <= 0.01 ~ "Up-regulated",
        logFC < -(2*sd(tab$logFC[tab$logFC > 0], na.rm = T)) & adj.P.Val <= 0.01 ~ "Down-regulated",
        TRUE ~ "Unchanged"
      )
    )
  tabs[[i]] <- tab
}

limma <- do.call(rbind, tabs)

write.table(limma, "limma.txt", sep = "\t", quote = F, row.names = F)

limma <- fread("limma.txt")
limma <- as.data.frame(limma)

#density plot fc
ggplot(limma, aes(x = logFC, color = contrast)) +
  geom_density(alpha = 0.5) +
  theme_minimal(base_size = 12)

# Count regulated ids ----------------------------------------------------
count <- limma %>%
  group_by(contrast, Regulation) %>%
  summarise(Count = n_distinct(PTM.CollapseKey))
write.table(count, "regulated_id_count_cutoff.txt", sep = "\t", quote = F, row.names = F)

count <- limma %>%
  mutate(
    Regulation = case_when(
      logFC > 0 & adj.P.Val <= 0.05 ~ "Up-regulated",
      logFC < 0 & adj.P.Val <= 0.05 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  ) %>%
  group_by(contrast, Regulation) %>%
  summarise(Count = n_distinct(PTM.CollapseKey))
write.table(count, "regulated_id_count.txt", sep = "\t", quote = F, row.names = F)

count$Regulation <- factor(count$Regulation, levels = c("Up-regulated", "Down-regulated", "Unchanged"))

count1 <- count %>% filter(
  contrast == "WT_5min" |
    contrast == "WT_15min" |
    contrast == "KO_5min" |
    contrast == "KO_15min" |
    contrast == "ND_5min" | 
    contrast == "ND_15min"
)
count1$contrast <- factor(
  count1$contrast,
  levels = c(
    "WT_5min",
    "WT_15min",
    "KO_5min",
    "KO_15min",
    "ND_5min",
    "ND_15min"
  )
)

#bar plot
ggplot(count1 %>% filter(Regulation != "Unchanged"), aes(x = contrast, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width = 0.9), size=2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Number of regulated phosphosites", fill = "Regulation") +
  facet_wrap(~ Regulation)

count2 <- count %>% filter(
  contrast == "KO_WT_5min" |
    contrast == "KO_WT_15min" |
    contrast == "ND_WT_5min" | 
    contrast == "ND_WT_15min" |
    contrast == "KO_WT_TNF_5min" |
    contrast == "KO_WT_TNF_15min" |
    contrast == "ND_WT_TNF_5min" |
    contrast == "ND_WT_TNF_15min"
)
count2$contrast <- factor(
  count2$contrast,
  levels = c(
    "KO_WT_5min",
    "KO_WT_15min",
    "ND_WT_5min",
    "ND_WT_15min",
    "KO_WT_TNF_5min",
    "KO_WT_TNF_15min",
    "ND_WT_TNF_5min",
    "ND_WT_TNF_15min"
  )
)

ggplot(count2 %>% filter(Regulation != "Unchanged"), aes(x = contrast, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width = 0.9), size=2) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Number of regulated phosphosites", fill = "Regulation") +
  facet_wrap(~ Regulation)

# Eulerr plot ------------------------------------------------------------
library(eulerr)

#5 minutes
WT <- limma %>%
  filter(contrast == "WT_5min") %>%
  #filter(adj.P.Val <= 0.05) %>%
  filter(Regulation != "Unchanged") %>%
  select(PTM.CollapseKey) %>%
  distinct() %>%
  pull()

KO <- limma %>%
  filter(contrast == "KO_5min") %>%
  #filter(adj.P.Val <= 0.05) %>%
  filter(Regulation != "Unchanged") %>%
  select(PTM.CollapseKey) %>%
  distinct() %>%
  pull()

ND <- limma %>%
  filter(contrast == "ND_5min") %>%
  #filter(adj.P.Val <= 0.05) %>%
  filter(Regulation != "Unchanged") %>%
  select(PTM.CollapseKey) %>%
  distinct() %>%
  pull()

overlap <- list(
  WT = WT,
  KO = KO,
  ND = ND
)

#plot(euler(overlap), quantities = T, fills = c("lightblue", "lightpink", "lightgreen"), main = "Regulated phosphosites 5 minutes wo cutoff")
plot(euler(overlap), quantities = T, fills = c("lightblue", "lightpink", "lightgreen"), main = "Regulated phosphosites 5 minutes with cutoff")

#15 minutes
WT <- limma %>%
  filter(contrast == "WT_15min") %>%
  #filter(adj.P.Val <= 0.05) %>%
  filter(Regulation != "Unchanged") %>%
  select(PTM.CollapseKey) %>%
  distinct() %>%
  pull()

KO <- limma %>%
  filter(contrast == "KO_15min") %>%
  #filter(adj.P.Val <= 0.05) %>%
  filter(Regulation != "Unchanged") %>%
  select(PTM.CollapseKey) %>%
  distinct() %>%
  pull()

ND <- limma %>%
  filter(contrast == "ND_15min") %>%
  #filter(adj.P.Val <= 0.05) %>%
  filter(Regulation != "Unchanged") %>%
  select(PTM.CollapseKey) %>%
  distinct() %>%
  pull()

overlap <- list(
  WT = WT,
  KO = KO,
  ND = ND
)

#plot(euler(overlap), quantities = T, fills = c("lightblue", "lightpink", "lightgreen"), main = "Regulated phosphosites 15 minutes wo cutoff")
plot(euler(overlap), quantities = T, fills = c("lightblue", "lightpink", "lightgreen"), main = "Regulated phosphosites 15 minutes with cutoff")

# Heatmaps ---------------------------------------------------------------

#5 minutes
#with cutoff
limma_reg_5min <- limma %>% 
  filter(contrast %in% c("WT_5min", "KO_5min", "ND_5min") & Regulation != "Unchanged")

#without cutoff
limma_reg_5min <- limma %>%
  filter(contrast %in% c("WT_5min", "KO_5min", "ND_5min") & adj.P.Val <= 0.05)

#different between TNF samples
limma2_reg_5min <- limma %>%
  filter(contrast %in% c("KO_WT_TNF_5min", "ND_WT_TNF_5min ") & P.Value <= 0.05)

length(unique(limma_reg_5min$PTM.CollapseKey))
length(unique(limma2_reg_5min$PTM.CollapseKey))

#z score
long_zscore_5min <- long %>%
  filter(PTM.CollapseKey %in% limma_reg_5min$PTM.CollapseKey) %>%
  #filter(PTM.CollapseKey %in% limma2_reg_5min$PTM.CollapseKey) %>%
  filter(tp == "5min") %>%
  mutate(id = paste0(PG.Genes, "_", PTM.SiteAA, PTM.SiteLocation, "_", PTM.Multiplicity)) %>%
  group_by(PTM.CollapseKey, id) %>%
  mutate(z = scale(quantity))

#mean
long_zscore_5min_mean <- long_zscore_5min %>%
  group_by(PTM.CollapseKey, id, genotype, treatment, rep) %>%
  summarise(z = mean(z, na.rm = T)) %>%
  mutate(condition = paste0(genotype, "_", treatment, "_", rep))

#wide format
wide_z_5min_mean <- reshape2::dcast(
  long_zscore_5min_mean,
  id ~ condition,
  value.var = "z"
)
rownames(wide_z_5min_mean) <- wide_z_5min_mean$id
wide_z_5min_mean <- wide_z_5min_mean %>% select(contains("WT"), contains("KO"), contains("ND"))

#column groups
condition <- ifelse(grepl("Ctrl", colnames(wide_z_5min_mean)), "Ctrl", "TNF")
genotype <- ifelse(grepl("WT", colnames(wide_z_5min_mean)), "WT", ifelse(grepl("KO", colnames(wide_z_5min_mean)), "KO", "ND"))
cols <- as.data.frame(cbind(condition, genotype))
cols$condition <- factor(cols$condition, levels = c("Ctrl", "TNF"))
cols$genotype <- factor(cols$genotype, levels = c("WT", "KO", "ND"))
rownames(cols) <- colnames(wide_z_5min_mean)

min(wide_z_5min_mean, na.rm = T)
max(wide_z_5min_mean, na.rm = T)

library(pheatmap)
library(RColorBrewer)
breaks <- seq(-2, 2, length.out = 101)
breaks <- seq(-2.5, 2.5, length.out = 101)
pheatmap(
  wide_z_5min_mean,
  color = colorRampPalette(rev(brewer.pal(n=9, name ="RdBu")))(100),
  show_rownames = T,
  fontsize_row = 7,
  cluster_cols = F,
  breaks = breaks,
  annotation_col = cols,
  #border_color = NA
)

#15 minutes
#with cutoff
limma_reg_15min <- limma %>%
  filter(contrast %in% c("WT_15min", "KO_15min", "ND_15min") & Regulation != "Unchanged")

#without cutoff
limma_reg_15min <- limma %>%
  filter(contrast %in% c("WT_15min", "KO_15min", "ND_15min") & adj.P.Val <= 0.05)

#different between TNF samples
limma2_reg_15min <- limma %>%
  filter(contrast %in% c("KO_WT_TNF_15min", "ND_WT_TNF_15min") & P.Value <= 0.05)

length(unique(limma_reg_15min$PTM.CollapseKey))
length(unique(limma2_reg_15min$PTM.CollapseKey))

#z score
long_zscore_15min <- long %>%
  filter(PTM.CollapseKey %in% limma_reg_15min$PTM.CollapseKey) %>%
  #filter(PTM.CollapseKey %in% limma2_reg_15min$PTM.CollapseKey) %>%
  filter(tp == "15min") %>%
  mutate(id = paste0(PG.Genes, "_", PTM.SiteAA, PTM.SiteLocation, "_", PTM.Multiplicity)) %>%
  group_by(PTM.CollapseKey, id) %>%
  mutate(z = scale(quantity))

#mean
long_zscore_15min_mean <- long_zscore_15min %>%
  group_by(PTM.CollapseKey, id, genotype, treatment, rep) %>%
  summarise(z = mean(z, na.rm = T)) %>%
  mutate(condition = paste0(genotype, "_", treatment, "_", rep))

#wide format
wide_z_15min_mean <- reshape2::dcast(
  long_zscore_15min_mean,
  id ~ condition,
  value.var = "z"
)
rownames(wide_z_15min_mean) <- wide_z_15min_mean$id
wide_z_15min_mean <- wide_z_15min_mean %>% select(contains("WT"), contains("KO"), contains("ND"))

#column groups
condition <- ifelse(grepl("Ctrl", colnames(wide_z_15min_mean)), "Ctrl", "TNF")
genotype <- ifelse(grepl("WT", colnames(wide_z_15min_mean)), "WT", ifelse(grepl("KO", colnames(wide_z_15min_mean)), "KO", "ND"))
cols <- as.data.frame(cbind(condition, genotype))
cols$condition <- factor(cols$condition, levels = c("Ctrl", "TNF"))
cols$genotype <- factor(cols$genotype, levels = c("WT", "KO", "ND"))
rownames(cols) <- colnames(wide_z_5min_mean)

min(wide_z_15min_mean, na.rm = T)
max(wide_z_15min_mean, na.rm = T)

breaks <- seq(-2.5, 2.5, length.out = 101)
pheatmap(
  wide_z_15min_mean,
  color = colorRampPalette(rev(brewer.pal(n=9, name ="RdBu")))(100),
  show_rownames = F,
  #fontsize_row = 7,
  cluster_cols = F,
  breaks = breaks,
  #border_color = NA,
  annotation_col = cols
)

# Plot proteins belonging to TNF and NFkB pathways -------------------------
library(KEGGREST)

# Load KEGG pathways
pathway_tnf <- keggGet("hsa04668") # TNF pathway
pathway_nfb <- keggGet("hsa04064") # NF-kB pathway

# Extract genes for TNF pathway
genes_tnf <- pathway_tnf[[1]]$GENE
descriptions_tnf <- genes_tnf[seq(2, length(genes_tnf), by = 2)]
gene_names_tnf <- sapply(descriptions_tnf, function(x) {
  strsplit(x, ";")[[1]][1]  # Split by ";" and take the first part
})

# Extract genes for NF-kB pathway
genes_nfb <- pathway_nfb[[1]]$GENE
descriptions_nfb <- genes_nfb[seq(2, length(genes_nfb), by = 2)]
gene_names_nfb <- sapply(descriptions_nfb, function(x) {
  strsplit(x, ";")[[1]][1]  # Split by ";" and take the first part
})

# Combine gene names from both pathways
gene_names <- unique(c(gene_names_tnf, gene_names_nfb))

#are there any ";" in the column PG.Genes in norm?
any(grepl(";", norm$PG.Genes))

#Remove text after ";" in the column PG.Genes in norm
norm$Gene <- gsub(";.*", "", norm$PG.Genes)

#regulated sites
#without cutoff
limma_reg <- limma %>% filter(adj.P.Val <= 0.05)
length(unique(limma_reg$PTM.CollapseKey))

#regulated sites
#with p value cutoff
limma_reg <- limma %>% filter(P.Value <= 0.01)
length(unique(limma_reg$PTM.CollapseKey))

#z score & mean by replicate
long_zscore_mean <- long %>%
  filter(PG.Genes %in% gene_names) %>%
  filter(PTM.CollapseKey %in% limma_reg$PTM.CollapseKey) %>%
  mutate(id = paste0(PG.Genes, "_", PTM.SiteAA, PTM.SiteLocation, "_", PTM.Multiplicity)) %>%
  group_by(PTM.CollapseKey, id) %>%
  mutate(z = scale(quantity)) %>%
  group_by(PTM.CollapseKey, id, genotype, treatment, tp, rep) %>%
  summarise(z = mean(z, na.rm = T)) %>%
  mutate(condition = paste0(genotype, "_", treatment, "_", tp, "_", rep))

#z score & mean by condition
long_zscore_mean <- long %>%
  filter(PG.Genes %in% gene_names) %>%
  filter(PTM.CollapseKey %in% limma_reg$PTM.CollapseKey) %>%
  mutate(id = paste0(PG.Genes, "_", PTM.SiteAA, PTM.SiteLocation, "_", PTM.Multiplicity)) %>%
  group_by(PTM.CollapseKey, id) %>%
  mutate(z = scale(quantity)) %>%
  group_by(PTM.CollapseKey, id, genotype, treatment, tp) %>%
  summarise(z = mean(z, na.rm = T)) %>%
  mutate(condition = paste0(genotype, "_", treatment, "_", tp, "_"))

#wide format
wide_z_mean <- reshape2::dcast(
  long_zscore_mean,
  id ~ condition,
  value.var = "z"
)
rownames(wide_z_mean) <- wide_z_mean$id
wide_z_mean <- wide_z_mean %>% select(contains("_5min"), contains("_15min"))
wide_z_mean <- wide_z_mean %>% select(contains("WT"), contains("KO"), contains("ND"))

#column groups
condition <- ifelse(grepl("Ctrl", colnames(wide_z_mean)), "Ctrl", "TNF")
genotype <- ifelse(grepl("WT", colnames(wide_z_mean)), "WT", ifelse(grepl("KO", colnames(wide_z_mean)), "KO", "ND"))
time.point <- ifelse(grepl("_5min", colnames(wide_z_mean)), "5min", "15min")
cols <- as.data.frame(cbind(condition, genotype, time.point))
cols$condition <- factor(cols$condition, levels = c("Ctrl", "TNF"))
cols$genotype <- factor(cols$genotype, levels = c("WT", "KO", "ND"))
cols$time.point <- factor(cols$time.point, levels = c("5min", "15min"))
rownames(cols) <- colnames(wide_z_mean)

min(wide_z_mean, na.rm = T)
max(wide_z_mean, na.rm = T)

library(pheatmap)
library(RColorBrewer)
breaks <- seq(-2, 2, length.out = 101)
pheatmap(
  wide_z_mean,
  color = colorRampPalette(rev(brewer.pal(n=9, name ="RdBu")))(100),
  show_rownames = T,
  fontsize_row = 6,
  cluster_cols = F,
  breaks = breaks,
  annotation_col = cols,
  display_numbers = T,
  fontsize_number = 5,
  main = "Regulated phosphosites in the KEGG TNF and NFkB signaling pathways"
)

#regulated sites
#with p value cutoff
#15 minutes
limma_reg_15min <- limma %>% filter(contrast %in% c("WT_15min", "KO_15min", "ND_15min") & P.Value <= 0.01)
length(unique(limma_reg_15min$PTM.CollapseKey))

#z score & mean by replicate
long_zscore_15min_mean <- long %>%
  filter(PG.Genes %in% gene_names) %>%
  filter(PTM.CollapseKey %in% limma_reg_15min$PTM.CollapseKey) %>%
  filter(tp == "15min") %>%
  mutate(id = paste0(PG.Genes, "_", PTM.SiteAA, PTM.SiteLocation, "_", PTM.Multiplicity)) %>%
  group_by(PTM.CollapseKey, id) %>%
  mutate(z = scale(quantity)) %>%
  group_by(PTM.CollapseKey, id, genotype, treatment, tp, rep) %>%
  summarise(z = mean(z, na.rm = T)) %>%
  mutate(condition = paste0(genotype, "_", treatment, "_", tp, "_", rep))

#wide format
wide_z_15min_mean <- reshape2::dcast(
  long_zscore_15min_mean,
  id ~ condition,
  value.var = "z"
)
rownames(wide_z_15min_mean) <- wide_z_15min_mean$id
wide_z_15min_mean <- wide_z_15min_mean %>% select(contains("WT"), contains("KO"), contains("ND"))

#column groups
condition <- ifelse(grepl("Ctrl", colnames(wide_z_15min_mean)), "Ctrl", "TNF")
genotype <- ifelse(grepl("WT", colnames(wide_z_15min_mean)), "WT", ifelse(grepl("KO", colnames(wide_z_15min_mean)), "KO", "ND"))
cols <- as.data.frame(cbind(condition, genotype))
cols$condition <- factor(cols$condition, levels = c("Ctrl", "TNF"))
cols$genotype <- factor(cols$genotype, levels = c("WT", "KO", "ND"))
rownames(cols) <- colnames(wide_z_15min_mean)

min(wide_z_15min_mean, na.rm = T)
max(wide_z_15min_mean, na.rm = T)

library(pheatmap)
library(RColorBrewer)
breaks <- seq(-2, 2, length.out = 101)
pheatmap(
  wide_z_15min_mean,
  color = colorRampPalette(rev(brewer.pal(n=9, name ="RdBu")))(100),
  show_rownames = T,
  fontsize_row = 6,
  cluster_cols = F,
  breaks = breaks,
  annotation_col = cols,
  display_numbers = T,
  fontsize_number = 4,
  main = "Regulated phosphosites in the KEGG TNF and NFkB signaling pathways"
)





