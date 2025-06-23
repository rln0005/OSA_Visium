#!/usr/bin/Rscript

# Load libraries
source("./custom_functions.R")
library(spatstat.explore)
library(Seurat)
library(SpatialExperiment)
library(tidyverse)
library(STdeconvolve)
library(spacexr)
library(corrplot)
library(ComplexHeatmap)
library(RColorBrewer)
library(patchwork)

# Helper function for colorization
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


################################################################################
### Purpose
# This script uses reference-based and reference free deconvolution approaches
# to help annotate integrated clsuters and also evaluate the data for cell type
# colocalization.

# Run RCTD deconcolution (reference-based method) - after running use 
# estimated predictions to assess integrated cluster composition. Then use
# Pearson correlation to evaluate the data for cell type co-localization. Lastly,
# extract the largest cell type proportion in each spot and use those values to
# annotate indiviudal tissues to compare with histo annotations.

### Input
# "merge.rds"  - starting point for all analysis
# "canine_naive_n6_annotated.rds" - annotated canine OS scRNA data

### Analysis note
# The lapply run RCTD step takes ~ 1.5 hours to run when completed on a compute
# node with 16 cores.

################################################################################


################################################################################
## Run RCTD deconcolution - reference based deconvolution ## <<<
outName <- "decon"

# Load in the integrated dataset and split by sample
seu.obj <- readRDS("../input/merge.rds")
seu.obj <- loadMeta(
    seu.obj = seu.obj, 
    metaFile = "./metaData/decoder.csv", 
    groupBy = "sample", metaAdd = "short_name"
)
seu.obj.list <- SplitObject(seu.obj, split.by = "sample")

# Prep the reference
sc_ref <- readRDS("../external_data/canine_naive_n6_annotated.rds")
counts <- sc_ref@assays$RNA@counts
sc_ref$celltype.l2 <- droplevels(sc_ref$celltype.l2)
Idents(sc_ref) <- "celltype.l2"
sc_ref <- RenameIdents(
  sc_ref,
  c(
    `B cell` = "B cell", `Plasma cell` = "Plasma cell", CD320_OC = "OC", 
    Cycling_OC = "OC", cDC1 = "DC", cDC2 = "DC", mregDC = "mregDC", 
    pDC = "pDC", preDC = "DC", `Endothelial cell` = "Endothelial cell", 
    Fibroblast = "Fibroblast", Hypoxic_osteoblast = "Hypoxic_osteoblast", 
    `IFN-TAM` = "TAM", `IFN-osteoblast` = "Osteoblast", `Mast cell` = "Mast cell", 
    Mature_OC = "OC", NK = "T_NK", Neutrophil = "Neutrophil", 
    Osteoblast_1 = "Osteoblast", Osteoblast_2 = "Osteoblast", Osteoblast_3 = "Osteoblast", 
    Osteoblast_cycling = "Osteoblast_cycling", `CD4 T cell` = "T_NK", 
    `CD8 T cell` = "T_NK", ANGIO_TAM = "TAM", `LA-TAM_C1QC_hi` = "TAM", 
    `LA-TAM_SPP2_hi` = "TAM", TAM_ACT = "TAM", TAM_INT = "TAM", `CD4+_TIM` = "TIM", 
    `CD4-_TIM` = "TIM", T_IFN = "T_NK", T_cycling = "T_NK"
  )
)
sc_ref$celltype.l2 <- Idents(sc_ref)
cell_types <- setNames(sc_ref$celltype.l2, colnames(sc_ref))
nUMI <- setNames(sc_ref$nCount_RNA, colnames(sc_ref))
reference <- Reference(counts, cell_types, nUMI)

colz.df <- data.frame(
  cell_type = sort(levels(cell_types)),
  cell_grp = c(
    "Immune", "Immune", "Stroma", "Stroma", "Tumor", "Immune", "Immune",
    "Stroma", "Tumor", "Tumor", "Immune", "Immune", "Immune", "Immune", 
    "Immune", "Immune"
  ),
  newCol = c(
    "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
    "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", 
    "darkturquoise", "darkorange4", "brown"
  )
) %>%
  arrange(cell_type) %>%
  mutate(
    cell_shape = case_when(
      cell_grp == "Immune" ~ 21,
      cell_grp == "Stroma" ~ 22,
      cell_grp == "Tumor" ~ 23
    )
  )

sample_decoder <- read.csv("./metaData/decoder.csv")

# Loop samples to run RCTD and save .rds files
lapply(seq_along(seu.obj.list)[12:(length(seu.obj.list))], function(i){
  #focus on one sample
  sampleName <- names(seu.obj.list)[i]
  seu.sub <- seu.obj.list[[i]]
  #extract required data
  counts <- seu.sub@assays$SCT$counts
  coords <- GetTissueCoordinates(seu.sub) %>%
    select(-cell, x, y)
  nUMI <- colSums(counts)
  puck <- SpatialRNA(coords, counts, nUMI)
  barcodes <- colnames(puck@counts)
  #Run RCTD
  myRCTD <- create.RCTD(puck, reference, max_cores = 16)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  saveRDS(myRCTD, file.path(paste0("../output/", outName, "/midres/", sampleName, "_RCTD.rds")))
}) # previously ran and saved output to "../output/decon/highres"

# Loop RCTD results to extract cell type predictions -- can resume here
dat.list <- lapply(seq_along(seu.obj.list), function(i){
  #focus on one sample
  sampleName <- names(seu.obj.list)[i]
  #load data
  myRCTD <- readRDS(paste0("../output/", outName, "/midres/", sampleName, "_RCTD.rds"))
  weights <- myRCTD@results$weights
  norm_weights <- normalize_weights(weights)
  #return req data
  dat <- rownames_to_column(as.data.frame(norm_weights))
  dat$name <- sample_decoder$short_name[match(sampleName, sample_decoder$sample)]
  return(dat)
})
dat <- do.call(rbind, dat.list)

#get predictions in long data format by cluster
dat.sum <- dat %>% 
  select(-name) %>%
    left_join(
      select(
        rownames_to_column(seu.obj@meta.data), 
        c(rowname, seurat_clusters, sample)), 
      by = "rowname"
    ) %>%
    group_by(seurat_clusters) %>%
    summarise_if(is.numeric, mean) %>%
    pivot_longer(cols = colnames(.)[2:ncol(.)]) 

# Plot pie charts for each cluster
cell_types <- sort(levels(sc_ref$celltype.l2))
cell_colz <- with(colz.df, setNames(newCol, cell_type))
#loop through each cluster to generate a pie chart for each
pis <- lapply(0:11, function(i){ # leave out ganglionic tissue
  #calc label position - will need to create "other" category
  cluster_freq.table <- dat.sum %>% 
    filter(seurat_clusters == i) %>%
    mutate(
      value = round(value, 2),
      pos = 1-(cumsum(value) - (0.5 * value))
    ) %>% 
    arrange(desc(value))
  cluster_freq.table$pos <- rev(cluster_freq.table$pos)

  df2 <- cluster_freq.table %>% 
    mutate(
      csum = rev(cumsum(rev(value))), 
      pos = value/2 + lead(csum, 1),
      pos = if_else(is.na(pos), value/2, pos),
      lab = ifelse(value >= 0.05,
                   paste0(name, " (", value*100, "%)"),
                   "")
    )
  #create the plot
  cluster_freq.table$name <- factor(cluster_freq.table$name, 
                                    levels = cluster_freq.table$name)
  p <- ggplot(cluster_freq.table, aes(x = "" , y = value, fill = name)) +
    geom_col(width = 1) +
    coord_polar(theta = "y", clip = "off") +
    ggrepel::geom_text_repel(
      data = df2, aes(y = pos, label = lab),size = 3, 
      nudge_x = rep(0.9, 5), show.legend = FALSE
    ) +
    guides(fill = guide_legend(title = "Cell type", nrow = 10)) +
    theme_void() +
    theme(
      plot.margin = unit(c(-7, -7, -7, -7), "pt"),
      plot.title = element_text(size = 14, vjust = -2, hjust = 0.5)
    ) +
    scale_fill_manual(values = cell_colz) +
    ggtitle(paste0("Cluster ", i)) + NoLegend()
  return(p)
})

p <- ggplot(dat.sum, aes(x = "" , y = value, fill = name)) +
    geom_col(width = 1) +
    coord_polar(theta = "y", clip = "off") +
    guides(fill = guide_legend(title = "Cell type", nrow = 3)) +
    theme_void() +
    theme(
      legend.text = element_text(size=12, margin = margin(7, 14, 7, 3)),
      legend.title = element_blank(),
      legend.key.spacing = unit(0.5, 'cm'),
      legend.key.size = unit(0.75, 'cm'),
      plot.margin = unit(c(-21, 21, -7, 21), "pt")
    ) +
    scale_fill_manual(values = cell_colz)
legg <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
#patchwork together the final plot
pp <- Reduce(`+`, pis) 
leg <- ggpubr::as_ggplot(legg[[1]]) + theme(legend.box.margin = unit(c(-21, 21, -7, 21), "pt"))
pp + leg + plot_layout(
  design = "ABCD
            EFGH
            IJKL
            MMMM")
ggsave(paste0("../output/", outName, "/", outName, "_pieChart.png"), width = 5, height = 5, scale = 2)
ggsave(plot = legg[[1]], file = paste0("../output/", outName, "/", outName, "_pieChart_legend.png"), width = 5, height = 1, scale = 2)

# Assess changes in the fibroblast abundance between primary and met samples
dat.sum <- dat %>% 
  select(-name) %>%
    left_join(
      select(
        rownames_to_column(seu.obj@meta.data), 
        c(rowname, seurat_clusters, sample)), 
      by = "rowname"
    ) %>%
    group_by(sample) %>%
    summarise_if(is.numeric, mean) %>%
    mutate(
      tissueSource = factor(
        ifelse(grepl("met", sample), "Metastatic", "Primary"),
        levels = c("Primary", "Metastatic")
      )
    ) %>% 
    select(sample, Fibroblast, tissueSource) %>%
    rename(Sample = "sample")

p <- ggplot(dat.sum, aes(x = tissueSource, y = Fibroblast, colour = tissueSource, fill = tissueSource)) + 
  geom_boxplot(linewidth = 1, alpha = 0.5) +
  guides(fill = "none") +
  guides(colour = "none") +
  ggnewscale::new_scale_colour() +
  geom_point(aes(colour = Sample), size = 2, position = position_jitter(width = 0.25)) +
  ggpubr::stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)), label.x.npc = "right", 
    label.y.npc = 1, vjust = 0, size = 3
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  labs(y = "Estimated fibroblast fraction")
ggsave(paste0("../output/", outName, "/fib_comp.png"), width = 5, height = 5)

# UMAP of clusters
colz <- gg_color_hue(13)
colz[13] <- "grey"

p <- DimPlot(
  seu.obj, reduction = "umap", group.by = "seurat_clusters", label.box = T, label = T,
  cols = colz
)
cusLabels(p, smallAxes = T, size = 10, textSize = 5, alpha = 0.8) #+ NoLegend()
ggsave(paste0("../output/", outName, "/", "umap_intClus.png"), height = 7, width = 7)

# UMAP of samples
colz <- gg_color_hue(14)

p <- DimPlot(
  seu.obj, reduction = "umap", group.by = "short_name", label.box = F, label = F,
  cols = colz
)
p <- formatUMAP(p)
p <- p + 
  theme(
    axis.title = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "pt")
  )
ggsave(paste0("../output/", outName, "/", "umap_intSample.png"), height = 7, width = 8)

# Spatial plots to find integrated cluster 5
seu.obj$highlight <- ifelse(seu.obj$seurat_clusters == 5, "OC", "Other")
p <- SpatialPlot(
  seu.obj, group.by = "highlight", pt.size.factor = 1.6, ncol = 7, crop = TRUE, 
  label = F
) & 
  theme(legend.title = element_blank()) & 
  NoLegend()
ggsave(paste0("../output/", outName, "/slide_visium.png"), height = 7, width = 12)    


# Extract barcodes for regions to exclude
inNames <- list.files(path = "../input/path_anno/", full.names = F)
inNames <- inNames[-c(6, 10)]
anno.list <- lapply(inNames, function(name){
  dat <- read.csv(paste0("../input/path_anno/", name))
  dat <- dat %>%
    mutate(
      sample_name = gsub(".csv", "", name),
      Barcode = paste0(sample_name, "_", Barcode)
    )
  colnames(dat) <- c("barcode", "path_anno", "name")
  return(dat)
})
anno <- do.call(rbind, anno.list)
exclude <- anno %>%
  filter(path_anno %in% c("Skeletal muscle", "Ganglionic tissue")) %>%
  pull(barcode)

# Assess cell types for colocalization - as a whole
datWhole <- dat %>% 
  filter(! rowname %in% exclude) %>%
  column_to_rownames() %>%
  select(-name)
png(file = paste0("../output/", outName, "/cor.png"), 
    width = 3000, height = 3000, res = 300)
par(mfcol = c(1, 1))         
corrplot(cor(datWhole), method = "number", type = "upper")
dev.off()

mat_melt <- function(mat){
  mat <- mat[rev(rownames(mat)), ]
  mat[row(mat) + col(mat) > nrow(mat) + 1] <- NA
  reshape2::melt(mat, na.rm = TRUE) %>%
    mutate(
      ct_cor = paste0(rowname, "-", variable)
    ) %>%
    select(-rowname, -variable)
}

# Assess cell types for colocalization - at individual sample level
dat_clean <- dat %>% 
  filter(! rowname %in% exclude)

datIndv.list <- dat_clean %>%
  select(-name) %>%
  split(., dat_clean$name)

cor_dat.list <- lapply(datIndv.list, function(x){
  sampleName <- substr(rownames(x)[1], 1, nchar(rownames(x)[1]) - 19)
  rownames(x) <- NULL
  x <- column_to_rownames(x)
  cor.mat <- rstatix::cor_mat(x)
  sig_mat <- rstatix::cor_get_pval(cor.mat)
  melted_cormat <- mat_melt(cor.mat)
  melted_pval <- mat_melt(sig_mat)
  melted_pval$value <- p.adjust(p = melted_pval$value)
  return(list(melted_cormat, melted_pval))
})
#compile cor results
cor_dat <- lapply(cor_dat.list, `[[`, 1) %>% 
  reduce(left_join, by = "ct_cor") %>%
  column_to_rownames("ct_cor")
cor_dat.mat <- as.matrix(cor_dat)
colnames(cor_dat.mat) <- names(cor_dat.list)
ind_dat <- cor_dat.mat
#compile pval results
pval_dat <- lapply(cor_dat.list, `[[`, 2) %>% 
  reduce(left_join, by = "ct_cor") %>%
  column_to_rownames("ct_cor")
pval_dat.mat <- as.matrix(pval_dat)
colnames(pval_dat.mat) <- names(cor_dat.list)
#remove the diagonal values
cor_dat.mat <- cor_dat.mat[rowSums(cor_dat.mat) != ncol(cor_dat.mat), ]
pval_dat.mat <- pval_dat[rowSums(cor_dat.mat) != ncol(cor_dat.mat), ]
#flag the sig values
sig_mat <- cor_dat.mat
sig_mat[abs(cor_dat.mat) > 0.3 & pval_dat < 0.05] <- "*"
sig_mat[sig_mat != "*"] <- ""
#only include interactions with at least one sig interaction across samples
cor_dat.mat <- cor_dat.mat[rowSums(sig_mat == "*") > 0, ]
sig_mat <- sig_mat[rowSums(sig_mat == "*") > 0, ]
#draw heatmap
ht <- Heatmap(
  name = "Pearson cor",
  mat = cor_dat.mat,
  cluster_rows = T,
  cluster_columns = T,
  border_gp = gpar(col = "black"),
  heatmap_legend_param = list(direction = "horizontal"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sig_mat[i, j], x, (y - convertHeight(grobHeight(textGrob("*")), "mm") * 0.45),
      gp = gpar(fontsize = 16, col = "grey80", vjust = 0))
  }
)
png(file = paste0("../output/", outName, "/cor_ht", ".png"), 
    width = 2750, height = 3500, res = 400)
par(mfcol = c(1, 1))
draw(ht, padding = unit(c(5, 5, 5, 15), "mm"), heatmap_legend_side = "top")
dev.off()

# Create cor plot for each sample
lapply(seq_along(seu.obj.list), function(i){
  sampleName <- sample_decoder$short_name[match(names(seu.obj.list)[i], sample_decoder$sample)]  
  dat_sub <- dat %>% 
    filter(name == sampleName) %>%
    column_to_rownames() %>%
    select(-name)
  png(file = paste0("../output/", outName, "/", sampleName, "_cor.png"), 
      width = 3000, height = 3000, res = 300)
  par(mfcol = c(1, 1))         
  corrplot(cor(dat_sub), method = "number", type = "upper", title = paste0("\n", sampleName))
  dev.off()
})

### New approach to find if cell type co-localize in the tumor or stroma
# Load in the SpaCET classifications
inFiles <- list.files("../output/spacet/", pattern = ".csv", full.names = T)
dat.list <- lapply(inFiles, function(inFile){
  dat <- read.csv(inFile)
  sampleName <- paste(gsub("../output/spacet//|_i.*\\.csv", "", inFile))
  dat %>% 
    mutate(
      bc = paste0(sampleName, "_", bc)
    ) %>%
    column_to_rownames(var = "bc") %>%
    rename(spacet_orig = "int")
})
spacet_meta <- do.call(rbind, dat.list)
seu.obj <- AddMetaData(seu.obj, spacet_meta)
seu.obj$spacet <- ifelse(seu.obj$spacet_orig != "Tumor", "Stroma", "Tumor")

# Loop the significant interactions
convert_df <- which(sig_mat == "*", arr.ind = TRUE)

cts <- str_split(rownames(convert_df), "-")
int_df <- data.frame(
  sample = colnames(sig_mat)[convert_df[ , 2]],
  ct1 = sapply(cts, `[[`, 1),
  ct2 = sapply(cts, `[[`, 2),
  r = cor_dat.mat[which(sig_mat == "*")]
) %>%
  filter(r > 0)

set.seed(1)
fisher_dat <- apply(int_df, MARGIN = 1, function(x){
  dat <- datIndv.list[[x[1]]]
  # Use permutation testing
  n_perm <- 1000
  fdr_cut <- 0.2
  perm_maxs <- replicate(
    n_perm, {
      perm_ct1 <- sample(dat[ , x[2]])
      perm_ct2 <- sample(dat[ , x[3]])
      perm_maxs <- sqrt(perm_ct1 * perm_ct2)
      return(perm_maxs)
    })

  actual_maxs <- dat[ , c("rowname", x[2], x[3])] %>%
    rename_at(vars(c(x[2], x[3])), ~c("ct1", "ct2")) %>%
    mutate(
      maxs =  sqrt(ct1 * ct2)
    )
  
  p_values <- sapply(dat$rowname, function(y) {
    actual_maxs <- actual_maxs %>%
      filter(rowname == y) %>%
      pull(maxs)

    sum(perm_maxs[which(dat$rowname == y), ] >= actual_maxs) / n_perm
  })

  # Filter for sig spots
  adjusted_p_values <- p.adjust(p_values, method = "fdr")
  co_loc <- names(adjusted_p_values[adjusted_p_values < fdr_cut])
  
  if(length(co_loc) > 0) {
    coloc_dat <- data.frame(
      coloc = rep("yes", length(co_loc)), 
      row.names = co_loc
    )

    seu.obj <- AddMetaData(seu.obj, coloc_dat)
    seu.obj$coloc <- factor(ifelse(is.na(seu.obj$coloc), "no", "yes"), levels = c("yes", "no"))

    # Extract significant cell type interactions by sample
    dat_sub <- filter(seu.obj@meta.data, short_name == x[1])
    mat <- table(dat_sub$coloc, dat_sub$spacet)

    # Conduct fisher's exact test to determine if co-localiations have a bias
    res.fisher <- fisher.test(mat)
    res.fisher <- data.frame(
      "Sample" = x[1], 
      "ct1" = x[2], 
      "ct2" = x[3],
      "OR" = unname(res.fisher$estimate), 
      "Pvalue" = res.fisher$p.value
    ) %>% 
      mutate(
        "SE_OR" = sqrt(1/mat[1, 1] + 1/mat[1, 2] + 1/mat[2, 1] + 1/mat[2, 2]),
        "lower_95CI" = exp(log(OR) - 1.96 * SE_OR),
        "upper_95CI" = exp(log(OR) + 1.96 * SE_OR)
      )

    # Plot the spots most likely to co-localize
    sampleName <- sample_decoder$sample[match(x[1], sample_decoder$short_name)]  
    seu.sub <- seu.obj.list[[sampleName]]
    sub_dat <- seu.obj@meta.data %>%
      rownames_to_column() %>%
      filter(sample == sampleName) %>%
      left_join(GetTissueCoordinates(seu.sub, scale = "lowres"), by = c("rowname" = "cell")) %>%
      select(-sample) %>%
      relocate(x, y)

    # Make the background histo image
    p1 <- SpatialPlot(
      seu.sub, images = sampleName, image.scale = "lowres", features = NULL,
      group.by = NULL, pt.size.factor = 0, crop = TRUE, alpha = c(1, 1)
    ) + NoLegend() + coord_cartesian(expand = TRUE)

    # Extract the and lift over the coordinate system
    base_plot <- p1[[1]]
    sub_dat <- sub_dat %>%
      select(-x, -y) %>%
      left_join(base_plot$data, by = c("rowname" = "cell"))

    # Calculate edges to crop and transpose the tissue to align with the pie charts
    crop_buffer <- 0
    ggbuild <- ggplot_build(base_plot)
    x_range <- range(sub_dat$y)
    y_range <- range(sub_dat$x)
    aspect_ratio <- diff(x_range) / diff(y_range)
    sub_dat$x <- (max(sub_dat$x) - sub_dat$x + min(sub_dat$x)) * aspect_ratio

    # Generate the plot
    p2 <- ggplot() + 
      geom_point(
        data = sub_dat, mapping = aes(x = y, y = x, fill = spacet_orig, color = coloc), shape = 21, size = 2
      ) +
      theme_void() +
      NoLegend() +
      scale_fill_manual(values = c("Tumor" = "gold", "Stroma" = "grey", "Interface" = "black")) +
      scale_color_manual(values = c("yes" = "red", "no" = NA)) +
      coord_equal(expand = TRUE) +
      theme(
        plot.margin = margin(6, 6, 8, 8, "pt")
      )
    p_final <- p1 + inset_element(
      p2, left = 0, bottom = 0, right = 1, top = 1, align_to = "full"
    )
    ggsave(plot = p_final, paste0("../output/", outName, "/", sampleName, "_", x[2], "_", x[3], "_high_int.png"), height = 7, width = 7)    

    return(res.fisher)
  } else {
    return(NULL)
  }
})
fisher.df <- do.call(rbind, fisher_dat)
fisher.df <- fisher.df[fisher.df$Pvalue < 0.05 & fisher.df$SE_OR != Inf, ] %>% # 
  mutate(
    label = paste0(ct1, "--", ct2, " (", Sample, ")"),
    `Co-localization\nregion` = ifelse(OR > 1, "Stroma/Interface", "Tumor"),
    Sample = paste0(Sample, " "),
#     Site = factor(ifelse(grepl("met", label), "Metastatic", "Primary"), levels = c("Primary", "Metastatic")),
    Interaction = as.factor(paste0(ct1, "--", ct2))
  )

p1 <- ggplot(fisher.df, aes(x = OR, y = label)) + 
  geom_vline(aes(xintercept = 1), size = 0.25, linetype = "dashed") + 
  geom_errorbarh(
    aes(xmax = upper_95CI, xmin = lower_95CI), linewidth = 0.5, height = 0.2, color = "gray50"
  ) +
  geom_point(aes(color = `Co-localization\nregion`), size = 3.5) +
  theme_bw() +
  coord_trans(x = "log10") +
  scale_x_continuous(breaks = c(0, 1, 10, 100, 1000, 5000)) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(7, 7, 7, -21, "pt")
  ) +
  ylab("") +
  xlab("Odds ratio") +
  scale_colour_manual(values = c("Tumor" = "gold", "Stroma/Interface" = "grey"))
ggsave(plot = p1, paste0("../output/", outName, "/or.png"), height = 4, width = 7)    

p2 <- ggplot(fisher.df, aes(y = label)) + 
  geom_tile(aes(fill = Interaction, x = 0), width = 1, show.legend = T) +
  theme_void() +
  theme(
    plot.margin = margin(7, -21, 7, 75, "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  geom_vline(aes(xintercept = 0.5), size = 0.25, linetype = "solid") +
  guides(fill = guide_legend(ncol = 1)) + 
  coord_cartesian(clip = 'off') + 
  geom_text(
    aes(x = -0.5, label = Sample), hjust = 1
  )

p3 <- ggplot(fisher.df, aes(y = label)) + 
  geom_tile(aes(x = 1), width = 1, show.legend = T) +
  theme_void() +
  theme(
    plot.margin = margin(7, -21, 7, 7, "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  geom_vline(aes(xintercept = 0.5), size = 0.25, linetype = "solid")

# Compile the plotting
main_plot <- p2 + plot_spacer() + p1 + plot_layout(widths = c(0.1, -0.1, 0.8), guides = "collect")
ggsave(paste0("../output/", outName, "/or.png"), height = 4, width = 5)    

# Annotate each spot a cell type based on max prediction
#compile the cell types
ct_dat <- data.frame(
  barcode = rownames(datWhole),
  cell_type = colnames(datWhole)[apply(datWhole, MARGIN = 1, which.max)]
)
#add the metadata
seu.obj <- AddMetaData(
  seu.obj, with(ct_dat, setNames(cell_type, barcode)), col.name = "ct_pred"
)
seu.obj$ct_pred <- as.factor(seu.obj$ct_pred)
Idents(seu.obj) <- "ct_pred"
seu.obj <- RenameIdents(seu.obj, cell_shape)
seu.obj$cell_shape <- Idents(seu.obj)
cell_shape <- setNames(seu.obj$cell_shape, colnames(seu.obj))
#set colors for uniform plotting and save plot for each individual sample
pis <- lapply(names(seu.obj@images), function(x){
  cell_shape <- unname(cell_shape[seu.obj$sample == x])
  p <- SpatialPlot(
    seu.obj, group.by = "ct_pred", crop = TRUE, images = x, 
    cols = cell_colz#, shape = as.numeric(cell_shape)
  ) + NoLegend()
  ggsave(paste0("../output/", outName, "/ct_pred_", x, ".png"), height = 7, width = 7)    
})

# Create pie chart overlay
lapply(seq_along(seu.obj.list), function(i){#  
  # Prep sample data
  sampleName <- names(seu.obj.list)[i]
  seu.sub <- seu.obj.list[[i]]
  sub_dat <- dat %>% 
    filter(name == sample_decoder$short_name[match(sampleName, sample_decoder$sample)]) %>%
    left_join(GetTissueCoordinates(seu.sub), by = c("rowname" = "cell")) %>%
    select(-name) %>%
    relocate(x, y)
  
  # Make the background histo image
  p1 <- SpatialPlot(
    seu.sub, images = sampleName, image.scale = "lowres", features = NULL,
    group.by = NULL, pt.size.factor = 0, crop = TRUE, image.alpha = 0.35
  ) + 
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "pt"),
      panel.spacing = unit(0, "pt")
    ) +
    NoLegend()
  
  # Extract the and lift over the coordinate system
  base_plot <- p1[[1]]
  sub_dat <- sub_dat %>%
    select(-x, -y) %>%
    left_join(base_plot$data, by = c("rowname" = "cell"))

  # Calculate edges to crop and transpose the tissue to align with the pie charts
  crop_buffer <- 0
  ggbuild <- ggplot_build(base_plot)
  x_range <- range(sub_dat$y)
  y_range <- range(sub_dat$x)
  aspect_ratio <- diff(x_range) / diff(y_range)
  sub_dat$x <- (max(sub_dat$x) - sub_dat$x + min(sub_dat$x)) * aspect_ratio
  img_size <- max(diff(x_range), diff(range(sub_dat$x)))

  # Generate the plot
  p2 <- ggplot() + 
    scatterpie::geom_scatterpie(
      data = sub_dat, mapping = aes(x = y, y = x, group = rowname, r = (0.0075 * img_size)),
      lwd = 0.01, 
      legend_name = "CellTypes",
      cols = colnames(sub_dat)[! colnames(sub_dat) %in% c("rowname", "x", "y", "ident")]
    ) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0, "pt")
    ) + 
    NoLegend() +
    scale_fill_manual(values = cell_colz) +
    coord_equal(expand = TRUE)
  p_final <- p1 + inset_element(
    p2, left = 0, bottom = 0, right = 1, top = 1, align_to = "full"
  )
  ggsave(plot = p_final, paste0("../output/", outName, "/", sampleName, "_ct_pie.png"), height = 7, width = 7)    
})
