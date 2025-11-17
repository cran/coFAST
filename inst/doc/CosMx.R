## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  set.seed(2024) # set a random seed for reproducibility.
#  library(coFAST) # load the package of coFAST method
#  library(Seurat)
#  data(CosMx_subset)
#  CosMx_subset

## ----eval = FALSE-------------------------------------------------------------
#  CosMx_subset <- NormalizeData(CosMx_subset)

## ----eval = FALSE-------------------------------------------------------------
#  CosMx_subset <- FindVariableFeatures(CosMx_subset)

## ----eval = FALSE, fig.width= 6, fig.height= 4.5------------------------------
#  dat_cor <- diagnostic.cor.eigs(CosMx_subset)
#  q_est <- attr(dat_cor, "q_est")
#  cat("q_est = ", q_est, '\n')

## ----eval = FALSE-------------------------------------------------------------
#  pos <- as.matrix(CosMx_subset@meta.data[,c("x", "y")]) # Extract the spatial coordinates
#  Adj_sp <- AddAdj(pos) ## calculate the adjacency matrix
#  CosMx_subset <- coFAST(CosMx_subset, Adj_sp = Adj_sp, q = q_est)
#  CosMx_subset

## ----eval = FALSE-------------------------------------------------------------
#  CosMx_subset <- AddCluster(CosMx_subset)
#  print(head(CosMx_subset))
#  table(Idents(CosMx_subset))

## ----eval = FALSE,  fig.height = 4, fig.width=10------------------------------
#  CosMx_subset <- Addcoord2embed(CosMx_subset, coord.name = c("x", "y"))
#  print(CosMx_subset)
#  cols_cluster <- PRECAST::chooseColors(palettes_name = 'Classic 20', n_colors=20, plot_colors = TRUE)
#  # DimPlot(CosMx_subset, reduction='Spatial', cols=cols_cluster, pt.size = 2)
#  DimPlot(CosMx_subset, reduction='Spatial', group.by = c("cofast.cluster", 'cell_type'), cols=cols_cluster, pt.size = 1.5)

## ----eval = FALSE-------------------------------------------------------------
#  dat.spa.score <- AggregationScore(CosMx_subset, reduction.name = 'Spatial')
#  print(dat.spa.score)

## ----eval = FALSE-------------------------------------------------------------
#  dat.embd.score <- AggregationScore(CosMx_subset, reduction.name = 'cofast')
#  print(dat.embd.score)

## ----eval = FALSE-------------------------------------------------------------
#  CosMx_subset <- pdistance(CosMx_subset, reduction = "cofast")

## ----eval = FALSE-------------------------------------------------------------
#  print(table(Idents(CosMx_subset)))
#  #Idents(CosMx_subset) <- CosMx_subset$cell_type
#  df_sig_list <- find.signature.genes(CosMx_subset)
#  str(df_sig_list)

## ----eval = FALSE-------------------------------------------------------------
#  dat <- get.top.signature.dat(df_sig_list, ntop = 2, expr.prop.cutoff = 0.1)
#  head(dat)

## ----eval = FALSE, fig.width=10,fig.height=7----------------------------------
#  CosMx_subset <- coembedding_umap(
#    CosMx_subset, reduction = "cofast", reduction.name = "UMAP",
#    gene.set = unique(dat$gene))

## ----eval = FALSE, fig.width=8,fig.height=8-----------------------------------
#  ## choose beutifual colors
#  cols_cluster2 <- c("black", cols_cluster)
#  p1 <- coembed_plot(
#     CosMx_subset, reduction = "UMAP",
#     gene_txtdata = subset(dat, label=='8'),
#     cols=cols_cluster2, pt_text_size = 3)
#  p1

## ----eval = FALSE, fig.width=9,fig.height=6-----------------------------------
#  p2 <- coembed_plot(
#     CosMx_subset, reduction = "UMAP",
#     gene_txtdata = dat, cols=cols_cluster2,
#     pt_text_size = 3, alpha=0.2)
#  p2

## ----eval = FALSE, fig.width=9,fig.height=6-----------------------------------
#  DimPlot(CosMx_subset, reduction = 'UMAP', cols=cols_cluster)

## ----eval = FALSE, fig.width=8,fig.height=3.6---------------------------------
#  FeaturePlot(CosMx_subset, reduction = 'UMAP', features = c("PSCA", "CEACAM6"))

## -----------------------------------------------------------------------------
sessionInfo()

