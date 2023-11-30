R Figure 1
================

``` r
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(BSgenome.Drerio.UCSC.danRer11)
  library(EnhancedVolcano)
  library(ggsci)
  library(patchwork)
})
```

# Read data

``` r
HB13 <- readRDS(file = "../data/HB13hpf_neural.RDS")
DefaultAssay(HB13) <- "SCT"
Idents(HB13) <- "Clusters"
```

``` r
HB16 <- readRDS(file = "../data/HB16hpf_neural.RDS")
DefaultAssay(HB16) <- "SCT"
Idents(HB16) <- "Clusters"
```

# Cluster names

## HB13hpf

Named clusters after clustering at resolution 8, needed to resolve all
rhombomeres into separate clusters.

``` r
DimPlot(HB13, reduction = "wnn.umap") + scale_color_igv()
```

![](Figure1_files/figure-gfm/umap_HB13-1.png)<!-- --> Combining multiple
clusters of same cell type.

``` r
Idents(HB13) <- "Clusters"
HB13 <- RenameIdents(HB13,
                       "CHB.1" = "CHB",
                       "CHB.2" = "CHB",
                       "CHB.3" = "CHB",
                       "SC.1" = "SC",
                       "SC.2" = "SC",
                       "SC.3" = "SC",
                       "r5.1" = "r5",
                       "r5.2" = "r5",
                       "MB.1" = "MB",
                       "MB.2" = "MB",
                       "MB.3" = "MB",
                       "MHB.1" = "MHB",
                       "MHB.2" = "MHB",
                       "MHB.3" = "MHB",
                       "MHB.4" = "MHB",
                       "MHB.5" = "MHB",
                       "FB.1" = "FB",
                       "FB.2" = "FB",
                       "FB.3" = "FB",
                       "FB.4" = "FB")
```

    ## Warning: Cannot find identity r5.2

    ## Warning: Cannot find identity r5.1

``` r
levels(HB13) <- c("FB","MB","MHB","r1","r1 & r2","r2","r3","r4","r5","r6",
                    "low_expression","unknown","Neuron","Ciliated","CHB","SC")
umap.HB13 <- DimPlot(HB13, reduction = "wnn.umap") + scale_color_igv()
umap.HB13
```

![](Figure1_files/figure-gfm/rename_HB13-1.png)<!-- -->

## HB16hpf

Named clusters after clustering at resolution 6, needed to resolve all
rhombomeres into separate clusters.

``` r
DimPlot(HB16, reduction = "wnn.umap") + scale_color_igv()
```

![](Figure1_files/figure-gfm/umap_HB16-1.png)<!-- --> Combining multiple
clusters of same cell type.

``` r
Idents(HB16) <- "Clusters"
HB16 <- RenameIdents(HB16,
                       "CHB.1" = "CHB",
                       "CHB.2" = "CHB",
                       "CHB.3" = "CHB",
                       "CHB.4" = "CHB",
                       "SC.1" = "SC",
                       "SC.2" = "SC",
                       "SC.3" = "SC",
                       "r5.1" = "r5",
                       "r5.2" = "r5",
                       "MB.1" = "MB",
                       "MB.2" = "MB",
                       "MHB.1" = "MHB",
                       "MHB.2" = "MHB",
                       "MHB.3" = "MHB")
levels(HB16) <- c("FB","MB","MHB","r1","r2","r3","r4","r5","r6",
                    "DorsNT & NC","Neuron","Ciliated","CHB","SC")
umap.HB16 <- DimPlot(HB16, reduction = "wnn.umap") + scale_color_igv()
umap.HB16
```

![](Figure1_files/figure-gfm/rename_HB16-1.png)<!-- --> \# Gene
expression and chromvar activity plots

``` r
GetChromVarPlots <- function(object, motifs, reduction){
  DefaultAssay(object) <- "chromvar"
  mylist <- list()
  for(motif in 1:length(motifs)){
    p <- FeaturePlot(
      object = object,
      reduction = reduction,
      features = motifs[[motif]],
      min.cutoff = 'q10',
      max.cutoff = 2,
      pt.size = 1)  +
      ggtitle(paste(motifs[[motif]],names(motifs[motif]), sep = " ")) +
      theme(plot.title = element_text(size = 15)) + NoLegend()
    mylist[[motifs[[motif]]]] <- p
  }
  
  return(mylist)
}
```

``` r
GetChromVarPlotsPlusLegend <- function(object, motifs, reduction){
  DefaultAssay(object) <- "chromvar"
  mylist <- list()
  for(motif in 1:length(motifs)){
    p <- FeaturePlot(
      object = object,
      reduction = reduction,
      features = motifs[[motif]],
      min.cutoff = 'q10',
      max.cutoff = 2,
      pt.size = 1)  +
      ggtitle(paste(motifs[[motif]],names(motifs[motif]), sep = " ")) +
      theme(plot.title = element_text(size = 15))
    mylist[[motifs[[motif]]]] <- p
  }
  
  return(mylist)
}
```

``` r
# HB13 <- RunChromVAR(
#   object = HB13,
#   genome = BSgenome.Drerio.UCSC.danRer11,
#   assay = "peaks"
# )
# saveRDS(HB13, file = "../data/HB13hpf_neural.RDS")
```

``` r
# HB16 <- RunChromVAR(
#   object = HB16,
#   genome = BSgenome.Drerio.UCSC.danRer11,
#   assay = "peaks"
# )
# saveRDS(HB16, file = "../data/HB16hpf_neural.RDS")
```

## HB13hpf

``` r
motif1.names <- c("CDX4","HOXB9","RARA::RXRG","HNF1B","EGR2","Mafb")
motifs1 <- c("MA1473.1","MA1503.1","MA1149.1","MA0153.2","MA0472.2","MA0117.2")
names(motifs1) <- motif1.names
```

``` r
HB13.cdx4 <- FeaturePlot(HB13, features = "cdx4", reduction = "wnn.umap", max.cutoff = 1.3)
HB13.hoxb7a <- FeaturePlot(HB13, features = "hoxb7a", reduction = "wnn.umap", max.cutoff = 1.3)
HB13.hoxb9a <- FeaturePlot(HB13, features = "hoxb9a", reduction = "wnn.umap", max.cutoff = 1.3)
```

``` r
HB13.exp <- HB13.cdx4 + HB13.hoxb7a + HB13.hoxb9a + plot_layout(ncol = 3, guides = "collect")
HB13.exp
```

![](Figure1_files/figure-gfm/HB13_exp-1.png)<!-- -->

``` r
HB13.cv.motif.legend <- GetChromVarPlotsPlusLegend(HB13, motifs1, "wnn.umap")
HB13.cv.list <- GetChromVarPlots(HB13, motifs1, "wnn.umap")
wrap_plots(HB13.cv.motif.legend)
```

![](Figure1_files/figure-gfm/getchromvarplots-1.png)<!-- -->

``` r
HB13.motif1 <- HB13.cv.motif.legend[[1]] + HB13.cv.list[[2]] + plot_layout(ncol = 2, guides = "collect")
HB13.motif1
```

![](Figure1_files/figure-gfm/HB13_motif1-1.png)<!-- -->

``` r
HB13.motif2 <- HB13.cv.motif.legend[[3]] + HB13.cv.list[[4]] + HB13.cv.list[[5]] + HB13.cv.list[[6]] + 
  plot_layout(ncol = 4, guides = "collect")
HB13.motif2
```

![](Figure1_files/figure-gfm/HB13_motif2-1.png)<!-- -->

``` r
HB13.crabp2a <- FeaturePlot(HB13, features = "crabp2a", reduction = "wnn.umap", max.cutoff = 1.3)
HB13.rorcb <- FeaturePlot(HB13, features = "rorcb", reduction = "wnn.umap", max.cutoff = 1.3)
HB13.egr2b <- FeaturePlot(HB13, features = "egr2b", reduction = "wnn.umap", max.cutoff = 1.3)
HB13.mafba <- FeaturePlot(HB13, features = "mafba", reduction = "wnn.umap", max.cutoff = 1.3)
```

``` r
HB13.exp2 <- HB13.crabp2a + HB13.rorcb + HB13.egr2b + HB13.mafba + plot_layout(ncol = 4, guides = "collect")
HB13.exp2
```

![](Figure1_files/figure-gfm/HB13_exp2-1.png)<!-- -->

## HB16hpf

``` r
HB16.cv.motif.legend <- GetChromVarPlotsPlusLegend(HB16, motifs1, "wnn.umap")
HB16.cv.list <- GetChromVarPlots(HB16, motifs1, "wnn.umap")
wrap_plots(HB16.cv.motif.legend)
```

![](Figure1_files/figure-gfm/HB16.getchromvarplots-1.png)<!-- -->

``` r
HB16.motif <- HB16.cv.motif.legend[[1]] + HB16.cv.list[[2]]  + 
  plot_layout(ncol = 2, guides = "collect")
HB16.motif
```

![](Figure1_files/figure-gfm/HB16_motif-1.png)<!-- -->

``` r
HB16.crabp2a <- FeaturePlot(HB16, features = "cdx4", reduction = "wnn.umap", max.cutoff = 1.3)
HB16.rorcb <- FeaturePlot(HB16, features = "hoxb9a", reduction = "wnn.umap", max.cutoff = 1.3)
HB16.egr2b <- FeaturePlot(HB16, features = "crabp2a", reduction = "wnn.umap", max.cutoff = 1.3)
HB16.mafba <- FeaturePlot(HB16, features = "rorcb", reduction = "wnn.umap", max.cutoff = 1.3)
```

``` r
HB16.exp <- HB16.crabp2a + HB16.rorcb + plot_layout(ncol = 2, guides = "collect")
HB16.exp
```

![](Figure1_files/figure-gfm/HB16_exp1-1.png)<!-- -->

``` r
HB16.exp2 <- HB16.egr2b + HB16.mafba + plot_layout(ncol = 2, guides = "collect")
HB16.exp2
```

![](Figure1_files/figure-gfm/HB16_exp2-1.png)<!-- --> \# Figure 1
combined plots

``` r
layout <- "
AABBB
AACCC
DDDDD
EEEEE
FFGGG
HHIII"

combined1 <- umap.HB13 + HB13.exp + HB13.motif1 + HB13.exp2 + HB13.motif2 + 
  umap.HB16 + HB16.exp + HB16.motif + HB16.exp2 + 
  plot_layout(design = layout)
combined1
```

![](Figure1_files/figure-gfm/combined-1.png)<!-- -->

``` r
layout <- "
AAAAAABBBB
CCCCCCDDDD
EEEEEEEEEE
FFFFFFFFFF
GGGGGGHHHH
IIIIIIHHHH"

combined2 <- umap.HB13 + umap.HB16 + 
  HB13.exp + HB13.motif1 + 
  HB13.exp2 + 
  HB13.motif2 + 
  HB16.exp + (HB16.motif + plot_layout(ncol = 1, guides = "collect")) + 
  HB16.exp2 + 
  plot_layout(design = layout, heights = c(2,1,1,1,1,1)) +
  theme(title = element_text(size = 3))
combined2
```

![](Figure1_files/figure-gfm/combined2-1.png)<!-- -->

``` r
ggsave(filename = "Plots/Figure1.png", width = 15, height = 20, plot = combined2)
```

``` r
sessionInfo()
```

    ## R version 4.2.3 (2023-03-15)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Monterey 12.6.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.1.2                     ggsci_3.0.0                        
    ##  [3] EnhancedVolcano_1.16.0              ggrepel_0.9.3                      
    ##  [5] ggplot2_3.4.2                       BSgenome.Drerio.UCSC.danRer11_1.4.2
    ##  [7] BSgenome_1.66.3                     rtracklayer_1.58.0                 
    ##  [9] Biostrings_2.66.0                   XVector_0.38.0                     
    ## [11] GenomicRanges_1.50.2                GenomeInfoDb_1.34.9                
    ## [13] IRanges_2.32.0                      S4Vectors_0.36.2                   
    ## [15] BiocGenerics_0.44.0                 Signac_1.10.0                      
    ## [17] SeuratObject_4.1.3                  Seurat_4.3.0.1                     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fastmatch_1.1-3             systemfonts_1.0.4          
    ##   [3] plyr_1.8.8                  igraph_1.4.2               
    ##   [5] lazyeval_0.2.2              sp_1.6-0                   
    ##   [7] splines_4.2.3               BiocParallel_1.32.6        
    ##   [9] listenv_0.9.0               scattermore_1.0            
    ##  [11] digest_0.6.31               htmltools_0.5.5            
    ##  [13] fansi_1.0.4                 magrittr_2.0.3             
    ##  [15] tensor_1.5                  cluster_2.1.4              
    ##  [17] ROCR_1.0-11                 globals_0.16.2             
    ##  [19] matrixStats_0.63.0          spatstat.sparse_3.0-1      
    ##  [21] colorspace_2.1-0            textshaping_0.3.6          
    ##  [23] xfun_0.39                   dplyr_1.1.2                
    ##  [25] crayon_1.5.2                RCurl_1.98-1.12            
    ##  [27] jsonlite_1.8.4              progressr_0.13.0           
    ##  [29] spatstat.data_3.0-1         survival_3.5-5             
    ##  [31] zoo_1.8-12                  glue_1.6.2                 
    ##  [33] polyclip_1.10-4             gtable_0.3.3               
    ##  [35] zlibbioc_1.44.0             leiden_0.4.3               
    ##  [37] DelayedArray_0.24.0         future.apply_1.10.0        
    ##  [39] abind_1.4-5                 scales_1.2.1               
    ##  [41] spatstat.random_3.1-4       miniUI_0.1.1.1             
    ##  [43] Rcpp_1.0.10                 viridisLite_0.4.2          
    ##  [45] xtable_1.8-4                reticulate_1.28            
    ##  [47] htmlwidgets_1.6.2           httr_1.4.6                 
    ##  [49] RColorBrewer_1.1-3          ellipsis_0.3.2             
    ##  [51] ica_1.0-3                   farver_2.1.1               
    ##  [53] pkgconfig_2.0.3             XML_3.99-0.14              
    ##  [55] uwot_0.1.14                 deldir_1.0-6               
    ##  [57] utf8_1.2.3                  labeling_0.4.2             
    ##  [59] tidyselect_1.2.0            rlang_1.1.1                
    ##  [61] reshape2_1.4.4              later_1.3.1                
    ##  [63] munsell_0.5.0               tools_4.2.3                
    ##  [65] cli_3.6.1                   generics_0.1.3             
    ##  [67] ggridges_0.5.4              evaluate_0.21              
    ##  [69] stringr_1.5.0               fastmap_1.1.1              
    ##  [71] ragg_1.2.5                  yaml_2.3.7                 
    ##  [73] goftest_1.2-3               knitr_1.42                 
    ##  [75] fitdistrplus_1.1-11         purrr_1.0.1                
    ##  [77] RANN_2.6.1                  pbapply_1.7-0              
    ##  [79] future_1.32.0               nlme_3.1-162               
    ##  [81] mime_0.12                   RcppRoll_0.3.0             
    ##  [83] compiler_4.2.3              rstudioapi_0.14            
    ##  [85] plotly_4.10.1               png_0.1-8                  
    ##  [87] spatstat.utils_3.0-2        tibble_3.2.1               
    ##  [89] stringi_1.7.12              highr_0.10                 
    ##  [91] lattice_0.21-8              Matrix_1.6-1.1             
    ##  [93] vctrs_0.6.2                 pillar_1.9.0               
    ##  [95] lifecycle_1.0.3             spatstat.geom_3.1-0        
    ##  [97] lmtest_0.9-40               RcppAnnoy_0.0.20           
    ##  [99] data.table_1.14.8           cowplot_1.1.1              
    ## [101] bitops_1.0-7                irlba_2.3.5.1              
    ## [103] httpuv_1.6.9                R6_2.5.1                   
    ## [105] BiocIO_1.8.0                promises_1.2.0.1           
    ## [107] KernSmooth_2.23-21          gridExtra_2.3              
    ## [109] parallelly_1.35.0           codetools_0.2-19           
    ## [111] MASS_7.3-60                 SummarizedExperiment_1.28.0
    ## [113] rjson_0.2.21                withr_2.5.0                
    ## [115] GenomicAlignments_1.34.1    sctransform_0.3.5          
    ## [117] Rsamtools_2.14.0            GenomeInfoDbData_1.2.9     
    ## [119] parallel_4.2.3              grid_4.2.3                 
    ## [121] tidyr_1.3.0                 rmarkdown_2.21             
    ## [123] MatrixGenerics_1.10.0       Rtsne_0.16                 
    ## [125] spatstat.explore_3.1-0      Biobase_2.58.0             
    ## [127] shiny_1.7.4                 restfulr_0.0.15
