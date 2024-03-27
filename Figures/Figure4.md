R Figure 4
================

# 1 libraries

``` r
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  library(dittoSeq)
  library(ComplexHeatmap)
  library(openxlsx)
})
```

    ## Warning: package 'dittoSeq' was built under R version 4.3.3

    ## Warning: package 'ComplexHeatmap' was built under R version 4.3.1

``` r
options(future.globals.maxSize = 4000 * 1024^2)
```

# 2 read data

``` r
HB16hpf <- readRDS("~/Documents/Projects/Sagerstrom/scRNA-seq_ATAC-seq_reanalysis_2021_07_12/Integratation/IntegrationOfIndependentNeuralSubsets/RDSfiles/HB16hpf.neural.int_peaks.assay.RDS")
HB13hpf <- readRDS("~/Documents/Projects/Sagerstrom/scRNA-seq_ATAC-seq_reanalysis_2021_07_12/Integratation/IntegrationOfIndependentNeuralSubsets/RDSfiles/HB13hpf.neural.int_peaks.assay.RDS")
```

``` r
DefaultAssay(HB13hpf) <- "chromvar"
DefaultAssay(HB16hpf) <- "chromvar"
HB13hpf@assays[["ATAC"]] <- NULL
HB13hpf@assays[["int_peaks"]] <- NULL
HB16hpf@assays[["ATAC"]] <- NULL
HB16hpf@assays[["int_peaks"]] <- NULL
```

combine CaudHB clusters, SC clusters and HB16hpf r5.1 and r5.2 clusters

``` r
HB13hpf$Clusters <- as.character(HB13hpf$Clusters)
HB13hpf$Clusters[HB13hpf$Clusters %in% c("SC.1","SC.2","SC.3")] <- "SC"
HB13hpf$Clusters[HB13hpf$Clusters %in% c("CaudHB.1","CaudHB.2","CaudHB.3")] <- "CHB"
HB13hpf$Clusters <- as.factor(HB13hpf$Clusters)

HB16hpf$Clusters <- as.character(HB16hpf$Clusters)
HB16hpf$Clusters[HB16hpf$Clusters %in% c("r5.1","r5.2")] <- "r5"
HB16hpf$Clusters[HB16hpf$Clusters %in% c("SC.1","SC.2","SC.3")] <- "SC"
HB16hpf$Clusters[HB16hpf$Clusters %in% c("CaudHB.1","CaudHB.2","CaudHB.3","CaudHB.4")] <- "CHB"
HB16hpf$Clusters <- as.factor(HB16hpf$Clusters)
```

``` r
p1 <- DimPlot(HB13hpf, reduction = "wnn.umap", group.by = "Clusters") + scale_color_igv() + ggtitle("HB13hpf")
p2 <- DimPlot(HB16hpf, reduction = "wnn.umap", group.by = "Clusters") + scale_color_igv() + ggtitle("HB16hpf")
p1 + p2 
```

![](Figure4_files/figure-gfm/dimplot-1.png)<!-- -->

``` r
HB.comb <- merge(HB13hpf, HB16hpf, add.cell.ids = c("HB13hpf","HB16hpf"), project = "Merge_13_16")
```

    ## 
    ## Binding matrix rows
    ## 
    ## Binding matrix rows
    ## 
    ## Binding matrix rows
    ## 
    ## Binding matrix rows

``` r
Idents(HB.comb) <- "Clusters"
HB.comb <- subset(HB.comb, idents = c("r6","SC","CHB"))
HB.comb$Clusters <- factor(HB.comb$Clusters, levels = c("r6","CHB","SC"))
HB.comb$Sample <- HB.comb$orig.ident
cells.rhb = WhichCells(HB.comb, idents = c("r6","CHB","SC"))
```

# 3 Get DE gene lists for various clusters

## 3.1 SC from HB13hpf

``` r
DefaultAssay(HB13hpf) <- "SCT"
Idents(HB13hpf) <- "Clusters"
SC.13.markers <- FindMarkers(HB13hpf, ident.1 = "SC", ident.2 = c("CHB","r6"), only.pos = TRUE, verbose = F)  
SC.13.markers$gene <- rownames(SC.13.markers) 
sprintf("Total number of DE genes = %d",nrow(SC.13.markers))
```

    ## [1] "Total number of DE genes = 344"

``` r
sprintf("Number of DE genes with p_val_adj < 0.05 = %d",nrow(SC.13.markers[SC.13.markers$p_val_adj < 0.05,]))
```

    ## [1] "Number of DE genes with p_val_adj < 0.05 = 89"

## 3.2 CHB from HB13hpf

``` r
CHB.13.markers <- FindMarkers(HB13hpf, ident.1 = "CHB", ident.2 = c("SC","r6"), only.pos = TRUE, verbose = F)  
CHB.13.markers$gene <- rownames(CHB.13.markers) 
sprintf("Total number of DE genes = %d",nrow(CHB.13.markers))
```

    ## [1] "Total number of DE genes = 134"

``` r
sprintf("Number of DE genes with p_val_adj < 0.05 = %d",nrow(CHB.13.markers[CHB.13.markers$p_val_adj < 0.05,]))
```

    ## [1] "Number of DE genes with p_val_adj < 0.05 = 7"

## 3.3 r6 from HB13hpf

``` r
r6.13.markers <- FindMarkers(HB13hpf, ident.1 = "r6", ident.2 = c("SC","CHB"), only.pos = TRUE, verbose = F)   
r6.13.markers$gene <- rownames(r6.13.markers) 
sprintf("Total number of DE genes = %d",nrow(r6.13.markers))
```

    ## [1] "Total number of DE genes = 344"

``` r
sprintf("Number of DE genes with p_val_adj < 0.05 = %d",nrow(r6.13.markers[r6.13.markers$p_val_adj < 0.05,]))
```

    ## [1] "Number of DE genes with p_val_adj < 0.05 = 25"

## 3.4 SC from HB16hpf

``` r
DefaultAssay(HB16hpf) <- "SCT"
Idents(HB16hpf) <- "Clusters"
SC.16.markers <- FindMarkers(HB16hpf, ident.1 = "SC", ident.2 = c("CHB","r6"), only.pos = TRUE, verbose = F)   
SC.16.markers$gene <- rownames(SC.16.markers) 
sprintf("Total number of DE genes = %d",nrow(SC.16.markers))
```

    ## [1] "Total number of DE genes = 254"

``` r
sprintf("Number of DE genes with p_val_adj < 0.05 = %d",nrow(SC.16.markers[SC.16.markers$p_val_adj < 0.05,]))
```

    ## [1] "Number of DE genes with p_val_adj < 0.05 = 105"

## 3.5 CHB from HB16hpf

``` r
CHB.16.markers <- FindMarkers(HB16hpf, ident.1 = "CHB", ident.2 = c("SC","r6"), only.pos = TRUE, verbose = F)   
CHB.16.markers$gene <- rownames(CHB.16.markers) 
sprintf("Total number of DE genes = %d",nrow(CHB.16.markers))
```

    ## [1] "Total number of DE genes = 124"

``` r
sprintf("Number of DE genes with p_val_adj < 0.05 = %d",nrow(CHB.16.markers[CHB.16.markers$p_val_adj < 0.05,]))
```

    ## [1] "Number of DE genes with p_val_adj < 0.05 = 25"

## 3.6 r6 from HB16hpf

``` r
r6.16.markers <- FindMarkers(HB16hpf, ident.1 = "r6", ident.2 = c("SC","CHB"), only.pos = TRUE, verbose = F)   
r6.16.markers$gene <- rownames(r6.16.markers) 
sprintf("Total number of DE genes = %d",nrow(r6.16.markers))
```

    ## [1] "Total number of DE genes = 335"

``` r
sprintf("Number of DE genes with p_val_adj < 0.05 = %d",nrow(r6.16.markers[r6.16.markers$p_val_adj < 0.05,]))
```

    ## [1] "Number of DE genes with p_val_adj < 0.05 = 50"

# 4. Heatmaps

``` r
hm_list <- list()
```

## 4.1 SC DE genes

``` r
SCgenes <- unique(c(SC.13.markers$gene[SC.13.markers$p_val_adj <0.05], SC.16.markers$gene[SC.16.markers$p_val_adj <0.05]))
DefaultAssay(HB.comb) <- "SCT"
hm_SC <- dittoHeatmap(HB.comb, 
                             genes = SCgenes, 
                             cells = cells.rhb, 
                             assay = "SCT", #slot = "data", 
                             annot.by = c("Sample","Clusters"), 
                             order.by = c("Sample","Clusters"), 
                             scaled.to.max = T, 
                             show_rownames = T,
                             cluster_rows = F,
                             drop_levels = T,
                      complex = T)
hm_SC@row_names_param[["gp"]][["fontsize"]] = 8
hm_SC
```

![](Figure4_files/figure-gfm/hmSC-1.png)<!-- -->

``` r
hm_list[["SC"]] <- hm_SC %>%
  draw(show_heatmap_legend = F, show_annotation_legend = F) %>%
  grid.grabExpr()
```

``` r
png(filename = "Plots/SC_heatmap_HB13_HB16.png", width=7,height=12.03,units="in",res=1200)
hm_SC <- dittoHeatmap(HB.comb, 
                      genes = SCgenes, 
                      cells = cells.rhb, 
                      assay = "SCT", #slot = "data", 
                      annot.by = c("Sample","Clusters"), 
                      order.by = c("Sample","Clusters"), 
                      scaled.to.max = T, 
                      show_rownames = T,
                      cluster_rows = F,
                      drop_levels = T,
                      complex = T)
hm_SC@row_names_param[["gp"]][["fontsize"]] = 8
hm_SC
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## 4.2 CHB DE genes

``` r
CHBgenes <- unique(c(CHB.13.markers$gene[CHB.13.markers$p_val_adj <0.05],
                     CHB.16.markers$gene[CHB.16.markers$p_val_adj <0.05]))
DefaultAssay(HB.comb) <- "SCT"
hm_CHB <- dittoHeatmap(HB.comb, 
                             genes = CHBgenes, 
                             cells = cells.rhb, 
                             assay = "SCT", #slot = "data", 
                             annot.by = c("Sample","Clusters"), 
                             order.by = c("Sample","Clusters"), 
                             scaled.to.max = T, 
                             show_rownames = T,
                             cluster_rows = F,
                             drop_levels = T,
                       complex = T)
hm_CHB@row_names_param[["gp"]][["fontsize"]] = 8
hm_CHB
```

![](Figure4_files/figure-gfm/hmCHB-1.png)<!-- -->

``` r
hm_list[["CHB"]] <- hm_CHB %>%
  draw(show_heatmap_legend = F, show_annotation_legend = F) %>%
  grid.grabExpr()
```

``` r
png(filename = "Plots/CHB_heatmap_HB13_HB16.png", width=7,height=2.68,units="in",res=1200)
hm_CHB <- dittoHeatmap(HB.comb, 
                       genes = CHBgenes, 
                       cells = cells.rhb, 
                       assay = "SCT", #slot = "data", 
                       annot.by = c("Sample","Clusters"), 
                       order.by = c("Sample","Clusters"), 
                       scaled.to.max = T, 
                       show_rownames = T,
                       cluster_rows = F,
                       drop_levels = T,
                       complex = T)
hm_CHB@row_names_param[["gp"]][["fontsize"]] = 8
hm_CHB
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## 4.3 r6 DE genes

``` r
r6genes <- unique(c(r6.13.markers$gene[r6.13.markers$p_val_adj <0.05], r6.16.markers$gene[r6.16.markers$p_val_adj <0.05]))
DefaultAssay(HB.comb) <- "SCT"
hm_r6 <- dittoHeatmap(HB.comb,
                             genes = r6genes,
                             cells = cells.rhb,
                             assay = "SCT", #slot = "data",
                             annot.by = c("Sample","Clusters"),
                             order.by = c("Sample","Clusters"),
                             scaled.to.max = T,
                             show_rownames = T,
                             cluster_rows = F,
                      complex = T)
hm_r6@row_names_param[["gp"]][["fontsize"]] = 8
hm_r6
```

![](Figure4_files/figure-gfm/hmr6-1.png)<!-- -->

``` r
hm_list[["r6"]] <- hm_r6 %>%
  draw(show_heatmap_legend = F, show_annotation_legend = F) %>%
  grid.grabExpr()
```

``` r
png(filename = "Plots/r6_heatmap_HB13_HB16.png", width=7,height=5.71,units="in",res=1200)
hm_r6 <- dittoHeatmap(HB.comb,
                      genes = r6genes,
                      cells = cells.rhb,
                      assay = "SCT", #slot = "data",
                      annot.by = c("Sample","Clusters"),
                      order.by = c("Sample","Clusters"),
                      scaled.to.max = T,
                      show_rownames = T,
                      cluster_rows = F,
                      complex = T)
hm_r6@row_names_param[["gp"]][["fontsize"]] = 8
hm_r6
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## 4.4 all SC, CHB, r6 DE genes together

``` r
allgenes <- unique(c(SC.13.markers$gene[SC.13.markers$p_val_adj <0.05], 
                     SC.16.markers$gene[SC.16.markers$p_val_adj <0.05],
                     CHB.13.markers$gene[CHB.13.markers$p_val_adj <0.05],
                     CHB.16.markers$gene[CHB.16.markers$p_val_adj <0.05],
                     r6.13.markers$gene[r6.13.markers$p_val_adj <0.05],
                     r6.16.markers$gene[r6.16.markers$p_val_adj <0.05]))
```

``` r
DefaultAssay(HB.comb) <- "SCT"
comb_heatmap <- dittoHeatmap(HB.comb, 
                             genes = allgenes, 
                             cells = cells.rhb, 
                             assay = "SCT", #slot = "data", 
                             annot.by = c("Sample","Clusters"), 
                             order.by = c("Sample","Clusters"), 
                             scaled.to.max = T, 
                             show_rownames = F,
                             clustering_method = "ward.D2")
```

![](Figure4_files/figure-gfm/hm_comb-1.png)<!-- -->

Use h=14 clustering of genes

``` r
as.data.frame(sort(cutree(comb_heatmap$tree_row, h=14)))
```

    ##                  sort(cutree(comb_heatmap$tree_row, h = 14))
    ## hoxc6b                                                     1
    ## cdx4                                                       1
    ## tanc2a                                                     1
    ## XKR4                                                       1
    ## ldlrad2                                                    1
    ## mllt3                                                      1
    ## col5a2a                                                    1
    ## apcdd1l                                                    1
    ## aopep                                                      1
    ## dnmt3ab                                                    1
    ## RIMBP2                                                     1
    ## nfia                                                       1
    ## nradd                                                      1
    ## CU639469.1                                                 1
    ## adamts3                                                    1
    ## arhgap29a                                                  1
    ## col28a2a                                                   1
    ## hoxd3a                                                     1
    ## spred2b                                                    1
    ## arhgef10                                                   1
    ## fam171a2a                                                  1
    ## pdzd2                                                      1
    ## ece2b                                                      1
    ## nrp2b                                                      1
    ## bcar3                                                      1
    ## madd                                                       1
    ## tnfrsf19                                                   1
    ## col14a1a                                                   1
    ## dachb                                                      1
    ## sox5                                                       1
    ## epha4a                                                     1
    ## rasgrp3                                                    1
    ## oc90                                                       1
    ## si:dkey-91m11.5                                            1
    ## bicc2                                                      1
    ## ptprz1b                                                    1
    ## pcdh2g28                                                   1
    ## col11a1a                                                   1
    ## ptprfa                                                     1
    ## abr                                                        1
    ## p3h2                                                       1
    ## map2                                                       1
    ## stard13b                                                   1
    ## tinagl1                                                    1
    ## hpca                                                       1
    ## itgb4                                                      1
    ## pax6b                                                      1
    ## adgrl3.1                                                   1
    ## hoxc3a                                                     2
    ## raraa                                                      2
    ## pkd1b                                                      2
    ## zgc:110158                                                 2
    ## wnt4                                                       2
    ## phip                                                       2
    ## evx1                                                       2
    ## prickle2b                                                  2
    ## igsf9a                                                     2
    ## ehbp1                                                      2
    ## mdka                                                       2
    ## serinc5                                                    2
    ## robo3                                                      2
    ## kaznb                                                      2
    ## phldb1b                                                    2
    ## notch3                                                     2
    ## antxr1c                                                    2
    ## hspg2                                                      2
    ## ssbp4                                                      2
    ## nectin1b                                                   2
    ## hoxa4a                                                     2
    ## yap1                                                       2
    ## tspan17                                                    2
    ## bmpr1ab                                                    2
    ## bcar1                                                      2
    ## nhsl1b                                                     2
    ## brsk2b                                                     2
    ## hoxb7a                                                     3
    ## kif26ab                                                    3
    ## hoxb9a                                                     3
    ## stm                                                        3
    ## plod2                                                      3
    ## arid3c                                                     3
    ## hoxb6b                                                     3
    ## hoxa9a                                                     3
    ## tcirg1b                                                    3
    ## phc2a                                                      3
    ## megf6b                                                     3
    ## hoxc8a                                                     3
    ## hoxb6a                                                     3
    ## kcnn3                                                      3
    ## esrrga                                                     3
    ## ly6pge                                                     3
    ## ppfia3                                                     3
    ## BX936418.1                                                 3
    ## rdh8a                                                      3
    ## lypd6b                                                     3
    ## cdkl1                                                      3
    ## LO018340.1                                                 3
    ## hoxa10b                                                    3
    ## angptl2b                                                   3
    ## itih6                                                      3
    ## rhoj                                                       3
    ## hoxb5b                                                     3
    ## pxdn                                                       3
    ## gdf6b                                                      3
    ## etv4                                                       3
    ## chrd                                                       3
    ## CABZ01060373.1                                             3
    ## BX640512.3                                                 3
    ## hoxb1b                                                     3
    ## nr6a1b                                                     3
    ## ralgps1                                                    3
    ## CR382300.2                                                 3
    ## skap1                                                      3
    ## si:dkey-25g12.4                                            3
    ## npr1a                                                      3
    ## necab2                                                     3
    ## dip2a                                                      3
    ## il17rd                                                     3
    ## mecom                                                      3
    ## stx1b                                                      3
    ## adgrd2                                                     3
    ## hoxb5a                                                     3
    ## prex2                                                      3
    ## adamts9                                                    3
    ## col5a3a                                                    3
    ## kdrl                                                       3
    ## phldb1a                                                    3
    ## kif13bb                                                    3
    ## thbs3b                                                     3
    ## sema3ga                                                    3
    ## rarga                                                      3
    ## hoxb3a                                                     4
    ## plxna3                                                     4
    ## col18a1a                                                   4
    ## macf1a                                                     4
    ## notch1b                                                    4
    ## jmjd1cb                                                    4
    ## foxp4                                                      4
    ## dag1                                                       4
    ## greb1l                                                     4
    ## zbtb16b                                                    4
    ## nr2f5                                                      4
    ## epha4b                                                     4
    ## zbtb16a                                                    4
    ## fgfr2                                                      4
    ## cdon                                                       5
    ## rorcb                                                      5
    ## dusp4                                                      5
    ## adgrl1a                                                    5
    ## mycn                                                       5
    ## crabp2a                                                    5
    ## nxph3                                                      5
    ## rxrgb                                                      5
    ## prkcdb                                                     5
    ## bicd1a                                                     5
    ## igdcc3                                                     5
    ## fam49al                                                    5
    ## meis3                                                      5
    ## SLC39A11                                                   5
    ## mafba                                                      5
    ## cyp26c1                                                    5
    ## cxcl12a                                                    5
    ## fgfrl1a                                                    5
    ## epha7                                                      5
    ## abcc8                                                      5
    ## csmd2                                                      5
    ## ptn                                                        5
    ## magi2a                                                     5
    ## prox2                                                      5
    ## palm1b                                                     5
    ## prdx6                                                      5
    ## samd10b                                                    5
    ## si:dkey-157g16.6                                           5
    ## si:ch73-252p3.1                                            5
    ## hoxb2a                                                     5
    ## fscn1a                                                     5
    ## mmp2                                                       5
    ## nhsl1a                                                     5
    ## mpz                                                        5
    ## zgc:158328                                                 5
    ## adgrg6                                                     5
    ## col15a1b                                                   5
    ## bmpr1ba                                                    5
    ## emid1                                                      5
    ## CT990561.1                                                 5
    ## smoc1                                                      5
    ## sfxn1                                                      5
    ## dchs1a                                                     5
    ## kif21b                                                     5
    ## igf2b                                                      5
    ## efnb1                                                      5
    ## lmf2b                                                      5
    ## arhgap39                                                   5
    ## zeb2b                                                      5
    ## mgat4c                                                     5
    ## olfm1a                                                     5
    ## wnt7bb                                                     5
    ## nppc                                                       5
    ## si:ch211-214c7.4                                           5
    ## slc4a4a                                                    5
    ## hcn2b                                                      5
    ## tgfbr3                                                     5
    ## ror1                                                       6
    ## hs6st1a                                                    6
    ## zfhx4                                                      6
    ## rxraa                                                      6
    ## tenm3                                                      6
    ## tenm4                                                      6
    ## adgrl2a                                                    6
    ## bahcc1b                                                    6
    ## ncoa3                                                      6
    ## fgfr3                                                      6
    ## efnb3b                                                     6
    ## jarid2b                                                    6
    ## jag2b                                                      6
    ## boc                                                        6
    ## tmem131l                                                   6
    ## si:ch73-386h18.1                                           6
    ## rasal2                                                     6
    ## ptprna                                                     7
    ## dbn1                                                       7
    ## skia                                                       7
    ## shroom4                                                    7
    ## nck2a                                                      7
    ## nr2f2                                                      7
    ## kdm5bb                                                     7
    ## sema4c                                                     7
    ## tox3                                                       7
    ## dtnba                                                      7
    ## kirrel3l                                                   7
    ## bach2b                                                     7

``` r
comb.split = data.frame(cutree(comb_heatmap$tree_row, h=14)) %>%
  dplyr::rename("hm_cluster" = 1) 
```

``` r
hm_all <- dittoHeatmap(HB.comb, genes = allgenes,
             cells = cells.rhb, assay = "SCT", slot = "data",
             annot.by = c("Sample","Clusters"), order.by = c("Sample","Clusters"), scaled.to.max = T, show_rownames = T,
             row_split = comb.split,complex = T,name = "chromVar\n scaled",row_title_gp = gpar(fontsize = 8),
             height = 39*unit(5, "mm"))
```

    ## Warning: argument `height` is not supported in pheatmap -> Heatmap translation,
    ## skip it.

``` r
hm_all@row_names_param[["gp"]][["fontsize"]] = 8
hm_all
```

![](Figure4_files/figure-gfm/hmall-1.png)<!-- -->

``` r
hm_list[["all"]] <- hm_all %>%
  draw(align_heatmap_legend = "heatmap_top") %>%
  grid.grabExpr()
```

``` r
png(filename = "Plots/SC_CHB_r6_heatmap_HB13_HB16.png", width=7,height=20,units="in",res=1200)
hm_all <- dittoHeatmap(HB.comb, genes = allgenes,
             cells = cells.rhb, assay = "SCT", slot = "data",
             annot.by = c("Sample","Clusters"), order.by = c("Sample","Clusters"), scaled.to.max = T, show_rownames = T,
             row_split = comb.split,complex = T,name = "chromVar\n scaled",row_title_gp = gpar(fontsize = 6),
             height = 39*unit(5, "mm"))
```

    ## Warning: argument `height` is not supported in pheatmap -> Heatmap translation,
    ## skip it.

``` r
hm_all@row_names_param[["gp"]][["fontsize"]] = 8
hm_all
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
#wrap_plots(hm_list, ncol = 4, heights = c(4.48,1,2.13,7.45))
p1 <- wrap_plots(hm_list[["r6"]]) / plot_spacer() + plot_layout(heights = c(2.13,5.32))
p2 <- wrap_plots(hm_list[["CHB"]]) / plot_spacer() + plot_layout(heights = c(1,6.45))
p3 <- wrap_plots(hm_list[["SC"]]) / plot_spacer() + plot_layout(heights = c(4.48,2.97))
p4 <- wrap_plots(hm_list[["all"]]) / plot_spacer() + plot_layout(heights = c(7.45,0.01))
plot_spacer() +  p1 + p2 + p3 + p4 + plot_layout(ncol = 5, widths = c(0.01,1,1,1,1.5), guides = 'collect')
```

![](Figure4_files/figure-gfm/hm_comb_plots-1.png)<!-- -->

``` r
#p2 + p3 + plot_layout(ncol = 2)
```

``` r
HB13.vp <- readRDS(file = "../VolcanoPlots/Plots/HB13_volcanoPlots.RDS")
HB16.vp <- readRDS(file = "../VolcanoPlots/Plots/HB16_volcanoPlots.RDS")
```

``` r
p1 <- wrap_plots(hm_list[["r6"]]) / plot_spacer() + plot_layout(heights = c(2.13,2.35))
p2 <- wrap_plots(hm_list[["CHB"]]) / plot_spacer() + plot_layout(heights = c(1,3.48))
p3 <- wrap_plots(hm_list[["SC"]]) / plot_spacer() + plot_layout(heights = c(4.48,0.01))
p4 <- wrap_plots(hm_list[["all"]]) / plot_spacer() + plot_layout(heights = c(7.45,0.01))

layout = "
ABCDE
AFFFE
AGGGE"

combined_plot <- plot_spacer() + p1 + p2 + p3 + p4 + HB13.vp + HB16.vp + plot_layout(design = layout, heights = c(4.48,1.48,1.48), widths = c(0.01,1,1.1,0.7,1.6))
combined_plot
```

![](Figure4_files/figure-gfm/figure4-1.png)<!-- -->

``` r
ggsave(filename = "Plots/Figure4.png", width = 28, height = 30, plot = combined_plot)
```

``` r
DEgenelist <- list()
DEgenelist[["HB13_SCvsCHBandr6"]] <- SC.13.markers[SC.13.markers$p_val_adj <0.05,]
DEgenelist[["HB13_CHBvsSCandr6"]] <- CHB.13.markers[CHB.13.markers$p_val_adj <0.05,]
DEgenelist[["HB13_r6vsSCandCHB"]] <- r6.13.markers[r6.13.markers$p_val_adj <0.05,]
DEgenelist[["HB16_SCvsCHBandr6"]] <- SC.16.markers[SC.16.markers$p_val_adj <0.05,]
DEgenelist[["HB16_CHBvsSCandr6"]] <- CHB.16.markers[CHB.16.markers$p_val_adj <0.05,]
DEgenelist[["HB16_r6vsSCandCHB"]] <- r6.16.markers[r6.16.markers$p_val_adj <0.05,]
write.xlsx(DEgenelist, file = "Markers/Fig4_SuppFig1_DEgenesTables_forHeatmaps.xlsx", rowNames = T)
```

``` r
sessionInfo()
```

    ## R version 4.3.0 (2023-04-21)
    ## Platform: x86_64-apple-darwin20 (64-bit)
    ## Running under: macOS Monterey 12.6.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Denver
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] openxlsx_4.2.5.2      ComplexHeatmap_2.18.0 dittoSeq_1.14.3      
    ##  [4] patchwork_1.2.0       ggsci_3.0.0           ggplot2_3.4.4        
    ##  [7] dplyr_1.1.4           Signac_1.10.0         SeuratObject_4.1.3   
    ## [10] Seurat_4.3.0.1       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppAnnoy_0.0.22            splines_4.3.0              
    ##   [3] later_1.3.2                 bitops_1.0-7               
    ##   [5] tibble_3.2.1                polyclip_1.10-6            
    ##   [7] lifecycle_1.0.4             doParallel_1.0.17          
    ##   [9] globals_0.16.2              lattice_0.22-5             
    ##  [11] MASS_7.3-60.0.1             magrittr_2.0.3             
    ##  [13] limma_3.58.1                plotly_4.10.4              
    ##  [15] rmarkdown_2.25              yaml_2.3.8                 
    ##  [17] httpuv_1.6.13               sctransform_0.4.1          
    ##  [19] zip_2.3.1                   sp_2.1-2                   
    ##  [21] spatstat.sparse_3.0-3       EnhancedVolcano_1.20.0     
    ##  [23] reticulate_1.34.0           cowplot_1.1.3              
    ##  [25] pbapply_1.7-2               RColorBrewer_1.1-3         
    ##  [27] abind_1.4-5                 zlibbioc_1.48.0            
    ##  [29] Rtsne_0.17                  GenomicRanges_1.54.1       
    ##  [31] purrr_1.0.2                 BiocGenerics_0.48.1        
    ##  [33] RCurl_1.98-1.14             circlize_0.4.16            
    ##  [35] GenomeInfoDbData_1.2.11     IRanges_2.36.0             
    ##  [37] S4Vectors_0.40.2            ggrepel_0.9.5              
    ##  [39] irlba_2.3.5.1               listenv_0.9.0              
    ##  [41] spatstat.utils_3.0-4        pheatmap_1.0.12            
    ##  [43] goftest_1.2-3               spatstat.random_3.2-2      
    ##  [45] fitdistrplus_1.1-11         parallelly_1.36.0          
    ##  [47] leiden_0.4.3.1              codetools_0.2-19           
    ##  [49] DelayedArray_0.28.0         RcppRoll_0.3.0             
    ##  [51] tidyselect_1.2.0            shape_1.4.6.1              
    ##  [53] farver_2.1.1                matrixStats_1.2.0          
    ##  [55] stats4_4.3.0                spatstat.explore_3.2-5     
    ##  [57] jsonlite_1.8.8              GetoptLong_1.0.5           
    ##  [59] ellipsis_0.3.2              progressr_0.14.0           
    ##  [61] ggridges_0.5.6              survival_3.5-7             
    ##  [63] iterators_1.0.14            systemfonts_1.0.6          
    ##  [65] foreach_1.5.2               tools_4.3.0                
    ##  [67] ragg_1.3.0                  ica_1.0-3                  
    ##  [69] Rcpp_1.0.12                 glue_1.7.0                 
    ##  [71] gridExtra_2.3               SparseArray_1.2.3          
    ##  [73] xfun_0.41                   MatrixGenerics_1.14.0      
    ##  [75] GenomeInfoDb_1.38.5         withr_3.0.0                
    ##  [77] fastmap_1.1.1               fansi_1.0.6                
    ##  [79] digest_0.6.34               R6_2.5.1                   
    ##  [81] mime_0.12                   textshaping_0.3.7          
    ##  [83] colorspace_2.1-0            scattermore_1.2            
    ##  [85] tensor_1.5                  spatstat.data_3.0-4        
    ##  [87] utf8_1.2.4                  tidyr_1.3.1                
    ##  [89] generics_0.1.3              data.table_1.14.10         
    ##  [91] httr_1.4.7                  htmlwidgets_1.6.4          
    ##  [93] S4Arrays_1.2.0              uwot_0.1.16                
    ##  [95] pkgconfig_2.0.3             gtable_0.3.4               
    ##  [97] lmtest_0.9-40               SingleCellExperiment_1.24.0
    ##  [99] XVector_0.42.0              htmltools_0.5.7            
    ## [101] clue_0.3-65                 scales_1.3.0               
    ## [103] Biobase_2.62.0              png_0.1-8                  
    ## [105] knitr_1.45                  rstudioapi_0.15.0          
    ## [107] reshape2_1.4.4              rjson_0.2.21               
    ## [109] nlme_3.1-164                zoo_1.8-12                 
    ## [111] GlobalOptions_0.1.2         stringr_1.5.1              
    ## [113] KernSmooth_2.23-22          parallel_4.3.0             
    ## [115] miniUI_0.1.1.1              pillar_1.9.0               
    ## [117] vctrs_0.6.5                 RANN_2.6.1                 
    ## [119] promises_1.2.1              xtable_1.8-4               
    ## [121] cluster_2.1.6               evaluate_0.23              
    ## [123] cli_3.6.2                   compiler_4.3.0             
    ## [125] Rsamtools_2.18.0            rlang_1.1.3                
    ## [127] crayon_1.5.2                future.apply_1.11.1        
    ## [129] labeling_0.4.3              plyr_1.8.9                 
    ## [131] stringi_1.8.3               viridisLite_0.4.2          
    ## [133] deldir_2.0-2                BiocParallel_1.36.0        
    ## [135] munsell_0.5.0               Biostrings_2.70.1          
    ## [137] lazyeval_0.2.2              spatstat.geom_3.2-7        
    ## [139] Matrix_1.6-5                future_1.33.1              
    ## [141] statmod_1.5.0               shiny_1.8.0                
    ## [143] highr_0.10                  SummarizedExperiment_1.32.0
    ## [145] ROCR_1.0-11                 igraph_1.6.0               
    ## [147] fastmatch_1.1-4
