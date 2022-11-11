source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/atac.R')
source('~/scripts/perturbator/de.R')

library(Pando)

filter <- dplyr::filter
select <- dplyr::select

setwd('~/projects/cutntag/')


#### Read data ####
marks <- read_rds('data/all_marks_list_v3.3motifs.rds')
rna <- read_rds('data/RNA/RNA_all_srt_v2.2matched.rds')


#### UMAPS ####
marks_meta <- marks %>% map(function(srt){
    srt@meta.data %>% 
        as_tibble(rownames='cell') %>% 
        inner_join(as_tibble(srt[['umap']]@cell.embeddings, rownames='cell')) %>% 
        return()
})

rna_meta <- rna@meta.data %>% 
    as_tibble(rownames='cell') %>% 
    inner_join(as_tibble(rna[['cssumap']]@cell.embeddings, rownames='cell')) 

rna_meta$age <- factor(rna_meta$age, levels=names(colours_timescale))

marks_meta$H3K27ac$age <- ifelse(marks_meta$H3K27ac$age=='128d', '4mo', marks_meta$H3K27ac$age)
marks_meta$H3K27me3$age <- ifelse(marks_meta$H3K27me3$age=='128d', '4mo', marks_meta$H3K27me3$age)
marks_meta$H3K4me3$age <- ifelse(marks_meta$H3K4me3$age=='128d', '4mo', marks_meta$H3K4me3$age)

marks_meta$H3K27ac$age <- factor(str_replace(marks_meta$H3K27ac$age, '-', '_'), levels=names(colours_timescale))
marks_meta$H3K27me3$age <- factor(str_replace(marks_meta$H3K27me3$age, '-', '_'), levels=names(colours_timescale))
marks_meta$H3K4me3$age <- factor(str_replace(marks_meta$H3K4me3$age, '-', '_'), levels=names(colours_timescale))


#### Celltype ####
ggplot(arrange(rna_meta, desc(celltype_jf)), aes(-CSSUMAP_1, CSSUMAP_2, fill=celltype_jf)) +
    geom_point(shape=21, size=2, stroke=0.05) +
    scale_fill_manual(values=pantone_celltype) +
    theme_void() + no_legend()

ggsave('plots/paper/fig1/fig1_RNA_celltype_umap.png', width=15, height=10)    


ggplot(arrange(marks_meta$H3K27me3, desc(celltype_jf)), aes(-UMAP_1, UMAP_2, fill=celltype_jf)) +
    geom_point(shape=21, size=2, stroke=0.05) +
    scale_fill_manual(values=pantone_celltype) +
    theme_void() + no_legend()
     
ggsave('plots/paper/fig1/fig1_H3K27me3_celltype_umap.png', width=15, height=10)   


ggplot(arrange(marks_meta$H3K27ac, desc(celltype_jf)), aes(-UMAP_1, UMAP_2, fill=celltype_jf)) +
    geom_point(shape=21, size=2, stroke=0.05) +
    scale_fill_manual(values=pantone_celltype) +
    theme_void() + no_legend()
     
ggsave('plots/paper/fig1/fig1_H3K27ac_celltype_umap.png', width=15, height=10)   


ggplot(arrange(marks_meta$H3K4me3, celltype_jf!='other'), aes(UMAP_1, UMAP_2, fill=celltype_jf)) +
    geom_point(shape=21, size=2, stroke=0.05) +
    scale_fill_manual(values=pantone_celltype) +
    theme_void() + no_legend()
     
ggsave('plots/paper/fig1/fig1_H3K4me3_celltype_umap.png', width=15, height=10)   


#### Age ####
ggplot(arrange(rna_meta, desc(age)), aes(-CSSUMAP_1, CSSUMAP_2, fill=age)) +
    geom_point(shape=21, size=2, stroke=0.05) +
    scale_fill_manual(values=colours_timescale) +
    theme_void() + no_legend()

ggsave('plots/paper/fig1/fig1_RNA_age_umap.png', width=15, height=10)    


ggplot(arrange(marks_meta$H3K27me3, desc(age)), aes(-UMAP_1, UMAP_2, fill=age)) +
    geom_point(shape=21, size=2, stroke=0.05) +
    scale_fill_manual(values=colours_timescale) +
    theme_void() + no_legend()
     
ggsave('plots/paper/fig1/fig1_H3K27me3_age_umap.png', width=15, height=10)   


ggplot(arrange(marks_meta$H3K27ac, desc(age)), aes(-UMAP_1, UMAP_2, fill=age)) +
    geom_point(shape=21, size=2, stroke=0.05) +
    scale_fill_manual(values=colours_timescale) +
    theme_void() + no_legend()
     
ggsave('plots/paper/fig1/fig1_H3K27ac_age_umap.png', width=15, height=10)   


ggplot(arrange(marks_meta$H3K4me3, desc(age)), aes(UMAP_1, UMAP_2, fill=age)) +
    geom_point(shape=21, size=2, stroke=0.05) +
    scale_fill_manual(values=colours_timescale) +
    theme_void() + no_legend()
     
ggsave('plots/paper/fig1/fig1_H3K4me3_age_umap.png', width=15, height=10)   


#### Barplots ####

rna_meta$age_plot <- fct_recode(rna_meta$age, '16'='4mo', '32'='8mo', '6'='ret_6w', '12'='ret_12w', '2'='15d', '8'='60d', '4'='35d')


p1 <- ggplot(rna_meta, aes(age_plot, fill=age)) +
    geom_bar() +
    scale_fill_manual(values=colours_timescale) +
    scale_y_continuous(expand=c(0,0), breaks=c(0,4000,8000)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    no_legend() +
    no_label() +
    article_text() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.y = element_line(size=0.2, color='lightgrey'),
        plot.margin = unit(c(0,0,0,0), 'lines')
    ) 


p2 <- ggplot(rna_meta, aes(age_plot, fill=celltype_jf)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=pantone_celltype) +
    scale_y_continuous(expand=c(0,0), breaks=c(0,1)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    no_legend() +
    no_label() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    ) 

p1 / p2 + plot_layout(heights=c(1,3))

ggsave('plots/paper/fig1/fig1_RNA_celltype_age_barplot.pdf', width=3, height=3, units='cm')  



marks_meta$H3K4me3$age_plot <- fct_recode(marks_meta$H3K4me3$age, '16'='4mo', '32'='8mo', '6'='ret_6w', '12'='ret_12w', '2'='15d', '8'='60d', '4'='35d')
p1 <- ggplot(marks_meta$H3K4me3, aes(age_plot, fill=age)) +
    geom_bar() +
    scale_fill_manual(values=colours_timescale) +
    scale_y_continuous(expand=c(0,0), breaks=c(0,4000,8000)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    no_legend() +
    no_label() +
    article_text() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.y = element_line(size=0.2, color='lightgrey'),
        plot.margin = unit(c(0,0,0,0), 'lines')
    ) 


p2 <- ggplot(marks_meta$H3K4me3, aes(age_plot, fill=celltype_jf)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=pantone_celltype) +
    scale_y_continuous(expand=c(0,0), breaks=c(0,1)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    no_legend() +
    no_label() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    ) 

p1 / p2 + plot_layout(heights=c(1,3))

ggsave('plots/paper/fig1/fig1_H3K4me3_celltype_age_barplot.pdf', width=3, height=3, units='cm')  



marks_meta$H3K27ac$age_plot <- fct_recode(marks_meta$H3K27ac$age, '16'='4mo', '32'='8mo', '6'='ret_6w', '12'='ret_12w', '2'='15d', '8'='60d', '4'='35d')
p1 <- ggplot(marks_meta$H3K27ac, aes(age_plot, fill=age)) +
    geom_bar() +
    scale_fill_manual(values=colours_timescale) +
    scale_y_continuous(expand=c(0,0), breaks=c(0,4000,8000)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    no_legend() +
    no_label() +
    article_text() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.y = element_line(size=0.2, color='lightgrey'),
        plot.margin = unit(c(0,0,0,0), 'lines')
    ) 


p2 <- ggplot(marks_meta$H3K27ac, aes(age_plot, fill=celltype_jf)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=pantone_celltype) +
    scale_y_continuous(expand=c(0,0), breaks=c(0,1)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    no_legend() +
    no_label() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    ) 

p1 / p2 + plot_layout(heights=c(1,3))

ggsave('plots/paper/fig1/fig1_H3K27ac_celltype_age_barplot.pdf', width=3, height=3, units='cm')  



marks_meta$H3K27me3$age_plot <- fct_recode(marks_meta$H3K27me3$age, '16'='4mo', '32'='8mo', '6'='ret_6w', '12'='ret_12w', '2'='15d', '8'='60d', '4'='35d')
p1 <- ggplot(marks_meta$H3K27me3, aes(age_plot, fill=age)) +
    geom_bar() +
    scale_fill_manual(values=colours_timescale) +
    scale_y_continuous(expand=c(0,0), breaks=c(0,4000,8000)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    no_legend() +
    no_label() +
    article_text() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.y = element_line(size=0.2, color='lightgrey'),
        plot.margin = unit(c(0,0,0,0), 'lines')
    ) 


p2 <- ggplot(marks_meta$H3K27me3, aes(age_plot, fill=celltype_jf)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=pantone_celltype) +
    scale_y_continuous(expand=c(0,0), breaks=c(0,1)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    no_legend() +
    no_label() +
    theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    ) 

p1 / p2 + plot_layout(heights=c(1,3))

ggsave('plots/paper/fig1/fig1_H3K27me3_celltype_age_barplot.pdf', width=3, height=3, units='cm')  



#### Tracks ####
ct_order <- c('psc', 'non_nect', 'nect', 'ctx_npc', 'ctx_ip', 'ctx_ex', 'dien_npc', 'dien_ex', 
              'nt_npc', 'mesen_ex', 'rhom_ex', 'RPC', 'RGC', 'astrocytes', 'OPC', 'choroid_plexus')

H3K27me3_cov <- subset(marks$H3K27me3, celltype_jf!='other' & celltype_jf!='OPC')
H3K27me3_cov$celltype_jf <- case_when(
    H3K27me3_cov$celltype_jf %in% c('psc', 'non_nect', 'nect', 'astrocytes', 'choroid_plexus') ~ H3K27me3_cov$celltype_jf,
    H3K27me3_cov$celltype_jf %in% c('ctx_npc', 'ctx_ip', 'ctx_ex') ~ 'ctx_npc',
    H3K27me3_cov$celltype_jf %in% c('dien_npc', 'dien_ex') ~ 'dien_npc',
    H3K27me3_cov$celltype_jf %in% c('nt_npc', 'mesen_ex', 'rhom_ex') ~ 'nt_npc',
    H3K27me3_cov$celltype_jf %in% c('RPC', 'RGC') ~ 'RPC'
    ) %>% factor(levels=ct_order)

H3K27ac_cov <- subset(marks$H3K27ac, celltype_jf!='other')
H3K27ac_cov$celltype_jf <- factor(H3K27ac_cov$celltype_jf, levels=ct_order)

H3K4me3_cov <- subset(marks$H3K4me3, celltype_jf!='other')
H3K4me3_cov$celltype_jf <- factor(H3K4me3_cov$celltype_jf, levels=ct_order)


#### FOXG1 ####
region <- FindRegion(H3K27me3_cov, 'chr14-28750000-28800000')

p1 <- CoveragePlot(H3K27me3_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K27me3_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_FOXG1_H3K27me3_tracks.pdf', width=7, height=25, units='cm')


p1 <- CoveragePlot(H3K27ac_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K27ac_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_FOXG1_H3K27ac_tracks.pdf', width=7, height=25, units='cm')


p1 <- CoveragePlot(H3K4me3_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K4me3_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_FOXG1_H3K4me3_tracks.pdf', width=7, height=25, units='cm')


#### POU5F1 ####
region <- FindRegion(H3K27me3_cov, 'chr6-31155000-31180000')

p1 <- CoveragePlot(H3K27me3_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K27me3_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_POU5F1_H3K27me3_tracks.pdf', width=7, height=25, units='cm')


p1 <- CoveragePlot(H3K27ac_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K27ac_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_POU5F1_H3K27ac_tracks.pdf', width=7, height=25, units='cm')


p1 <- CoveragePlot(H3K4me3_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K4me3_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_POU5F1_H3K4me3_tracks.pdf', width=7, height=25, units='cm')



#### SIX6 ####
region <- FindRegion(H3K27me3_cov, 'chr14-60490000-60530000')

p1 <- CoveragePlot(H3K27me3_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K27me3_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_SIX6_H3K27me3_tracks.pdf', width=7, height=25, units='cm')


p1 <- CoveragePlot(H3K27ac_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K27ac_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_SIX6_H3K27ac_tracks.pdf', width=7, height=25, units='cm')


p1 <- CoveragePlot(H3K4me3_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K4me3_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_SIX6_H3K4me3_tracks.pdf', width=7, height=25, units='cm')




#### HOX ####
region <- FindRegion(H3K27me3_cov, 'chr17-48500000-48700000')

p1 <- CoveragePlot(H3K27me3_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K27me3_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_HOX_H3K27me3_tracks.pdf', width=7, height=25, units='cm')


p1 <- CoveragePlot(H3K27ac_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K27ac_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_HOX_H3K27ac_tracks.pdf', width=7, height=25, units='cm')


p1 <- CoveragePlot(H3K4me3_cov, region, annotation=F, peaks=F, group.by='celltype_jf', window=500) +
    scale_fill_manual(values=pantone_celltype)
p2 <- AnnotationPlot(H3K4me3_cov, region)
(p1 / p2 & article_text() & theme_void() & no_legend()) + plot_layout(heights=c(10,1))
ggsave('plots/paper/fig1/fig1_HOX_H3K4me3_tracks.pdf', width=7, height=25, units='cm')












