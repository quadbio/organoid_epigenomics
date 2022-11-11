source('~/scripts/single_cell/graphs.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/perturbator/de.R')

library(Pando)
library(destiny)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd('~/projects/cutntag/')

data('motifs')

select <- dplyr::select

#### Read data ###
marks <- read_rds('data/all_marks_list_v3.2annot.rds')
motifs <- read_rds('~/resources/DB/JASPAR/JASPAR2020_hs_all_collections_motifs.rds')

marks <- map(marks, function(srt){
    srt@active.assay <- 'peaks'
    print('Adding motifs')
    srt <- AddMotifs(srt, genome=BSgenome.Hsapiens.UCSC.hg38, pfm=motifs, assay='peaks')
    print('Running chromvar')
    srt <- RunChromVAR(srt, genome=BSgenome.Hsapiens.UCSC.hg38, assay='peaks')
    return(srt)
})

marks %>% write_rds('data/all_marks_list_v3.3motifs.rds')
