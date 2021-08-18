# for windows version R 3.5.1
# run_SPRING.filter_cells.R
# filter cells for running SPRING
# AUTHOR: Tuo Zhang
# DATE: 12/19/2018
# VERSION: 1.0
# NEW: filter cells based on percent.mito and nFeature_RNA (for SPRING)
# NOTE: first version
# 

library(Seurat)
library(dplyr)
library(Matrix)
#library(Rmagic)
library(SingleR)
library(Rmagic)
#library("gplots")
#library("RColorBrewer")

basedir <- "/Users/novawang/testdir/scRNA-seq" 
#basedir <- "/Users/freshtuo/Work/SingleCell/DropSeq" "/Users/novawang/testdir/Bioinformatics"
workdir <- paste(basedir, "DropSeq_MJ_171005", sep="\\")
cellcyclegenes <- "D:\\Tools\\Seurat\\CellCycles.v2\\regev_lab_cell_cycle_genes.uppercase.txt"
countfile <- countfile <- file.path(basedir, "counts.merged.clean.with_MT.no_zero_counts_genes.txt.gz")
infodir <- paste(workdir, "info", sep="\\")
figdir <- file.path(basedir,"Plots")
countsfile.forSPRING <- file.path(basedir,"Spring","forSPRING.expression_data.tsv.gz")
genefile.forSPRING <- file.path(basedir,"Spring","forSPRING.gene_list.txt.gz")
cellfile.forSPRING <- file.path(basedir,"Spring","forSPRING.cell_group.csv.gz")
coordfile.fromSPRING <- paste(infodir, "SPRING.coordinates.txt", sep="\\")
colorfile.fromPaper <- paste(infodir, "cell.colors.txt", sep="\\")
orderfile.fromPaper <- paste(infodir, "cell.plotting.order.txt", sep="\\")
normscaler <- 1e4
rseed <- 98
project <- "Ming"

# read in raw counts data
matt.data <- read.table("counts.merged.clean.with_MT.no_zero_counts_genes.txt.gz", sep="\t", header=TRUE, row.names=1)
dim(matt.data)
# 23844  5000

##### convert to sparsed matrix (cause error in calling colSums)
####matt.data.matrix <- Matrix(as.matrix(matt.data), sparse=T)
####class(matt.data.matrix)
####dim(matt.data.matrix)
##### 19198  4725

# set a random seed
set.seed(rseed)

# Initialize the Seurat object with the raw (non-normalized data)
# Do not do any normalizations at this step
# genes expressed in >= 1 cell; keep all cells with at least 100 detected genes
matt <- CreateSeuratObject(counts = matt.data, raw.data=matt.data, min.cells =1, min.features = 100, is.expr=0, names.field=1, names.delim="_", do.scale=F, do.center=F, project=project)
matt
# An object of class seurat in project Ming 
#  23844 genes across 4976 samples.

# calculate the percentage of mitochondrial genes and ribosome genes
# and store the results in percent.mito and percent.ribo using AddMetaData.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric
# !!! needs to modify: 
# ddSeq data, no mito genes
# use '^mt-' for mouse mm10; and "^MT-" for human hg19
mito.genes <- grep(pattern="^MT-", x=rownames(matt@assays[["RNA"]]@data), value=1)
percent.mito <- colSums(matt@assays[["RNA"]]@data[mito.genes,])/colSums(matt@assays[["RNA"]]@data)
length(mito.genes)
# 35
ribo.genes <- grep(pattern="^RPL|^RPS", x=rownames(matt@assays[["RNA"]]@data), value=1)
percent.ribo <- colSums(matt@assays[["RNA"]]@data[ribo.genes,])/colSums(matt@assays[["RNA"]]@data)
length(ribo.genes)
# 128

# take a quick look at the percent.mito and percent.ribo
summary(percent.mito)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001339 0.057369 0.082088 0.097237 0.113798 0.888681 
length(percent.mito)
# 4976
sum(percent.mito>0.1)
# 1703
sum(percent.mito>0.05)
# 4088
sum(percent.mito>0.08)
# 2586
sum(percent.mito>0.15)
# 537
sum(percent.mito>0.2)
# 210

summary(percent.ribo)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.004073 0.154976 0.210816 0.208987 0.266454 0.498458 
length(percent.ribo)
# 4976
sum(percent.ribo>0.1)
# 4482
sum(percent.ribo>0.05)
# 4729
sum(percent.ribo>0.08)
# 4622
sum(percent.ribo>0.15)
# 3827
sum(percent.ribo>0.2)
# 2732

# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
matt <- AddMetaData(object=matt, metadata=percent.mito, col.name="percent.mito")
matt <- AddMetaData(object=matt, metadata=percent.ribo, col.name="percent.ribo")

head(matt@meta.data, n=3)

# fix orig.ident
sample.rename <- c("D0_unstimulated","D0_unstimulated","D8_stimulated","D8_stimulated","D8_other","D8_other","D8_CD42bCD41","D8_CD42bCD41")
names(sample.rename) <- c("R1M1","R1M2","R1M3","R1M4","R1M5","R1M6","R1M7","R1M8")
matt@meta.data$orig.ident <- as.vector(sample.rename[substr(rownames(matt@meta.data),1,4)])
#matt@meta.data$orig.ident <- factor(matt@meta.data$orig.ident, levels=c("D0_unstimulated","D8_other","D8_stimulated","D8_CD42bCD41"))

# Basic QC and selecting cells for further analysis
# nFeature_RNA and nCount_RNA are automatically calculated for every object by Seurat.
# used to be called nGene and nUMI
# Draw violin plot for number of genes/UMIs per gene in each batch of samples
png(file=file.path('Plots', "VlnPlot.QC.png"), width=12800, height=4800, units="px", pointsize=6, res=600)
VlnPlot(object=matt, features=c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"))
dev.off()
# FeatureScatter is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
png(file=file.path('Plots', "VlnPlot.QC.png"), width=7200, height=2400, units="px", pointsize=6, res=600)
par(mfrow=c(1,3))
FeatureScatter(object=matt, feature1="nCount_RNA", feature2="nFeature_RNA")
FeatureScatter(matt, "nCount_RNA", "percent.mito")
FeatureScatter(matt, "nCount_RNA", "percent.ribo")
par(mfrow=c(1,1))
dev.off()
####png(file=file.path("Plots", "FeatureScatter.eg.png"), width = 6400, height = 6400, units = "px", pointsize = 6, res=600)
####par(mfrow=c(2,2))
####FeatureScatter(matt, "Shox2", "Gapdh")
####FeatureScatter(matt, "Shox2", "Gapdh", cell.ids = WhichCells(matt,"R1M1"))
####FeatureScatter(matt, "Shox2", "Gapdh", cell.ids = WhichCells(matt,"R1M2"))
####FeatureScatter(matt, "Shox2", "Gapdh", cell.ids = WhichCells(matt,"R1M3"))
####dev.off()
####par(mfrow=c(1,1))

# refined per-gene violin plot
# nFeature_RNA 
g <- ggplot(matt@meta.data, aes(x=orig.ident, y=nFeature_RNA, fill=orig.ident))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
g <- g + geom_hline(yintercept=4000, linetype="dashed", color = "red")
ggsave(file.path('Plots', "Violin.nFeature_RNA.png"), width=8, height=6, dpi=600)
# nCount_RNA
g <- ggplot(matt@meta.data, aes(x=orig.ident, y=nCount_RNA, fill=orig.ident))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
ggsave(file.path('Plots',"Violin.nCount_RNA.png"), width=8, height=6, dpi=600)
# percent.mito
g <- ggplot(matt@meta.data, aes(x=orig.ident, y=percent.mito, fill=orig.ident))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
g <- g + geom_hline(yintercept=0.2, linetype="dashed", color = "red")
ggsave(file.path('Plots', "Violin.percent.mito.png"), width=8, height=6, dpi=600)
# percent.ribo
g <- ggplot(matt@meta.data, aes(x=orig.ident, y=percent.ribo, fill=orig.ident))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
ggsave(file.path('Plots', "Violin.percent.ribo.png"), width=8, height=6, dpi=600)

# # run NormalizeData first
# # Normalizing the data
# # employ a global-scaling normalization method "LogNormalize" that
# # normalizes the gene expression measurements for each cell by the total expression,
# # multiplies this by a scale factor, and log-transform the result.
matt <- NormalizeData(object=matt, normalization.method="LogNormalize", scale.factor=normscaler)
# 
# # correlation between PRMT1 expression and percent.mito/#genes per cell
cell.info <- matt@meta.data
cell.info$PRMT1 <- as.vector(t(matt@assays[["RNA"]]@data["PRMT1",rownames(matt@meta.data)]))
# # PRMT1 expression and percent.mito
g <- ggplot(cell.info, aes(x=percent.mito, y=PRMT1, color=orig.ident))
g <- g + geom_point(shape=19, size=2, alpha=0.4)
g <- g + geom_vline(xintercept=0.2, linetype="dashed", color = "red")
g <- g + theme_classic()
ggsave(file.path('Plots', "correlation.percent_mito.PRMT1.png"), width=8, height=6, dpi=600)

# # PRMT1 expression and #genes per cell
g <- ggplot(cell.info, aes(x=nFeature_RNA, y=PRMT1, color=orig.ident))
g <- g + geom_point(shape=19, size=2, alpha=0.4)
g <- g + geom_vline(xintercept=4000, linetype="dashed", color = "red")
g <- g + theme_classic()
ggsave(file.path('Plots', "correlation.gene_per_cell.PRMT1.png"), width=8, height=6, dpi=600)
# 
# # ACTB expression and #genes per cell
cell.info$ACTB <- as.vector(t(matt@assays[["RNA"]]@data["ACTB", rownames(matt@meta.data)]))
g <- ggplot(cell.info, aes(x=nFeature_RNA, y=ACTB, color=orig.ident))
g <- g + geom_point(shape=19, size=2, alpha=0.4)
g <- g + geom_vline(xintercept=4000, linetype="dashed", color = "red")
g <- g + theme_classic()
ggsave(file.path('Plots', "correlation.gene_per_cell.ACTB.png"), width=8, height=6, dpi=600)
# 
# # GAPDH expression and #genes per cell
cell.info$GAPDH <- as.vector(t(matt@assays[["RNA"]]@data["GAPDH", rownames(matt@meta.data)]))
g <- ggplot(cell.info, aes(x=nFeature_RNA, y=GAPDH, color=orig.ident))
g <- g + geom_point(shape=19, size=2, alpha=0.4)
g <- g + geom_vline(xintercept=4000, linetype="dashed", color = "red")
g <- g + theme_classic()
ggsave(file.path('Plots', "correlation.gene_per_cell.GAPDH.png"), width=8, height=6, dpi=600)

# percent.mito 0.2
# nFeature_RNA 4000

# not used in new version Seurat v2.0, replaced by subset
##### We can select a subset of cells (i.e. remove cells that are potentially outliers)
##### Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nFeature_RNA but on anything in the object (as in FeatureScatter above)
##### we filter out:
##### cells that have unique gene counts over 2500 or below 200
####matt <- subset(matt, subset = nFeature_RNA <= 2500)
####matt <- subset(matt, subset= nFeature_RNA >= 200)
##### cells that have unique UMI counts over 50000
####matt <- subset(matt, subset = nCount_RNA <= 50000)
##### cells that have mitochondrial genes percentage over 0.1
####matt <- subset(matt, subset= percent.mito <= 0.1)

# Filter out cells that have mitochondrial genes greater than 0.1 (10%),
# and cells that have unique gene counts over 4000
matt <- subset(matt, subset= nFeature_RNA <= 4000)
matt <- subset(matt, subset= percent.mito <= 0.2)

# redraw nFeature_RNA/nCount_RNA/percent.mito/percent.ribo distribution plot
png(file=file.path('Plots', "VlnPlot.QC.clean.png"), width=12800, height=4800, units="px", pointsize=6, res=600)
VlnPlot(matt, c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"))
dev.off()
# redraw nFeature_RNA - nCount_RNA scatter plot
png(file=file.path('Plots', "FeatureScatter.QC.clean.png"), width=7500, height=2400, units="px", pointsize=6, res=600)
par(mfrow=c(1,3))
FeatureScatter(matt, "nCount_RNA", "nFeature_RNA")
FeatureScatter(matt, "nCount_RNA", "percent.mito")
FeatureScatter(matt, "nCount_RNA", "percent.ribo")
par(mfrow=c(1,1))
dev.off()

# refined per-gene violin plot
# nFeature_RNA
g <- ggplot(matt@meta.data, aes(x=orig.ident, y=nFeature_RNA, fill=orig.ident))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
g <- g + geom_hline(yintercept=4000, linetype="dashed", color = "red")
ggsave(file.path('Plots', "Violin.nFeature_RNA.clean.png"), width=8, height=6, dpi=600)
# nCount_RNA
g <- ggplot(matt@meta.data, aes(x=orig.ident, y=nCount_RNA, fill=orig.ident))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
ggsave(file.path('Plots', "Violin.nCount_RNA.clean.png"), width=8, height=6, dpi=600)
# percent.mito
g <- ggplot(matt@meta.data, aes(x=orig.ident, y=percent.mito, fill=orig.ident))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
g <- g + geom_hline(yintercept=0.2, linetype="dashed", color = "red")
ggsave(file.path('Plots', "Violin.percent.mito.clean.png"), width=8, height=6, dpi=600)
# percent.ribo
g <- ggplot(matt@meta.data, aes(x=orig.ident, y=percent.ribo, fill=orig.ident))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
ggsave(file.path('Plots', "Violin.percent.ribo.clean.png"), width=8, height=6, dpi=600)

# after removing unwanted cells, how many cells/genes left?
matt
# An object of class seurat in project Ming 
#  23844 genes across 4678 samples.

# run NormalizeData first
# Normalizing the data
# employ a global-scaling normalization method "LogNormalize" that
# normalizes the gene expression measurements for each cell by the total expression,
# multiplies this by a scale factor, and log-transform the result.
matt <- NormalizeData(object=matt, normalization.method="LogNormalize", scale.factor=normscaler)

# save seurat object
save(matt, file=file.path('Original Data', 'matt.before_tSNE.Robj'))

# output log-scale expression data
write.table(as.data.frame(as.matrix(matt@assays[["RNA"]]@data)), file = file.path("Spring","expression.table.log-scale.txt.gz"), quote=FALSE, na="", sep="\t", col.names=NA)

# write info to file
write.table(cell.info[rownames(matt@meta.data),], file=file.path("Spring", "info.table.txt"), quote=FALSE, na="", sep="\t", col.names=NA)

# prepare inputs for SPRING
# prepare input files for SPRING
counts.forSPRING <- matt@assays[["RNA"]]@counts[,rownames(matt@meta.data)]
counts.forSPRING <- counts.forSPRING[rowSums(counts.forSPRING) > 0,]
dim(counts.forSPRING)
# 23489  4678
# gene list file
gene.file <- gzfile(file.path("Spring","forSPRING.gene_list.txt.gz"), 'w')
write.table(rownames(counts.forSPRING), file=gene.file, row.names=FALSE, col.names=FALSE, quote=FALSE, na="")
# cell info file
cell.file <- gzfile(file.path("Spring","forSPRING.cell_group.csv.gz"), 'w')
write.table(t(matt@meta.data[,c("orig.ident"),drop=F]), file=cell.file, row.names=T, col.names=F, quote=F, na="",sep=",")
# expression data
exp.file <- gzfile(file.path("Spring","forSPRING.expression_data.tsv.gz"), 'w')
write.table(as.data.frame(as.matrix(counts.forSPRING)), file=exp.file, row.names=T, col.names=NA, quote=F, na="", sep="\t")

# collect info for each cell
cell.info <- matt@meta.data
cell.info$PRMT1 <- as.vector(t(matt@assays[["RNA"]]@data["PRMT1",rownames(matt@meta.data)]))

# load coordinates result from SPRING
coord <- read.table(coordfile.fromSPRING, header=F, check.names=F, sep=",", row.names=1)
dim(coord)
# 4678    2
colnames(coord) <- c("X","Y")
rownames(coord) <- rownames(cell.info)
head(coord)

# incorporate coordinates to cell.info
cell.info.merged <- merge(cell.info, coord, by=0, all=T)
rownames(cell.info.merged) <- cell.info.merged$Row.names
cell.info.merged <- cell.info.merged[,-1]
dim(cell.info.merged)
# 4678    8

# check PRMT1 expressionfocus on cells in five blocks
# -Y: -400 ~ -600
# -X: -1250~-1000; -1000~-750; -750~-500; -500~-250
my.cate <- function(x, y){
  if (-y > -600 & -y <= -400){
    if (-x > -1250 & -x <= -1000){
      return("G0")
    } else if (-x > -1000 & -x <= -750){
      return("G1")
    } else if (-x > -750 & -x <= -500){
      return("G2")
    } else if (-x > -500 & -x <= 250){
      return("G3")
    } else {
      return("Other")
    }
  } else {
    return("Other")
  }
}

# SPRING plot
###################################################### all cells ####################################################
# mirror x and y-axis
g <- ggplot(cell.info.merged, aes(x=-X, y=-Y, color=orig.ident))
g <- g + geom_point(shape=19, size=4, alpha=0.6)
g <- g + scale_color_manual(values=c("D0_unstimulated"="gray80","D8_CD42bCD41"="mediumorchid4","D8_other"="skyblue4","D8_stimulated"="mistyrose4"))
g <- g + theme_classic()
g <- g + theme(axis.title=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
ggsave(paste(figdir,"SPRING.filtered.by_condition.png",sep="\\"), width=11, height=8, dpi=300)
ggsave(paste(figdir,"SPRING.filtered.by_condition.svg",sep="\\"), width=11, height=8, dpi=300)

# SPRING plot, separated by condition
# switch x and y-axis
g <- ggplot(cell.info.merged, aes(x=-X, y=-Y, color=orig.ident))
g <- g + geom_point(shape=19, size=4, alpha=0.6) + facet_wrap(~orig.ident, ncol=2)
g <- g + scale_color_manual(values=c("D0_unstimulated"="gray80","D8_CD42bCD41"="mediumorchid4","D8_other"="skyblue4","D8_stimulated"="mistyrose4"))
g <- g + theme_classic()
g <- g + theme(axis.title=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
ggsave(paste(figdir,"SPRING.filtered.by_condition.separated.png",sep="\\"), width=21, height=16, dpi=300)
ggsave(paste(figdir,"SPRING.filtered.by_condition.separated.svg",sep="\\"), width=21, height=16, dpi=300)

# SPRING plot, expression of PRMT1
g <- ggplot(cell.info.merged[with(cell.info.merged, order(PRMT1)),], aes(x=-X, y=-Y, color=PRMT1))
g <- g + geom_point(shape=19, size=2, alpha=0.6)
g <- g + scale_color_gradient(low="gray", high="red")
g <- g + theme_classic()
ggsave(paste(figdir,"SPRING.filtered.exp_PRMT1.png",sep="\\"), width=10, height=8, dpi=300)

# SPRING plot, expression of PRMT1, separated by condition
g <- ggplot(cell.info.merged[with(cell.info.merged, order(PRMT1)),], aes(x=-X, y=-Y, color=PRMT1))
g <- g + geom_point(shape=19, size=2, alpha=0.6) + facet_wrap(~orig.ident, ncol=2)
g <- g + scale_color_gradient(low="gray", high="red")
g <- g + theme_classic()
ggsave(paste(figdir,"SPRING.filtered.exp_PRMT1.separated_conditions.png",sep="\\"), width=20, height=16, dpi=300)

# save cell annotations to file
write.table(cell.info.merged, file=paste(infodir,"cell.filtered.annotations.SPRING.txt",sep="\\"), quote=FALSE, na="", sep="\t", col.names=NA)

# assign block
cell.info.merged$my.category <- mapply(my.cate, cell.info.merged$X, cell.info.merged$Y)
head(cell.info.merged)
table(cell.info.merged$my.category)
#   G0    G1    G2    G3 Other 
# 1212   276   589   370  2231 

# violin plot of PRMT1 expression (origin) per category
g <- ggplot(subset(cell.info.merged, !my.category %in% c("Other")), aes(x=my.category, y=PRMT1, color=my.category))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
g <- g + theme_classic()
ggsave(paste(figdir, "violin.exp_PRMT1.origin.all_conditions.per_group.png", sep="\\"), width=8, height=6, dpi=300)

aggregate(PRMT1 ~ my.category, data=cell.info.merged, FUN=mean)
#   my.category     PRMT1
# 1          G0 0.3156121
# 2          G1 0.2952740
# 3          G2 0.3186058
# 4          G3 0.3262956
# 5       Other 0.1946917
######################################################################################################################

################################################ Minkui's request ####################################################
# check PRMT1 expression on cells in seven regions
#1 (X <-1000; -625 <Y< -375);  
#2 (X <-1000; -375 <Y< -125); 
#3 (-1000< X <-750; -625 <Y< -375); 
#4 (-750< X <-500; -625 <Y< -375); 
#5 (X > -500; -625 <Y< -375); 
#6 (-500< X; -375 <Y); 
#7 (Y<-625). 
# Please extract the PRMT1 expression of the seven regions for the four conditions (D0, D8-CD41CD42, D8-others; D8_stimulated). 
# Meanwhile, can you also extract the cell numbers of of the seven regions for the the four conditions
my.cate.check <- function(x, y){
  if (-y > -625){
    if (-y < -375){
      if (-x < -1000){
        return("R1")
      } else if (-x < -750){
        return("R3")
      } else if (-x < -500){
        return("R4")
      } else{
        return("R5")
      }
    } else if (-x > -500){
      return("R6")
    } else if (-y < -125){
      if (-x < -1000){
        return("R2")
      } else{
        return("Other")
      }
    } else {
      return("Other")
    }
  } else {
    return("R7")
  }
}

# assign block
cell.info.merged$my.category.check <- mapply(my.cate.check, cell.info.merged$X, cell.info.merged$Y)
head(cell.info.merged)
table(cell.info.merged$my.category.check)
# Other    R1    R2    R3    R4    R5    R6    R7 
#   247  1387  1119   294   590   436   467   138 

# violin plot of PRMT1 expression (origin) per category
g <- ggplot(subset(cell.info.merged, !my.category.check %in% c("Other")), aes(x=my.category.check, y=PRMT1, color=my.category.check))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
g <- g + theme_classic()
ggsave(paste(figdir, "violin.exp_PRMT1.origin.all_conditions.minkui.png", sep="\\"), width=8, height=6, dpi=300)

# mean expression of PRMT1 per region
aggregate(PRMT1 ~ my.category.check, data=cell.info.merged, FUN=mean)
#   my.category.check     PRMT1
# 1             Other 0.1012146
# 2                R1 0.3503965
# 3                R2 0.1367292
# 4                R3 0.2619048
# 5                R4 0.3525424
# 6                R5 0.3256881
# 7                R6 0.1477516
# 8                R7 0.1304348

# cell numbers per condition
as.data.frame.matrix(table(cell.info.merged[,c("my.category.check","orig.ident")]))
#       D0_unstimulated D8_CD42bCD41 D8_other D8_stimulated
# Other             216            3       12            16
# R1               1014            0      217           156
# R2                525            0      330           264
# R3                 41            7      143           103
# R4                 14          346       79           151
# R5                  2          332       22            80
# R6                  1          384        0            82
# R7                  0            5       61            72

# mean expression of PRMT1 per region in each condition
aggregate(PRMT1 ~ my.category.check, data=subset(cell.info.merged, orig.ident %in% c("D0_unstimulated")), FUN=mean)
# my.category.check     PRMT1
# 1             Other 0.1018519
# 2                R1 0.3500986
# 3                R2 0.1733333
# 4                R3 0.2195122
# 5                R4 0.3571429
# 6                R5 0.0000000
# 7                R6 0.0000000
aggregate(PRMT1 ~ my.category.check, data=subset(cell.info.merged, orig.ident %in% c("D8_CD42bCD41")), FUN=mean)
# my.category.check     PRMT1
# 1             Other 0.0000000
# 2                R3 0.0000000
# 3                R4 0.2485549
# 4                R5 0.2409639
# 5                R6 0.1276042
# 6                R7 0.0000000
aggregate(PRMT1 ~ my.category.check, data=subset(cell.info.merged, orig.ident %in% c("D8_other")), FUN=mean)
# my.category.check      PRMT1
# 1             Other 0.00000000
# 2                R1 0.12903226
# 3                R2 0.04242424
# 4                R3 0.16083916
# 5                R4 0.17721519
# 6                R5 0.00000000
# 7                R7 0.08196721
aggregate(PRMT1 ~ my.category.check, data=subset(cell.info.merged, orig.ident %in% c("D8_stimulated")), FUN=mean)
# my.category.check     PRMT1
# 1             Other 0.1875000
# 2                R1 0.6602564
# 3                R2 0.1818182
# 4                R3 0.4368932
# 5                R4 0.6821192
# 6                R5 0.7750000
# 7                R6 0.2439024
# 8                R7 0.1805556

######################################################################################################################

##################################### below work on D8 stimulated condition only #####################################
# subset cell info on D8 stimulated cells
cell.info.merged.sti <- subset(cell.info.merged, orig.ident %in% c("D8_stimulated"))

# SPRING plot, expression of PRMT1 (original), D8_stimulated cells 
g <- ggplot(cell.info.merged.sti[with(cell.info.merged.sti, order(PRMT1)),], aes(x=-X, y=-Y, color=PRMT1))
g <- g + geom_point(shape=19, size=2, alpha=0.6)
g <- g + scale_color_gradient(low="gray", high="red")
g <- g + theme_classic()
ggsave(paste(figdir,"SPRING.filtered.exp_PRMT1.orig.D8_stimulated.png",sep="\\"), width=10, height=8, dpi=300)

# assign block
cell.info.merged.sti$my.category <- mapply(my.cate, cell.info.merged.sti$X, cell.info.merged.sti$Y)
head(cell.info.merged.sti)
table(cell.info.merged.sti$my.category)
# G0    G1    G2    G3 Other 
# 127    98   150    69   480 

# violin plot of PRMT1 expression (origin) per category
g <- ggplot(subset(cell.info.merged.sti, !my.category %in% c("Other")), aes(x=my.category, y=PRMT1, color=my.category))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6)
g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
g <- g + theme_classic()
ggsave(paste(figdir, "violin.exp_PRMT1.origin.D8_stimulated.per_group.png", sep="\\"), width=8, height=6, dpi=300)

aggregate(PRMT1 ~ my.category, data=cell.info.merged.sti, FUN=mean)
#   my.category     PRMT1
# 1          G0 0.5607689
# 2          G1 0.3993660
# 3          G2 0.4962380
# 4          G3 0.6035621
# 5       Other 0.2538624

# # run magic on D8_stimulate sample
# # read in raw counts data
# counts <- read.table(countsfile.forSPRING, header=T, check.names=F, sep="\t", row.names=1)
# dim(counts)
# # 23489  4678
# 
# table(cell.info$orig.ident)
# # D0_unstimulated    D8_CD42bCD41        D8_other   D8_stimulated 
# #            1813            1077             864             924 
# 
# # select D8_stimulated and run MAGIC
# counts.sti <- counts[,rownames(subset(cell.info, orig.ident %in% c("D8_stimulated")))]
# # remove all-zero genes
# counts.sti <- counts.sti[rowSums(counts.sti) > 0,]
# dim(counts.sti)
# # 19041   924
# 
# # switch two dimensions
# counts.sti.t <- t(as.matrix(counts.sti))
# 
# # filter data
# # check how many cells expressing PRMT1 and DUSP4
# sum(counts.sti.t[,"PRMT1"]>0)
# # 261
# sum(counts.sti.t[,"DUSP4"]>0)
# # 3
# # # keep genes expressed in at least 10 cells
# # keep_cols <- colSums(counts.sti.t > 0) > 10
# # sum(keep_cols)
# # # 11778
# # counts.sti.t.clean <- counts.sti.t[, keep_cols]
# # dim(counts.sti.t.clean)
# # # 939 11778
# # NO cleaning
# counts.sti.t.clean <- counts.sti.t
# 
# # look at the distribution of library sizes
# g <- ggplot() + geom_histogram(aes(x=rowSums(counts.sti.t.clean)), bins=50)
# g <- g + geom_vline(xintercept=1000, color="red")
# g <- g + theme_classic()
# ggsave(paste(figdir, "histogram.nCount_RNA.per.cell.png", sep="\\"), width=8, height=6, dpi=300)
# 
# # use all cells
# # # keep cells with at least 1000 UMIs and at most 40000 UMIs
# # keep_rows <- rowSums(counts.sti.t.clean) > 1000 & rowSums(counts.sti.t.clean) < 40000
# # sum(keep_rows)
# # # 801
# # counts.sti.t.clean <- counts.sti.t.clean[keep_rows,]
# # dim(counts.sti.t.clean)
# # # 801 11778
# 
# # normalize data
# counts.norm <- library.size.normalize(counts.sti.t.clean)
# counts.norm <- sqrt(counts.norm)
# 
# # run MAGIC
# # on PRMT1 gene
# counts.magic.PRMT1 <- magic(counts.norm, genes=c("PRMT1","VIM","GAPDH"))
# # before
# g <- ggplot(as.data.frame(counts.norm)) + geom_point(aes(x=VIM, y=GAPDH, color=PRMT1))
# g <- g + scale_color_viridis_c(option="B")
# g <- g + theme_classic()
# ggsave(paste(figdir,"exp.PRMT1.before_magic.png",sep="\\"), width=8, height=6, dpi=300)
# # after
# g <- ggplot(counts.magic.PRMT1) + geom_point(aes(x=VIM, y=GAPDH, color=PRMT1))
# g <- g + scale_color_viridis_c(option="B")
# g <- g + theme_classic()
# ggsave(paste(figdir,"exp.PRMT1.after_magic.png",sep="\\"), width=8, height=6, dpi=300)
# # on DUSP4 gene
# counts.magic.DUSP4 <- magic(counts.norm, genes=c("DUSP4","VIM","GAPDH"))
# # before
# g <- ggplot(as.data.frame(counts.norm)) + geom_point(aes(x=VIM, y=GAPDH, color=DUSP4))
# g <- g + scale_color_viridis_c(option="B")
# g <- g + theme_classic()
# ggsave(paste(figdir,"exp.DUSP4.before_magic.png",sep="\\"), width=8, height=6, dpi=300)
# # after
# g <- ggplot(counts.magic.DUSP4) + geom_point(aes(x=VIM, y=GAPDH, color=DUSP4))
# g <- g + scale_color_viridis_c(option="B")
# g <- g + theme_classic()
# ggsave(paste(figdir,"exp.DUSP4.after_magic.png",sep="\\"), width=8, height=6, dpi=300)
# 
# # all_genes
# counts.magic.all <- magic(counts.norm, genes="all_genes", init=counts.magic.PRMT1)
# # output
# write.table(t(as.matrix(counts.magic.all$result)), file=gzfile(paste(infodir,"exp.sqrt.after_magic.tsv.gz",sep="\\")), quote=FALSE, na="", sep="\t", col.names=NA)
# #write.table(t(as.matrix((counts.magic.all$result)^2)), file=gzfile(paste(infodir,"exp.pseudo_counts.after_magic.tsv.gz",sep="\\")), quote=FALSE, na="", sep="\t", col.names=NA)
# 
# # load MAGIC-corrected expression 
# exp.magic <- t(as.matrix(counts.magic.all$result))
# # add DUSP4
# cell.info.merged$DUSP4 <- as.vector(t(matt@assays[["RNA"]]@data["DUSP4", rownames(matt@meta.data)]))
# # add ACTB
# cell.info.merged$ACTB <- as.vector(t(matt@assays[["RNA"]]@data["ACTB", rownames(matt@meta.data)]))
# # add GAPDH
# cell.info.merged$GAPDH <- as.vector(t(matt@assays[["RNA"]]@data["GAPDH", rownames(matt@meta.data)]))
# # subset cell info on D8 stimulated cells
# cell.info.merged.sti <- subset(cell.info.merged, orig.ident %in% c("D8_stimulated"))
# # add PRMT1 magic-corrected expression
# cell.info.merged.sti$PRMT1.magic <- as.vector(t(exp.magic["PRMT1", rownames(cell.info.merged.sti)]))
# # add DUSP4 magic-corrected expression
# cell.info.merged.sti$DUSP4.magic <- as.vector(t(exp.magic["DUSP4", rownames(cell.info.merged.sti)]))
# # add ACTB magic-corrected expression
# cell.info.merged.sti$ACTB.magic <- as.vector(t(exp.magic["ACTB", rownames(cell.info.merged.sti)]))
# # add GAPDH magic-corrected expression
# cell.info.merged.sti$GAPDH.magic <- as.vector(t(exp.magic["GAPDH", rownames(cell.info.merged.sti)]))
# 
# head(cell.info.merged.sti)
# 
# # assign block
# cell.info.merged.sti$my.category <- mapply(my.cate, cell.info.merged.sti$X, cell.info.merged.sti$Y)
# head(cell.info.merged.sti)
# table(cell.info.merged.sti$my.category)
# # G0    G1    G2    G3 Other 
# # 127    98   150    69   480 
# 
# # violin plot of magic-corrected PRMT1 expression per category
# g <- ggplot(subset(cell.info.merged.sti, !my.category %in% c("Other")), aes(x=my.category, y=PRMT1.magic, color=my.category))
# g <- g + geom_violin()
# g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6)
# g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
# g <- g + theme_classic()
# ggsave(paste(figdir, "violin.exp_PRMT1.magic.D8_stimulated.per_group.png", sep="\\"), width=8, height=6, dpi=300)
# # mean expression per group
# aggregate(PRMT1.magic ~ my.category, data=cell.info.merged.sti, FUN=mean)
# #   my.category PRMT1.magic
# # 1          G0   0.2888453
# # 2          G1   0.2968668
# # 3          G2   0.3334369
# # 4          G3   0.3276475
# # 5       Other   0.2508511
# 
# # violin plot of magic-corrected DUSP4 expression per category
# g <- ggplot(subset(cell.info.merged.sti, !my.category %in% c("Other")), aes(x=my.category, y=DUSP4.magic, color=my.category))
# g <- g + geom_violin()
# g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6)
# g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
# g <- g + theme_classic()
# ggsave(paste(figdir, "violin.exp_DUSP4.magic.D8_stimulated.per_group.png", sep="\\"), width=8, height=6, dpi=300)
# # mean expression per group
# aggregate(DUSP4.magic ~ my.category, data=cell.info.merged.sti, FUN=mean)
# #   my.category  DUSP4.magic
# # 1          G0 0.0008297119
# # 2          G1 0.0031883925
# # 3          G2 0.0012762344
# # 4          G3 0.0014516636
# # 5       Other 0.0053571517
# 
# # violin plot of magic-corrected ACTB expression per category
# g <- ggplot(subset(cell.info.merged.sti, !my.category %in% c("Other")), aes(x=my.category, y=ACTB.magic, color=my.category))
# g <- g + geom_violin()
# g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6)
# g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
# g <- g + theme_classic()
# ggsave(paste(figdir, "violin.exp_ACTB.magic.D8_stimulated.per_group.png", sep="\\"), width=8, height=6, dpi=300)
# # mean expression per group
# aggregate(ACTB.magic ~ my.category, data=cell.info.merged.sti, FUN=mean)
# #   my.category ACTB.magic
# # 1          G0   2.029044
# # 2          G1   2.300338
# # 3          G2   2.723813
# # 4          G3   2.989455
# # 5       Other   2.325374
# 
# # violin plot of magic-corrected GAPDH expression per category
# g <- ggplot(subset(cell.info.merged.sti, !my.category %in% c("Other")), aes(x=my.category, y=GAPDH.magic, color=my.category))
# g <- g + geom_violin()
# g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.6)
# g <- g + stat_summary(fun.data=mean_sdl, geom="pointrange", color="gray75")
# g <- g + theme_classic()
# ggsave(paste(figdir, "violin.exp_GAPDH.magic.D8_stimulated.per_group.png", sep="\\"), width=8, height=6, dpi=300)
# # mean expression per group
# aggregate(GAPDH.magic ~ my.category, data=cell.info.merged.sti, FUN=mean)
# #   my.category GAPDH.magic
# # 1          G0    2.322502
# # 2          G1    2.338436
# # 3          G2    2.384483
# # 4          G3    2.388863
# # 5       Other    2.292924
######################################################################################################################

########################################### check expression on mark genes ###########################################
# add expression of several markers
cell.info.merged$CD34 <- as.vector(t(matt@assays[["RNA"]]@data["CD34",rownames(matt@meta.data)]))
cell.info.merged$CD38 <- as.vector(t(matt@assays[["RNA"]]@data["CD38",rownames(matt@meta.data)]))

# how many cells express a given marker
round(sum(cell.info.merged$CD34 > 0)/nrow(cell.info.merged) * 100, digits=2)
# 17.4%
round(sum(cell.info.merged$CD38 > 0)/nrow(cell.info.merged) * 100, digits=2)
# 3.95%

# how many cell express both CD34 and CD38
nrow(subset(cell.info.merged, CD34 > 0 & CD38 > 0))
# 57 cells (1.22%)

# locate these cells co-expressing CD34&CD38 in SPRING plot
cell.info.merged$target <- with(cell.info.merged, ifelse(CD34>0 & CD38>0, "yes", "no"))
cell.info.merged$target <- factor(cell.info.merged$target, levels=c("no","yes"))
cell.info.merged <- cell.info.merged[with(cell.info.merged, order(target)),]
g <- ggplot(cell.info.merged, aes(x=-X, y=-Y, color=target))
g <- g + geom_point(shape=19, size=2, alpha=0.6)
#g <- g + scale_color_manual(values=c("yes"="red4","no"="royalblue1"))
g <- g + theme_classic()
ggsave(paste(figdir, "locate.CD34.CD38.positive.cells.png", sep="\\"), width=8, height=6, dpi=300)

# violin plot showing expression of CD38
g <- ggplot(cell.info.merged, aes(x=factor(1), y=CD38))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=median_hilow, geom="pointrange", color="gray75")
ggsave(paste(figdir, "violin.expression.CD38.png" ,sep="\\"), width=4, height=6, dpi=300)

# violin plot showing expression of CD38, on cells that co-express CD34 and CD38
g <- ggplot(subset(cell.info.merged, target %in% c("yes")), aes(x=factor(1), y=CD38))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=median_hilow, geom="pointrange", color="gray75")
ggsave(paste(figdir, "violin.expression.CD38.target.png" ,sep="\\"), width=4, height=6, dpi=300)

# locate cells co-expressing CD34&CD38 + CD38 in c(0.75, 1.25) in SPRING plot
cell.info.merged$target2 <- with(cell.info.merged, ifelse(target %in% c("yes") & CD38 <= 1.25 & CD38 > 0.75, "yes", "no"))
cell.info.merged$target2 <- factor(cell.info.merged$target2, levels=c("no","yes"))
cell.info.merged <- cell.info.merged[with(cell.info.merged, order(target2)),]
table(cell.info.merged$target2)
#   no  yes 
# 4654   24
g <- ggplot(cell.info.merged, aes(x=-X, y=-Y, color=target2))
g <- g + geom_point(shape=19, size=2, alpha=0.6)
#g <- g + scale_color_manual(values=c("yes"="red4","no"="royalblue1"))
g <- g + theme_classic()
ggsave(paste(figdir, "locate.CD34.CD38.positive.restrict.cells.png", sep="\\"), width=8, height=6, dpi=300)

######################################################################################################################

########################################### assign cell type using SingleR ###########################################
# build input files for SingleR
counts.forSingleR <- matt@assays[["RNA"]]@counts[,rownames(matt@meta.data)]
counts.forSingleR <- counts.forSingleR[rowSums(counts.forSingleR) > 0,]
dim(counts.forSingleR)
# 23489  4678
# cell annotation file
singler.cellann.file <- paste(infodir,"annotations.for.SingleR.txt.gz",sep="\\")
singler.cellann <- data.frame(orig.ident=matt@meta.data[,"orig.ident"])
rownames(singler.cellann) <- rownames(matt@meta.data)
write.table(singler.cellann, file=gzfile(singler.cellann.file), quote=FALSE, na="", sep="\t", col.names=NA)

# run SingleR
singler <- CreateSinglerObject(counts=counts.forSingleR, annot=singler.cellann.file,
                               project.name=project, min.genes=0,
                               technology="10X", species="Human", citation="",
                               ref.list=list(), normalize.gene.length=F, variable.genes="de",
                               fine.tune=F, do.main.types=T, reduce.file.size=F,
                               numCores=SingleR.numCores)

#------------------------------------- data based on HPCA database --------------------------------------------------#
# draw heatmap by the main cell types (before refine-tuning)
png(filename=paste(figdir,"singler.heatmap.HPCA.main_cell_type.png",sep="\\"), width=10, height=8, units="in", res=300)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n = Inf,
                    clusters=singler$meta.data$orig.ident)
dev.off()

# prepare SPRING coordinates
spring.n <- data.frame(X=-cell.info.merged$X, Y=-cell.info.merged$Y)
rownames(spring.n) <- rownames(cell.info.merged)

# present the annotations in a t-SNE plot (before refine-tuning):
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single.main,
                       spring.n, do.label=FALSE,
                       do.letters=T, labels=singler$singler[[1]]$SingleR.single.main$labels,
                       label.size=8, dot.size=5, alpha=0.4)
ggsave(filename=paste(figdir,"singler.SPRING.HPCA.main_cell_type.png", sep="\\"), plot=out$p, width=10, height=8, units="in", dpi=300)

# manually plot SPRING - cell type
g <- ggplot(out$df, aes(x=x, y=y, color=ident))
g <- g + geom_point(shape=19, size=2, alpha=0.6)
g <- g + theme_classic()
ggsave(filename=paste(figdir,"singler.SPRING.HPCA.main_cell_type.manual.png", sep="\\"), width=10, height=8, units="in", dpi=300)
g <- g + facet_wrap(~ident, ncol=4)
ggsave(filename=paste(figdir,"singler.SPRING.HPCA.main_cell_type.manual.split.png", sep="\\"), width=20, height=16, units="in", dpi=300)

#--------------------------------------------------------------------------------------------------------------------#

######################################################################################################################

############################################## assign cell type manually #############################################
# infer cell types based on known marker genes
# Ba (basophilic or mast cell)
ba.markers <- intersect(rownames(matt@assays[["RNA"]]@data), c("LMO4","IFITM1","LY6E","TPSAB1"))
length(ba.markers)
# 4
# Meg (megakaryocytic)
meg.markers <- intersect(rownames(matt@assays[["RNA"]]@data), c("PF4","ITGA2B","VWF","PBX1","MEF2C"))
length(meg.markers)
# 5
# MPP (multipotential progenitors)
mpp.markers <- intersect(rownames(matt@assays[["RNA"]]@data), c("HLF","GCNT2","PROCR"))
length(mpp.markers)
# 3
# Ly (lymphocytic)
ly.markers <- intersect(rownames(matt@assays[["RNA"]]@data), c("CD79A","IGLL1","VPREB3","VPREB1","LEF1"))
length(ly.markers)
# 5
# M (monocytic)
m.markers <- intersect(rownames(matt@assays[["RNA"]]@data), c("CSF1R","CCR2"))
length(m.markers)
# 2
# GN (granulocytic neutrophil)
gn.markers <- intersect(rownames(matt@assays[["RNA"]]@data), c("MPO","ELANE","PRTN3"))
length(gn.markers)
# 3

#------------------------------------------- run AddModuleScore -----------------------------------------------------#
# use pre-set seed
enrich.name <- "CellType"
genes.list <- list("Ba"=ba.markers, "Meg"=meg.markers, "MPP"=mpp.markers, "Ly"=ly.markers, "M"=m.markers, "GN"=gn.markers)
my.ct <- AddModuleScore(object=matt, genes.list=genes.list, enrich.name=enrich.name, seed.use=rseed, ctrl.size=min(vapply(X = genes.list, FUN = length, FUN.VALUE = numeric(1))))
# which columns are cell type scores?
cc.columns <- grep(pattern=enrich.name, x=colnames(my.ct@meta.data))
# extract cell type scores
cc.scores <- my.ct@meta.data[,cc.columns]
# clean intermediate object
#rm(my.ct)
#gc(verbose=F)
# assign each cell a cell type
assignments <- apply(
  X = cc.scores,
  MARGIN = 1,
  FUN = function(scores, first="Ba", second="Meg", third="MPP", fourth="Ly", fifth="M", sixth="GN", null="Unknown"){
    if (all(scores <= 0)){
      return(null)
    } else {
      if (sum(scores==max(scores)) > 1){# identical score in multiple cell types
        return(null)
      } else {
        return(c(first,second,third,fourth,fifth,sixth)[which(scores==max(scores))])
      }
    }
  }
)
table(assignments)
# assignments
#   Ba      GN      Ly       M     Meg     MPP Unknown 
# 1417     907     413      61    1132     127     621 

# check for entries with identical scores
# for (k in c(1:4678)){
#   if (length(assignments[[k]]) != 1){
#     print(k)
#   }
# }
# [1] 59
# [1] 1235
# [1] 1431

# merge cell type score with assigned cell type
cc.scores <- merge(cc.scores, data.frame(assignments),by=0)
colnames(cc.scores) <- c("rownames","Ba.Score","Meg.Score","MPP.Score","Ly.Score","M.Score","GN.Score","CellType")
rownames(cc.scores) <- cc.scores$rownames
cc.scores <- cc.scores[, -1]
matt <- AddMetaData(object=matt, metadata=cc.scores)

table(matt@meta.data$CellType)
#   Ba      GN      Ly       M     Meg     MPP Unknown 
# 1417     907     413      61    1132     127     621 

save(matt, file=paste(infodir, "matt.celltype.Robj", sep="\\"))

# fetch SPRING coordinates
spring.n <- data.frame(X=-cell.info.merged$X, Y=-cell.info.merged$Y)
rownames(spring.n) <- rownames(cell.info.merged)

# add cell type assignment
spring.n.celltype <- merge(spring.n, cc.scores, by=0, all=T)
rownames(spring.n.celltype) <- spring.n.celltype$Row.names
spring.n.celltype <- spring.n.celltype[,-1]
dim(spring.n.celltype)
# 4678   9
write.table(spring.n.celltype, file=paste(infodir, "info.table.celltype.txt", sep="\\"), quote=FALSE, na="", sep="\t", col.names=NA)

# check cell type on SPRING plot
g <- ggplot(spring.n.celltype, aes(x=X,y=Y,color=CellType))
g <- g + geom_point(shape=19, size=2, alpha=.6)
g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
g <- g + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=18)) + theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))
g <- g + ggtitle("SPRING plot by cell type")
ggsave(paste(figdir, "SPRING.Plot.by_cell_type.png", sep="\\"), width=11, height=8, dpi=300)
g <- g + facet_wrap(~CellType, ncol=3)
ggsave(paste(figdir, "SPRING.Plot.by_cell_type.separated.png", sep="\\"), width=22, height=16, dpi=300)

# violin plot of scores for each assigned cell type
scores.list <- list()
k <- 1
for (ct in c("Ba","Meg","MPP","M","Ly","GN")){
  tsubset <- subset(spring.n.celltype, CellType %in% c(ct))
  tscores <- tsubset[,c(paste(ct,"Score",sep="."),"CellType")]
  colnames(tscores) <- c("Score","CellType")
  scores.list[[k]] <- tscores
  k <- k + 1
}
# merge scores
my.scores <- rbind(scores.list[[1]], scores.list[[2]], scores.list[[3]], scores.list[[4]], scores.list[[5]], scores.list[[6]])
# plotting
g <- ggplot(my.scores, aes(x=CellType, y=Score, fill=CellType))
g <- g + geom_violin()
g <- g + geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.3)
g <- g + stat_summary(fun.data=median_hilow, geom="pointrange", color="gray75")
ggsave(paste(figdir, "Violin.celltype.scores.png", sep="\\"), width=8, height=6, dpi=600)

# locate points with highest Ba.Score
ba.scores <- subset(my.scores, CellType %in% c("Ba"))
ba.top.cells <- rownames(subset(ba.scores, Score > quantile(ba.scores$Score, 0.8)))
# mark cells
my.spring.ba.top.cells <- spring.n.celltype[, c("X","Y")]
my.spring.ba.top.cells$cellid <- rownames(my.spring.ba.top.cells)
my.spring.ba.top.cells$Mark <- with(my.spring.ba.top.cells, ifelse(cellid %in% ba.top.cells, "yes", "no"))
# plotting
g <- ggplot(my.spring.ba.top.cells, aes(x=X, y=Y, color=Mark))
g <- g + geom_point(shape=19, size=2, alpha=0.6)
g <- g + theme_classic()
ggsave(paste(figdir, "SPRING.top.ba.cells.png", sep="\\"), width=10, height=8, dpi=600)
#--------------------------------------------------------------------------------------------------------------------#

#------------------------------------- create own function (all markers) --------------------------------------------#
# infer cell types based on selected marker genes
# ALL marker genes should express otherwise a cell will not be considered as the given cell type
# expression matrix
my.exp <- as.matrix(matt@assays[["RNA"]]@data)
my.find.cells <- function(texp, tmarkers){
  # find cells with all markers expressed
  # count how many marker genes were expressed in a cell
  tcounts <- colSums(texp[tmarkers,] > 0)
  # return cells with all markers expressed
  return(names(tcounts[tcounts == length(tmarkers)]))
}

# apply the function to six cell types
ba.cells <- my.find.cells(my.exp, ba.markers)
meg.cells <- my.find.cells(my.exp, meg.markers)
mpp.cells <- my.find.cells(my.exp, mpp.markers)
ly.cells <- my.find.cells(my.exp, ly.markers)
m.cells <- my.find.cells(my.exp, m.markers)
gn.cells <- my.find.cells(my.exp, gn.markers)

length(ba.cells)
# 105
length(meg.cells)
# 372
length(mpp.cells)
# 0
length(ly.cells)
# 10
length(m.cells)
# 5
length(gn.cells)
# 370

# check if any overlapping between different inferred cell sets
all.cell.hits <- c(ba.cells, meg.cells, mpp.cells, ly.cells, m.cells, gn.cells)
all.cell.hits.counts <- as.data.frame.array(table(all.cell.hits))
colnames(all.cell.hits.counts) <- "counts"
duplicate.cells <- rownames(subset(all.cell.hits.counts, counts > 1))
length(duplicate.cells)
# 14 cells
# which cell type population do those 14 cells in
for (cid in duplicate.cells){
  print(cid)
  if (cid %in% ba.cells){
    print("ba.cells")
  }
  if (cid %in% meg.cells){
    print("meg.cells")
  }
  if (cid %in% mpp.cells){
    print("mpp.cells")
  }
  if (cid %in% ly.cells){
    print("ly.cells")
  }
  if (cid %in% m.cells){
    print("m.cells")
  }
  if (cid %in% gn.cells){
    print("gn.cells")
  }
}
# all duplicates fall into below two cell type populations
# [1] "ba.cells"
# [1] "meg.cells"

# merge cell annotations for plotting
# five target cell types, (MPP no hits)
cohort <- c("Ba", "Meg", "Ly", "M", "GN")
my.cell.type.list <- list()
k <- 1
for (y in list(ba.cells, meg.cells, ly.cells, m.cells, gn.cells)){
  my.cell.type.list[[k]] <- data.frame(cellid=y[! y %in% duplicate.cells])
  my.cell.type.list[[k]]$cohort <- cohort[k]
  k <- k + 1
}
# unknown cell types
all.cells <- rownames(spring.n)
unknown.cells <- data.frame(cellid=c(all.cells[! all.cells %in% all.cell.hits], duplicate.cells), cohort="Unknown")
# merge
my.cell.type <- rbind(my.cell.type.list[[1]],my.cell.type.list[[2]],my.cell.type.list[[3]],my.cell.type.list[[4]],my.cell.type.list[[5]],unknown.cells)
rownames(my.cell.type) <- my.cell.type$cellid
my.cell.type <- my.cell.type[,-1,drop=F]

# mark cells on SPRING plot
my.cell.type.forPlot <- merge(spring.n, my.cell.type, by=0, all=T)
rownames(my.cell.type.forPlot) <- my.cell.type.forPlot$Row.names
my.cell.type.forPlot <- my.cell.type.forPlot[,-1]
dim(my.cell.type.forPlot)
# 4678   3

# sort cells
my.cell.type.forPlot$cohort <- factor(my.cell.type.forPlot$cohort, levels=c("Unknown","Ba","GN","Ly","M","Meg"))
# plot
my.cols <- c("Unknown"="lightgray", "Ba"="tan3", "GN"="royalblue1", "Ly"="lightskyblue2", "M"="Olivedrab3", "Meg"="red4")
g <- ggplot(my.cell.type.forPlot.sorted, aes(x=X, y=Y, color=cohort))
g <- g + geom_point(shape=19, size=2, alpha=0.6)
g <- g + theme_classic()
g <- g + scale_color_manual(values=my.cols)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.markers.png", sep="\\"), width=10, height=8, dpi=600)
g <- g + facet_wrap(~cohort, ncol=3)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.markers.separated.png", sep="\\"), width=15, height=8, dpi=600)

# rename cell type
rename.cell.type <- c("Others", "Ba", "GN", "Ly", "M", "Mk")
names(rename.cell.type) <- c("Unknown", "Ba", "GN", "Ly", "M", "Meg")
my.cell.type.forPlot.sorted$cohort.rename <- rename.cell.type[my.cell.type.forPlot.sorted$cohort]

# add PRMT1 expression
my.cell.type.forPlot.sorted.exp <- merge(my.cell.type.forPlot.sorted, cell.info.merged[,c("nFeature_RNA","nCount_RNA","orig.ident","percent.mito","percent.ribo","PRMT1")], by=0, all=T)
rownames(my.cell.type.forPlot.sorted.exp) <- my.cell.type.forPlot.sorted.exp$Row.names
my.cell.type.forPlot.sorted.exp <- my.cell.type.forPlot.sorted.exp[,-1]
dim(my.cell.type.forPlot.sorted.exp)

# change color scheme
my.cell.type.forPlot.sorted.exp$cohort.rename <- factor(my.cell.type.forPlot.sorted.exp$cohort.rename, levels=c("Others","Ba","GN","Ly","M","Mk"))
my.cell.type.forPlot.sorted.exp <- my.cell.type.forPlot.sorted.exp[with(my.cell.type.forPlot.sorted.exp, order(cohort.rename)),]

my.cols.v2 <- c("Others"="lightgray", "Ba"="steelblue1", "GN"="slategray3", "Ly"="navy", "M"="red", "Mk"="sienna1")
g <- ggplot(my.cell.type.forPlot.sorted.exp, aes(x=X, y=Y, color=cohort.rename))
g <- g + geom_point(shape=19, size=3, alpha=0.6)
g <- g + theme_classic()
g <- g + theme(legend.title=element_blank(), legend.text=element_text(size=18))
g <- g + scale_color_manual(values=my.cols.v2)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.markers.v2.png", sep="\\"), width=10, height=8, dpi=600)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.markers.v2.svg", sep="\\"), width=10, height=8)
g <- g + facet_wrap(~orig.ident, ncol=2)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.markers.v2.split.png", sep="\\"), width=10, height=8, dpi=600)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.markers.v2.split.svg", sep="\\"), width=10, height=8)

# draw figure showing PRMT1 expression
my.cell.type.forPlot.sorted.exp <- my.cell.type.forPlot.sorted.exp[with(my.cell.type.forPlot.sorted.exp, order(PRMT1)),]
g <- ggplot(my.cell.type.forPlot.sorted.exp, aes(x=X, y=Y, color=PRMT1))
g <- g + geom_point(shape=19, size=3, alpha=0.6)
g <- g + theme_classic()
g <- g + scale_color_gradient(low="lightgray", high="red4")
g <- g + theme(legend.title=element_blank(), legend.text=element_text(size=18))
ggsave(paste(figdir, "SPRING.exp_PRMT1.v2.png", sep="\\"), width=10, height=8, dpi=600)
ggsave(paste(figdir, "SPRING.exp_PRMT1.v2.svg", sep="\\"), width=10, height=8)
g <- g + facet_wrap(~orig.ident, ncol=2)
ggsave(paste(figdir, "SPRING.exp_PRMT1.v2.split.png", sep="\\"), width=10, height=8, dpi=600)
ggsave(paste(figdir, "SPRING.exp_PRMT1.v2.split.svg", sep="\\"), width=10, height=8)

# save cell type results
write.table(my.cell.type.forPlot.sorted.exp, file=paste(infodir, "info.table.manual_celltype_all.txt", sep="\\"), quote=FALSE, na="", sep="\t", col.names=NA)

# add 4 more genes
for(tgene in c('ITGA2B','GP1BA','PF4','IFNAR1')){
  print(tgene)
  cell.info.merged[,'Gene'] <- as.vector(matt@assays[["RNA"]]@data[tgene,rownames(cell.info.merged)])
  cell.info.merged <- cell.info.merged[with(cell.info.merged, order(Gene)),]
  g <- ggplot(cell.info.merged, aes(x=-X, y=-Y, color=Gene))
  g <- g + geom_point(shape=19, size=3, alpha=0.6)
  g <- g + theme_classic()
  g <- g + scale_color_gradient(low="lightgray", high="red4")
  g <- g + theme(legend.title=element_blank(), legend.text=element_text(size=18))
  ggsave(paste(figdir, paste("SPRING.exp_",tgene,"v2.png",sep='.'), sep="\\"), width=10, height=8, dpi=600)
  ggsave(paste(figdir, paste("SPRING.exp_",tgene,"v2.svg",sep='.'), sep="\\"), width=10, height=8)
  g <- g + facet_wrap(~orig.ident, ncol=2)
  ggsave(paste(figdir, paste("SPRING.exp_",tgene,"v2.split.png",sep='.'), sep="\\"), width=10, height=8, dpi=600)
  ggsave(paste(figdir, paste("SPRING.exp_",tgene,"v2.split.svg",sep='.'), sep="\\"), width=10, height=8)
}

#--------------------------------------------------------------------------------------------------------------------#

#------------------------------------- create own function (n-1 markers) --------------------------------------------#
# infer cell types based on selected marker genes
# n-1 marker genes should express otherwise a cell will not be considered as the given cell type
# expression matrix
my.exp <- as.matrix(matt@assays[["RNA"]]@data)
my.find.cells.v2 <- function(texp, tmarkers){
  # find cells with all markers expressed
  # count how many marker genes were expressed in a cell
  tcounts <- colSums(texp[tmarkers,] > 0)
  # return cells with all markers expressed
  return(names(tcounts[tcounts == length(tmarkers)-1]))
}

# apply the function to six cell types
ba.cells.v2 <- my.find.cells.v2(my.exp, ba.markers)
meg.cells.v2 <- my.find.cells.v2(my.exp, meg.markers)
mpp.cells.v2 <- my.find.cells.v2(my.exp, mpp.markers)
ly.cells.v2 <- my.find.cells.v2(my.exp, ly.markers)
m.cells.v2 <- my.find.cells.v2(my.exp, m.markers)
gn.cells.v2 <- my.find.cells.v2(my.exp, gn.markers)

length(ba.cells.v2)
# 439
length(meg.cells.v2)
# 446
length(mpp.cells.v2)
# 11
length(ly.cells.v2)
# 23
length(m.cells.v2)
# 99
length(gn.cells.v2)
# 276

# check if any overlapping between different inferred cell sets
all.cell.hits.v2 <- c(ba.cells.v2, meg.cells.v2, mpp.cells.v2, ly.cells.v2, m.cells.v2, gn.cells.v2)
all.cell.hits.counts.v2 <- as.data.frame.array(table(all.cell.hits.v2))
colnames(all.cell.hits.counts.v2) <- "counts"
duplicate.cells.v2 <- rownames(subset(all.cell.hits.counts.v2, counts > 1))
length(duplicate.cells.v2)
# 109 cells
# which cell type population do those 109 cells in
for (cid in duplicate.cells.v2){
  print(cid)
  if (cid %in% ba.cells.v2){
    print("ba.cells")
  }
  if (cid %in% meg.cells.v2){
    print("meg.cells")
  }
  if (cid %in% mpp.cells.v2){
    print("mpp.cells")
  }
  if (cid %in% ly.cells.v2){
    print("ly.cells")
  }
  if (cid %in% m.cells.v2){
    print("m.cells")
  }
  if (cid %in% gn.cells.v2){
    print("gn.cells")
  }
}

# merge cell annotations for plotting
# six target cell types
cohort.v2 <- c("Ba", "Meg", "MPP", "Ly", "M", "GN")
my.cell.type.list.v2 <- list()
k <- 1
for (y in list(ba.cells.v2, meg.cells.v2, mpp.cells.v2, ly.cells.v2, m.cells.v2, gn.cells.v2)){
  my.cell.type.list.v2[[k]] <- data.frame(cellid=y[! y %in% duplicate.cells.v2])
  my.cell.type.list.v2[[k]]$cohort <- cohort.v2[k]
  k <- k + 1
}
# unknown cell types
all.cells <- rownames(spring.n.celltype)
unknown.cells.v2 <- data.frame(cellid=c(all.cells[! all.cells %in% all.cell.hits.v2], duplicate.cells.v2), cohort="Unknown")
# merge
my.cell.type.v2 <- rbind(my.cell.type.list.v2[[1]],my.cell.type.list.v2[[2]],my.cell.type.list.v2[[3]],my.cell.type.list.v2[[4]],my.cell.type.list.v2[[5]],my.cell.type.list.v2[[6]],unknown.cells.v2)
rownames(my.cell.type.v2) <- my.cell.type.v2$cellid
my.cell.type.v2 <- my.cell.type.v2[,-1,drop=F]
dim(my.cell.type.v2)
# 4678    1
table(my.cell.type.v2$cohort)
#  Ba      GN      Ly       M     Meg     MPP Unknown 
# 355     246      22      72     372       7    3604 

# mark cells on SPRING plot
my.cell.type.forPlot.v2 <- merge(spring.n.celltype, my.cell.type.v2, by=0, all=T)
dim(my.cell.type.forPlot.v2)
# 4678   11
# sort cells
my.cell.type.forPlot.v2$cohort <- factor(my.cell.type.forPlot.v2$cohort, levels=c("Unknown","Ba","GN","Ly","MPP","M","Meg"))
my.cell.type.forPlot.v2.sorted <- my.cell.type.forPlot.v2[with(my.cell.type.forPlot.v2, order(cohort)),]
# plot
my.cols.v2 <- c("Unknown"="lightgray", "Ba"="tan3", "GN"="royalblue1", "Ly"="lightskyblue2", "MPP"="lightsalmon2", "M"="Olivedrab3", "Meg"="red4")
g <- ggplot(my.cell.type.forPlot.v2.sorted, aes(x=X, y=Y, color=cohort))
g <- g + geom_point(shape=19, size=2, alpha=0.6)
g <- g + theme_classic()
g <- g + scale_color_manual(values=my.cols.v2)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.n-1.markers.png", sep="\\"), width=10, height=8, dpi=600)
g <- g + facet_wrap(~cohort, ncol=3)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.n-1.markers.separated.png", sep="\\"), width=15, height=8, dpi=600)
#--------------------------------------------------------------------------------------------------------------------#

#------------------------------------ create own function (n-1/n markers) -------------------------------------------#
# infer cell types based on selected marker genes
# what's new?
# - add two new markers for M
# - for Ly and M: n-1 marker genes should express otherwise a cell will not be considered as the given cell type
# - for Ba, GN, MPP, Meg: ALL(coord) marker genes should express otherwise a cell will not be considered as the given cell type
# expression matrix
my.exp <- as.matrix(matt@assays[["RNA"]]@data)

# all(coord) markers should express
my.find.cells <- function(texp, tmarkers){
  # find cells with all markers expressed
  # count how many marker genes were expressed in a cell
  tcounts <- colSums(texp[tmarkers,] > 0)
  # return cells with all markers expressed
  return(names(tcounts[tcounts == length(tmarkers)]))
}

# n-1 markers should express
my.find.cells.v2 <- function(texp, tmarkers){
  # find cells with all markers expressed
  # count how many marker genes were expressed in a cell
  tcounts <- colSums(texp[tmarkers,] > 0)
  # return cells with all markers expressed
  return(names(tcounts[tcounts == length(tmarkers)-1]))
}

# updated markers for M
m.markers.v2 <- intersect(rownames(matt@assays[["RNA"]]@data), c("CSF1R","CCR2","CTSS","IRF8"))
length(m.markers.v2)
# 4

# apply the function to six cell types
# all (coord)
ba.cells.v3 <- my.find.cells(my.exp, ba.markers)
meg.cells.v3 <- my.find.cells(my.exp, meg.markers)
mpp.cells.v3 <- my.find.cells(my.exp, mpp.markers)
gn.cells.v3 <- my.find.cells(my.exp, gn.markers)
# n-1
ly.cells.v3 <- my.find.cells.v2(my.exp, ly.markers)
m.cells.v3 <- my.find.cells.v2(my.exp, m.markers.v2)

length(ba.cells.v3)
# 105
length(meg.cells.v3)
# 372
length(mpp.cells.v3)
# 0
length(ly.cells.v3)
# 23
length(m.cells.v3)
# 20
length(gn.cells.v3)
# 370

# check if any overlapping between different inferred cell sets
all.cell.hits.v3 <- c(ba.cells.v3, meg.cells.v3, ly.cells.v3, m.cells.v3, gn.cells.v3)
all.cell.hits.counts.v3 <- as.data.frame.array(table(all.cell.hits.v3))
colnames(all.cell.hits.counts.v3) <- "counts"
duplicate.cells.v3 <- rownames(subset(all.cell.hits.counts.v3, counts > 1))
length(duplicate.cells.v3)
# 15 cells
# which cell type population do those 109 cells in
for (cid in duplicate.cells.v3){
  print(cid)
  if (cid %in% ba.cells.v3){
    print("ba.cells")
  }
  if (cid %in% meg.cells.v3){
    print("meg.cells")
  }
  if (cid %in% mpp.cells.v3){
    print("mpp.cells")
  }
  if (cid %in% ly.cells.v3){
    print("ly.cells")
  }
  if (cid %in% m.cells.v3){
    print("m.cells")
  }
  if (cid %in% gn.cells.v3){
    print("gn.cells")
  }
}

# 1    m.cells & gn.cells
# 14   ba.cells & meg.cells

# merge cell annotations for plotting
# five target cell types, (MPP no hits)
cohort.v3 <- c("Ba", "Meg", "Ly", "M", "GN")
my.cell.type.list.v3 <- list()
k <- 1
for (y in list(ba.cells.v3, meg.cells.v3, ly.cells.v3, m.cells.v3, gn.cells.v3)){
  my.cell.type.list.v3[[k]] <- data.frame(cellid=y[! y %in% duplicate.cells.v3])
  my.cell.type.list.v3[[k]]$cohort <- cohort.v3[k]
  k <- k + 1
}
# unknown cell types
all.cells <- rownames(cell.info.merged)
unknown.cells.v3 <- data.frame(cellid=c(all.cells[! all.cells %in% all.cell.hits.v3], duplicate.cells.v3), cohort="Unknown")
# merge
my.cell.type.v3 <- rbind(my.cell.type.list.v3[[1]],my.cell.type.list.v3[[2]],my.cell.type.list.v3[[3]],my.cell.type.list.v3[[4]],my.cell.type.list.v3[[5]],unknown.cells.v3)
rownames(my.cell.type.v3) <- my.cell.type.v3$cellid
my.cell.type.v3 <- my.cell.type.v3[,-1,drop=F]
dim(my.cell.type.v3)
# 4678    1
table(my.cell.type.v3$cohort)
# Ba      GN      Ly       M     Meg Unknown 
# 91     369      23      19     358    3818 

# mark cells on SPRING plot
my.cell.type.forPlot.v3 <- merge(spring.n, my.cell.type.v3, by=0, all=T)
dim(my.cell.type.forPlot.v3)
# 4678   4
# sort cells
my.cell.type.forPlot.v3$cohort <- factor(my.cell.type.forPlot.v3$cohort, levels=c("Unknown","Ba","GN","Ly","M","Meg"))
my.cell.type.forPlot.v3.sorted <- my.cell.type.forPlot.v3[with(my.cell.type.forPlot.v3, order(cohort)),]
# plot
my.cols.v3 <- c("Unknown"="lightgray", "Ba"="tan3", "GN"="royalblue1", "Ly"="lightskyblue2", "M"="Olivedrab3", "Meg"="red4")
g <- ggplot(my.cell.type.forPlot.v3.sorted, aes(x=X, y=Y, color=cohort))
g <- g + geom_point(shape=19, size=2, alpha=0.6)
g <- g + theme_classic()
g <- g + scale_color_manual(values=my.cols.v3)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.markers.v3.png", sep="\\"), width=10, height=8, dpi=600)
g <- g + facet_wrap(~cohort, ncol=3)
ggsave(paste(figdir, "SPRING.cell.types.manual.by.markers.v3.separated.png", sep="\\"), width=15, height=8, dpi=600)
#--------------------------------------------------------------------------------------------------------------------#

######################################################################################################################

#-------------------------------------- expression on D8_CD42bCD41 cells --------------------------------------------#
# get D8_CD42bCD41 cells
tcells <- rownames(subset(matt@meta.data, orig.ident %in% c("D8_CD42bCD41")))
length(tcells)
# 1077 cells

# get expression
texp <- as.matrix(matt@assays[["RNA"]]@data)
dim(texp)
# 23844  4678

# draw scatter plot highlighting a pair of genes:
my.scatter.plot <- function(geneA, geneB, expMat, tlabel=NULL){
  # check if gene is available
  if (! geneA %in% rownames(expMat)){
    print(paste("gene",geneA,"not","available","due","to","limit","number","of","cells","expressing","it",sep=" "))
    return()
  }
  if (! geneB %in% rownames(expMat)){
    print(paste("gene",geneB,"not","available","due","to","limit","number","of","cells","expressing","it",sep=" "))
    return()
  }
  # extract expression
  tdata <- as.data.frame(t(expMat[c(geneA,geneB),tcells]))
  colnames(tdata) <- c("geneA", "geneB")
  g <- ggplot(tdata, aes(x=geneA, y=geneB))
  g <- g + geom_point(shape=19, size=2, color="royalblue2", alpha=0.6)
  g <- g + theme_classic()
  g <- g + xlab(geneA) + ylab(geneB)
  toutfile <- paste(figdir, paste("scatter","plot","exp",geneA,geneB,"png",sep="."), sep="\\")
  if (! is.null(tlabel)){
    toutfile <- paste(figdir, paste("scatter","plot","exp",geneA,geneB,tlabel,"png",sep="."), sep="\\")
  }
  ggsave(toutfile, width=6, height=5, dpi=300)
  #return(NULL)
}

my.scatter.plot("PRMT1","ITGA2B", texp, "D8_CD42bCD41")
my.scatter.plot("PRMT1","GP1BA", texp, "D8_CD42bCD41")
my.scatter.plot("PRMT1","PF4", texp, "D8_CD42bCD41")
my.scatter.plot("PRMT1","IFNAR1", texp, "D8_CD42bCD41")
my.scatter.plot("ITGA2B","GP1BA", texp, "D8_CD42bCD41")

# try MAGIC correction on D8_CD42bCD41 cells
tcounts <- matt@assays[["RNA"]]@counts[,tcells]
# remove all-zero genes
tcounts <- tcounts[rowSums(tcounts)>0,]
dim(tcounts)
# 17880  1077

# switch two dimensions
tcounts.t <- t(as.matrix(tcounts))

# filter data
# keep genes expressed in at least 10 cells
keep_cols <- colSums(tcounts.t > 0) > 10
sum(keep_cols)
# 10656
tcounts.t.clean <- tcounts.t[, keep_cols]
dim(tcounts.t.clean)
# 1077 10656

# look at the distribution of library sizes
g <- ggplot() + geom_histogram(aes(x=rowSums(tcounts.t.clean)), bins=50)
g <- g + geom_vline(xintercept=1000, color="red")
g <- g + theme_classic()
ggsave(paste(figdir, "histogram.nCount_RNA.per.cell.D8_CD42bCD41.png", sep="\\"), width=8, height=6, dpi=300)

# use all cells
# # keep cells with at least 1000 UMIs and at most 40000 UMIs
# keep_rows <- rowSums(tcounts.t.clean) > 1000 & rowSums(tcounts.t.clean) < 40000
# sum(keep_rows)
# # 801
# tcounts.t.clean <- tcounts.t.clean[keep_rows,]
# dim(tcounts.t.clean)
# # 801 11778

# normalize data
tcounts.norm <- library.size.normalize(tcounts.t.clean)
tcounts.norm <- sqrt(tcounts.norm)

# run MAGIC
# on PRMT1 gene
tcounts.magic <- magic(tcounts.norm, genes=c("PRMT1","VIM","GAPDH"))
# before
g <- ggplot(as.data.frame(tcounts.norm)) + geom_point(aes(x=VIM, y=GAPDH, color=PRMT1))
g <- g + scale_color_viridis_c(option="B")
g <- g + theme_classic()
ggsave(paste(figdir,"exp.PRMT1.before_magic.D8_CD42bCD41.png",sep="\\"), width=8, height=6, dpi=300)
# after
g <- ggplot(tcounts.magic) + geom_point(aes(x=VIM, y=GAPDH, color=PRMT1))
g <- g + scale_color_viridis_c(option="B")
g <- g + theme_classic()
ggsave(paste(figdir,"exp.PRMT1.after_magic.D8_CD42bCD41.png",sep="\\"), width=8, height=6, dpi=300)

# all_genes
tcounts.magic.all <- magic(tcounts.norm, genes="all_genes", init=tcounts.magic)
# output
texp.magic <- t(as.matrix(tcounts.magic.all$result))
dim(texp.magic)
# 10656  1077
write.table(texp.magic, file=gzfile(paste(infodir,"exp.sqrt.after_magic.D8_CD42bCD41.tsv.gz",sep="\\")), quote=FALSE, na="", sep="\t", col.names=NA)
###write.table(t(as.matrix((counts.magic.all$result)^2)), file=gzfile(paste(infodir,"exp.pseudo_counts.after_magic.tsv.gz",sep="/")), quote=FALSE, na="", sep="\t", col.names=NA)

# redraw scatter plot, using MAGIC corrected expression
my.scatter.plot("PRMT1","ITGA2B", texp.magic, "after_MAGIC.D8_CD42bCD41")
my.scatter.plot("PRMT1","GP1BA", texp.magic, "after_MAGIC.D8_CD42bCD41")
my.scatter.plot("PRMT1","PF4", texp.magic, "after_MAGIC.D8_CD42bCD41")
my.scatter.plot("PRMT1","IFNAR1", texp.magic, "after_MAGIC.D8_CD42bCD41")
my.scatter.plot("ITGA2B","GP1BA", texp.magic, "after_MAGIC.D8_CD42bCD41")

# save magic object
save(tcounts.magic.all, file=paste(infodir, "tcounts.magic.D8_CD42bCD41.Robj", sep="\\"))

# calculate Spearman correlation between each pair of genes
# after magic
scc <- cor(x=t(texp.magic), method="spearman")
dim(scc)
# 10656 10656
# rank genes by their correlation with PRMT1
write.table(data.frame(PRMT1 = sort(scc["PRMT1",], decreasing=T)), file=paste(infodir, "spearman.correlation.D8_CD42bCD41.after_magic.PRMT1.txt",sep="\\"), quote=FALSE, na="", sep="\t", col.names=F, row.names=T)

# before magic
texp <- as.matrix(matt@assays[["RNA"]]@data)[,tcells]
texp <- texp[rowSums(texp)>0,]
dim(texp)
# 17880  1077
scc.raw <- cor(x=t(texp), method="spearman")
dim(scc.raw)
# 17880 17880
# rank genes by their correlation with PRMT1
write.table(data.frame(PRMT1 = sort(scc.raw["PRMT1",], decreasing=T)), file=paste(infodir, "spearman.correlation.D8_CD42bCD41.before_magic.PRMT1.txt",sep="\\"), quote=FALSE, na="", sep="\t", col.names=F, row.names=T)

######################################################################################################################

#--------------------------------------- MAGIC on other three conditions --------------------------------------------#
for (cdt in c("D0_unstimulated", "D8_other", "D8_stimulated")){
  print(paste("processing","condition",cdt,sep=" "))
  tcells <- rownames(subset(matt@meta.data, orig.ident %in% c(cdt)))
  print(paste(length(tcells),"cells",sep=" "))
  # try MAGIC correction on selected cells
  tcounts <- matt@assays[["RNA"]]@counts[,tcells]
  # remove all-zero genes
  tcounts <- tcounts[rowSums(tcounts)>0,]
  print(dim(tcounts))
  
  # switch two dimensions
  tcounts.t <- t(as.matrix(tcounts))
  
  # filter data
  # keep genes expressed in at least 10 cells
  keep_cols <- colSums(tcounts.t > 0) > 10
  sum(keep_cols)
  # 10656
  tcounts.t.clean <- tcounts.t[, keep_cols]
  print(dim(tcounts.t.clean))
  # 1077 10656
  
  # look at the distribution of library sizes
  g <- ggplot() + geom_histogram(aes(x=rowSums(tcounts.t.clean)), bins=50)
  g <- g + geom_vline(xintercept=1000, color="red")
  g <- g + theme_classic()
  ggsave(paste(figdir, paste("histogram","nCount_RNA","per","cell",cdt,"png", sep="."), sep="\\"), width=8, height=6, dpi=300)
  
  # use all cells
  # # keep cells with at least 1000 UMIs and at most 40000 UMIs
  # keep_rows <- rowSums(tcounts.t.clean) > 1000 & rowSums(tcounts.t.clean) < 40000
  # sum(keep_rows)
  # # 801
  # tcounts.t.clean <- tcounts.t.clean[keep_rows,]
  # dim(tcounts.t.clean)
  # # 801 11778
  
  # normalize data
  tcounts.norm <- library.size.normalize(tcounts.t.clean)
  tcounts.norm <- sqrt(tcounts.norm)
  
  # run MAGIC
  # on PRMT1 gene
  tcounts.magic <- magic(tcounts.norm, genes=c("PRMT1","VIM","GAPDH"))
  # before
  g <- ggplot(as.data.frame(tcounts.norm)) + geom_point(aes(x=VIM, y=GAPDH, color=PRMT1))
  g <- g + scale_color_viridis_c(option="B")
  g <- g + theme_classic()
  ggsave(paste(figdir,paste("exp","PRMT1","before_magic",cdt,"png",sep="."),sep="\\"), width=8, height=6, dpi=300)
  # after
  g <- ggplot(tcounts.magic) + geom_point(aes(x=VIM, y=GAPDH, color=PRMT1))
  g <- g + scale_color_viridis_c(option="B")
  g <- g + theme_classic()
  ggsave(paste(figdir,paste("exp","PRMT1","after_magic",cdt,"png",sep="."),sep="\\"), width=8, height=6, dpi=300)
  
  # all_genes
  tcounts.magic.all <- magic(tcounts.norm, genes="all_genes", init=tcounts.magic)
  # output
  texp.magic <- t(as.matrix(tcounts.magic.all$result))
  print(dim(texp.magic))
  # 10656  1077
  write.table(texp.magic, file=gzfile(paste(infodir,paste("exp","sqrt","after_magic",cdt,"tsv","gz",sep="."),sep="\\")), quote=FALSE, na="", sep="\t", col.names=NA)
  ###write.table(t(as.matrix((counts.magic.all$result)^2)), file=gzfile(paste(infodir,"exp.pseudo_counts.after_magic.tsv.gz",sep="/")), quote=FALSE, na="", sep="\t", col.names=NA)
  
  # redraw scatter plot, using MAGIC corrected expression
  my.scatter.plot("PRMT1","ITGA2B", texp.magic, paste("after_MAGIC",cdt,sep="."))
  my.scatter.plot("PRMT1","GP1BA", texp.magic, paste("after_MAGIC",cdt,sep="."))
  my.scatter.plot("PRMT1","PF4", texp.magic, paste("after_MAGIC",cdt,sep="."))
  my.scatter.plot("PRMT1","IFNAR1", texp.magic, paste("after_MAGIC",cdt,sep="."))
  my.scatter.plot("ITGA2B","GP1BA", texp.magic, paste("after_MAGIC",cdt,sep="."))
  
  # save magic object
  save(tcounts.magic.all, file=paste(infodir, paste("tcounts","magic","all",cdt,"Robj",sep="."), sep="\\"))
}

# D0_unstimulated: gene PF4 not available due to limit number of cells expressing it

######################################################################################################################

#--------------------------------------- plot SPRING with color by Paper --------------------------------------------#
# run Fig1b_functions to infer cell types (python scripts)
# load the colors
my.paper.colors <- read.table(colorfile.fromPaper, header=F, check.names=F, sep="\t")
rownames(my.paper.colors) <- rownames(coord)
colnames(my.paper.colors) <- c("R","G","B")
dim(my.paper.colors)

# convert RGB values (red, green and blue) to hex code
my.paper.colors.hex <- apply(my.paper.colors, 1, FUN=function(x) {rgb(x[1],x[2],x[3])})
length(my.paper.colors.hex)
# 4678

# load cell plotting order
# 0-based --> 1-based
my.cell.order <- read.table(orderfile.fromPaper, header=F, check.names=F, sep="\t")$V1 + 1
length(my.cell.order)
# 4678

# prepare data
tdata <- n
tdata$CellID <- rownames(tdata)
# add condition info
tdata.merged <- merge(tdata, matt@meta.data[rownames(tdata),c("orig.ident"),drop=F], by=0, all=T)
rownames(tdata.merged) <- tdata.merged$Row.names
tdata.merged <- tdata.merged[rownames(tdata),-1]# keep the original cell order

# make plots
# all cells
g <- ggplot(tdata.merged[my.cell.order,], aes(x=-X, y=-Y, color=CellID))
g <- g + geom_point(shape=19, size=4, alpha=0.6)
g <- g + scale_color_manual(values=my.paper.colors.hex[my.cell.order])
g <- g + theme(legend.position="none")
g <- g + theme(axis.title=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
ggsave(paste(figdir,"SPRING.Plot.use_paper_codes.ggplot.png",sep="\\"), width=10, height=8, dpi=600)
ggsave(paste(figdir,"SPRING.Plot.use_paper_codes.ggplot.svg",sep="\\"), width=10, height=8)
# separate cells from 4 conditions
g <- ggplot(tdata.merged[my.cell.order,], aes(x=-X, y=-Y, color=CellID))
g <- g + geom_point(shape=19, size=4, alpha=0.6) + facet_wrap(~orig.ident, ncol=2)
g <- g + scale_color_manual(values=my.paper.colors.hex[my.cell.order])
g <- g + theme(legend.position="none")
g <- g + theme(axis.title=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.ticks=element_blank())
ggsave(paste(figdir,"SPRING.Plot.use_paper_codes.ggplot.separated.png",sep="\\"), width=20, height=16, dpi=600)
ggsave(paste(figdir,"SPRING.Plot.use_paper_codes.ggplot.separated.svg",sep="\\"), width=20, height=16)

######################################################################################################################

sessionInfo()


