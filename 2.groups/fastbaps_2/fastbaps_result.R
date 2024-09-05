# Title     : TODO
# Objective : TODO
# Created by: A
# Created on: 2021/12/2

options(device="windows")


library(ape)
library(maps)
library(fastbaps)
library(ggtree)
library(phytools)
library(ggplot2)

# 从faalignment文件开始计算
#   导入数据
fasta.file.name <- "D:\\python\\3.MLST_SC\\5.MLST\\2.groups\\fastbaps\\all_core_genes.fa"
sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#   设定先验参数
#   type: 'symmetric', 'baps', 'optimise.symmetric', 'optimise.baps', 'hc'
sparse.data <- optimise_prior(sparse.data, type = 'hc')
#   进行fast_baps运算，设定初始的k值，期待的分组结果 , k.init=
baps.hc <- fast_baps(sparse.data)
best.partition <- best_baps_partition(sparse.data, baps.hc)

# 把计算结果放到进化树上显示
newick.file.name <- "D:\\python\\3.MLST_SC\\5.MLST\\2.groups\\fastbaps\\all_sc_mlst.treefile"
iqtree <- phytools::read.newick(newick.file.name)
plot.df <- data.frame(id = colnames(sparse.data$snp.matrix), fastbaps = best.partition,
    stringsAsFactors = FALSE)
# 保存分组结果为excel
library(xlsx)
path <- "D:\\python\\3.MLST_SC\\5.MLST\\2.groups\\fastbaps"
setwd(path)
write.xlsx(plot.df, "groups_fastbaps.xlsx")

# 可视化
gg <- ggtree(iqtree)
f2 <- facet_plot(gg, panel = "fastbaps", data = plot.df, geom = geom_tile, aes(x = fastbaps),
    color = "blue")

# 保存图片
library(Ipaper)
write_fig(f2, file = "groups_fastbaps.pdf", width = 10, height = 5, devices = NULL, res = 300, show = TRUE)



