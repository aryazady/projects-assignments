library(GEOquery)
library(limma)
library(Biobase)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(Rtsne)


series <- "GSE48558"
platform <- "GPL6244"

curDir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(sub(paste0("/", sub("(.+)/", "", curDir)), "", curDir))

gset <- getGEO(series, GSEMatrix=TRUE, AnnotGPL=TRUE, destdir = 'Data/')

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
 

gr1 <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
              "XXXXXXXXXXXXXXXXXX5X4XXX1X1225X4XX44XX44X5X4X5X4X3",
              "XXX3XXX3XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111004000",
              "55555554222214444444")
gr2 <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
              "XXXXXXXXXXXXXXXXXX1X1XXX1X1111X1XX11XX11X1X1X1X1X1",
              "XXX1XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111001000",
              "11111111111111111111")
gr_normal <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX2XXX2XXXXX",
                    "XXXXXXXXXXXXXXXXXX1X0XXX2X2331X0XX00XX00X1X0X1X0X4",
                    "XXX4XXX4XXXXXXXXXXXXXXXXXXXXXXXXXXXXX2222222XX0XXX",
                    "11111110333320000000")

gr1 <- strsplit(gr1, split="")[[1]]
gr2 <- strsplit(gr2, split="")[[1]]
gr_normal <- strsplit(gr_normal, split="")[[1]]


sel <- which(gr1 != "X")
gr1 <- gr1[sel]
gset1 <- gset[ ,sel]

sel <- which(gr2 != "X")
gr2 <- gr2[sel]
gset2 <- gset[ ,sel]

sel <- which(gr_normal != "X")
gr_normal <- gr_normal[sel]
gset_normal <- gset[ ,sel]

ex1 = exprs(gset1)
ex2 = exprs(gset2)
ex_normal = exprs(gset_normal)

### checking if data is logarithmic and are normalized
pdf("Result/Plots/boxplot.pdf", width =  19.375, height = 27.40625)
boxplot(ex2, main="Samples Boxplot", horizontal=TRUE)
dev.off()

gset1$group <- gr1
gset2$group <- gr2
gset_normal$group <- gr_normal
gr1 <- factor(gr1)
gr2 <- factor(gr2)
gr_normal <- factor(gr_normal)
groups1 <- make.names(c("AML","Granul","Mono","CD34+","T Cell","B Cell"))
groups2 <- make.names(c("AML","Normal"))
groups_normal <- make.names(c("T Cell","B Cell","Granulocytes","Monocytes","CD34+"))
levels(gr1) <- groups1
levels(gr2) <- groups2
levels(gr_normal) <- groups_normal

pdf("Result/Plots/PHeatMap.pdf")
pheatmap(cor(ex2), labels_row = gr2, labels_col = gr2, border_color = NA, color = bluered(500), fontsize=7, main = "AML-Normal")
pheatmap(cor(ex1), labels_row = gr1, labels_col = gr1, border_color = NA, color = bluered(500), fontsize=7, main = "AML-Detailed")
pheatmap(cor(ex_normal), labels_row = gr_normal, labels_col = gr_normal, border_color = NA, color = bluered(500), fontsize=7, main = "Normal")
dev.off()

### unscaled pca
pc1 <- prcomp(ex1)
pc2 <- prcomp(ex2)

pdf("Result/Plots/pc_not_scaled.pdf")
plot(pc1, main = "PC Vectors Variance")
plot(pc1$x[,1:2], main = "PC Gene Expression")
dev.off()

### scaling gene expression
ex1.scale <- t(scale(t(ex1), scale = F))
ex2.scale <- t(scale(t(ex2), scale = F))
ex_normal.scale <- t(scale(t(ex_normal), scale = F))


pc1 <- prcomp(ex1.scale)
pc2 <- prcomp(ex2.scale)

pc_normal <- prcomp(ex_normal.scale)
tsne_normal_r <- Rtsne(t(ex_normal.scale),dims=3, perplexity=16)

# tsne_x1 <- Rtsne(ex1.scale, dims=2, perplexity=50)
tsne_r1 <- Rtsne(t(ex1.scale),dims=3, perplexity=22, theta = 0.0)

# tsne_x2 <- Rtsne(ex2.scale, dims=2, perplexity=50)
tsne_r2 <- Rtsne(t(ex2.scale),dims=3, perplexity=22, theta = 0.0)

pdf("Result/Plots/pHeatMap_2D.pdf")
pheatmap(cor(t(tsne_normal_r$Y)), labels_row = gr_normal, labels_col = gr_normal,
         border_color = NA, color = bluered(500), fontsize=7,
         main = "Normal Samples Correlation After Dimentionality Reduction")
pheatmap(cor(t(tsne_r1$Y)), labels_row = gr1, labels_col = gr1,
         border_color = NA, color = bluered(500), fontsize=7,
         main = "AML-Detailed Correlation After Dimentionality Reduction")
dev.off()

pdf("Result/Plots/pc_scaled.pdf")
plot(pc2, main = "PC Vectors Variance (Scaled)")
plot(pc2$x[,1:2], main = "PC Gene Expression (Scaled)")
# plot(pc1, main = "AML-Detailed-prcomp")
# plot(pc1$x[,1:2], main = "AML-Detailed-prcomp")
# 
# plot(tsne_x2$Y, main = "AML-Normal-tSNE")
# plot(tsne_x1$Y, main = "AML-Detailed-tSNE")
dev.off()

pcr1 <- data.frame(pc1$rotation[,1:3], Group = gr1)
pcr2 <- data.frame(pc2$rotation[,1:3], Group = gr2)

tsne1 <- data.frame(tsne_r1$Y, Group = gr1)
tsne2 <- data.frame(tsne_r2$Y, Group = gr2)

pdf("Result/Plots/clustered_ samples.pdf")
ggplot(pcr2, aes(PC1, PC2, color=Group)) + geom_point(size = 4) + theme_bw() + ggtitle("AML-Normal-prcomp")
ggplot(pcr1, aes(PC1, PC2, color=Group)) + geom_point(size = 4) + theme_bw() + ggtitle("AML-Detailed-prcomp")

ggplot(tsne2, aes(X1, X2, color=Group)) + geom_point(size = 4) + theme_bw() + ggtitle("AML-Normal-tSNE")
ggplot(tsne1, aes(X1, X2, color=Group)) + geom_point(size = 4) + theme_bw() + ggtitle("AML-Detailed-tSNE")
dev.off()

design <- model.matrix(~group + 0, gset1)

colnames(design) <- levels(gr1)

fit <- lmFit(gset1, design)

cont.matrix_mono <- makeContrasts(AML-Mono, levels=design)
fit2_mono <- contrasts.fit(fit, cont.matrix_mono)
fit2_mono <- eBayes(fit2_mono, 0.01)

cont.matrix_cd34 <- makeContrasts(AML-CD34., levels=design)
fit2_cd34 <- contrasts.fit(fit, cont.matrix_cd34)
fit2_cd34 <- eBayes(fit2_cd34, 0.01)

# cont.matrix_b.cell <- makeContrasts(AML-B.Cell, levels=design)
# fit2_b.cell <- contrasts.fit(fit, cont.matrix_b.cell)
# fit2_b.cell <- eBayes(fit2_b.cell, 0.01)

# cont.matrix_cd34_mono <- makeContrasts(CD34.-Mono, levels=design)
# fit2_cd34_mono <- contrasts.fit(fit, cont.matrix_cd34_mono)
# fit2_cd34_mono <- eBayes(fit2_cd34_mono, 0.01)
# tT_c_m <- topTable(fit2_cd34_mono, adjust="fdr", sort.by = "B", number = Inf)
# tT_c_m <- subset(tT_c_m, select=c("Gene.ID","Gene.symbol","adj.P.Val","logFC"))
# up.gene <- subset(tT_c_m, logFC > 1 & adj.P.Val < 0.05)
# up.gene <- unique(as.character(strsplit2(up.gene$Gene.symbol, "///")))
# write.table(up.gene, "Result/Significant Changes in Gene Expression/test_.txt", row.names=F, col.names = F, quote = F)


tT_mono <- topTable(fit2_mono, adjust="fdr", sort.by="B", number = Inf)
tT_cd34 <- topTable(fit2_cd34, adjust="fdr", sort.by="B", number = Inf)
# tT_b.cell <- topTable(fit2_b.cell, adjust="fdr", sort.by="B", number = Inf)

tT_mono <- subset(tT_mono, select=c("Gene.ID","Gene.symbol","adj.P.Val","logFC"))
tT_cd34 <- subset(tT_cd34, select=c("Gene.ID","Gene.symbol","adj.P.Val","logFC"))
# tT_b.cell <- subset(tT_b.cell, select=c("Gene.ID","Gene.symbol","adj.P.Val","logFC"))

write.table(tT_mono, "Result/Significant Changes in Gene Expression//AML_Mono.txt", row.names=F, sep="\t", quote = F)
write.table(tT_cd34, "Result/Significant Changes in Gene Expression//AML_CD34+.txt", row.names=F, sep="\t", quote = F)
# write.table(tT_b.cell, "Result/Significant Changes in Gene Expression//AML_B-Cell.txt", row.names=F, sep="\t", quote = F)

aml.up <- subset(tT_mono, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, "Result/Significant Changes in Gene Expression//AML_Mono_Up_Genes.txt", row.names=F, col.names = F, quote = F)

aml.down <- subset(tT_mono, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, "Result/Significant Changes in Gene Expression//AML_Mono_Down_Genes.txt", row.names=F, col.names = F, quote = F)

aml.up <- subset(tT_cd34, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, "Result/Significant Changes in Gene Expression//AML_CD34+_Up_Genes.txt", row.names=F, col.names = F, quote = F)

aml.down <- subset(tT_cd34, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, "Result/Significant Changes in Gene Expression//AML_CD34+_Down_Genes.txt", row.names=F, col.names = F, quote = F)

# aml.up <- subset(tT_b.cell, logFC > 1 & adj.P.Val < 0.05)
# aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
# write.table(aml.up.genes, "Result/Significant Changes in Gene Expression//AML_B-Cell_Up_Genes.txt", row.names=F, col.names = F, quote = F)
# 
# aml.down <- subset(tT_b.cell, logFC < -1 & adj.P.Val < 0.05)
# aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
# write.table(aml.down.genes, "Result/Significant Changes in Gene Expression//AML_B-Cell_Down_Genes.txt", row.names=F, col.names = F, quote = F)
