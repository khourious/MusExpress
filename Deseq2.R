# carregando pacotes
library(sva)
library(DESeq2)
library(readr)
library(dplyr)
library(calibrate)
library(genefilter)
library(gplots)
library(tidyr)
library(PCAtools)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(cowplot)
library(factoextra)
library(BiocParallel)
library(pheatmap)

getwd()
setwd("C:/Users/.../")
# Daddos processados seguindo pipeline SARA para dados single-end 

# 1. --- carregando arquivos ----
gene_counts <- read.csv("counts/gene_count_matrix.csv", row.names = 1)
colnames(gene_counts) <- gsub("^X", "", colnames(gene_counts))
new_row_names <- sub("^[^|]*\\|", "", row.names(gene_counts))

# ---- Verificar se há duplicatas ----
if (any(duplicated(new_row_names))) {
  # Adicionar um sufixo numérico para duplicatas
  new_row_names <- make.unique(new_row_names)}

# --- Ajustar os nomes nas colunas
row.names(gene_counts) <- new_row_names

sample_ids <- gsub("\\.\\d+$", "", colnames(gene_counts))
keep <- !duplicated(sample_ids)
gene_counts_clean <- gene_counts[, keep]
colnames(gene_counts_clean) <- sample_ids[keep]

colnames(gene_counts_clean) <- gsub("_S\\d+$", "", colnames(gene_counts_clean))

# ---- Conferir resultado
colnames(gene_counts_clean)

# ---- Carregar metadados ----
info_Balb <- read_delim("info_BALBC.txt", delim = "\t", 
                   escape_double = FALSE, trim_ws = TRUE)

row.names(info_Balb) <- info_Balb$Sample

# ---- Filtrar apenas as amostras "BALBC" ----
gene_counts_Balbc <- gene_counts_clean[, colnames(gene_counts_clean) %in% info_Balb$Sample]

info_Balb <- info_Balb[match(colnames(gene_counts_Balbc), info_Balb$Sample), ]

# 2. --- DeSeq2 ---
# --- criar o DESeqDataSet com contagens brutas ---
dds <- DESeqDataSetFromMatrix(countData = gene_counts_Balbc, 
                              colData = info_Balb, 
                              design = ~ Time)

# --- Not required ---
dds <- dds[rowSums(counts(dds)) >= 5, ]
dds <- dds[rowSums(counts(dds) > 50) > 3, ] 


# --- Wald test ----
dds <- DESeq(dds)
resultsNames(dds)

# 3. --- PCA ---

# --- Variance stabilizing transformation ---
vsd <- vst(dds, blind=FALSE)

# --- Plotando PCA ---
plotPCA(vsd, intgroup="Time")
plotPCA(vsd, intgroup="Time", pc = c(1,3))

# --- Organizando as informações do PCA
pcaobj <- prcomp(t(assay(vsd)))
pcaobjto <- prcomp(assay(vsd))

# --- Plotando eigene values dos PCx ---
png("eig.balbc.png", w=660, h=480, pointsize=40)
fviz_eig(pcaobj, 
         addlabels = TRUE, 
         ylim = c(0, 50))
dev.off()

# --- Transformando para novo plot ---
pca <- prcomp(t(assay(vsd)))
af <- cbind(pca$x[,1:3]) %>% as.data.frame()
af$PC1 <- as.numeric(af$PC1) / (pca$sdev[1] * sqrt(nrow(info_Balb)))
af$PC2 <- as.numeric(af$PC2) / (pca$sdev[2] * sqrt(nrow(info_Balb)))
af$PC3 <- as.numeric(af$PC3) / (pca$sdev[3] * sqrt(nrow(info_Balb)))
info_Balb$Time <- factor(info_Balb$Time,
  levels = c("D0", "D3", "D7", "D21","D35", "D70"))
af[,"Time"] <- as.factor(info_Balb$Time)
af[,"sampleId"]  <- as.factor(info_Balb$Sample)

# --- Plot PC1 + PC2 ---
pca <- ggplot(af, aes(PC1, PC2, colour = `Time`)) +
  geom_point(size = 4.5) +
  scale_color_manual(values = c("D0" = "#252525", "D3" = "#ff7f00", 
                                "D7" = "#b2182b", "D21" = "#6a51a3",
                                "D35" = "#41ab5d", "D70" = "#2171b5")) +
  geom_text_repel(aes(label = sampleId ), size = 2.5, box.padding = 0.5, point.padding = 0.3) +
  theme_bw () +
  labs(title = "PCA",
       x = "PC1: 54% variance",
       y = "PC2: 12% variance") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

pca #650, 500

# --- Add density curves to y and x axis ---
xdens <- axis_canvas(pca, axis = "x") + 
  geom_density(data = af, aes(x = PC1, fill = Type, colour = Type), alpha = 0.3) +
  scale_fill_manual(values=c("#4F94CD", "#EE6363")) +
  scale_color_manual(values = c("#104E8B", "#CD2626"))
ydens <- axis_canvas(pca, axis = "y", coord_flip = TRUE) + 
  geom_density(data = af, aes(x = PC2, fill = Type, colour = Type), alpha = 0.3) +
  scale_fill_manual(values=c("#4F94CD", "#EE6363")) +
  scale_color_manual(values = c("#104E8B", "#CD2626")) +
  coord_flip()
pdf("PCA.zikv.ctrl.pdf", w=8, h=6, pointsize=20)
png("PCA.zikv.ctrl.pdf", w=825, h=650, pointsize=20)
# W = 825, H = 650
pca %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()
dev.off()

# --- Plot PC1 + PC3 ---
pca <- ggplot(af, aes(PC1, PC3, colour = Time)) +
  geom_point(size = 4.5) +
  scale_color_manual(values = c("D0" = "#252525", "D3" = "#ff7f00", 
                                "D7" = "#b2182b", "D21" = "#6a51a3",
                                "D35" = "#41ab5d", "D70" = "#2171b5")) +
  geom_text_repel(aes(label = sampleId ), size = 2.5, box.padding = 0.5, point.padding = 0.3) +
  theme_bw() +
  labs(title = "PCA",
       x = "PC1: 54% variance",
       y = "PC3: 10% variance")  +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
# See
pca

# --- Add density curves to y and x axis ---
xdens <- axis_canvas(pca, axis = "x") + 
  geom_density(data = af, aes(x = PC1, fill = `Type`, colour = `Type`), alpha = 0.3) +
  scale_fill_manual(values=c("#4F94CD", "#EE6363")) +
  scale_color_manual(values = c("#104E8B", "#CD2626"))
ydens <- axis_canvas(pca, axis = "y", coord_flip = TRUE) + 
  geom_density(data = af, aes(x = PC3, fill = `Type`, colour = `Type`), alpha = 0.3) +
  scale_fill_manual(values=c("#4F94CD", "#EE6363")) +
  scale_color_manual(values = c("#104E8B", "#CD2626")) +
  coord_flip()
pdf("PCA.zikv.ctrl.2.pdf", w=8, h=6, pointsize=20)
png("PCA.zikv.ctrl.pdf", w=825, h=650, pointsize=20)
# W = 825, H = 650
pca %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()
dev.off()


# 4. --- Plotar QCs ---
# --- Plot de dispersão ---
pdf("qc-dispersions.pdf", 50, 50, pointsize=100)
plotDispEsts(dds, main="dispersion plot")
dev.off()

# --- sample distance heatmap ---
distsRL <- as.matrix(dist(t(assay(vsd))))
hmcol <- colorRampPalette(c("darkgreen","yellow","orange", "red","red3"))(100)

# --- Definir as anotações do heatmap ---
annotation_col <- data.frame(
  Time = info_Balb$Time)
rownames(annotation_col) <- colnames(distsRL)

ann_colors <- list(Time = c("D0" = "#252525", "D3" = "#ff7f00", 
                   "D7" = "#b2182b", "D21" = "#6a51a3",
                   "D35" = "#41ab5d", "D70" = "#2171b5"))  


# --- Gerar o heatmap com anotação ----
pheatmap(distsRL, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = hmcol,  
         breaks = seq(0, max(distsRL), length.out = 101),
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         main = "Sample Distance Matrix")

# --- Dot Plot ---
par(mfrow = c( 1, 2))
dds <- estimateSizeFactors(dds)
par(mfrow=c(1,2))
plot(log2(1+counts(dds)[,1:2]),col=rgb(0,0,0,.2),pch=16,cex=1.0,main="log2")
plot(assay(vsd)[,1:2],col=rgb(0,0,0,.2),pch=16,cex=1.0,main="vsd")

# 5. ---- Capturar os resultados ----
resultsNames(dds)
# --- aplicar coeficiente nas comparações ----
resShrink <- lfcShrink(dds, coef="Time_D3_vs_D0", type = "apeglm")

# --- Pegar os resultados de DEGs ----
table(resShrink$padj<0.05)

# --- Ordenar por p-adj ---
resShrink <- resShrink[order(resShrink$padj), ]
resShrink <-na.omit(resShrink)

# --- Juntar com a contagem normalizada ---
resShrinkdata <- merge(as.data.frame(resShrink),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(resShrinkdata)[1] <- "GeneSymbol"
row.names(resShrinkdata) <- resShrinkdata$GeneSymbol
head(resShrinkdata)

# --- Volcano plot ---
EnhancedVolcano(resShrink, lab=rownames(resShrink), 
                x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=2, pointSize = 3.0, labSize = 5.0)


# --- Volcano plot mais bonito ---- 
data <- resShrinkdata
D3.down <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1),]
D3.up <-data[(data$log2FoldChange >= 1 & data$padj <= 0.05),]
D3.degs <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1) | 
                    (data$log2FoldChange >= 1 & data$padj <= 0.05),]

keyvals <- ifelse(data$log2FoldChange < -1 & data$padj <= 0.05, '#253494',
                  ifelse(data$log2FoldChange > 1 & data$padj <= 0.05, '#cb181d', 'grey70'))
keyvals[is.na(keyvals)] <- 'grey70'
names(keyvals)[keyvals == '#cb181d'] <- 'Up-regulated'
names(keyvals)[keyvals == '#253494'] <- 'Down-regulated'
names(keyvals)[keyvals == 'grey70'] <- 'Not significant'

data_sub <- data[(data$padj <= 0.00001 & abs(data$log2FoldChange) >= 6), ]

EnhancedVolcano(data,
                lab = rownames(data),
                title = 'CHIKV BALBC: D3 vs D0',
                subtitle = " ",
                selectLab = rownames(data_sub),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3,
                labSize = 4,
                colCustom = keyvals,
                xlim = c(-12, 12),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


# --- Salvar resultados ---
write.csv(D3.down, file="D3-down.csv")
write.csv(D3.up, file="D3-up.csv")
write.csv(D3.degs, file="D3-degs.csv")
write.csv(data, file="D3-results.csv")

# ---- D7 vs D0

# decrease the fold change noise with shrinkage function
resultsNames(dds)
resShrink <- lfcShrink(dds, coef="Time_D7_vs_D0", type = "apeglm")

# get differential expression results
table(resShrink$padj<0.05)

# order by adjusted p-value
resShrink <- resShrink[order(resShrink$padj), ]
resShrink <-na.omit(resShrink)

# merge with normalized count data
resShrinkdata <- merge(as.data.frame(resShrink),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(resShrinkdata)[1] <- "GeneSymbol"
row.names(resShrinkdata) <- resShrinkdata$GeneSymbol
head(resShrinkdata)

# Volcano plot
EnhancedVolcano(resShrink, lab=rownames(resShrink), 
                x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=2, pointSize = 3.0, labSize = 5.0)


# DEGs (volcano plot)
data <- resShrinkdata
D7.down <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1),]
D7.up <-data[(data$log2FoldChange >= 1 & data$padj <= 0.05),]
D7.degs <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1) | 
                  (data$log2FoldChange >= 1 & data$padj <= 0.05),]

keyvals <- ifelse(data$log2FoldChange < -1 & data$padj <= 0.05, '#253494',
                  ifelse(data$log2FoldChange > 1 & data$padj <= 0.05, '#cb181d', 'grey70'))
keyvals[is.na(keyvals)] <- 'grey70'
names(keyvals)[keyvals == '#cb181d'] <- 'Up-regulated'
names(keyvals)[keyvals == '#253494'] <- 'Down-regulated'
names(keyvals)[keyvals == 'grey70'] <- 'Not significant'

data_sub <- data[(data$padj <= 0.000000000000001 & data$log2FoldChange >= 4 | 
                    data$padj <= 0.001 & data$log2FoldChange <= -4), ]

EnhancedVolcano(data,
                lab = rownames(data),
                title = 'CHIKV BALBC: D7 vs D0',
                subtitle = " ",
                selectLab = rownames(data_sub),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3,
                labSize = 4,
                colCustom = keyvals,
                xlim = c(-12, 12),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


# write results
write.csv(D7.down, file="D7-down.csv")
write.csv(D7.up, file="D7-up.csv")
write.csv(D7.degs, file="D7-degs.csv")
write.csv(data, file="D7-results.csv")


# ---- D21 vs D0

# decrease the fold change noise with shrinkage function
resultsNames(dds)
resShrink <- lfcShrink(dds, coef="Time_D21_vs_D0", type = "apeglm")

# get differential expression results
table(resShrink$padj<0.05)

# order by adjusted p-value
resShrink <- resShrink[order(resShrink$padj), ]
resShrink <-na.omit(resShrink)

# merge with normalized count data
resShrinkdata <- merge(as.data.frame(resShrink),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(resShrinkdata)[1] <- "GeneSymbol"
row.names(resShrinkdata) <- resShrinkdata$GeneSymbol
head(resShrinkdata)

# Volcano plot
EnhancedVolcano(resShrink, lab=rownames(resShrink), 
                x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=2, pointSize = 3.0, labSize = 5.0)


# DEGs (volcano plot)
data <- resShrinkdata
D21.down <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1),]
D21.up <-data[(data$log2FoldChange >= 1 & data$padj <= 0.05),]
D21.degs <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1) | 
                  (data$log2FoldChange >= 1 & data$padj <= 0.05),]

keyvals <- ifelse(data$log2FoldChange < -1 & data$padj <= 0.05, '#253494',
                  ifelse(data$log2FoldChange > 1 & data$padj <= 0.05, '#cb181d', 'grey70'))
keyvals[is.na(keyvals)] <- 'grey70'
names(keyvals)[keyvals == '#cb181d'] <- 'Up-regulated'
names(keyvals)[keyvals == '#253494'] <- 'Down-regulated'
names(keyvals)[keyvals == 'grey70'] <- 'Not significant'

data_sub <- data[(data$padj <= 0.0000000001 & data$log2FoldChange >= 2 | 
                    data$padj <= 0.00001 & data$log2FoldChange <= -2), ]

EnhancedVolcano(data,
                lab = rownames(data),
                title = 'CHIKV BALBC: D21 vs D0',
                subtitle = " ",
                selectLab = rownames(data_sub),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3,
                labSize = 4,
                colCustom = keyvals,
                xlim = c(-12, 12),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


# write results
write.csv(D21.down, file="D21-down.csv")
write.csv(D21.up, file="D21-up.csv")
write.csv(D21.degs, file="D21-degs.csv")
write.csv(data, file="D21-results.csv")


# ---- D35 vs D0

# decrease the fold change noise with shrinkage function
resultsNames(dds)
resShrink <- lfcShrink(dds, coef="Time_D35_vs_D0", type = "apeglm")

# get differential expression results
table(resShrink$padj<0.05)

# order by adjusted p-value
resShrink <- resShrink[order(resShrink$padj), ]
resShrink <-na.omit(resShrink)

# merge with normalized count data
resShrinkdata <- merge(as.data.frame(resShrink),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(resShrinkdata)[1] <- "GeneSymbol"
row.names(resShrinkdata) <- resShrinkdata$GeneSymbol
head(resShrinkdata)

# Volcano plot
EnhancedVolcano(resShrink, lab=rownames(resShrink), 
                x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=2, pointSize = 3.0, labSize = 5.0)


# DEGs (volcano plot)
data <- resShrinkdata
D35.down <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1),]
D35.up <-data[(data$log2FoldChange >= 1 & data$padj <= 0.05),]
D35.degs <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1) | 
                   (data$log2FoldChange >= 1 & data$padj <= 0.05),]

keyvals <- ifelse(data$log2FoldChange < -1 & data$padj <= 0.05, '#253494',
                  ifelse(data$log2FoldChange > 1 & data$padj <= 0.05, '#cb181d', 'grey70'))
keyvals[is.na(keyvals)] <- 'grey70'
names(keyvals)[keyvals == '#cb181d'] <- 'Up-regulated'
names(keyvals)[keyvals == '#253494'] <- 'Down-regulated'
names(keyvals)[keyvals == 'grey70'] <- 'Not significant'

data_sub <- data[(abs(data$log2FoldChange) >= 3), ]

EnhancedVolcano(data,
                lab = rownames(data),
                title = 'CHIKV BALBC: D35 vs D0',
                subtitle = " ",
                selectLab = rownames(data_sub),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3,
                labSize = 4,
                colCustom = keyvals,
                xlim = c(-12, 12),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


# write results
write.csv(D35.down, file="D35-down.csv")
write.csv(D35.up, file="D35-up.csv")
write.csv(D35.degs, file="D35-degs.csv")
write.csv(data, file="D35-results.csv")

# ---- D70 vs D0

# decrease the fold change noise with shrinkage function
resultsNames(dds)
resShrink <- lfcShrink(dds, coef="Time_D70_vs_D0", type = "apeglm")

# get differential expression results
table(resShrink$padj<0.05)

# order by adjusted p-value
resShrink <- resShrink[order(resShrink$padj), ]
resShrink <-na.omit(resShrink)

# merge with normalized count data
resShrinkdata <- merge(as.data.frame(resShrink),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(resShrinkdata)[1] <- "GeneSymbol"
row.names(resShrinkdata) <- resShrinkdata$GeneSymbol
head(resShrinkdata)

# Volcano plot
EnhancedVolcano(resShrink, lab=rownames(resShrink), 
                x='log2FoldChange', y='padj', pCutoff=0.05, 
                FCcutoff=2, pointSize = 3.0, labSize = 5.0)


# DEGs (volcano plot)
data <- resShrinkdata
D70.down <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1),]
D70.up <-data[(data$log2FoldChange >= 1 & data$padj <= 0.05),]
D70.degs <- data[(data$padj <= 0.05 & data$log2FoldChange <= -1) | 
                   (data$log2FoldChange >= 1 & data$padj <= 0.05),]

keyvals <- ifelse(data$log2FoldChange < -1 & data$padj <= 0.05, '#253494',
                  ifelse(data$log2FoldChange > 1 & data$padj <= 0.05, '#cb181d', 'grey70'))
keyvals[is.na(keyvals)] <- 'grey70'
names(keyvals)[keyvals == '#cb181d'] <- 'Up-regulated'
names(keyvals)[keyvals == '#253494'] <- 'Down-regulated'
names(keyvals)[keyvals == 'grey70'] <- 'Not significant'

data_sub <- data[(abs(data$log2FoldChange) >= 3), ]

EnhancedVolcano(data,
                lab = rownames(data),
                title = 'CHIKV BALBC: D70 vs D0',
                subtitle = " ",
                selectLab = rownames(data_sub),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3,
                labSize = 4,
                colCustom = keyvals,
                xlim = c(-12, 12),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


# write results
write.csv(D70.down, file="D70-down.csv")
write.csv(D70.up, file="D70-up.csv")
write.csv(D70.degs, file="D70-degs.csv")
write.csv(data, file="D70-results.csv")
