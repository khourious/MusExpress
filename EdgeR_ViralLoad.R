library(dplyr)
library(tidyr)
library(ggplot2)
library(edgeR)

# 1. --- Carregar arquivos ----
# --- Contagem do genoma de CHIKV nas amostras ---
chikv_samples_counts <- read.delim("featurecounts/chikv_samples_counts.txt", 
                                   comment.char="#")

# --- Pegar os nomes das amostras ---
colnames(chikv_samples_counts)

# --- Remover prefixo longo ---
colnames(chikv_samples_counts) <- gsub(
  "X.mnt.HDD.1.CHIKV.", "", 
  colnames(chikv_samples_counts)
)

# --- Remover extensão .hisat2.sorted.bam ---
colnames(chikv_samples_counts) <- gsub(
  ".hisat2.sorted.bam", "", 
  colnames(chikv_samples_counts)
)

# --- Remover sufixo _S seguido de número ---
colnames(chikv_samples_counts) <- gsub(
  "_S[0-9]+", "", 
  colnames(chikv_samples_counts)
)

# --- Conferir resultado ---
colnames(chikv_samples_counts)

# --- Remover colunas de anotação ---
counts <- chikv_samples_counts[, -(1:6)]
linhagem <- ifelse(grepl("Balbc", colnames(counts), ignore.case = TRUE), "Balbc",
                   ifelse(grepl("Black", colnames(counts), ignore.case = TRUE), "Black", NA))

# --- Verificar se todas foram classificadas ---
table(linhagem)

# 2. --- Separar as colunas em duas matrizes
counts_balbc <- counts[, linhagem == "Balbc", drop = FALSE]
counts_black <- counts[, linhagem == "Black", drop = FALSE]

# 3. --- Normalizar ---

# --- BALBC ---
# --- Converter em matriz ---
counts_balbc <- as.matrix(counts_balbc)

# --- Identificar amostras com soma zero ---
col_sums_balbc <- colSums(counts_balbc)
zero_samples_balbc <- names(col_sums_balbc[col_sums_balbc == 0])

# ---Filtrar apenas amostras válidas ---
counts_filtered_balbc <- counts_balbc[, col_sums_balbc > 0]

# --- Criar objeto DGEList e calcular CPM ----
dge_balbc <- DGEList(counts = counts_filtered_balbc)
cpm_filtered_balbc <- cpm(dge_balbc, log = FALSE)

# --- Inserir contagens zero ---
# --- Agora criar uma matriz CPM completa, com todas as amostras ---
# --- Inicialmente, uma matriz de zeros com mesma dimensão ---
cpm_counts_balbc <- matrix(0, nrow = nrow(counts_balbc), ncol = ncol(counts_balbc))
rownames(cpm_counts_balbc) <- rownames(counts_balbc)
colnames(cpm_counts_balbc) <- colnames(counts_balbc)

# --- Preencher CPM apenas para amostras válidas ---
cpm_counts_balbc[, col_sums_balbc > 0] <- cpm_filtered_balbc

# --- checar ---
head(cpm_counts_balbc)

# --- BLACK ---
# --- Converter em matriz ---
counts_black <- as.matrix(counts_black)

# --- Identificar amostras com soma zero ---
col_sums_black <- colSums(counts_black)
zero_samples_black <- names(col_sums_black[col_sums_black == 0])

# ---Filtrar apenas amostras válidas ---
counts_filtered_black <- counts_black[, col_sums_black > 0]

# --- Criar objeto DGEList e calcular CPM ----
dge_black <- DGEList(counts = counts_filtered_black)
cpm_filtered_black <- cpm(dge_black, log = FALSE)

# --- Inserir contagens zero ---
# --- Agora criar uma matriz CPM completa, com todas as amostras ---
# --- Inicialmente, uma matriz de zeros com mesma dimensão ---
cpm_counts_black <- matrix(0, nrow = nrow(counts_black), ncol = ncol(counts_black))
rownames(cpm_counts_black) <- rownames(counts_black)
colnames(cpm_counts_black) <- colnames(counts_black)

# --- Preencher CPM apenas para amostras válidas ---
cpm_counts_black[, col_sums_black > 0] <- cpm_filtered_black

# --- checar ---
head(cpm_counts_black)


# 4. --- VIsualizar os dados ---
library(reshape2)
library(ggplot2)
library(stringr)

# --- BALBC ---
# --- Transformar CPM em formato longo ---
cpm_long_balbc <- melt(cpm_counts_balbc)
colnames(cpm_long_balbc) <- c("gene", "sample", "cpm")

# --- Extrair grupo e tempo dos nomes das amostras ---
cpm_long_balbc$grupo <- str_extract(cpm_long_balbc$sample, "BALBC|BLACK")
cpm_long_balbc$tempo <- str_extract(cpm_long_balbc$sample, "D[0-9]+")

# --- Definir ordem dos tempos ---
cpm_long_balbc$tempo <- factor(cpm_long_balbc$tempo, levels = c("D0","D3","D7","D21","D35","D70"))

# --- BLACK ---
# --- Transformar CPM em formato longo ---
cpm_long_black <- melt(cpm_counts_black)
colnames(cpm_long_black) <- c("gene", "sample", "cpm")

# --- Extrair grupo e tempo dos nomes das amostras ---
cpm_long_black$grupo <- str_extract(cpm_long_black$sample, "BALBC|BLACK")
cpm_long_black$tempo <- str_extract(cpm_long_black$sample, "D[0-9]+")

# --- Definir ordem dos tempos ---
cpm_long_black$tempo <- factor(cpm_long_black$tempo, levels = c("D0","D3","D7","D21","D35","D70"))

# 5. --- Criar gráficos ---
library(ggplot2)

# --- BALBC em rosa ---
ggplot(cpm_long_balbc, aes(x = tempo, y = cpm)) +
  geom_violin(fill = "#FF69B4", color = "white", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  theme_linedraw(base_size = 14) +
  labs(title = "CHIKV viral load - BALBC (CPM)",
       x = " ",
       y = "Counts") +
  scale_y_continuous(limits = c(1, max(cpm_long_balbc$cpm))) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )


ggplot(subset(cpm_long, grupo == "BALBC"), aes(x = tempo, y = cpm)) +
  geom_violin(fill = "#FF69B4", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  theme_linedraw(base_size = 14) +
  labs(title = "CHIKV viral load - BALBC (CPM)",
       x = "Days post infection",
       y = "Counts") +
  scale_y_continuous(limits = c(0, max(cpm_long$cpm))) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )

ggplot(subset(cpm_long, grupo == "BALBC"), aes(x = tempo, y = cpm)) +
  geom_violin(fill = "#FF69B4", color = "white", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  theme_linedraw(base_size = 14) +
  labs(title = "CHIKV viral load - BALBC (CPM)",
       x = "Days post infection",
       y = "Counts") +
  scale_y_continuous(limits = c(0, max(cpm_long$cpm))) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )


# --- BLACK em verde ---
ggplot(cpm_long_black, aes(x = tempo, y = cpm)) +
  geom_violin(fill = "#32CD32", color = "white", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  theme_linedraw(base_size = 14) +
  labs(title = "CHIKV viral load - BLACK (CPM)",
       x = " ",
       y = "Counts") +
  scale_y_continuous(limits = c(1, max(cpm_long_black$cpm))) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )

ggplot(subset(cpm_long, grupo == "BLACK"), aes(x = tempo, y = cpm)) +
  geom_violin(fill = "#32CD32", color = "black", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  theme_linedraw(base_size = 14) +
  labs(title = "CHIKV viral load - BLACK (CPM)",
       x = "Days post infection",
       y = "Counts") +
  scale_y_continuous(limits = c(0, max(cpm_long$cpm))) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )


ggplot(subset(cpm_long, grupo == "BLACK"), aes(x = tempo, y = cpm)) +
  geom_violin(fill = "#32CD32", color = "white", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  theme_linedraw(base_size = 14) +
  labs(title = "CHIKV viral load - BLACK (CPM)",
       x = "Days post infection",
       y = "Counts") +
  scale_y_continuous(limits = c(0, max(cpm_long$cpm))) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )
