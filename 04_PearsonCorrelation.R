# --- Carregar pacotes ---
library(dplyr)
library(tidyr)
library(PerformanceAnalytics)

# --- Carregue os dados ---
# --- Dados de quantificação viral ----
dados_virais <- read.delim("featurecounts/chikv_samples_counts.txt", 
                           comment.char="#", row.names = 1)

# --- Remover colunas técnicas ---
dados_virais <- dados_virais[, -(1:5)]

# ---Renomear colunas de amostras ---
colnames(dados_virais) <- colnames(dados_virais) %>%
  gsub("X.mnt.HDD.1.CHIKV.", "", .) %>%   # remove prefixo longo
  gsub(".hisat2.sorted.bam", "", .) %>%   # remove sufixo final
  gsub("_S[0-9]+", "", .)

# --- Dados de expressão normalizado das amostras de BALBc ---
expr_balbc <- read.delim("BALBC/balbc-cemitool.csv", row.names = 1)

# --- Dados de expressão normalizado das amostras de BLACK  ---
expr_black <- read.delim("BLACK/black-cemitool.csv", row.names = 1)

# ---- Medidas das patas ----
patas <- read.delim("patas.txt", 
                           comment.char="#")

# --- Organizar os dados ---

# --- Dados BALBC ---
#  --- Separar a carga viral  ---
viral_balbc <- dados_virais[, grep("BALBC", colnames(dados_virais))]

#  --- Separar a medida das patas  ---
patas_balbc <- patas[, grep("BALBC", colnames(patas))]

# --- Dados BLACK ---
#  --- Separar a carga viral  ---
viral_black <- dados_virais[, grep("BLACK", colnames(dados_virais))]

#  --- Separar a medida das patas  ---
patas_black <- patas[, grep("BLACK", colnames(patas))]


# --- Fluxo da correlação de Pearson ---
# --- Interseção dos nomes das amostras presentes em três objetos ---
common_samples <- Reduce(intersect, list(colnames(expr_balbc),
                                         colnames(patas_balbc),
                                         colnames(viral_balbc)))

# --- Seleciona apenas as colunas correspondentes ---
expr <- expr_balbc[, common_samples, drop = FALSE] # Garante que continue sendo uma matriz

# --- Extrai da matriz e converte em numérico ---
viral <- as.numeric(viral_balbc ["CHIKV", common_samples])  
patas <- as.numeric(patas_balbc[1, common_samples]) 


# --- Interseção dos nomes das amostras presentes em três objetos ---
common_samples <- Reduce(intersect, list(colnames(expr_black),
                                         colnames(patas_black),
                                         colnames(viral_black)))

# --- Seleciona apenas as colunas correspondentes ---
expr <- expr_black[, common_samples, drop = FALSE]# Garante que continue sendo uma matriz

# --- Extrai da matriz e converte em numérico ---
viral <- as.numeric(viral_black ["CHIKV", common_samples])  
patas <- as.numeric(patas_black[1, common_samples]) 



# --- Transformar expressão e carga viral em log10 (+1 para evitar log(0)) ---
expr_log <- log10(expr + 1)
viral_log <- log10(as.numeric(viral) + 1)
patas_log <- log10(as.numeric(patas) + 1)

# --- Calcular correlação gene x carga viral (ambos em log10) ---
res_list <- apply(expr_log, 1, function(gene_expr) {
  ct_viral <- cor.test(as.numeric(gene_expr), viral_log, method = "pearson")
  ct_patas <- cor.test(as.numeric(gene_expr), patas_log, method = "pearson")

  # --- Cria um vetor com os resultados principais ---
  c(cor_viral = ct_viral$estimate, pval_viral = ct_viral$p.value,
    cor_patas = ct_patas$estimate, pval_patas = ct_patas$p.value)
})

# --- Transforma lista em data frame ---
res_df <- as.data.frame(t(res_list))
rownames(res_df) <- rownames(expr)

# --- Ajusta os p-valores das correlações gene vs viral usando o método BH ---
res_df$padj_viral <- p.adjust(res_df$pval_viral, method = "BH")
res_df$padj_patas <- p.adjust(res_df$pval_patas, method = "BH")

# --- Calcula o R² ---
res_df$R2_viral <- res_df$cor_viral.cor^2
res_df$R2_patas <- res_df$cor_patas.cor^2

# --- Filtra de acordo com cutoff ---
genes_fortes <- rownames(res_df)[abs(res_df$R2_viral) > 0.64 ]#  & abs(res_df$pval_viral) < 0.05]
genes_fortes <- rownames(res_df)[abs(res_df$R2_patas) > 0.64 ]# & abs(res_df$pval_patas) < 0.05]


# --- Visualização ---
# --- Selecionar os top 10 nomes de genes ---
top10 <- rownames(res_df)[order(abs(res_df$R2_viral), decreasing = TRUE)][1:10]
top10 <- rownames(res_df)[order(abs(res_df$R2_patas), decreasing = TRUE)][1:10]

# --- filtrar os top 10  da carga viral ---
my_data <- data.frame(
  viral_load = viral,
  t(expr[top10, ])  
)

# --- filtrar os top 10  das patas ---
my_data <- data.frame(
  paw = patas,
  t(expr[top10, ])   # genes depois
)

# --- Transformar em log10 ---
my_data_log <- log10(my_data + 1) # consistência da escala

# --- Plotar gráfico ---
chart.Correlation(my_data_log, histogram = TRUE, pch = 19)



# --- Filtrar os genes com maior correlação nas tabelas --- 
# --- Tabelas serão aplicadas no IPA ---

# --- Carregar tabelas ---
D3.results_balbc <- read.csv("/BALBC/D3-results.csv", row.names = 1)
D7.results_balbc <- read.csv("/BALBC/D7-results.csv", row.names = 1)
D21.results_balbc <- read.csv("/BALBC/D21-results.csv", row.names = 1)
D35.results_balbc <- read.csv("/BALBC/D35-results.csv", row.names = 1)
D70.results_balbc <- read.csv("/BALBC/D70-results.csv", row.names = 1)


D3.results_black <- read.csv("/BLACK/D3-results.csv", row.names = 1)
D7.results_black <- read.csv("/BLACK/D7-results.csv", row.names = 1)
D21.results_black <- read.csv("/BLACK/D21-results.csv", row.names = 1)
D35.results_black <- read.csv("/BLACK/D35-results.csv", row.names = 1)
D70.results_black <- read.csv("/BLACK/D70-results.csv", row.names = 1)


# --- Preparar função ---
prepara_tabela <- function(df, genes_fortes, timepoint) {
  # Seleciona colunas principais
  df_sub <- df[, c("GeneSymbol", "log2FoldChange", "pvalue", "padj")] # selecionar colunas de interesse
  # Renomeia para incluir o timepoint
  colnames(df_sub)[2:4] <- paste0(colnames(df_sub)[2:4], "_", timepoint)
  
  # Adiciona fake_pval com sufixo do timepoint
  df_sub[[paste0("fake_pval_", timepoint)]] <- ifelse(df_sub$GeneSymbol %in% genes_fortes, 0.05, 1)
  
  return(df_sub)
}

# --- Aplicar função ----
D3  <- prepara_tabela(D3.results_balbc,  genes_fortes, "D3")
D7  <- prepara_tabela(D7.results_balbc,  genes_fortes, "D7")
D21 <- prepara_tabela(D21.results_balbc, genes_fortes, "D21")
D35 <- prepara_tabela(D35.results_balbc, genes_fortes, "D35")
D70 <- prepara_tabela(D70.results_balbc, genes_fortes, "D70")

D3  <- prepara_tabela(D3.results_black,  genes_fortes, "D3")
D7  <- prepara_tabela(D7.results_black,  genes_fortes, "D7")
D21 <- prepara_tabela(D21.results_black, genes_fortes, "D21")
D35 <- prepara_tabela(D35.results_black, genes_fortes, "D35")
D70 <- prepara_tabela(D70.results_black, genes_fortes, "D70")

# --- Combina por gene (rownames) ---
tabela_final <- Reduce(function(x, y) merge(x, y, by = "GeneSymbol", all = TRUE),
                       list(D3, D7, D21, D35, D70))

# --- Adiciona coluna GeneSymbol ---
rownames(tabela_final)  <- tabela_final$GeneSymbol 

# --- Reorganiza para GeneSymbol primeiro ---
tabela_final <- tabela_final[, c("GeneSymbol", setdiff(colnames(tabela_final), "GeneSymbol"))]

# --- Exportar ---
write.csv(tabela_final, "BALBC_cargaviral_R8.csv", row.names = FALSE)
write.csv(tabela_final, "BALBC_patas_R8.csv", row.names = FALSE)

write.csv(tabela_final, "BLACK_cargaviral_R8.csv", row.names = FALSE)
write.csv(tabela_final, "BLACK_patas_R8.csv", row.names = FALSE)


# --- Visualizaçãos ---
# --- Heatmaps 
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)


# --- Balbc ---
BALBC_cargaviral_R8 <- read.csv("/correlações/BALBC_cargaviral_R8.csv")
BALBC_patas_R8 <- read.csv("/correlações/BALBC_patas_R8.csv")

# --- Colunas de interesse ---
lfc_cols <- c("log2FoldChange_D3", "log2FoldChange_D7", 
              "log2FoldChange_D21", "log2FoldChange_D35", "log2FoldChange_D70")

# --- Genes detectados na correlação com carga viral ---
genes_viral <- BALBC_cargaviral_R8 %>%
  filter(fake_pval_D3 <= 0.05 |
           fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 |
           fake_pval_D35 <= 0.05 |
           fake_pval_D70 <= 0.05) %>%
  select(GeneSymbol, all_of(lfc_cols))

# --- Genes detectados na correlação com patas ---
genes_patas <- BALBC_patas_R8 %>%
  filter(fake_pval_D3 <= 0.05 |
           fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 |
           fake_pval_D35 <= 0.05 |
           fake_pval_D70 <= 0.05) %>%
  select(GeneSymbol, all_of(lfc_cols))

# --- Unir os dois conjuntos e remover duplicados ---
genes_union <- bind_rows(genes_viral, genes_patas) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  arrange(GeneSymbol)   # ordenar alfabeticamente

# --- Converter em tabela longa ---
df_long <- genes_union %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "Condition",
               values_to = "log2FC")

# --- Ajustar tabela ---
df_long$Condition <- gsub("log2FoldChange_", "", df_long$Condition)
df_long$log2FC <- as.numeric(df_long$log2FC)

# --- Garantir que log2FC é numérico ---
df_long$log2FC <- as.numeric(df_long$log2FC)

# --- Converter Condition em fator com ordem ---
df_long$Condition <- factor(df_long$Condition,
                            levels = c("D3", "D7", "D21", "D35", "D70"))

# --- Plot único heatmap Balbc ---
ggplot(df_long, aes(x = Condition, y = GeneSymbol, fill = log2FC)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradientn(
    colours = c("#225ea8", "#3690c0", "white", "#f16913", "#ce1256", "#ae017e"),
    values = scales::rescale(c(min(df_long$log2FC, na.rm = TRUE),
                               -2, 0, 2,
                               max(df_long$log2FC, na.rm = TRUE))),
    limits = range(df_long$log2FC, na.rm = TRUE),
    na.value = "grey10",
    name = "log2FC"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 9, angle = 0, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),  # mantém nomes dos genes
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Heatmap BALB/c")

# --- Genes detectados na correlação com carga viral ---
genes_viral <- BALBC_cargaviral_R8 %>%
  filter(fake_pval_D3 <= 0.05 |
           fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 |
           fake_pval_D35 <= 0.05 |
           fake_pval_D70 <= 0.05) %>%
  select(GeneSymbol, all_of(lfc_cols))

# --- Genes detectados na correlação com patas ---
genes_patas <- BALBC_patas_R8 %>%
  filter(fake_pval_D3 <= 0.05 |
           fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 |
           fake_pval_D35 <= 0.05 |
           fake_pval_D70 <= 0.05) %>%
  select(GeneSymbol, all_of(lfc_cols))

# --- Unir os dois conjuntos e remover duplicados ---
genes_union <- bind_rows(genes_viral, genes_patas) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  arrange(GeneSymbol)   # ordenar alfabeticamente

# --- Criar listas de genes detectados em cada correlação ---
genes_viral <- BALBC_cargaviral_R8 %>%
  filter(fake_pval_D3 <= 0.05 | fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 | fake_pval_D35 <= 0.05 | fake_pval_D70 <= 0.05) %>%
  pull(GeneSymbol)

genes_patas <- BALBC_patas_R8 %>%
  filter(fake_pval_D3 <= 0.05 | fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 | fake_pval_D35 <= 0.05 | fake_pval_D70 <= 0.05) %>%
  pull(GeneSymbol)

# --- União dos genes ---
genes_union <- sort(unique(c(genes_viral, genes_patas)))

# --- Criar tabela para o heatmap de origem ---
df_origem <- expand.grid(GeneSymbol = genes_union,
                         Correlacao = c("Viral Load", "Paw"),
                         stringsAsFactors = FALSE)

df_origem <- df_origem %>%
  mutate(presente = case_when(
    Correlacao == "Viral Load" & GeneSymbol %in% genes_viral ~ 1,
    Correlacao == "Paw" & GeneSymbol %in% genes_patas ~ 1,
    TRUE ~ 0
  ))


# --- Plot Heatmap Balbc ---
ggplot(df_origem, aes(x = Correlacao, y = GeneSymbol, fill = factor(presente))) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_manual(values = c("0" = "grey90", "1" = "#ce1256"),
                    name = "Correlação",
                    labels = c("Não det", "Det")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 6),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right") +
  labs(title = "Heatmap BALB/c")






# --- Black ---
BLACK_cargaviral_R8 <- read.csv("/correlações/BLACK_cargaviral_R8.csv")
BLACK_patas_R8 <- read.csv("/correlações/BLACK_patas_R8.csv")

# --- Colunas de interesse ---
lfc_cols <- c("log2FoldChange_D3", "log2FoldChange_D7", 
              "log2FoldChange_D21", "log2FoldChange_D35", "log2FoldChange_D70")

# --- Genes detectados na correlação com carga viral ---
genes_viral <- BLACK_cargaviral_R8 %>%
  filter(fake_pval_D3 <= 0.05 |
           fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 |
           fake_pval_D35 <= 0.05 |
           fake_pval_D70 <= 0.05) %>%
  select(GeneSymbol, all_of(lfc_cols))

# --- Genes detectados na correlação com patas ---
genes_patas <- BLACK_patas_R8 %>%
  filter(fake_pval_D3 <= 0.05 |
           fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 |
           fake_pval_D35 <= 0.05 |
           fake_pval_D70 <= 0.05) %>%
  select(GeneSymbol, all_of(lfc_cols))

# --- Unir os dois conjuntos e remover duplicados ---
genes_union <- bind_rows(genes_viral, genes_patas) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  arrange(GeneSymbol)   # ordenar alfabeticamente


df_long <- genes_union %>%
  pivot_longer(cols = starts_with("log2FoldChange"),
               names_to = "Condition",
               values_to = "log2FC")

df_long$Condition <- gsub("log2FoldChange_", "", df_long$Condition)
df_long$log2FC <- as.numeric(df_long$log2FC)

# --- Garantir que log2FC é numérico ---
df_long$log2FC <- as.numeric(df_long$log2FC)

# --- Converter Condition em fator com ordem explícita ---
df_long$Condition <- factor(df_long$Condition,
                            levels = c("D3", "D7", "D21", "D35", "D70"))

# --- Plot heatmap ---
ggplot(df_long, aes(x = Condition, y = GeneSymbol, fill = log2FC)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradientn(
    colours = c("#225ea8", "#3690c0", "white", "#31a354", "#216F63FF"),
    values = scales::rescale(c(min(df_long$log2FC, na.rm = TRUE),
                               -12, 0, 16,
                               max(df_long$log2FC, na.rm = TRUE))),
    limits = range(df_long$log2FC, na.rm = TRUE),
    na.value = "grey10",
    name = "log2FC"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 9, angle = 0, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),  # mantém nomes dos genes
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Heatmap BLACK")



# --- Criar listas de genes detectados em cada correlação ---
genes_viral <- BLACK_cargaviral_R8 %>%
  filter(fake_pval_D3 <= 0.05 | fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 | fake_pval_D35 <= 0.05 | fake_pval_D70 <= 0.05) %>%
  pull(GeneSymbol)

genes_patas <- BLACK_patas_R8 %>%
  filter(fake_pval_D3 <= 0.05 | fake_pval_D7 <= 0.05 |
           fake_pval_D21 <= 0.05 | fake_pval_D35 <= 0.05 | fake_pval_D70 <= 0.05) %>%
  pull(GeneSymbol)

# --- União dos genes ---
genes_union <- sort(unique(c(genes_viral, genes_patas)))

# --- Criar tabela para o heatmap de origem ---
df_origem <- expand.grid(GeneSymbol = genes_union,
                         Correlacao = c("Viral Load", "Paw"),
                         stringsAsFactors = FALSE)

df_origem <- df_origem %>%
  mutate(presente = case_when(
    Correlacao == "Viral Load" & GeneSymbol %in% genes_viral ~ 1,
    Correlacao == "Paw" & GeneSymbol %in% genes_patas ~ 1,
    TRUE ~ 0
  ))


# --- Plot Heatmap Black ---
ggplot(df_origem, aes(x = Correlacao, y = GeneSymbol, fill = factor(presente))) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_manual(values = c("0" = "grey90", "1" = "#31a354"),
                    name = "Correlação",
                    labels = c("Não det", "Det")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 6),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right") +
  labs(title = "Heatmap BLACK")

