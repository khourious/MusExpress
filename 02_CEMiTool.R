# --- Instalando pacotes ---
BiocManager::install("CEMiTool")
library("CEMiTool")
# --- carregando pacotes ---
library(CEMiTool)
library(ggplot2)
library(DESeq2)
library(readr)
library(dplyr)
library(gplots)
library(stringr)

# --- Configurando local de trabalho ---
getwd()
setwd("C:/Users/.../")

# --- carregando arquivos --- 
sample_phenotypes_BALBC <- read.delim("sample_phenotypes_BALBC.txt")
gmt_file <- "m5.go.mf.v2026.1.Mm.symbols.gmt"
gmt_sets <- read_gmt(gmt_file)

# --- Usar dados do 01_DeSeq2 ----
rld <- rlog(dds, blind = FALSE)

# --- Transformar para o formato apropriado
dados_rlog <- assay(rld) 
dados_rlog <- as.data.frame(dados_rlog)

# ---  criar o objeto sem gráficos automáticos
cem <- cemitool(expr = dados_rlog,
                annot = sample_phenotypes_BALBC,
                filter = TRUE,
                filter_pval = 0.1,
                force_beta = TRUE,
                min_ngen = 30,
                plot = FALSE)

# --- gerar gráficos manualmente
plot_mean_var(cem)
table(sapply(module_genes(cem), length))  # tamanho dos módulos
cem <- plot_profile(cem)
nmodules(cem)
plots <- show_plot(cem, "profile")
cem <- plot_gsea(cem)
plots <- show_plot(cem, "gsea")

# --- Visualizar plots
plots[1]
plots[2]
plots[3]

# --- Realizar analise ORA ---

# --- Converter o nome das vias para maior legebilidade- --
gmt_sets <- gmt_sets %>%
  mutate(
    # tira os underlines
    term = str_replace_all(term, "_", " "),
    # troca DN/UP por Down/Up
    term = str_replace(term, "DN$", "Down"),
    term = str_replace(term, "UP$", "Up"),
    # coloca o primeiro "palavra" (autor) entre parênteses no fim
    term = str_replace(term, "^([A-Za-z0-9]+) (.*)", "\\2 (\\1)")
  )

# --- Análise ORA ---
cem <- mod_ora(cem, gmt_sets, verbose = TRUE)
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
plots[1]
plots[2]
plots[3]

# --- Análise de interações ---
# --- Carregar biblioteca STRINGdb ---
library(STRINGdb)

# --- Aumentar o tempo limite de download para 10 minutos ---
options(timeout = 6000)

# ---  Especificar a espécie (Mus musculus) ---
species_id <- 10090

# --- Conectar ao banco de dados STRING (versão 11.5 é de 2026) ---
string_db <- STRINGdb$new(
  version = "11.5",
  species = species_id,
  score_threshold = 400,  # Filtra interações com escore combinado >= 400
  input_directory = getwd() # Pasta para download dos dados
)

# --- Converter as interações ---
interacoes_df <- data.frame(
  from = dict$gene[match(interacoes_raw$from, dict$STRING_id)],
  to   = dict$gene[match(interacoes_raw$to,   dict$STRING_id)]
)

# --- Remover quaisquer NAs ---
interacoes_df <- na.omit(interacoes_df)

# --- Remover duplicatas ---
interacoes_df <- unique(interacoes_df)

# --- Tabela de interações ---
head(interacoes_df)
message(paste("Número de interações únicas:", nrow(interacoes_df)))

# --- Organizando para o CEMiTool
colnames(interacoes_df) <- c("Gene1", "Gene2")
head(interacoes_df)


# --- Colocar genes em caixa baixa
interacoes_df$Gene1 <- tools::toTitleCase(tolower(interacoes_df$Gene1))
interacoes_df$Gene2 <- tools::toTitleCase(tolower(interacoes_df$Gene2))


# --- Incluir tabela de interações no objeto CEMiTool ---
interactions_data(cem) <- interacoes_df

# --- Plotar redes ---
cem <- plot_interactions(cem)

# --- Visualizar
plots <- show_plot(cem, "interaction")
plots[1]
plots[2]
plots[3]


# --- Salvar resultados ---
generate_report(cem, directory = "cemitool_balbc/MF", force = TRUE)

# --- Salvar plots ---
save_plots(
  cem,
  value = c("all", "profile", "gsea", "ora", "interaction", "beta_r2", "mean_k",
            "sample_tree", "mean_var", "hist", "qq"),
  force = FALSE,
  directory = "cemitool_balbc_test/MF/Plots"
)

# ---  Executar e salvar diagnóstico ---
diagnostic_report(
  cem,
  title = "Diagnostics",
  directory = "./cemitool_balbc_test/MF/Diagnostics",
  force = FALSE
)

write_files(cem, directory = "./cemitool_balbc_test/MF/Tables", force = FALSE)

