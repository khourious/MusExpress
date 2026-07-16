# MusExpress

## Descrição
Este repositório tem como objetivo documentar e armazenar o método de processamento feito para análises downstream no R de dados de RNA-seq em tecido muscular de camundongos.  
Não contém dados brutos, apenas scripts e fluxos de análise utilizados após o pré-processamento.

---

## Fluxo de trabalho
1. **Entrada de dados**  
   - Matrizes de expressão gênica já processadas pelo pipeline [SARA](https://github.com/khourious/sara)
   - Dados de RNA-seq previamente passaram pelo QC, foram alinhados e quantificados

2. **Processamento no R**
   - Normalização das matrizes de expressão (via DESeq2)
   - Identificação de genes diferencialmente expressos 
   - Correlações com variáveis experimentais  
   - Visualizações gráficas (heatmaps, PCA, volcano plots) 
   - Análise de co-expressão (via CEMiTool)
        - enriquecimento com Gene Ontoloty
        - Interactoma feito com base dados do StringDB, filtrado para _mus musculus_
   - Estimativa de sequências virais:  
     - As leituras foram alinhadas contra o genoma viral para estimar a quantidade relativa de vírus sequenciado 
     - Esses valores foram contador com o featureCounts e normalizados com edgeR
       
3. **Deconvolução celular**  
   - Ferramenta: [xCIBERSORT]([ca://s?q=ImmuCC_deconvolucao](https://cibersortx.stanford.edu/)).  
   - Assinatura específica para tecido muscular:  
     [`Muscle.sig.matrix.csv`](https://github.com/wuaipinglab/ImmuCC/blob/master/tissue_immucc/SignatureMatrix/Muscle.sig.matrix.csv).  
   - Objetivo: estimar proporções de células imunes presentes no tecido muscular.

4. **Enriquecimento de vias**:
   - Processado pelo IPA
   - Extraídas informações de upstream, canonical pathways e disease and biological functions
