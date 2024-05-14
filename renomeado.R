#GENEPLAST - DIFERENCIAÇÃO SEXUAL

# ============ INSTALAÇÃO DO ANNOTATIONHUB (APENAS 1ª VEZ) ===========

#install.packages("AnnotationHub")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("AnnotationHub")

# ============ INSTALAÇÃO DO GENEPLAST/GENEBRIDGE (APENAS 1ª VEZ) ===========

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("geneplast.data")

# Orientações do README

#install.packages("knitr")
#install.packages("rmarkdown")
#install.packages("BiocManager")
#BiocManager::install("BiocStyle")

#install.packages("remotes")
#remotes::install_github("sysbiolab/GeneBridge", build_vignettes=TRUE)

# ==========================================================
# BIBLIOTECAS ESSENCIAIS

#library(geneplast.data)
#install.packages(c("curl", "httr"))

library(GeneBridge)
vignette("GeneBridge")
library(tibble)
library(stringi)
library(AnnotationHub)
library(dplyr)

# ===== CUSTOMIZÁVEIS =====
library(vroom)
library(igraph)
library(ggraph)
library(purrr)
library(tidyr)
library(here)

# ===== Informações sobre ortólogos =====
data("ogdata")
# ================= CITAÇÕES =================
# citation("pkgname") para citar cada pacote utilizado ao final do TCC

# =============== CITAÇÃO DO GENEPLAST==============

# Dalmolin RJS, Castro MAA. Geneplast: evolutionary rooting and plasticity inference. R package, 2015.
# ============================================

set.seed(1024)

# Define uma função para receber uma lista de identificadores
# que será usada para fazer uma 
# busca no String-DB (banco de dados de interações de proteínas)
# Os parâmetros da consulta incluem os identificadores formatados e 
# a espécie (default "9606" para humano). 
# O resultado é convertido de JSON para um objeto R usando a função 'fromJSON'.

# Get identifiers from STRINGdb
get_string_ids <- function(ids, species = "9606") {
  ids_collapsed <- paste0(ids, collapse = "%0d")

  jsonlite::fromJSON(
    RCurl::postForm(
      "https://string-db.org/api/json/get_string_ids",
      identifiers = ids_collapsed,
      echo_query  = "1",
      species = species
    ),
  )
}

# ===================================================
# PROCESSAMENTO DE GENES ANALISADOS
# ===================================================

genes_interesse <- read.csv("./genes/genes.csv")

# 451: TABELADOS

# Verifica e remove termos duplicados
genes_interesse <- genes_interesse %>%
  distinct()

# 331: PÓS RETIRADA DE DUPLICAÇÕES

# Salvar a tabela com os dados após a remoção dos termos duplicados
write.csv(genes_interesse, file = "./genes/genes_sem_duplicatas.csv", row.names = FALSE)

# eliminação dos espaços
genes_interesse <- gsub(" ", "", genes_interesse$GENES)

# ===================================================

# Get interaction network from STRINGdb
get_string_network <-
  function(ids,
           species = "9606",
           required_score = 0) {
    ids_collapsed <- paste0(ids, collapse = "%0d")

    jsonlite::fromJSON(
      RCurl::postForm(
        "https://string-db.org/api/json/network",
        identifiers = ids_collapsed,
        echo_query  = "1",
        required_score = as.character(required_score),
        species = species
      ),
    )
  }


# Function to combine scores according to the STRINGdb algorithm
combinescores <- function(dat,
                          evidences = "all",
                          confLevel = 0.4) {
  if (evidences[1] == "all") {
    edat <- dat[,-c(1, 2, ncol(dat))]
  } else {
    if (!all(evidences %in% colnames(dat))) {
      stop("NOTE: one or more 'evidences' not listed in 'dat' colnames!")
    }
    edat <- dat[, evidences]
  }
  if (any(edat > 1)) {
    edat <- edat / 1000
  }
  edat <- 1 - edat
  sc <- apply(
    X = edat,
    MARGIN = 1,
    FUN = function(x)
      1 - prod(x)
  )
  dat <- cbind(dat[, c(1, 2)], combined_score = sc)
  idx <- dat$combined_score >= confLevel
  dat <- dat[idx,]
  return(dat)
}

subset_graph_by_root <-
  function(geneplast_result, root_number, graph) {
    filtered <- geneplast_result %>%
      filter(root >= root_number) %>%
      pull(protein_id)

    induced_subgraph(graph, which(V(graph)$name %in% filtered))
  }

# =================== ENRAIZAMENTO ===================

gene <- genes_interesse

# genes_de_interesse <- vroom("tabela_genes.csv")
# genes <- genes_interesse$GENES

# gene <- c("FGFR3", "ALDH1L1", "S100B")

#  ======== !!!! Só funciona no linux !!!! ========

string_id <- get_string_ids(gene) %>%
  janitor::clean_names() %>%
  tidyr::separate(string_id,
                  into = c("ssp_id", "string_id"),
                  sep = "\\.")
# =============================================

# Construindo ids_collapsed
ids_collapsed <- paste0(gene, collapse = "%0d")

# Fazendo a solicitação à API do STRINGdb
string_id <- RCurl::postForm("https://string-db.org/api/json/get_string_ids",
                             identifiers = ids_collapsed,
                             echo_query = "1",
                             species = species,
                             ssl.verifypeer = FALSE) %>%
  janitor::clean_names() %>%
  tidyr::separate(string_id,
                  into = c("ssp_id", "string_id"),
                  sep = "\\.")

# ============= A parte de cima de problema de validação de certificado

## Get geneplast databases

ah <- AnnotationHub()
meta <- query(ah, "geneplast")
head(meta)

# Carregue os objetos na sessão usando o ID do conjunto de dados escolhido do STRING database v11.0
load(meta[["AH83116"]])

# Takes a lot of RAM to do:
# Cog Mappings (cogdata.tsv) generated by bin/create_cog_map.sh
cogs_of_interest <- cogdata %>%
  filter(ssp_id %in% string_id$ssp_id) %>%
  filter(protein_id %in% string_id$string_id) %>%
  select(-ssp_id) %>%
  left_join(string_id, by = c("protein_id" = "string_id"))

gc()

# ================ ERRO NESSE PONTO  ================
#-Preprocessing input data...
#Erro: objeto 'cogs_of_interest' não encontrado
#> ogr <- groot(ogr, nPermutations = 1000, verbose = FALSE)
#Error in h(simpleError(msg, call)) : 
#  erro na avaliação do argumento 'object' na seleção do método para a função 'groot': 'objeto 'ogr' não encontrado'

# =======================================

ogr <-
  newBridge(
    ogdata = cogdata,
    phyloTree = phyloTree,
    ogids = cogs_of_interest$cog_id,
    refsp = "9606"
  )

ogr <- runBridge(ogr, threshold = 0.3)
                 
                 
                 nPermutations = 1000, verbose = FALSE)
ogr <- runPermutation(ogr, nPermutations = 1000)

res <- getBridge(ogr, what = "results")

## Nomeando os clados enraizados e obtendo a tabela de resultados finais

CLADE_NAMES <-
  "https://raw.githubusercontent.com/dalmolingroup/neurotransmissionevolution/ctenophora_before_porifera/analysis/geneplast_clade_names.tsv"

lca_names <- vroom(CLADE_NAMES)

groot_df <- res %>%
  tibble::rownames_to_column("cog_id") %>%
  select(cog_id, root = Root) %>%
  left_join(lca_names) %>%
  left_join(cogs_of_interest)

# -=================

# SALVAR O ARQUIVO PRA CONSULTA POSTERIOR
write.csv(groot_df, file = "~/UFRN - BIOMEDICINA/TCC/genes/groot_df.csv", row.names = FALSE)

View(groot_df)
# -=================

# EXEMPLO: ÁRVORE PARA A KISSPEPTINA

# groot.plot(ogr, whichOG = "NOG26751")
# -=================

##### Network assembly ----------
#  ======== !!!! Só funciona no linux !!!! ========
network <- get_string_network(string_id$string_id)

network_separated <-  network %>%
  separate(stringId_A,
           into = c("ncbi_taxon_id", "stringId_A"),
           sep = "\\.") %>%
  separate(stringId_B,
           into = c("ncbi_taxon_id", "stringId_B"),
           sep = "\\.")

nodelist <-
  data.frame(node = unique(
    c(network_separated$stringId_A, network_separated$stringId_B)
  )) %>%
  left_join(string_id, by = c("node" = "string_id"))

network_filtered <- network |>
    separate(stringId_A,
             into = c("ncbi_taxon_id", "stringId_A"),
             sep = "\\.") %>%
    separate(stringId_B,
             into = c("ncbi_taxon_id", "stringId_B"),
             sep = "\\.") %>%
    dplyr::select(stringId_A, stringId_B)

# REALIZAR COM O DADO COMPLETO \/
# Filtra a rede por canais de evidência e confiabilidade das ligações

# network_filtered <- network %>%
#   combinescores(.,
#                 evidences = c("ascore", "escore", "dscore"),
#                 confLevel = 0.4) %>%
#   separate(stringId_A,
#            into = c("ncbi_taxon_id", "stringId_A"),
#            sep = "\\.") %>%
#   separate(stringId_B,
#            into = c("ncbi_taxon_id", "stringId_B"),
#            sep = "\\.") %>%
#   dplyr::select(stringId_A, stringId_B)

graph <-
  graph_from_data_frame(network_filtered, directed = FALSE, vertices = nodelist)


roots <- unique(groot_df$root) %>%
  set_names(unique(groot_df$clade_name))

subsets <-
  map(roots, ~ subset_graph_by_root(groot_df, .x, graph))

# Plotando a rede mais nova

ggraph(subsets[[1]], "kk") +
  geom_edge_link(color = "#90909020") +
  geom_node_point() +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom")

# Plotando a rede mais antiga

ggraph(subsets[[length(subsets)]], "kk") +
  geom_edge_link(color = "#90909020") +
  geom_node_point() +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom")
