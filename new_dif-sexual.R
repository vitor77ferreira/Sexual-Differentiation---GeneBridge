#GENEPLAST - DIFERENCIAÇÃO SEXUAL

# ============ INSTALAÇÃO DO ANNOTATIONHUB (APENAS 1ª VEZ) ===========

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("AnnotationHub")

# =============== MEDIDAS CORRETIVAS =============== 

# Se aparecer a mensagem de que o AnnotationHub não está insrtalado, dê o seguinte comando no terminal do Linux:

# sudo apt install libcurl4-openssl-dev

# Depois um install.packages("RCurl") no RStudio.

# ============ INSTALAÇÃO DO GENEPLAST/GENEBRIDGE (APENAS 1ª VEZ) ===========

# ANTIGO. SÓ EXECUTAR SE USAR OS COMANDOS DO GENEPLAST

#BiocManager::install("geneplast.data")
# ==============================================

# # ============ Orientações do README # ============ 

#install.packages("knitr")
#install.packages("rmarkdown")
#BiocManager::install("BiocStyle")

#install.packages("remotes")
#remotes::install_github("sysbiolab/GeneBridge", build_vignettes=TRUE)

#install.packages("RCurl")

# ==========================================================
# BIBLIOTECAS ESSENCIAIS

library("geneplast")
# Algumas funções do GeneBridge ainda não funcionam como no geneplast

library(GeneBridge)
library(RCurl)
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

# POODE CUSTOMIZAR PARA O DIRETÓRIO DO SEU COMPUTADOR
genes_interesse <- read.csv("./genes/genes.csv")

# 451: TABELADOS

# Verifica e remove termos duplicados
genes_interesse <- genes_interesse %>%
  distinct()

# 331: PÓS RETIRADA DE DUPLICAÇÕES

# Salvar a tabela com os dados após a remoção dos termos duplicados
#write.csv(genes_interesse, file = "./genes/genes_sem_duplicatas.csv", row.names = FALSE)

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

# ======= SÓ FUNCIONA NO LINUX =========
string_id <- get_string_ids(gene) %>%
  janitor::clean_names() %>%
  tidyr::separate(string_id,
                  into = c("ssp_id", "string_id"),
                  sep = "\\.")
# =======================================

# Construindo ids_collapsed
ids_collapsed <- paste0(gene, collapse = "%0d")

## Get geneplast databases

ah <- AnnotationHub()
meta <- query(ah, "geneplast")

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

#ogr <-
#  newBridge(
#    ogdata = cogdata,
#    phyloTree = phyloTree,
#    ogids = cogs_of_interest$cog_id,
#    refsp = "9606"
#  )

ogr <-
  groot.preprocess(
    cogdata = cogdata,
    phyloTree = phyloTree,
    cogids = cogs_of_interest$cog,
    spid = "9606"
  )

# ogr <- runBridge(ogr, threshold = 0.3)

# ogr <- runPermutation(ogr, nPermutations = 1000, verbose = FALSE)

ogr <- groot(ogr, nPermutations = 1000, verbose = FALSE)


#res <- getBridge(ogr, what = "results")

res <- groot.get(ogr, what = "results")

## Nomeando os clados enraizados e obtendo a tabela de resultados finais

CLADE_NAMES <-
  "https://raw.githubusercontent.com/dalmolingroup/neurotransmissionevolution/ctenophora_before_porifera/analysis/geneplast_clade_names.tsv"

lca_names <- vroom(CLADE_NAMES)

groot_df <- res %>%
  tibble::rownames_to_column("cog_id") %>%
  select(cog_id, root = Root) %>%
  left_join(lca_names) %>%
  left_join(cogs_of_interest)

# SALVAR O ARQUIVO PRA CONSULTA POSTERIOR
#write.csv(groot_df, file = "~/Documentos/TCC/groot_df.csv", row.names = FALSE) #linux
write.csv(groot_df, file = "~/UFRN - BIOMEDICINA/TCC/groot_df.csv", row.names = FALSE) #Windows

# Extraição de colunas  da matriz groot_df
tabela_genes <- groot_df %>%
  select(clade_name, preferred_name, annotation)
View(groot_df)

# LEMBRETE: salvar a tabela e incluir nos resultados do TCC
write.csv(tabela_genes, file = "~/UFRN - BIOMEDICINA/TCC/tabela_de_genes.csv", row.names = FALSE)

# -================= EXEMPLO: ÁRVORE PARA A SRARP # -================= 
#(Steroid receptor-associated and regulated protein; 
#May regulate the transcriptional function of androgen and estrogen receptors.)

#groot.plot(ogr, whichOG = "NOG89338")

# Processando ogr com groot
ogr <- groot(ogr)

# Checando se o status de Rooting é correto agora
print(ogr@status["Rooting"])

# Se o status for adequado, prossiga com o plot
if (ogr@status["Rooting"] == "[x]") {
  # Envolva a função de plotagem com print() para garantir que seja exibida
  print(
    groot.plot(ogr, whichOG = "NOG89338", fname="gproot", width=4.5, height=6.5, cex.lab=0.3, 
               cex.nodes=0.6, adj.tips=c(1, 0.5), lab.offset=1.5, col.tips=c("green2","grey"), 
               col.edges=c("black","grey"), col.root="red", plot.sspnames=TRUE, 
               plot.subtree=FALSE, plot.lcas=FALSE)
  )
} else {
  stop("ogr não está adequadamente enraizado.")
}

if (is.na(ogr@status["Rooting"])) {
  ogr <- groot(ogr)  # Processa novamente se necessário
  # Após reprocessar, você pode querer visualizar o gráfico novamente
  if (ogr@status["Rooting"] == "[x]") {
    print(
      groot.plot(ogr, whichOG = "NOG89338", fname="gproot", width=4.5, height=6.5, cex.lab=0.3, 
                 cex.nodes=0.6, adj.tips=c(1, 0.5), lab.offset=1.5, col.tips=c("green2","grey"), 
                 col.edges=c("black","grey"), col.root="red", plot.sspnames=TRUE, 
                 plot.subtree=FALSE, plot.lcas=FALSE)
    )
  } else {
    stop("ogr ainda não está adequadamente enraizado após reprocessamento.")
  }
}

# -=================

##### Network assembly ----------

# Só funciona no LINUX
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

# =============== REALIZAR COM O DADO COMPLETO \/
# Filtra a rede por canais de evidência e confiabilidade das ligações

network_filtered <- network %>%
   combinescores(.,
                 evidences = c("ascore", "escore", "dscore"),
                 confLevel = 0.4) %>%
   separate(stringId_A,
            into = c("ncbi_taxon_id", "stringId_A"),
            sep = "\\.") %>%
   separate(stringId_B,
            into = c("ncbi_taxon_id", "stringId_B"),
            sep = "\\.") %>%
   dplyr::select(stringId_A, stringId_B)

# ==============================================================

graph <-
  graph_from_data_frame(network_filtered, directed = FALSE, vertices = nodelist)


roots <- unique(groot_df$root) %>%
  set_names(unique(groot_df$clade_name))

subsets <-
  map(roots, ~ subset_graph_by_root(groot_df, .x, graph))

# Plotando a rede mais nova e salvando como arquivo
g1 <- ggraph(subsets[[1]], "kk") +
  geom_edge_link(color = "#90909020") +
  geom_node_point() +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom")

#ggsave("~/Documentos/TCC/plots/rede_nova.png", plot = g1, width = 10, height = 8) #linux
ggsave("~/Documentos/TCC/plots/rede_nova.png", plot = g1, width = 10, height = 8) #windows

# Plotando a rede mais antiga e salvando como arquivo
g2 <- ggraph(subsets[[length(subsets)]], "kk") +
  geom_edge_link(color = "#90909020") +
  geom_node_point() +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom")

ggsave("~/Documentos/TCC/plots/rede_antiga.png", plot = g2, width = 10, height = 8)

# ======= GRÁFICOS =========
plot(x = groot_df$clade_name, y = groot_df$preferred_name) # deu errado
