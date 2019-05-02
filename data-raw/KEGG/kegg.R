# curl -s  "http://rest.kegg.jp/list/organism" | grep "Homo sapiens"

# T01001  hsa     Homo sapiens (human) Eukaryotes;Animals;Vertebrates;Mammals

# cd ~/Documents/Batcave/GEO/crossmeta/data-raw/KEGG/hsa
# curl "http://rest.kegg.jp/list/pathway/T01001" | cut -f 1 | while read A; do  curl -o "${A}.xml" "http://rest.kegg.jp/get/${A}/kgml" ; done

# log date downloaded
# touch "$('%Y-%m-%d')"

library(data.table)

# KEGG database
mapkKGML <- list.files('~/Documents/Batcave/GEO/crossmeta/data-raw/KEGG/hsa', pattern = '^path:hsa', full.names = TRUE)
gslist <- lapply(mapkKGML, KEGGgraph::parseKGML)

names(gslist) <- lapply(gslist, function(path) path@pathwayInfo@number)
gs.names <- sapply(gslist, function(path) path@pathwayInfo@title)

gslist <- lapply(gslist, function(path) {
    # nodes and node types
    nodes <- path@nodes
    types <- sapply(nodes, function(node) node@type)

    # get genes in pathway
    unlist(lapply(nodes[types == 'gene'], function(node) gsub('^hsa:', '', node@name)), use.names = FALSE)
})

# get symbols
load("~/Documents/Batcave/GEO/crossmeta/R/sysdata.rda")
enids <- unique(unlist(gslist))
enids <- enids[enids %in% hs$ENTREZID]
syms <- toupper(hs[enids, SYMBOL_9606])
names(syms) <- enids

# add symbols names
gslist <- lapply(gslist, function(path) {
    names(path) <- syms[path]
    return(path)
})

# pathways in both only
inboth <- intersect(names(gslist), names(gs.names))
gslist <- gslist[inboth]
gs.names <- gs.names[inboth]

# save
usethis::use_data(gslist, gs.names, overwrite = TRUE)
