# curl -s  "http://rest.kegg.jp/list/organism" | grep "Homo sapiens"

# T01001  hsa     Homo sapiens (human) Eukaryotes;Animals;Vertebrates;Mammals

# $ curl "http://rest.kegg.jp/list/pathway/T00041" | cut -f 1 | while read A; do  curl -o "${A}.xml" "http://rest.kegg.jp/get/${A}/kgml" ; done

# put files in /home/alex/R/library/SPIA/extdata/keggxml/hsa

# open R

# SPIA::makeSPIAdata('/home/alex/R/library/SPIA/extdata/keggxml/hsa', out.path = NULL)




# KEGG database
process_kegg <- function() {

    mapkKGML <- list.files('/home/alex/R/library/SPIA/extdata/keggxml/hsa/', full.names = TRUE)
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
    return(list(gslist = gslist, gs.names=gs.names))
}

