# homologene.data available from http://1.usa.gov/1TGoIy7

# Get homologene data frame.
#
# Function loads and, if necessary, sets up homologene data frame.
#
# Setup results in a dataframe with all entrez ids in one column ('ENTREZID') and
# the homologous human entrez ids in another column ('ENTREZID_HS').
#
# @param homologene_path String specifying full path to homologene file.
#
# @return data.frame with column ENTREZID and ENTREZID_HS.

get_homologene <- function(homologene_path) {

    homologene <- utils::read.delim(homologene_path, header = FALSE)

    if (ncol(homologene) == 6) {

        message("Setting up homologene data. Will take a while (one time only).")

        # get homologous human (9606) entrez ids for all entrez ids (V3)
        entrez_HS <- annotationTools::getHOMOLOG(homologene$V3,
                                                 9606,
                                                 homologene)

        # remove entrez ids without homologous human entrez id
        homologene <- homologene[!is.na(entrez_HS), ]
        entrez_HS  <- entrez_HS[!is.na(entrez_HS)]

        # expand homologene where multiple homologous human entrez ids
        homologene$entrez_HS <- sapply(entrez_HS,
                                       function(x) paste(x, collapse = ","))

        homologene <- data.table::data.table(homologene)
        homologene <- homologene[,
                                 list(entrez_HS = unlist(strsplit(entrez_HS, ","))),
                                 by = V3]

        homologene <- as.data.frame(homologene)

        # remove rows where entrez id is human (creates 1:many issue)
        #
        homologene <- homologene[!homologene$ENTREZID %in% homologene$ENTREZID_HS, ]

        # save to disc
        utils::write.table(homologene, file = homologene_path,
                           sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    colnames(homologene) <- c("ENTREZID", "ENTREZID_HS")

    # add single entry for human entrez id (used to lookup entrez fData column)
    #
    hs <- unique(homologene$ENTREZID_HS)
    hs <- data.frame(ENTREZID = hs, ENTREZID_HS = hs)
    homologene <- rbind(homologene, hs)

    # character class
    homologene[] <- lapply(homologene, as.character)

    return(data.table(homologene, key='ENTREZID'))
}

devtools::use_data(homologene, gpl_bioc, hs, sources, org_pkg, org_taxid, token, internal = TRUE, overwrite = TRUE)

