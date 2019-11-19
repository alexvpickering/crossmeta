library(data.table)
library(DBI)

# download from ftp://ftp.ncbi.nih.gov/gene/DATA/
dt <- fread('data-raw/entrezdt/gene2accession', select = c('#tax_id', 'GeneID', 'Symbol'))
colnames(dt) <- c('TAXID', 'ENTREZID', 'SYMBOL')
dt[] <- lapply(dt, as.character)
dt <- unique(dt)
setkey(dt, ENTREZID)

# homo sapiens
hs <- dt[TAXID == '9606']
hs[,TAXID := NULL]
colnames(hs) <- c('ENTREZID', 'SYMBOL_9606')

# save human for mapping
saveRDS(hs, 'data-raw/entrezdt/hs.rds')

# remove taxid
dt <- dt[, TAXID := NULL]
saveRDS(dt, 'data-raw/entrezdt/endt.rds')

# make sqlite
if (!file.exists('inst/extdata')) dir.create('inst/extdata', recursive = TRUE)

s <- sprintf("create table %s(%s, primary key(%s))", "ensql",
             paste(names(dt), collapse = ", "),
             names(dt)[1])

db <- dbConnect(RSQLite::SQLite(), "inst/extdata/ensql.sqlite")
dbGetQuery(db, s)
dbWriteTable(db, "ensql", dt, append = TRUE, row.names = FALSE)






