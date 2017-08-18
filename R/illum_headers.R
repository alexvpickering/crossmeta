

fix_illum_headers <- function(data_paths, eset) {

    # functions
    get_hline <- function (file, columns, sep = "\t") {
        if (missing(columns) || !length(columns))
            stop("must specify column headings to find")
        if (!length(columns))
            stop("column headings must be specified")
        con <- file(file, "r")
        on.exit(close(con))
        out <- list()
        Found <- FALSE
        i <- 0
        repeat {
            i <- i + 1
            txt <- readLines(con, n = 1)
            if (!length(txt))
                stop("Specified column headings not found in file")
            Found <- TRUE
            for (a in columns) Found <- Found && length(grep(a, txt, TRUE))
            if (Found)
                break
        }
        return(i)
    }

    for (path in data_paths) {
        # fixed path
        fpath <- gsub(".txt", "_fixed.txt", path, fixed = TRUE)

        # read raw file
        rawf <- readLines(path)

        # make tab seperated if currently not
        delim <- reader::get.delim(path, n=50, skip=100)
        if (delim != '\t') rawf <- gsub(delim, '\t', rawf)

        # exclude lines starting with hashtag or tab
        exclude <- grepl('^.?#|^\t', rawf)
        rawf <- rawf[!exclude]

        # remove trailing tabs
        rawf <- gsub('\t*$', '', rawf)

        # save as will read from
        writeLines(rawf, fpath)

        # fread first 1000 rows as example
        hline <- get_hline(fpath, 'signal|pval|detection|id')
        ex <- data.table::fread(fpath, skip = hline-1, nrows = 1000, fill = TRUE)
        ex <- as.data.frame(ex)

        # fix annotation columns ----

        # look for column with ILMN entries
        nilmn  <- apply(ex, 2, function(col) sum(grepl('ILMN_', col)))
        idcol  <- which(nilmn > 950)

        # fix if idcol is not ID_REF
        if (length(idcol) && names(idcol) != 'ID_REF') {
            names(ex)[idcol] <- names(idcol)
            rawf[hline] <- paste0(names(ex), collapse = '\t')
        }

        # identify other annotation columns
        isnum   <- sapply(ex, class) == 'numeric'
        anncols <- setdiff(names(ex)[!isnum], 'ID_REF')

        if (!length(anncols)) anncols <- NULL

        # rename Signal and Pvalue identifiers ----
        pcols <- grep('pval|detection', colnames(ex), TRUE)
        scols <- grep('signal', colnames(ex), TRUE)

        # rename pvalue columns
        if (length(pcols)) {

            # longest common prefix or suffix in pvalue columns
            pcol_prefix <- lcPrefix(colnames(ex)[pcols])
            pcol_sufix  <- lcSuffix(colnames(ex)[pcols])

            if (grepl('pval|detection', pcol_prefix, TRUE)) {
                colnames(ex)[pcols] <- gsub(pcol_prefix, '', colnames(ex)[pcols])
                colnames(ex)[pcols] <- paste0('Detection-', colnames(ex)[pcols])

                rawf[hline] <- paste0(colnames(ex), collapse = '\t')

            } else if (grepl('pval|detection', pcol_sufix, TRUE)) {
                colnames(ex)[pcols] <- gsub(pcol_sufix, '', colnames(ex)[pcols])
                colnames(ex)[pcols] <- paste0('Detection-', colnames(ex)[pcols])

                rawf[hline] <- paste0(colnames(ex), collapse = '\t')
            }
        }

        # rename signal columns
        if (length(scols)) {

            # longest common prefix or suffix in pvalue columns
            scol_prefix <- lcPrefix(colnames(ex)[scols])
            scol_sufix  <- lcSuffix(colnames(ex)[scols])

            if (grepl('signal', scol_prefix, TRUE)) {
                colnames(ex)[scols] <- gsub(scol_prefix, '', colnames(ex)[scols])
                colnames(ex)[scols] <- paste0('AVG_Signal-', colnames(ex)[scols])

                rawf[hline] <- paste0(colnames(ex), collapse = '\t')

            } else if (grepl('signal', scol_sufix, TRUE)) {
                colnames(ex)[scols] <- gsub(scol_sufix, '', colnames(ex)[scols])
                colnames(ex)[scols] <- paste0('AVG_Signal-', colnames(ex)[scols])

                rawf[hline] <- paste0(colnames(ex), collapse = '\t')
            }
        }


        # if Pvalue every second column, set Signal to every first ----
        if (!length(pcols))
            pcols <- unname(which(sapply(ex, function(col) ifelse(is.numeric(col), max(col) < 1.01, FALSE))))

        if (length(pcols))
            p2nd <- all.equal(seq(min(pcols), max(pcols), 2), pcols)

        if (length(pcols) && p2nd) {
            # make Signal columns to left of Pvalue columns
            # (if couldn't detect)
            if (!length(scols))
                scols <- pcols-1

            # use Signal column for sample names
            sample_names <- gsub('AVG_Signal-', '', colnames(ex)[scols])

            colnames(ex)[scols] <- paste0('AVG_Signal-', sample_names)
            colnames(ex)[pcols] <- paste0('Detection-', sample_names)

            rawf[hline] <- paste(colnames(ex), collapse = '\t')
        }


        # if num numeric columns is n samples, set Signal to numeric columns ----
        nsamp <- ncol(eset)

        if (sum(isnum) == nsamp) {
            # add Signal identifier to all numeric columns
            colnames(ex)[isnum] <- gsub('AVG_Signal-', '', colnames(ex)[isnum])
            colnames(ex)[isnum] <- paste0('AVG_Signal-', colnames(ex)[isnum])
        }
    }
    writeLines(rawf, fpath)
    return(anncols)
}
