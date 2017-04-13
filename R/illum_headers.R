


# used by load_raw to fix illumina headers
fix_illum_headers <- function(gse_names, data_dir) {

    # functions
    can_load  <- function(path) {
        tryCatch({limma::read.ilmn(path, probeid = "ID_REF", verbose = FALSE); TRUE},
                 error = function(e) FALSE)
    }
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


    fixed <- NULL
    cat('Trying to fix headers for Illumina data:', paste(gse_names, collapse = ', '), '\n')

    for (gse_name in gse_names) {
        gse_dir <- file.path(data_dir, gse_name)
        paths   <- list.files(gse_dir, pattern = "non.norm.*txt$|raw.*txt$|nonorm.*txt$",
                              full.names = TRUE, ignore.case = TRUE)


        for (path in paths) {

            # get header line
            col1 <- 'ID_REF\\s*\t'
            col2 <- 'TARGET_ID\\s*\t|TargetID\\s*\t'
            col3 <- 'PROBE_ID\\s*\t|ProbeID\\s*\t'

            hline <- tryCatch(get_hline(path, col1), error = function(e) NULL)

            if (is.null(hline))
                hline <- tryCatch(get_hline(path, col2), error = function(e) NULL)

            if (is.null(hline))
                hline <- tryCatch(get_hline(path, col3), error = function(e) NULL)

            if (is.null(hline))
                next()

            # read raw file
            rawf <- readLines(path)

            # fixed path
            fpath <- gsub(".txt", "_fixed.txt", path, fixed = TRUE)


            # Target ID issue?
            pat0 <- "\\bTARGET_ID\\b|\\bTargetID\\b"

            if (grepl(pat0, rawf[hline], TRUE)) {
                # cat('pat0\n')
                rawf[hline] <- gsub(pat0, "ID_REF", rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }

            } else {

                # Probe ID issue?
                pat1 <- "\\bPROBE_ID\\b|\\bProbeID\\b"

                if (grepl(pat1, rawf[hline], TRUE)) {
                    # cat('pat1\n')
                    rawf[hline] <- gsub(pat1, "ID_REF", rawf[hline], TRUE)
                    writeLines(rawf, fpath)

                    # check if fixed
                    if (can_load(fpath)) {
                        fixed <- c(fixed, gse_name)
                        next()
                    }
                }
            }

            # not REF_ID issue if here
            # not header issue if can load original
            if (can_load(path)) next()

            #----

            # Pattern 1.1 issue?
            pat1_1 <- "\\(*.p.Value\\)*"

            if (grepl(pat1_1, rawf[hline], TRUE)) {
                # cat('pat1_1\n')
                rawf[hline] <- gsub(pat1_1, "-Detection", rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            # Pattern 1.2 issue?
            pat1_2 <- "\\(*.AVERAGE.Signal\\)*"

            if (grepl(pat1_2, rawf[hline], TRUE)) {
                # cat('pat1_2\n')
                rawf[hline] <- gsub(pat1_2, "-AVG_Signal", rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            # Pattern 1.3 issue?
            pat1_3 <- "([^\t]+).Signal\t\\1.Detection[^\t]*"
            rep1_3 <- "AVG_Signal-\\1\tDetection-\\1"

            if (grepl(pat1_3, rawf[hline], TRUE)) {
                # cat('pat1_3\n')
                rawf[hline] <- gsub(pat1_3, rep1_3, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            #----

            # Pattern 2 issue?
            pat2 <- "([^\t]+)\tDetection[^\t]*"
            rep2 <- "AVG_Signal-\\1\tDetection-\\1"


            if (grepl(pat2, rawf[hline], TRUE)) {
                # cat('pat2\n')
                rawf[hline] <- gsub(pat2, rep2, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            #----

            # Pattern 3 issue?
            pat3 <- "([^\t]+)\t\\1.Avg_NBEADS"
            rep3 <- "AVG_Signal-\\1\t\\1.Avg_NBEADS"

            if (grepl(pat3, rawf[hline], TRUE)) {
                # cat('pat3\n')
                rawf[hline] <- gsub(pat3, rep3, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            #----

            # Pattern 4 issue?
            pat4 <- "([^\t]+) \\(Average signal\\)\t\\1 \\(p-Value\\)"
            rep4 <- "AVG_Signal-\\1\tDetection-\\1"

            if (grepl(pat4, rawf[hline], TRUE)) {
                # cat('pat4\n')
                rawf[hline] <- gsub(pat4, rep4, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            #----
            # Pattern 5 issue?
            pat5 <- "([^\t]+)\t\\1[^\t]*Detection[^\t]*"
            rep5 <- "AVG_Signal-\\1\tDetection-\\1"

            if (grepl(pat5, rawf[hline], TRUE)) {
                # cat('pat5\n')
                rawf[hline] <- gsub(pat5, rep5, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

        }
    }
    return(unique(fixed))
}
