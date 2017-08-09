


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

        paths   <- list.files(gse_dir, "non.norm.*txt$|raw.*txt$|nonorm.*txt$",
                              full.names = TRUE, ignore.case = TRUE)



        for (path in paths) {

            # get header line
            col1 <- 'ID_REF\\s*\t'
            col2 <- 'TARGET_ID\\s*\t|TargetID\\s*\t'
            col3 <- 'PROBE_ID\\s*\t|ProbeID\\s*\t'

            hline <- tryCatch(get_hline(path, col1), error = function(e) NULL)

            if (!is.null(hline)) {
                idref <- FALSE

            } else {
                idref <- TRUE
                hline <- tryCatch(get_hline(path, col2), error = function(e) NULL)
            }

            if (is.null(hline))
                hline <- tryCatch(get_hline(path, col3), error = function(e) NULL)

            if (is.null(hline))
                next()

            # read raw file
            rawf <- readLines(path)

            # fixed path
            fpath <- gsub(".txt", "_fixed.txt", path, fixed = TRUE)


            # ID_REF issue? ----

            if (idref) {


                pat <- "\\bTARGET_ID\\b|\\bTargetID\\b"

                if (grepl(pat, rawf[hline], TRUE)) {
                    # cat('target\n')
                    rawf[hline] <- gsub(pat, "ID_REF", rawf[hline], TRUE)
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
            }

            # not ID_REF issue if here
            # not header issue if can load original
            if (can_load(path)) next()

            # Pattern 1 issue? ----

            pat <- ".\\(p.Value\\)"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat1\n')
                rawf[hline] <- gsub(pat, "-Detection", rawf[hline], TRUE)
                writeLines(rawf, fpath)

                pat <- ".\\(AVERAGE.Signal\\)"

                if (grepl(pat, rawf[hline], TRUE)) {
                    # cat('pat1\n')
                    rawf[hline] <- gsub(pat, "-AVG_Signal", rawf[hline], TRUE)
                    writeLines(rawf, fpath)

                    # check if fixed
                    if (can_load(fpath)) {
                        fixed <- c(fixed, gse_name)
                        next()
                    }
                }
            }

            # Pattern 2 issue? ----

            pat <- ".p.Value"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat2\n')
                rawf[hline] <- gsub(pat, "-Detection", rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            pat <- ".AVERAGE.Signal"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat2\n')
                rawf[hline] <- gsub(pat, "-AVG_Signal", rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            # Pattern 3 issue? ----

            pat <- "([^\t]+).Signal\t\\1.Detection[^\t]*"
            rep <- "AVG_Signal-\\1\tDetection-\\1"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat3\n')
                rawf[hline] <- gsub(pat, rep, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            # Pattern 4 issue? ----

            pat <- "([^\t]+)\tDetection[^\t]*"
            rep <- "AVG_Signal-\\1\tDetection-\\1"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat4\n')
                rawf[hline] <- gsub(pat, rep, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            # Pattern 5 issue? ----

            pat <- "([^\t]+)\t\\1.Avg_NBEADS"
            rep <- "AVG_Signal-\\1\t\\1.Avg_NBEADS"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat5\n')
                rawf[hline] <- gsub(pat, rep, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            # Pattern 6 issue? ----

            pat <- "([^\t]+) \\(Average signal\\)\t\\1 \\(p-Value\\)"
            rep <- "AVG_Signal-\\1\tDetection-\\1"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat6\n')
                rawf[hline] <- gsub(pat, rep, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            # Pattern 7 issue? ----

            pat <- "([^\t]+)\t\\1[^\t]*Detection[^\t]*"
            rep <- "AVG_Signal-\\1\tDetection-\\1"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat7\n')
                rawf[hline] <- gsub(pat, rep, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            # Pattern 8 ----

            pat <- "(GSM[0-9]+\t*)"
            rep <- "AVG_Signal-\\1"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat8\n')
                rawf[hline] <- gsub(pat, rep, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }

            # Pattern 9 ----

            pat <- "(SAMPLE[0-9]+\t*)"
            rep <- "AVG_Signal-\\1"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat9\n')
                rawf[hline] <- gsub(pat, rep, rawf[hline], TRUE)
                writeLines(rawf, fpath)

                # check if fixed
                if (can_load(fpath)) {
                    fixed <- c(fixed, gse_name)
                    next()
                }
            }


            # Pattern 10 ----

            pat <- "\t([^\t]+)"
            rep <- "\tAVG_Signal-\\1"

            if (grepl(pat, rawf[hline], TRUE)) {
                # cat('pat10\n')
                rawf[hline] <- gsub(pat, rep, rawf[hline], TRUE)
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
