# Illumina loader utility for load_plat.
#
# Used by load_plat to load an eset.
#
# @param eset Expression set obtained by \code{getGEO}.
# @param gse_name String specifying GSE name.
# @param gse_dir String specifying path to GSE folder.
#
# @seealso \code{\link{load_plat}}.
# @return Annotated eset.

load_illum_plat <- function(eset, gse_name, gse_dir, ensql) {
  
  try(fData(eset)[fData(eset) == ""] <- NA)
  try(fData(eset)[] <- lapply(fData(eset), as.character))
  
  # convert .xls to txt
  xls_paths <- list.files(gse_dir, pattern = "_supplementary_.*.xlsx?$", full.names = TRUE, ignore.case = TRUE)
  xls_to_txt(xls_paths)
  
  # fix header issues
  illum_pat <- "^GSM[0-9]+.*txt$|non.*norm.*txt$|raw.*txt$|nonorm.*txt$|_supplementary_.*.txt$"
  elist_paths <- list.files(gse_dir, pattern = illum_pat, full.names = TRUE, ignore.case = TRUE)
  elist_paths <- elist_paths[!grepl('fixed[.]txt$', elist_paths)]
  annotation  <- fix_illum_headers(elist_paths, eset)
  
  # load fixed elist paths
  elist_paths <- gsub(".txt", "_fixed.txt", elist_paths, fixed = TRUE)
  elist <- limma::read.ilmn(elist_paths, probeid = "ID_REF", annotation = annotation)
  
  
  # don't correct if already log transformed (already corrected?)
  logd <- max(elist$E, na.rm = TRUE) < 100
  if (!logd) {
    elist <- tryCatch (
      limma::neqc(elist),
      error = function(c) {
        # PMID:19068485 recommends mle and offset 50
        elist <- limma::backgroundCorrect(elist, method = "normexp",
                                          normexp.method = "mle",
                                          offset = 50)
        
        return(limma::normalizeBetweenArrays(elist, method = "quantile"))
      })
  }
  
  # merge eset and elist fdata
  elist <- merge_elist(eset, elist)
  
  if ('ID_REF' %in% colnames(elist$genes)) {
    row.names(elist$E) <- row.names(elist$genes) <- make.unique(elist$genes$ID_REF)
  } else {
    row.names(elist$E) <- row.names(elist$genes) <- NULL
  }
  
  # determine best sample matches
  res   <- match_samples(eset, elist)
  elist <- elist[, res$elist_order]
  eset  <- eset[, res$eset_order]
  warn  <- res$warn
  
  # keep gse matrix and raw elist title
  pData(eset)$title.gsemat <- pData(eset)$title
  pData(eset)$title.raw    <- colnames(elist)
  
  
  if (warn) {
    # use raw elist titles to ensure correct contrasts
    pData(eset)$title <- colnames(elist)
    
    # add illum colname to warn about pData
    pData(eset)$illum <- NA
  }
  colnames(elist) <- sampleNames(eset)
  
  # convert limma object to eset
  eset <- to_eset(elist, eset)
  
  # add SYMBOL annotation
  eset <- symbol_annot(eset, gse_name, ensql)
  return(eset)
}


#' Covert .xls files to .txt
#' 
#' For converting Illumina _Supplementary_.*.xls files to .txt for load_illum_plat.
#'
#' @param xls_paths Paths to .xls files
#'
#' @return NULL
#' @export
#'
xls_to_txt <- function(xls_paths) {
  for (xls_path in xls_paths) {
    d <- readxl::read_excel(xls_path)
    txt_path <- gsub('.xlsx?$', '.txt', xls_path)
    write.table(d, txt_path, sep='\t', quote = FALSE, row.names = FALSE)
  }
}


# like base pmatch
# partial match occurs if the whole of the element of x matches any part of the element of table
fuzzy_pmatch <- function(x, table) {
  x <- tolower(x)
  table <- tolower(table)
  
  # first look for perfect matches
  perfect <- match(x, table)
  
  # is every x has a perfect match in table, return
  if (sum(is.na(perfect)) == 0) return(perfect)
  
  # otherwise first grep
  tomatch <- x[is.na(perfect)]
  gmatch  <- sapply(tomatch, function(val) {
    res <- grep(val, table, fixed = TRUE)[1]
    if (!length(res)) return(NA_integer_)
    return(res)
  })
  
  # fill in grep result where NA in perfect
  perfect[is.na(perfect)] <- gmatch[is.na(perfect)]
  return(perfect)
}

match_samples <- function(eset, elist) {
  
  # determine if elist has fewer samples
  data_fewer <- ncol(elist) < ncol(eset)
  
  # check if colnames match ----
  if (data_fewer) {
    # check if all elist colnames in eset colnames
    if (all(colnames(elist) %in% colnames(eset))) {
      cat('Illumina samples matched by column names.\n')
      return(list(elist_order = colnames(elist), eset_order = colnames(elist), warn = FALSE))
    }
    
  } else {
    # check if all eset colnames in elist colnames
    if (all(colnames(eset) %in% colnames(elist))) {
      cat('Illumina samples matched by column names.\n')
      return(list(elist_order = colnames(eset), eset_order = colnames(eset), warn = FALSE))
    }
  }
  
  # check if eset pdata col matches elist colnames ----
  
  if (!is.null(colnames(elist))) {
    
    # matrix of positions of matches for elist colnames among those for each pdata column
    matches <- sapply(Biobase::pData(eset), function(col) {
      fuzzy_pmatch(colnames(elist), col)
    })
    
    # number of unique non NA matches for each pdata column
    nunique <- apply(matches, 2, function(match) length(unique(match[!is.na(match)])))
    
    # number of unique non NA matches should be the min of number of eset or pdata samples
    nmin <- min(ncol(eset), ncol(elist))
    if (any(nunique == nmin)) {
      cat('Illumina samples matched by pdata column.\n')
      
      # matches where satisfied
      bestcol <- names(which(nunique == nmin))[1]
      matches <- matches[, bestcol]
      
      # elist_order is positions where matches are not NA
      elist_order <- which(!is.na(matches))
      
      # eset_order is non NA matches
      eset_order <- matches[!is.na(matches)]
      
      return(list(elist_order = elist_order, eset_order = eset_order, warn = FALSE))
    }
  }
  
  # check if similarity offers unique match ----
  
  
  # make sure eset is log2 transformed
  logd <- max(exprs(eset), na.rm = TRUE) < 1000
  if (!logd) {
    exprs(eset) <- log2(exprs(eset) + abs(min(exprs(eset), na.rm = TRUE)) + 16)
  }
  
  # row names are the best match columns for elist and eset
  elist <- elist[!is.na(elist$genes[[elist$elistcol]]), ]
  row.names(elist) <- make.unique(elist$genes[[elist$elistcol]])
  
  eset <- eset[!is.na(fData(eset)[[elist$esetcol]]), ]
  row.names(eset)  <- make.unique(fData(eset)[[elist$esetcol]])
  
  # only include rows without missing values
  eset   <- eset[stats::complete.cases(exprs(eset)), ]
  elist  <- elist[stats::complete.cases(elist$E), ]
  
  qres   <- list()
  ngenes <- min(nrow(eset), nrow(elist))
  
  if (data_fewer) {
    # determine most similar eset sample for each sample in elist
    for (i in 1:ncol(elist)) {
      qsamp <- elist$E[, i]
      qres[[colnames(elist)[i]]] <- query_ref(qsamp, exprs(eset), sorted = FALSE, ngenes = ngenes)
    }
    
  } else {
    # determine most similar elist sample for each sample in eset
    for (i in 1:ncol(eset)) {
      qsamp <- exprs(eset)[, i]
      qres[[colnames(eset)[i]]] <- query_ref(qsamp, elist$E, sorted = FALSE, ngenes = ngenes)
    }
  }
  
  # eset sample to most similar elist sample
  qres <- as.data.frame(qres)
  best <- sapply(qres, which.max)
  
  if (length(best) == length(unique(best))) {
    cat('Illumina samples matched by similarity.\n')
    
    if (data_fewer) {
      elist_order <- colnames(elist)
      eset_order <- best
    } else {
      elist_order <- best
      eset_order <- colnames(eset)
    }
    
    return(list(elist_order = elist_order, eset_order = eset_order, warn = FALSE))
    
  } else {
    # look for misses in non-first query results
    dups   <- unique(best[duplicated(best)])
    misses <- setdiff(1:nrow(qres), unique(best))
    
    n <- nrow(qres)
    for (dup in dups) {
      # query results for duplicate
      i <- 1
      qres_dup  <- qres[, best == dup]
      
      while (dup %in% dups & i < n) {
        ibest_dup <- sapply(qres_dup, function(col) which(col == sort(col, partial=n-i)[n-i]))
        
        # for each miss
        for (miss in misses) {
          # check if one ibest is miss
          imiss <- ibest_dup == miss
          
          if (sum(imiss) == 1){
            # if so, replace best with ibest
            ibest_repl <- ibest_dup[imiss]
            best[names(ibest_repl)] <- ibest_repl
            
            # also update duplicates and misses
            dups   <- best[duplicated(best)]
            misses <- setdiff(1:nrow(qres), unique(best))
            
            # if no more misses, break
            if (!length(misses)) {
              break()
            }
          }
        }
        i <- i + 1
      }
    }
    
    if (!length(dups)) {
      cat('Illumina samples matched by similarity using non-first ranks.\n')
      if (data_fewer) {
        elist_order <- colnames(elist)
        eset_order <- best
      } else {
        elist_order <- best
        eset_order <- colnames(eset)
      }
      return(list(elist_order = elist_order, eset_order = eset_order, warn = FALSE))
      
    } else {
      cat('Illumina samples not matched.\n')
      return(list(elist_order = colnames(elist), eset_order = colnames(eset), warn = TRUE))
    }
  }
}


#' Get correlation between query and reference signatures.
#' 
#' Determines the pearson correlation between the query and each reference signature.
#'
#' @param query Named numeric vector of differentual expression values for
#'   query genes. Usually 'meta' slot of \code{get_dprimes} result.
#' @param ref A matrix of differential expression
#'   to query against (rows are genes, columns are samples).
#' @param sorted Would you like the results sorted by decreasing similarity?
#'   Default is TRUE.
#' @param ngenes The number of top differentially-regulated (up and down) query genes to use. 
#'
#' @return Vector of pearson correlations between query and reference signatures.
#'
query_ref <- function(query, ref, sorted = TRUE, ngenes = 200) {
  
  # use only common genes
  query <- query[names(query) %in% row.names(ref)]
  
  # top up/down ngenes
  top_ngenes  <- utils::head(names(sort(abs(query), TRUE)), ngenes)
  query <- query[top_ngenes]
  ref   <- ref[names(query), ,drop = FALSE]
  
  # pearson correlation
  sim <- stats::cor(query, ref, method="pearson")
  sim <- structure(c(sim), names=colnames(sim))
  
  if (sorted) {
    return(sort(sim, decreasing = TRUE))
  } else {
    return(sim)
  }
  
}



merge_elist <- function(eset, elist) {
  
  if (is.null(elist$genes)) stop('Raw elist lacks feature names.')
  
  # get eset and elist fdata columns
  esetcols  <- fData(eset)
  elistcols <- elist$genes
  
  
  # find eset fData column that best matches elist features
  best  <- c(esetcol=NA, elistcol=NA)
  bestf <- 0
  
  for (i in seq_along(elistcols)) {
    
    elistcol <- elistcols[[i]]
    
    # get fraction of fdata column that has a match
    matches <- sapply(names(esetcols), function(esetcol) {
      sum(elistcol %in% esetcols[, esetcol]) / length(elistcol)
    })
    
    # update best
    if (max(matches) >= bestf) {
      bestf <- max(matches)
      best['elistcol'] <- names(elistcols)[i]
      best['esetcol']  <- names(matches[which.max(matches)])
    }
  }
  
  if (bestf > 0.3) {
    # merge eset and elist fdata columns
    esetcols  <- esetcols[!duplicated(esetcols[best['esetcol']]),, drop = FALSE]
    elistcols <- merge(elistcols, esetcols, all.x = TRUE, by.x = best['elistcol'], by.y = best['esetcol'], sort = FALSE)
    elistcols[elistcols == ""] <- NA
    
    elist$genes <- elistcols
    
    # add best info for illumina sample matching
    elist$elistcol <- best[['elistcol']]
    elist$esetcol  <- best[['esetcol']]
  }
  elist$genes[] <- lapply(elist$genes, as.character)
  return(elist)
}

#' Open raw Illumina microarray files.
#'
#' Helper function to open raw Illumina microarray files in order to check that
#' they are formatted correctly. For details on correct format, please see
#' 'Checking Raw Illumina Data' in vignette.
#'
#' @param gse_names Character vector of Illumina GSE names to open.
#' @param data_dir String specifying directory with GSE folders.
#'
#' @return Character vector of successfully formated Illumina GSE names.
#' @export
#'
#' @examples
#' library(lydata)
#'
#' # Illumina GSE names
#' illum_names <- c("GSE50841", "GSE34817", "GSE29689")
#'
#' # location of raw data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # open raw data files with default text editor
#' # open_raw_illum(illum_names)

open_raw_illum <- function (gse_names, data_dir = getwd()) {
  
  out_names <- gse_names
  for (i in seq_along(gse_names)) {
    # get data paths
    gse_dir <- paste(data_dir, gse_names[i], sep = "/")
    data_paths <- list.files(gse_dir, pattern = "non.norm.*txt",
                             full.names = TRUE, ignore.case = TRUE)
    data_paths <- c(data_paths, list.files(gse_dir, pattern = ".xls",
                                           full.names = TRUE))
    # open data file
    for (j in seq_along(data_paths)) system2("xdg-open", data_paths[j])
    
    # check success
    success <- tcltk::tk_select.list(choices = c("Yes", "No"),
                                     title = paste(gse_names[i],
                                                   "formated successfully?"))
    # remove unsuccessful
    if (success == "No") out_names <- setdiff(out_names, gse_names[i])
  }
  return(out_names)
}

