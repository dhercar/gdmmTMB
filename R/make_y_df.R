#' Format community data
#' @description
#' Formats species x site table into a long data frame of community dissimilarities
#'
#' @param com community matrix (sites x species)
#' @param D Optional. Numeric matrix with two columns
#' @param id Optional. Unique id for each sample
#' @param method Method used to calculate community dissimilarity between pairs of sites. Either a distance metric supported by vegan, 'abcd' to obtain each component of the binary contingecy matrix, 'decomp1' or 'decomp2' for Pierre Legendre's or Andres Baselga's turnover and nestedness decomposition
#' @param num_den For some dissimilarity indices, whether function should return the numerator and denominator in two different columns
#' @param trans Transformation to be applied to raw data before calculating distance matrix. Either a function or 'binary' to transform to presence-absence)
#' @param drop_empty_rows Indicates whether empty rows (i.e., sites containing 0s only) should be removed
#' @param drop_empty_cols Indicates whether empty columns (i.e., species without observations) should be removed
#' @param binary  Indicates whether abundance data should be transform into presence absence
#' @param check_trans Check for non-finite values in transformation
#' @param na.rm Remove rows with NA?
#' @param na.replace Replace NAs with specific value
#'
#' @return A data frame with the following columns
#' \itemize{
#'   \item s1, s2  -  A combination of sites (rows).
#'   \item Additional columns indicating the dissimilarity between sites.
#' }
#'
#' @importFrom vegan vegdist designdist
#'
#'@keywords internal

make_y_df <-
  function(com = NULL,
           D = NULL,
           id = NULL,
           method = 'bray',
           num_den = FALSE,
           trans = NULL,
           drop_empty_rows = FALSE,
           drop_empty_cols = TRUE,
           binary = FALSE,
           check_trans = TRUE,
           na.rm = FALSE,
           na.replace = NULL){

    if(is.null(com)){
      stop("'com' must be supplied. You have not provided a community matrix")
    }
    # ---- Check data class ----
    if(!(is.matrix(com) | is.data.frame(com))) {
      stop(paste0("'com' must be a site x species matrix or a data frame. You have supplied a ", class(com)[1]))
    }

    #Convert to matrix
    try(com <- as.matrix(com, rownames = TRUE, colnames = TRUE))

    tryCatch({
      class(com) <- 'numeric'
    },
    # Error message
    error = function(cond) {
      stop("Community matrix 'com' could not be transformed into numeric format")
    })

    if (method == 'sorensen') {
      binary = TRUE
    }


    # ---- NA handling ----

    if (!is.null(na.replace)) {
      tryCatch({
        n_na <- sum(is.na(com))
        if (n_na > 0) {
          com[is.na(com)] <- na.replace
          warning(paste0("Values replaced! ", n_na, " NAs have been replaced with ", na.replace))
        }
      }, error = function(cond){
        stop(paste0("Could not substitute NAs with ", na.replace))
      })

    } else if (na.rm == TRUE) {

      # rows containing missing values
      NA_row <- apply(com, 1, function(x) any(is.na(x)))

      # remove samples from id and env
      id <- id[!NA_row]
      com <- com[!NA_row,]

      warning(paste0("Rows removed! ", sum(NA_row), " rows containing NAs have been removed"))
    }

    # ---- Transform matrix ----
    # Remove empty row (sites)
    if (drop_empty_rows) {
      row_drop <- rowSums(com) == 0 #Check for empty rows

      if (any(row_drop)) {
        com <- com[!row_drop,]

        if (is.null(id)) {
          site_char <- paste(which(row_drop), collapse = ",")
        } else {
          site_char <- paste(id[which(row_drop)], collapse = ",")
          id <- id[!row_drop]
        }
        warning(paste0("Sites removed! ", sum(row_drop), " empty rows (all 0s) have been removed. Rows: ", site_char))
      }
    }

    # Remove empty col (sp)
    if (drop_empty_cols) {
      col_drop <- colSums(com) == 0 #Check for empty cols

      if (any(col_drop)) {
        com <- com[,!col_drop]

        sp_char <- paste(names(col_drop)[(which(col_drop))], collapse = ",")
        warning(paste0("Species removed! ", sum(col_drop), " empty columns (all 0s) have been removed. Columns: ", sp_char))
      }
    }

    if (!is.null(trans)) {
      #Apply provided function
      if (is.function(trans)) {
        # Try function
        tryCatch({
          com <- trans(com)
        },
        # Error message
        error = function(cond) {
          stop("Unable to transform data with function provided")
        })

        if (all(!is.finite(com)) & check_trans) {
          stop("Non-infinite values! Transformation produced non-infinite values only")

        }else if (any(!is.finite(com)) & check_trans) {
          warning(paste0("Non-infinite values created! Transformation produced ", sum(!is.finite(com)), " non-infinite values"))
        }
      } else {
        stop("'trans' must be a function. 'trans' must be the name of a function such as log, log1p or sqrt. For instance: trans = sqrt, not trans = \"sqrt\" or trans = sqrt())")

      }
    }

    if (binary) {
      com = (com > 0)
    }

    # ---- Calculate distance matrix ----

    if (!is.null(id)) {
      if (anyDuplicated(id)) {
        stop("The 'id' variable must contain unique values for each site (no duplicates).")
      }
    }

    # handle D if provided
    if (is.null(D)) {
      # Generate all pairs by numeric index
      D <- t(combn(1:nrow(com), 2))
      colnames(D) <- c('s1', 's2')
    } else {
      if (!is.matrix(D) && !is.data.frame(D)) {
        stop("'D' must be a matrix or data frame with two columns")
      }
      if (ncol(D) != 2) {
        stop("'D' must have exactly two columns")
      }

      D <- as.matrix(D)

      if (is.character(D[,1]) || is.character(D[,2])) {
        # Require id to be provided for name matching
        if (is.null(id)) {
          stop("When 'D' contains site names, 'id' must be provided to map site names to indices")
        }

        # all D values are in id
        if (!all(D[,1] %in% id) || !all(D[,2] %in% id)) {
          stop("Some site names in 'D' are not found in 'id'")
        }

        idx1 <- match(D[,1], id)
        idx2 <- match(D[,2], id)

        D <- cbind(idx1, idx2)
        colnames(D) <- c('s1', 's2')
      }
    }


    # This could be parallelised in the future´
    if (num_den) {
      if (method == 'sorensen') {

        betas <- lapply(c('a', 'b', 'c'),
                        function(x) {
                          as.matrix(vegan::designdist(com,method = x,abcd = TRUE))
                        })

        a = betas[[1]][cbind(D[,1], D[,2])]
        b = betas[[2]][cbind(D[,1], D[,2])]
        c = betas[[3]][cbind(D[,1], D[,2])]

        dist.data <- data.frame(D,
                                num_sor = (b + c),
                                den_sor = (2*a + b + c))

      } else if (method == 'jaccard') {

        betas <- lapply(c('a', 'b', 'c'),
                        function(x) {
                          as.matrix(vegan::designdist(com,method = x,abcd = TRUE))
                        })

        a = betas[[1]][cbind(D[,1], D[,2])]
        b = betas[[2]][cbind(D[,1], D[,2])]
        c = betas[[3]][cbind(D[,1], D[,2])]

        dist.data <- data.frame(D,
                                num_jac = (b + c),
                                den_jac = (a + b + c))

      } else if (method == 'bray') {

        bray1 <- data.frame(t(apply(D,
                                    MARGIN = 1,
                                    function(x){
                                      bray1(com[x[1],], com[x[2],])
                                    })))

        names(bray1) <- c('num_bra','den_bra')
        dist.data <- data.frame(D,bray1)

      }else{
        stop(paste0("'num_den' not supported for method '", method, "'! 'num_den' only supported for methods \"bray\", \"sorensen\" and \"jaccard\". Method provided: ", method))
      }
    }else{
      if (method == 'abcd') {
        # List of a, b, c, d, components
        betas <- lapply(c('a', 'b', 'c', 'd'), function(x) {
          as.matrix(vegan::designdist(com, method = x, abcd = TRUE))
        })

        # Create data frame
        dist.data <- data.frame(
          D,
          a = betas[[1]][cbind(D[,1], D[,2])],
          b = betas[[2]][cbind(D[,1], D[,2])],
          c = betas[[3]][cbind(D[,1], D[,2])],
          d = betas[[4]][cbind(D[,1], D[,2])]
        )
      } else if (method %in% c('decomp1', 'decomp2')) {
        betas <- lapply(c('a', 'b', 'c'), function(x) {
          as.matrix(vegan::designdist(com, method = x, abcd = TRUE))
        })

        a = betas[[1]][cbind(D[,1], D[,2])]
        b = betas[[2]][cbind(D[,1], D[,2])]
        c = betas[[3]][cbind(D[,1], D[,2])]

        if (method == 'decomp1') {
          #Legendre's apprach
          dist.data <- data.frame(
            D,
            sim = 2 * pmin(b, c) / (2 * a + b + c),
            sne = abs(b - c) / (2 * a + b + c),
            sor = (b + c) / (2 * a + b + c)
          )

        } else if (method == 'decomp2') {
          #Baselga's approach
          dist.data <- data.frame(
            D,
            sim = pmin(b, c) / (a + pmin(b + c)),
            sor = (b + c) / (2 * a + b + c)
          )
          dist.data$sne <- dist.data$sor - dist.data$sim
        }
      } else {

        if (method == 'sorensen') {
          method = 'bray'
        }

        tryCatch({
          betas <- as.matrix(vegan::vegdist(com, method = method))
        }, error = function(cond) {
          stop("Unable to calculate community dissimilarity matrix! 'method' should be one of:\n- A method compatible with vegan `vegdist`\n- \"decomp1\" or \"decomp2\" for Sorensen-based decomposition of beta-diversity into species replacement and richness differences\n- \"abcd\" for individual elements of the binary contingency table")
        })

        dist.data <- data.frame(D, diss = betas[cbind(D[,1], D[,2])])
      }
    }
    return(dist.data)
  }
