# Parsing of random terms, data preparation and checks shared by all the code
# paths of pglmm() ----

# Whether a (symmetric) matrix is positive definite, numerically: its
# Cholesky decomposition exists and its smallest pivot is not negligible
pglmm_is_pd <- function(V, tol = 1e-10) {
  V <- as.matrix(V)
  R <- tryCatch(chol(V), error = function(e) NULL)
  !is.null(R) && min(diag(R))^2 > tol * max(abs(diag(V)))
}

#' Process the var-cov matrices of random terms
#'
#' @param x A named list of var-cov matrices of random terms. The names should be the
#' group variables that are used as random terms and must be presented. The actual object
#' can be a phylogeny with class "phylo" or a prepared var-cov matrix. If it is a phylogeny,
#' we will prune it and then convert it to a var-cov matrix assuming brownian motion evolution.
#' We will also standardize all var-cov matrices to have determinant of one
#' (unless \code{standardize = FALSE}).
#' @param df A data frame that includes the group variables, i.e., names of \code{x}.
#' @param standardize Whether to standardize the var-cov matrices.
#' @return A named list, which includes the processed var-cov matrices of random terms.
#' @noRd

parse_conv_ranef <- function(x, df, standardize = TRUE) {
  if (is.null(names(x))) stop("conv_ranef list must have names")
  x2 <- x
  out_list <- lapply(1:length(x), function(i) {
    xx <- x[[i]]
    spl <- levels(df[, names(x)[i]])
    if ("phylo" %in% class(xx)) {
      # phylogeny
      if (length(setdiff(spl, xx$tip.label))) {
        pglmm_stop_missing_levels(names(x)[i], spl, xx$tip.label)
      }
      if (length(setdiff(xx$tip.label, spl))) {
        warning(paste0("Drop species from the phylogeny that are not in the variable ", names(x)[i]),
          call. = FALSE, immediate. = TRUE
        )
        xx <- ape::drop.tip(xx, setdiff(xx$tip.label, spl))
      }
      Vphy <- ape::vcv(xx)
      if (standardize) {
        Vphy <- Vphy / max(Vphy)
        Vphy <- Vphy / exp(determinant(Vphy)$modulus[1] / ape::Ntip(xx))
      }
      Vphy <- Vphy[spl, spl] # same order as species levels
    }

    if (inherits(xx, "dist")) xx <- as.matrix(xx)

    if (inherits(xx, c("matrix", "Matrix"))) {
      # already a cov matrix
      if (length(setdiff(spl, row.names(xx)))) {
        pglmm_stop_missing_levels(names(x)[i], spl, row.names(xx), "covariance matrix")
      }
      # keep the levels in the data only, in the same order as the levels
      # (also when the matrix has the same names but in a different order)
      if (is.null(colnames(xx)) || !identical(rownames(xx), colnames(xx))) {
        stop(paste0("The row and column names of the matrix of variable ", names(x)[i],
                    " must be the same and in the same order"))
      }
      xx <- xx[spl, spl, drop = FALSE]
      if (!pglmm_is_pd(xx)) {
        stop("The covariance matrix of variable ", names(x)[i], " in cov_ranef is not positive ",
             "definite (e.g. two levels with the same rows, such as tips of a phylogeny joined by ",
             "branches of length zero); pglmm() needs positive definite covariance matrices.",
             call. = FALSE)
      }
      if (standardize && abs(det(xx) - 1) > 0.0001) {
        warning("The cov matrix is not standardized, we will do this now...", call. = FALSE, immediate. = TRUE)
        xx <- xx / max(xx)
        xx <- xx / exp(determinant(xx)$modulus[1] / nrow(xx))
        if (abs(det(xx) - 1) > 0.0001) warning("Failed to standardize the cov matrix", call. = FALSE, immediate. = TRUE)
      }
      Vphy <- xx
    }
    x2[[i]] <<- xx # update the phylo or cov matrix
    Vphy
  })
  names(out_list) <- names(x)
  list(updated_orgi_list = x2, cleaned_list = out_list)
}

# ---- Random terms of the formula ----
#
# pglmm_parse_re() is the only function that reads the random terms of the
# formula: pglmm(), pglmm_compare(), prep_dat_pglmm(), the sparse and Laplace
# code (pglmm_sparse_terms(), pglmm_glmm_data()) and the Bayesian code all use
# its output. It checks the terms (with the same rules and messages for all
# code paths), writes them in a standard form, and returns:
# - formula: the formula with the terms in their standard form;
# - specs: one element per random term of the formula (after the changes
#   below), with lhs ("1" or the covariate of a slope), group (the grouping
#   variables, without "__"), phy (whether each grouping variable has a
#   phylogeny or covariance matrix, i.e. ends with "__"), nested (two grouping
#   variables, e.g. sp@site), grp (e.g. "sp__@site") and name ("1|sp__@site");
# - terms: the terms of the model in the order of the formula, as fitted by
#   all code paths: a non-nested phylogenetic term such as 1|sp__ gives two
#   terms, 1|sp (phy = FALSE) and 1|sp__ (phy = TRUE), unless phylo_only = TRUE
#   (pglmm_compare(): the phylogenetic term only, with phylo_only = TRUE);
# - vars: the variables of the random terms (grouping variables and covariates
#   of slopes).
# Changes to the formula:
# - (0 + x | sp) is written (x | sp) (a slope without intercept in phyr);
# - interactions of grouping factors (lme4), e.g. 1|sp:site (also from
#   1|site/sp), are written as nested terms, 1|sp@site;
# - a term that a phylogenetic term already includes as its non-phylogenetic
#   partner (1|sp with 1|sp__) is dropped, with a message (fitting it twice
#   gives two variances that cannot be told apart), and so is the same term
#   twice.
# The checks: grouping variables and covariates of slopes are columns of data,
# covariates are numeric, nested terms have two grouping variables, and terms
# have one covariate (no x1 + x2 | g). Parsing a formula already in its
# standard form gives the same formula (no messages).
pglmm_parse_re <- function(formula, data, phylo_only = FALSE) {
  bars <- unique(reformulas::findbars(formula))
  if (!length(bars)) {
    return(list(formula = formula, specs = list(), terms = list(), vars = character(0)))
  }
  rewritten <- FALSE
  specs <- vector("list", length(bars))
  for (k in seq_along(bars)) {
    x2 <- as.character(bars[[k]])
    lhs <- gsub("^0 ?[+] ?", "", x2[2])
    grp <- x2[3]
    term <- paste0("(", x2[2], " | ", grp, ")")
    if (grepl("+", lhs, fixed = TRUE)) {
      stop("Random terms with several covariates, such as ", term, ", are not allowed yet. ",
           "Split them into terms with one covariate each (e.g. (1 | sp) + (x | sp), or ",
           "(x || sp)); the terms are independent.", call. = FALSE)
    }
    if (lhs != x2[2]) {
      bars[[k]][[2]] <- str2lang(lhs)
      rewritten <- TRUE
    }
    if (grepl(":", grp, fixed = TRUE)) {
      # the interaction of two grouping factors (lme4), e.g. 1|sp:site, also
      # from 1|site/sp, is the nested term 1|sp@site: one level per
      # combination, independent
      grp <- gsub(" ?: ?", "@", grp)
      bars[[k]][[3]] <- str2lang(grp)
      rewritten <- TRUE
    }
    gp <- strsplit(grp, "@", fixed = TRUE)[[1]]
    if (length(gp) > 2) {
      stop("Nested terms have two grouping variables (e.g. 1 | sp@site), not ", length(gp),
           " as in ", term, ".", call. = FALSE)
    }
    g <- sub("__$", "", gp)
    miss <- setdiff(g, names(data))
    if (length(miss)) {
      stop("The grouping variable", if (length(miss) > 1) "s", " ", paste(miss, collapse = " and "),
           " of the random term ", term, if (length(miss) > 1) " are" else " is",
           " not in the data.", call. = FALSE)
    }
    if (lhs != "1") {
      if (lhs %nin% names(data)) {
        stop("The covariate ", lhs, " of the random slope ", term, " must be a column of data",
             if (make.names(lhs) != lhs) " (e.g. compute it as a new column first)", ".",
             call. = FALSE)
      }
      v <- data[[lhs]]
      if (!is.numeric(v) && !is.logical(v)) {
        stop("The covariate ", lhs, " of the random slope ", term, " is not numeric (class ",
             class(v)[1], "). Random slopes need numeric covariates; for a factor, use numeric ",
             "(0/1) columns for its levels, each in its own random term.", call. = FALSE)
      }
    }
    specs[[k]] <- list(lhs = lhs, group = g, phy = grepl("__$", gp), nested = length(gp) == 2,
                       grp = grp, name = paste0(lhs, "|", grp))
  }
  nm <- vapply(specs, `[[`, "", "name")
  partners <- unlist(lapply(specs, function(s) {
    if (!s$nested && s$phy && !phylo_only) paste0(s$lhs, "|", s$group)
  }))
  dup <- nm %in% partners
  if (any(dup)) {
    message("The random term", if (sum(dup) > 1) "s", " ",
            paste0("(", sub("|", " | ", nm[dup], fixed = TRUE), ")", collapse = " and "),
            " ", if (sum(dup) > 1) "are" else "is", " dropped: ",
            "a phylogenetic term such as (1 | sp__) already includes its non-phylogenetic ",
            "term (1 | sp).")
  }
  # the same term twice (e.g. 1|sp:site and 1|sp@site): kept once
  dup <- dup | duplicated(nm)
  if (any(dup) || rewritten) {
    rhs <- reformulas::nobars(formula)
    rhs <- rhs[[length(rhs)]]
    for (b in bars[!dup]) rhs <- call("+", rhs, call("(", b))
    formula[[length(formula)]] <- rhs
  }
  specs <- specs[!dup]

  # the terms of the model
  terms <- list()
  for (s in specs) {
    if (!s$nested && s$phy) {
      if (!phylo_only) {
        terms[[length(terms) + 1]] <- c(s[c("lhs", "group")], list(
          phy = FALSE, nested = FALSE, phylo_only = FALSE, name = paste0(s$lhs, "|", s$group)))
      }
      terms[[length(terms) + 1]] <- c(s[c("lhs", "group", "phy", "nested", "name")],
                                      list(phylo_only = phylo_only))
    } else {
      terms[[length(terms) + 1]] <- c(s[c("lhs", "group", "phy", "nested", "name")],
                                      list(phylo_only = FALSE))
    }
  }
  vars <- unique(unlist(lapply(specs, function(s) c(s$group, if (s$lhs != "1") s$lhs))))
  list(formula = formula, specs = specs, terms = terms, vars = vars)
}

# Rows of data without missing values in the response, the fixed effects
# (offsets included) and the variables of the random terms (vars)
pglmm_complete_rows <- function(formula, data, vars) {
  mf <- tryCatch(model.frame(reformulas::nobars(formula), data, na.action = NULL),
                 error = function(e) NULL)
  ok <- stats::complete.cases(data[, intersect(vars, names(data)), drop = FALSE])
  if (!is.null(mf)) ok <- ok & stats::complete.cases(mf)
  ok
}

# The phylogenetic variables of the terms of pglmm_parse_re() (those with a
# phylogeny or covariance matrix in cov_ranef)
pglmm_re_phy_vars <- function(terms) {
  unique(unlist(lapply(terms, function(t) t$group[t$phy])))
}

# repulsion of each term of pglmm_parse_re(): one value for each phylogenetic
# variable of nested terms, in the order of the formula. A single value
# applies to the phylogenies of species only, i.e. the first variable of
# nested terms (e.g. sp in 1|sp__@site__); phylogenies of the second variable
# (sites) are not inverted. As in Ives & Helmus (2011, Appendix A), repulsion
# comes from competition among related species. Returns the terms, with
# `repulsion` (one logical per grouping variable) for nested terms.
pglmm_re_repulsion <- function(terms, repulsion) {
  nested <- which(vapply(terms, function(t) t$nested, TRUE))
  n.rep <- sum(vapply(terms[nested], function(t) sum(t$phy), 1L))
  if (length(repulsion) == 1) {
    repulsion <- unlist(lapply(terms[nested], function(t) (repulsion & c(TRUE, FALSE))[t$phy]))
  }
  if (length(repulsion) != n.rep) {
    stop("The number of repulsion terms specified is not correct: please double check")
  }
  i <- 0
  for (k in nested) {
    phy <- terms[[k]]$phy
    r <- c(FALSE, FALSE)
    r[phy] <- as.logical(repulsion[i + seq_len(sum(phy))])
    i <- i + sum(phy)
    terms[[k]]$repulsion <- r
  }
  terms
}

# Whether a nested term without phylogeny (e.g. 1|sp@site) has one level per
# observation, i.e. each combination of its two variables appears once: it is
# then an observation-level term
pglmm_one_obs_per_cell <- function(data, g) {
  !anyDuplicated(data.frame(data[[g[1]]], data[[g[2]]]))
}

# Whether the model has the observation-level term 1|obs (terms of
# pglmm_parse_re()): poisson models and binomial models with more than one
# trial (cbind(successes, failures)), with add.obs.re = TRUE, unless a nested
# intercept without phylogeny (e.g. 1|sp@site) already has one level per
# observation. Shows the messages of earlier versions.
pglmm_re_obs <- function(terms, data, formula, family, add.obs.re = TRUE) {
  resp <- tryCatch(model.response(model.frame(reformulas::nobars(formula), data, na.action = NULL)),
                   error = function(e) NULL)
  # binomial data with trials: cbind(successes, failures), or proportions with
  # prior weights (the numbers of trials, see pglmm_glmm_data())
  trials <- is.array(resp) || !is.null(pglmm_prior_weights(data))
  if (!(family == "poisson" || (family == "binomial" && trials))) return(FALSE)
  obs.term <- any(vapply(terms, function(t) {
    t$nested && !any(t$phy) && t$lhs == "1" && pglmm_one_obs_per_cell(data, t$group)
  }, TRUE))
  if (obs.term) {
    if (add.obs.re) {
      message("It seems that you specified an observation-level random term already e.g. (1|sp@site);
                       we will set 'add.obs.re' to FALSE.")
    }
    return(FALSE)
  }
  if (add.obs.re) {
    message("We add an observation-level random term '1|obs' for poisson and binomial data.")
  } else {
    message("For poisson and binomial data, it would be a good idea to add an observation-level random term (add.obs.re = TRUE).")
  }
  add.obs.re
}

# Stop if the suggested packages pkgs (needed by the function fun) are not
# installed
pglmm_need_pkgs <- function(pkgs, fun) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, TRUE, quietly = TRUE)]
  if (length(miss)) {
    stop(fun, "() needs the package", if (length(miss) > 1) "s", " ",
         paste(miss, collapse = ", "), ". Install ", if (length(miss) > 1) "them" else "it",
         " with install.packages(", paste0("c(", paste0('"', miss, '"', collapse = ", "), ")"),
         ").", call. = FALSE)
  }
  invisible(TRUE)
}

# The prior weights of pglmm() and pglmm_compare(): evaluated as lm() and
# glm() do (a column of data, given by name, or a vector, with one value per
# row of data: wexpr is evaluated in data, then in env, the environment of the
# formula; if that fails, e.g. through a wrapper function, value() forces the
# argument where it was given), and kept in data as the column .pglmm_weights,
# so that the code paths that drop the rows with missing values keep them
# aligned; see pglmm_prior_weights().
pglmm_add_weights <- function(data, wexpr, env, value) {
  if (is.null(wexpr)) return(data)
  if (is.null(data)) stop("weights need data.", call. = FALSE)
  w <- tryCatch(eval(wexpr, data, env), error = function(e) value())
  if (is.null(w)) return(data)
  if (is.character(w) && length(w) == 1 && nrow(data) != 1) {
    if (w %nin% names(data)) stop("column \"", w, "\" (weights) not found in data.", call. = FALSE)
    w <- data[[w]]
  }
  if (!is.numeric(w) || length(w) != nrow(data)) {
    stop("weights should be numeric, with one value per row of data.", call. = FALSE)
  }
  if (any(!is.na(w) & (!is.finite(w) | w <= 0))) {
    stop("weights should be positive (drop the observations with weight zero).", call. = FALSE)
  }
  data$.pglmm_weights <- as.numeric(w)
  data
}

# the environment of a formula (as model.frame()), or env
pglmm_formula_env <- function(formula, env) {
  fe <- if (inherits(formula, "formula")) environment(formula)
  if (is.null(fe)) env else fe
}

# The prior weights of the rows of data (after rows were dropped), or NULL
pglmm_prior_weights <- function(data) {
  if (is.null(data) || ".pglmm_weights" %nin% names(data)) NULL else data$.pglmm_weights
}

# The arguments tree and tree_site of pglmm() and pglmm_plot_ranef(),
# deprecated since 2019 and defunct since phyr 1.9.0: the phylogenies go in
# cov_ranef
pglmm_defunct_tree <- function(tree, tree_site) {
  if (is.null(tree) && is.null(tree_site)) return(invisible(NULL))
  stop("the arguments `tree` and `tree_site` were removed in phyr 1.9.0; give the phylogenies ",
       "in `cov_ranef` instead, e.g. `cov_ranef = list(sp = tree, site = tree_site)`.",
       call. = FALSE)
}

# Errors for the phylogenies and covariance matrices of cov_ranef ----

# Levels of the grouping variable g that are not in the phylogeny (or
# covariance matrix) of cov_ranef: which ones, and a hint when they differ from
# the tip labels only in spaces vs underscores or capitalization
pglmm_stop_missing_levels <- function(g, lev, labels, what = "phylogeny") {
  miss <- setdiff(lev, labels)
  stop(length(miss), " of the ", length(lev), " levels of ", g, if (length(miss) > 1) " are" else " is",
       " not in the ", what, " of cov_ranef$", g, ": ", ps_some(miss), if (length(miss) <= 5) ".",
       ps_label_hint(miss, labels),
       "\nDrop these rows from the data, or give a ", what, " that has them.", call. = FALSE)
}

# A phylogenetic term (e.g. 1|sp__) without its phylogeny in cov_ranef
pglmm_stop_no_cov <- function(g, cov_ranef) {
  nm <- names(cov_ranef)
  stop("The phylogenetic term(s) of ", g, " (", g, "__) need a phylogeny or a covariance ",
       "matrix named ", g, " in cov_ranef, e.g. cov_ranef = list(", g, " = tree); ",
       if (length(nm)) paste0("cov_ranef has ", paste(nm, collapse = ", "), ".") else "cov_ranef is empty.",
       call. = FALSE)
}

# Checks of the arguments of pglmm() before fitting ----

# The family: a name of the families of pglmm(), with a suggestion for a
# misspelled one
pglmm_check_family <- function(family, bayes) {
  fams <- c("gaussian", pglmm_laplace_families)
  if (bayes) fams <- c("gaussian", "binomial", "poisson", "nbinomial", "zeroinflated.binomial",
                       "zeroinflated.poisson", "zeroinflated.nbinomial")
  if (is.character(family) && length(family) == 1 && family %in% fams) return(invisible(NULL))
  if (bayes && identical(family, "Gamma")) {
    stop('family = "Gamma" is only available with bayes = FALSE.', call. = FALSE)
  }
  near <- if (is.character(family) && length(family) == 1) {
    fams[utils::adist(tolower(family), tolower(fams)) <= 2]
  }
  stop("family = \"", paste(family, collapse = " "), "\" is not available",
       if (bayes) " with bayes = TRUE", "; the families are ", paste0("\"", fams, "\"", collapse = ", "),
       ".", if (length(near)) paste0(" Did you mean \"", near[1], "\"?"), call. = FALSE)
}

# The variables of the fixed effects that are neither in data nor in the
# environment of the formula
pglmm_check_vars <- function(formula, data) {
  v <- all.vars(reformulas::nobars(formula))
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  miss <- v[!v %in% names(data) & !vapply(v, exists, TRUE, envir = env)]
  if (length(miss)) {
    stop("The variable", if (length(miss) > 1) "s", " ", paste(miss, collapse = ", "),
         " of the formula ", if (length(miss) > 1) "are" else "is", " not in data.", call. = FALSE)
  }
  invisible(NULL)
}

# Random terms whose grouping variable (or combination of variables) has a
# single level in the data: their variance cannot be estimated (as lme4)
pglmm_check_levels <- function(terms, data) {
  for (t in terms) {
    g <- t$group
    key <- unique(do.call(paste, c(lapply(g, function(v) as.character(data[[v]])), sep = ":")))
    if (length(key) < 2) {
      stop("The grouping variable", if (length(g) > 1) "s", " of the random term ", t$name,
           " ha", if (length(g) > 1) "ve" else "s", " a single level in the data (", key[1],
           "), so its variance cannot be estimated: remove the term.", call. = FALSE)
    }
  }
  invisible(NULL)
}

# Warn for a nested intercept without phylogeny (e.g. 1|sp@site) with one level
# per observation (one observation per combination of its grouping
# variables): an observation-level term, whose variance cannot be estimated
# with gaussian data (it is confounded with the residual variance) and binary
# data, and overlaps the dispersion parameter of the families that have one
# (it is the observation-level term 1|obs of poisson and binomial data with
# trials, see pglmm_re_obs()). Returns the names of the terms that cannot be
# estimated (gaussian and binary data), which pglmm_fit() drops.
pglmm_check_obs_terms <- function(terms, data, formula, family) {
  family <- as.character(family)
  resp <- tryCatch(model.response(model.frame(reformulas::nobars(formula), data, na.action = NULL)),
                   error = function(e) NULL)
  binary <- family == "binomial" && !is.array(resp) && is.null(pglmm_prior_weights(data))
  spec <- pglmm_family_spec(family)
  disp <- if (!is.null(spec$disp)) {
    lab <- spec$disp$name[spec$disp$name != "zi"]
    if (length(lab)) lab
  }
  why <- if (family == "gaussian") {
    "its variance is confounded with the residual variance, so it cannot be estimated"
  } else if (binary) {
    "its variance cannot be estimated with binary (0/1) data"
  } else if (!is.null(disp)) {
    paste0("it overlaps the dispersion parameter (", disp[1], ") of the ", family,
           " family, so that the two are hard to separate")
  }
  if (is.null(why)) return(invisible(character(0)))
  drop <- family == "gaussian" || binary
  out <- character(0)
  for (t in terms) {
    if (t$nested && !any(t$phy) && t$lhs == "1" && pglmm_one_obs_per_cell(data, t$group)) {
      warning("The random term ", t$name, " has one level per observation (one observation for each ",
              "combination of ", paste(t$group, collapse = " and "), "): it is an observation-level ",
              "term, and ", why, ". ",
              if (drop) "It is dropped from the model (with bayes = FALSE)." else "Consider removing it",
              if (!drop) " (with several observations for some combinations, it is not an observation-level term).",
              call. = FALSE)
      if (drop) out <- c(out, t$name)
    }
  }
  invisible(out)
}

# The formula without the random terms named drop (as the names of the terms
# of pglmm_parse_re(), e.g. "1|sp@site")
pglmm_drop_terms <- function(formula, drop) {
  bars <- reformulas::findbars(formula)
  nm <- vapply(bars, function(b) paste0(deparse1(b[[2]]), "|", deparse1(b[[3]])), "")
  nm <- gsub(" ", "", nm)
  keep <- bars[!nm %in% gsub(" ", "", drop)]
  fixed <- reformulas::nobars(formula)
  rhs <- c(deparse1(fixed[[length(fixed)]]), vapply(keep, function(b) paste0("(", deparse1(b), ")"), ""))
  f <- stats::as.formula(paste(deparse1(fixed[[2]]), "~", paste(rhs, collapse = " + ")))
  environment(f) <- environment(formula)
  f
}
