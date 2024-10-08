#' Solve Linear Mixed Models
#'
#' Solve Linear Mixed Models using REML.
#'
#' A Linear Mixed Model (LMM) has the form
#' \deqn{y = X \beta + Z u + e, u ~ N(0,G), e ~ N(0,R)} where
#' \eqn{y} is a vector of observations, \eqn{\beta} is a vector with the fixed
#' effects, \eqn{u} is a vector with the random effects, and \eqn{e} a vector of
#' random residuals. \eqn{X} and \eqn{Z} are design matrices.
#'
#' LMMsolve can fit models where the matrices \eqn{G^{-1}} and \eqn{R^{-1}} are
#' a linear combination of precision matrices \eqn{Q_{G,i}} and \eqn{Q_{R,i}}:
#' \deqn{G^{-1} = \sum_{i} \psi_i Q_{G,i} \;, R^{-1} = \sum_{i} \phi_i Q_{R,i}}
#' where the precision parameters \eqn{\psi_i} and \eqn{\phi_i} are estimated
#' using REML. For most standard mixed models \eqn{1/{\psi_i}} are the variance
#' components and \eqn{1/{\phi_i}} the residual variances. We use a formulation
#' in terms of precision parameters to allow for non-standard mixed models using
#' tensor product splines.
#'
#' @param fixed A formula for the fixed part of the model. Should be of the
#' form "response ~ pred"
#' @param random A formula for the random part of the model. Should be of the
#' form "~ pred".
#' @param spline A formula for the spline part of the model. Should be of the
#' form "~ spl1D()", ~ spl2D()" or "~spl3D()". Generalized Additive Models (GAMs) can
#' also be used, for example "~ spl1D() + spl2D()"
#' @param group A named list where each component is a numeric vector
#' specifying contiguous fields in data that are to be considered as a
#' single term.
#' @param ginverse A named list with each component a symmetric matrix, the
#' precision matrix of a corresponding random term in the model. The row and
#' column order of the precision matrices should match the order of the
#' levels of the corresponding factor in the data.
#' @param weights A character string identifying the column
#' of data to use as relative weights in the fit. Default value NULL, weights are
#' all equal to one.
#' @param data A data.frame containing the modeling data.
#' @param residual A formula for the residual part of the model. Should be of
#' the form "~ pred".
#' @param family An object of class family specifying the distribution
#' and link function.
#' @param offset An optional numerical vector containing an a priori
#' known component to be included in the linear predictor during fitting.
#' @param tolerance A numerical value. The convergence tolerance for the
#' modified Henderson algorithm to estimate the variance components.
#' @param trace Should the progress of the algorithm be printed? Default
#' \code{trace = FALSE}.
#' @param maxit A numerical value. The maximum number of iterations for the
#' algorithm. Default \code{maxit = 250}.
#' @param theta initial values for penalty or precision parameters. Default
#' \code{NULL}, all precision parameters set equal to 1.
#' @param grpTheta a vector to give components the same penalty. Default
#' \code{NULL}, all components have a separate penalty.
#'
#' @return An object of class \code{LMMsolve} representing the fitted model.
#' See \code{\link{LMMsolveObject}} for a full description of the components in
#' this object.
#'
#' @examples
#' ## Fit models on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Fit the same model with genotype as random effect.
#' LMM1_rand <- LMMsolve(fixed = yield ~ rep,
#'                      random = ~gen,
#'                      data = john.alpha)
#'
#' ## Fit the model with a 1-dimensional spline at the plot level.
#' LMM1_spline <- LMMsolve(fixed = yield ~ rep + gen,
#'                        spline = ~spl1D(x = plot, nseg = 20),
#'                        data = john.alpha)
#'
#' ## Fit models on multipop data included in the package.
#' data(multipop)
#'
#' ## The residual variances for the two populations can be different.
#' ## Allow for heterogeneous residual variances using the residual argument.
#' LMM2 <- LMMsolve(fixed = pheno ~ cross,
#'                 residual = ~cross,
#'                 data = multipop)
#'
#' ## QTL-probabilities are defined by the columns pA, pB, pC.
#' ## They can be included in the random part of the model by specifying the
#' ## group argument and using grp() in the random part.
#'
#' # Define groups by specifying columns in data corresponding to groups in a list.
#' # Name used in grp() should match names specified in list.
#' lGrp <- list(QTL = 3:5)
#' LMM2_group <- LMMsolve(fixed = pheno ~ cross,
#'                       group = lGrp,
#'                       random = ~grp(QTL),
#'                       residual = ~cross,
#'                       data = multipop)
#'
#' @seealso \code{\link{LMMsolveObject}}, \code{\link{spl1D}},
#' \code{\link{spl2D}}, \code{\link{spl3D}}
#'
#' @importFrom stats model.frame terms model.matrix contrasts as.formula
#' terms.formula aggregate model.response var formula gaussian
#'
#' @export
LMMsolve <- function(fixed,
                     random = NULL,
                     spline = NULL,
                     group = NULL,
                     ginverse = NULL,
                     weights = NULL,
                     data,
                     residual = NULL,
                     family = gaussian(),
                     offset = 0,
                     tolerance = 1.0e-6,
                     trace = FALSE,
                     maxit = 250,
                     theta = NULL,
                     grpTheta = NULL) {
  ## Input checks.
  if (!inherits(data, "data.frame")) {
    stop("data should be a data.frame.\n")
  }
  if (!inherits(fixed, "formula") || length(terms(fixed)) != 3) {
    stop("fixed should be a formula of the form \"resp ~ pred\".\n")
  }
  if (!is.null(random) &&
      (!inherits(random, "formula") || length(terms(random)) != 2)) {
    stop("random should be a formula of the form \" ~ pred\".\n")
  }
  splErr <- paste("spline should be a formula of form \"~ spl1D() + ... + ",
                  "spl1D()\", \"~ spl2D()\" or \"~spl3D()\"\n")
  if (!is.null(spline)) {
    if (!inherits(spline, "formula")) stop(splErr)
    spline <- formula(paste((gsub(pattern = "LMMsolver::",
                                  replacement = "",
                                  as.character(spline))), collapse = ""))
    splTrms <- terms(spline, specials = c("spl1D", "spl2D", "spl3D"))
    splSpec <- attr(splTrms, "specials")
    if (length(terms(splTrms)) != 2 ||
        ## Spline formula should consist of splxD() terms and nothing else.
        length(unlist(splSpec)) != length(labels(terms(spline))) ||
        length(splSpec$spl2D) > 1 ||
        length(splSpec$spl3D) > 1) {
      stop(splErr)
    }
  }
  if (!is.null(ginverse) &&
      (!is.list(ginverse) ||
       length(names(ginverse)) == 0 ||
       (!all(sapply(X = ginverse, FUN = function(x) {
         (is.matrix(x) || spam::is.spam(x)) && isSymmetric(x)}))))) {
    stop("ginverse should be a named list of symmetric matrices.\n")
  }
  if (!is.null(residual) &&
      (!inherits(residual, "formula") || length(terms(residual)) != 2)) {
    stop("residual should be a formula of the form \" ~ pred\".\n")
  }
  if (!is.numeric(tolerance) || length(tolerance) > 1 || tolerance < 0) {
    stop("tolerance should be a positive numerical value.")
  }
  if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 0) {
    stop("maxit should be a positive numerical value.")
  }
  if (!is.null(grpTheta) &&
      (!is.numeric(grpTheta) || !isTRUE(all.equal(round(grpTheta),grpTheta)) ||
       max(grpTheta) != length(unique(grpTheta)))) {
    stop("grpTheta should be integers 1,2,...nGrp")
  }
  ## Check that all variables used in fixed formula are in data.
  data <- checkFormVars(fixed, data, naAllowed = FALSE)
  ## Remove NA for response variable from data.
  respVar <- all.vars(fixed)[attr(terms(fixed), "response")]
  respVarNA <- is.na(data[[respVar]])
  if (sum(respVarNA) > 0) {
    warning(sum(respVarNA), " observations removed with missing value for ",
            respVar, ".\n", call. = FALSE)
    data <- data[!respVarNA, ]
  }
  ## Check for weights should be done after removing observations with NA for
  ## response var. This prevents error messages when both response var and
  ## weight are NA.
  if (!is.null(weights)) {
    if (!(weights %in% colnames(data))) {
      stop("weights not defined in dataframe data")
    }
    w <- data[[weights]]
    if (!is.numeric(w) || sum(is.na(w)) != 0 || min(w) < 0) {
      stop("weights should be a numeric vector with non-negative values")
    }
  } else {
    w <- rep(1, nrow(data))
  }
  ## Remove observations with zero weights
  weightsZero <- w == 0
  if (sum(weightsZero) > 0) {
    # warning(sum(weightsZero), " observations removed with zero weights \n ",
    #        call. = FALSE)
    data <- data[!weightsZero, ]
    ## remove missing values for weight (default w=1).
    w <- w[!weightsZero]
  }
  ## Drop unused factor levels from data.
  data <- droplevels(data)
  ## Check random term for conditional factor
  condFactor <- condFactor(random, data)
  if (!is.null(condFactor)) {
    random <- condFactor$random
  }
  ## Check that all variables used in formulas are in data.
  chkGroup <- checkGroup(random, group, data)
  random <- chkGroup$random
  group <- chkGroup$group
  data <- checkFormVars(random, data)
  data <- checkFormVars(residual, data)
  ## Make random part.
  if (!is.null(random)) {
    mf <- model.frame(random, data, drop.unused.levels = TRUE, na.action = NULL)
    mt <- terms(mf)
    f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
    Z1 <- Matrix::sparse.model.matrix(mt, data = mf,
                                      contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                                             FUN = contrasts,
                                                             contrasts = FALSE))
    dim1.r <- table(attr(Z1, "assign"))[-1]
    term1.labels.r <- attr(mt, "term.labels")
    scFactor1 <- rep(1, length(dim1.r))
    ## Number of variance parameters (see Gilmour 1995) for each variance component
    varPar1 <- rep(1, length(dim1.r))
    if (ncol(Z1) > 1) {
      Z1 <- Z1[, -1, drop = FALSE]
      Z1 <- spam::as.spam.dgCMatrix(Z1)
    }
    else {
      Z1 <- NULL
    }
  } else {
    scFactor1 <- NULL
    dim1.r <- NULL
    term1.labels.r <- NULL
    Z1 <- NULL
    varPar1 <- NULL
  }
  if (!is.null(group)) {
    ndx <- unlist(group)
    dim2.r <- sapply(X = group, FUN = length)
    term2.labels.r <- names(group)
    scFactor2 <- rep(1, length(dim2.r))
    varPar2 <- rep(1, length(dim2.r))
    Z2 <- spam::as.spam(as.matrix(data[, ndx]))
  } else {
    dim2.r <- NULL
    term2.labels.r <- NULL
    scFactor2 <- NULL
    Z2 <- NULL
    varPar2 <- NULL
  }
  if (!is.null(condFactor)) {
    dim3.r <- condFactor$dim.r
    term3.labels.r <- condFactor$term.labels.r
    scFactor3 <- rep(1, length(dim3.r))
    varPar3 <- rep(1, length(dim3.r))
    Z3 <- condFactor$Z
  } else {
    dim3.r <- NULL
    term3.labels.r <- NULL
    scFactor3 <- NULL
    Z3 <- NULL
    varPar3 <- NULL
  }
  if (!(is.null(random) && is.null(group) && is.null(condFactor))) {
    Z <- spam::cbind.spam(Z1, Z2, Z3)
    dim.r <- c(dim1.r, dim2.r, dim3.r)
    term.labels.r <- c(term1.labels.r, term2.labels.r, term3.labels.r)
    scFactor <- c(scFactor1, scFactor2, scFactor3)
    varPar <- c(varPar1, varPar2, varPar3)
    e <- cumsum(dim.r)
    s <- e - dim.r + 1
    lGinv <- list()
    for (i in 1:length(dim.r)) {
      tmp <- rep(0, sum(dim.r))
      tmp[s[i]:e[i]] <- 1
      lGinv[[i]] <- spam::cleanup(spam::diag.spam(tmp))
    }
    names(lGinv) <- term.labels.r
  } else {
    Z <- NULL
    lGinv <- NULL
    dim.r <- NULL
    term.labels.r <- NULL
    scFactor <- NULL
    varPar <- NULL
  }
  ## Replace identity matrix with ginverse for random terms.
  if (!is.null(ginverse)) {
    ## Convert to spam.
    e <- cumsum(dim.r)
    s <- e - dim.r + 1
    for (i in seq_along(ginverse)) {
      ginvVar <- names(ginverse)[i]
      ginvMat <- ginverse[[i]]
      k <- which(term.labels.r == ginvVar)
      if (length(k) == 0) {
        stop("ginverse element ", ginvVar, " not defined in random part.\n")
      }
      if (dim.r[k] != nrow(ginvMat)) {
        stop("Dimensions of ", ginvVar, " should match number of levels ",
             "for corresponding factor in data.\n")
      }
      ndx <- s[k]:e[k]
      ## as.spam drops row and column names, so only converting to spam
      ## after all checking is done.
      lGinv[[k]][ndx, ndx] <- spam::as.spam(ginvMat)
    }
  }
  ## Make fixed part.
  mf <- model.frame(fixed, data, drop.unused.levels = TRUE)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
  X <- Matrix::sparse.model.matrix(mt, data = mf,
                                   contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                                          FUN = contrasts, contrasts = TRUE))
  term.labels.f <- attr(mt, "term.labels")

  q <- qr(as.matrix(X))
  if (q$rank != ncol(X)) {
    remCols <- q$pivot[-seq(q$rank)]
    ## Compare terms before and after removing extra columns.
    ## If a complete term is removed, it also has to be removed from the labels.
    f.terms.orig <- as.numeric(names(table(attr(X, "assign"))))

    dim.f.tab <- table(attr(X, "assign")[-remCols])
    dim.f <- as.numeric(dim.f.tab)
    X <- X[ , -remCols, drop = FALSE]
    f.terms.new <- as.numeric(names(dim.f.tab))
    if (!setequal(f.terms.orig, f.terms.new)) {
      term.labels.f <- term.labels.f[-setdiff(f.terms.orig, f.terms.new)]
    }
  } else {
    dim.f <- as.numeric(table(attr(X, "assign")))
  }
  nNonSplinesRandom <- length(dim.r)

  ## Convert to spam matrix and cleanup
  #Xs <- spam::as.spam.dgCMatrix(X)
  #Xs <- spam::cleanup(Xs)

  ## Add spline part.
  splResList <- NULL
  NomEffDimRan <- NULL
  if (!is.null(spline)) {
    splTerms <- labels(splTrms)
    Nterms <- length(splTerms)
    splResList <- list()
    for (i in 1:Nterms) {
      splRes <- eval(parse(text = splTerms[i]), envir = data, enclos = parent.frame())
      ## Multiple 1D gam models should have unique x variables.
      if (!is.null(term.labels.f) && !is.null(splRes$term.labels.f) &&
          splRes$term.labels.f %in% term.labels.f) {
        stop("x variables in 1D splines should be unique.\n")
      }
      splResList[[i]] <- splRes
      ## Add to design matrix fixed effect X.
      X <- cbind(X, splRes$X)
      ## Add to design matrix random effect Z.
      Z <- spam::cbind.spam(Z, splRes$Z)
      ## Expand matrices Ginv to the updated Z.
      lGinv <- expandGinv(lGinv, splRes$lGinv)
      ## A splxD model has x parameters.
      varPar <- c(varPar, length(splRes$lGinv))
      ## Add dims.
      dim.f <- c(dim.f, splRes$dim.f)
      dim.r <- c(dim.r, splRes$dim.r)
      ## Add nominal ED.
      NomEffDimRan <- c(NomEffDimRan, splRes$EDnom)
      ## Add labels.
      term.labels.f <- c(term.labels.f, splRes$term.labels.f)
      term.labels.r <- c(term.labels.r, splRes$term.labels.r)
      scFactor <- c(scFactor, splRes$scaleFactor)
    }
    # if splines of different dimensions are combined:
    if (length(splSpec[!sapply(X = splSpec, FUN = is.null)]) > 1) {
       #check whether variables are uniquely defined
        nameVar <- unlist(lapply(splResList, function(z) {names(z$x)}))
        if (length(nameVar)!=length(unique(nameVar))) {
        stop("variables in splines should be unique.\n")
      }
    }
  }

  ## Convert to spam matrix and cleanup
  Xs <- spam::as.spam.dgCMatrix(X)
  Xs <- spam::cleanup(Xs)

  if (nNonSplinesRandom > 0) {
    ## calculate NomEff dimension for non-spline part
    NomEffDimNonSplines <- calcNomEffDim(Xs, Z, dim.r[c(1:nNonSplinesRandom)], term.labels.r)
    ## combine with splines part
    NomEffDimRan <- c(NomEffDimNonSplines, NomEffDimRan)
  }

  ## Add intercept.
  if (attr(mt, "intercept") == 1) {
    term.labels.f <- c("(Intercept)", term.labels.f)
  }
  ## construct inverse of residual matrix R.
  lRinv <- constructRinv(df = data, residual = residual, weights = w)
  nRes <- length(lRinv)
  scFactor <- c(scFactor, rep(1, nRes))
  y <- model.response(mf)
  ## check whether the variance for response is not zero.
  if (is.null(residual)) {
    if (var(y) < .Machine$double.eps / 2) {
      stop("Variance response variable zero or almost zero.\n")
    }
  } else {
    resVar <- all.vars(residual)
    varGrp <- tapply(X = y, INDEX = data[[resVar]], FUN = var)
    ndxVar <- which(varGrp < .Machine$double.eps / 2)
    if (length(ndxVar) > 0) {
      levels_f <- levels(data[[resVar]])
      levelsNoVar <- paste(levels_f[ndxVar], collapse = ", ")
      stop("Variance response variable zero or almost zero for levels:\n",
           levelsNoVar, "\n")
    }
  }
  ## Fit the model.
  if (!is.null(theta)) {
    if (length(theta) != length(scFactor)) {
      stop("Argument theta has wrong length \n")
    }
    theta <- theta / scFactor
  } else {
    theta <- 1 / scFactor
  }
  if (!is.null(grpTheta)) {
    if (length(grpTheta) != length(scFactor)) {
      stop("Argument grpTheta has wrong length \n")
    }
  } else {
    grpTheta <- c(1:length(scFactor))
  }

  if (family$family == "gaussian") {
    obj <- sparseMixedModels(y = y, X = Xs, Z = Z, lGinv = lGinv, lRinv = lRinv,
                             tolerance = tolerance, trace = trace, maxit = maxit,
                             theta = theta, grpTheta = grpTheta)
    dev.residuals <- family$dev.resids(y, obj$yhat, w)
    deviance <- sum(dev.residuals)
  } else {
    ## MB, 23 jan 2023
    ## binomial needs global weights
    weights <- w
    nobs <- length(y)
    mustart <- etastart <- NULL
    eval(family$initialize)
    mu <- mustart
    eta <- family$linkfun(mustart)
    nNonRes <- length(theta) - nRes
    fixedTheta <- c(rep(FALSE, nNonRes), rep(TRUE, nRes))
    theta[(nNonRes + 1):(nNonRes + nRes)] <- 1
    trace_GLMM <- NULL
    for (i in 1:maxit) {
      deriv <- family$mu.eta(eta)
      z <- (eta - offset) + (y - mu)/deriv
      wGLM <- as.vector(deriv^2 / family$variance(mu))
      wGLM <- wGLM*w
      lRinv <- constructRinv(df = data, residual = residual, weights = wGLM)
      obj <- sparseMixedModels(y = z, X = Xs, Z = Z, lGinv = lGinv, lRinv = lRinv,
                               tolerance = tolerance, trace = trace, maxit = maxit,
                               theta = theta, fixedTheta = fixedTheta,
                               grpTheta = grpTheta)
      eta.old <- eta
      eta <- obj$yhat + offset
      mu <- family$linkinv(eta)
      theta <- obj$theta
      tol <- sum((eta - eta.old)^2) / sum(eta^2)

      aux.df <- data.frame(itOuter = rep(i, obj$nIter),
                           tol=c(rep(NA,obj$nIter-1), tol))
      trace_GLMM <- rbind(trace_GLMM, cbind(aux.df, obj$trace))
      if (trace) {
        cat("Generalized Linear Mixed Model iteration", i, ", tol=", tol, "\n")
      }
      if (tol < tolerance) {
        dev.residuals <- family$dev.resids(y, mu, w)
        deviance <- sum(dev.residuals)
        break;
      }
    }
  }
  ## Add names to ndx of coefficients.
  ndxCf <- seq_along(obj$a)
  ## Fixed terms.
  ef <- cumsum(dim.f)
  sf <- ef - dim.f + 1
  ndxCoefF <- nameCoefs(coefs = ndxCf, desMat = X, termLabels = term.labels.f,
                        s = sf, e = ef, data = data, type = "fixed")
  ## Random terms.
  er <- sum(dim.f) + cumsum(dim.r)
  sr <- er - dim.r + 1
  ndxCoefR <- nameCoefs(coefs = ndxCf, termLabels = term.labels.r, s = sr,
                        e = er, data = data, group = group, type = "random")
  ## Combine result for fixed and random terms.
  ndxCoefTot <- c(ndxCoefF, ndxCoefR)
  names(ndxCoefTot) <- c(term.labels.f, term.labels.r)
  ## Extract effective dimensions from fitted model.
  EffDimRes <- attributes(lRinv)$cnt
  EffDimNamesRes <- attributes(lRinv)$names
  NomEffDim <- c(NomEffDimRan, EffDimRes)
  # Calc upper bound for nominal effective dimension:
  N <- nrow(X)
  p <- ncol(X)
  NomEffDim <- pmin(NomEffDim, N - p)
  ## Make ED table for fixed effects.
  EDdf1 <- data.frame(Term = term.labels.f,
                      Effective = dim.f,
                      Model = dim.f,
                      Nominal = dim.f,
                      Ratio = rep(1, length(dim.f)),
                      Penalty = rep(0, length(dim.f)),
                      VarComp = rep(NA, length(dim.f)))
  ## Make ED table for random effects.
  EDdf2 <- data.frame(Term = obj$EDnames,
                      Effective = obj$ED,
                      Model = c(rep(dim.r, varPar), EffDimRes),
                      Nominal = NomEffDim,
                      Ratio = obj$ED / NomEffDim,
                      Penalty = obj$theta,
                      VarComp = c(rep(term.labels.r, varPar), EffDimNamesRes))
  ## Make full ED table.
  EDdf <- rbind(EDdf1, EDdf2)
  EDdf <- EDdf[-which(colnames(EDdf) == "VarComp")]
  rownames(EDdf) <- NULL
  ## Make variance table.
  varComp <- factor(EDdf2[["VarComp"]], levels = unique(EDdf2[["VarComp"]]))
  VarDf <- aggregate(x = EDdf2$Penalty, by = list(VarComp = varComp), FUN = sum)
  VarDf[["Variance"]] = 1 / VarDf[["x"]]
  VarDf <- VarDf[c("VarComp", "Variance")]
  ## Compute REML constant.
  constantREML <- -0.5 * log(2 * pi) * (length(y) - sum(dim.f))
  dim <- c(dim.f, dim.r)
  if (family$family == "gaussian") {
    trace <- obj$trace
  } else {
    trace <- trace_GLMM
  }
  return(LMMsolveObject(logL = obj$logL,
                        sigma2e = obj$sigma2e,
                        tau2e = obj$tau2e,
                        EDdf = EDdf,
                        varPar = varPar,
                        VarDf = VarDf,
                        theta = obj$theta,
                        coefMME = obj$a,
                        ndxCoefficients = ndxCoefTot,
                        yhat = obj$yhat,
                        residuals = obj$residuals,
                        nIter = obj$nIter,
                        y = y,
                        X = Xs,  # sparse format
                        Z = Z,
                        lGinv = lGinv,
                        lRinv = lRinv,
                        C = obj$C,
                        cholC = obj$cholC,
                        constantREML = constantREML,
                        dim = dim,
                        Nres = length(lRinv),
                        term.labels.f = term.labels.f,
                        term.labels.r = term.labels.r,
                        splRes = splResList,
                        family = family,
                        deviance = deviance,
                        trace = trace))
}
