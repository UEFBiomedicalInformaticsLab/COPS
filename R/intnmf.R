
nmf.opt.k <- function (dat = dat, n.runs = 30, n.fold = 5, k.range = 2:8, 
                       result = TRUE, make.plot = TRUE, progress = TRUE, st.count = 10, 
                       maxiter = 100, wt = if (is.list(dat)) rep(1, length(dat)) else 1) 
{
  if (!is.list(dat) & length(wt) > 1) 
    stop("Weight is not applicable for single data")
  if (!is.list(dat)) 
    dat <- list(dat)
  n.dat <- length(dat)
  if (n.dat != length(wt)) 
    stop("Number of weights must match number of data")
  for (i in 1:n.dat) assign(paste("d", i, sep = ""), 
                            eval(parse(text = paste("dat[[", i, "]]", 
                                                    sep = ""))))
  for (i in 1:n.dat) {
    if (!all(eval(parse(text = paste("d", i, sep = ""))) >= 
             0)) 
      stop(paste("All values must be positive. There are -ve entries in dat", 
                 i, sep = ""))
  }
  set.seed(12345)
  n.sample <- nrow(dat[[1]])
  CPI <- matrix(NA, length(k.range), n.runs)
  dimnames(CPI) <- list(paste("k", k.range, sep = ""), 
                        paste("run", 1:n.runs, sep = ""))
  count <- 0
  for (i in 1:n.runs) {
    for (k in k.range) {
      R.ind <- NULL
      random.sample <- sample(seq(n.sample), n.sample)
      for (j in 1:n.fold) {
        test.sample <- random.sample[(round((j - 1) * 
                                              n.sample/n.fold) + 1):round(j * n.sample/n.fold)]
        train.sample <- setdiff(random.sample, test.sample)
        d.train <- lapply(dat, function(x) x[train.sample, 
                                             , drop = FALSE])
        d.test <- lapply(dat, function(x) x[test.sample, 
                                            , drop = FALSE])
        fit.train <- nmf.mnnals(dat = d.train, k = k, 
                                maxiter = maxiter, st.count = st.count, n.ini = 1, 
                                ini.nndsvd = FALSE, seed = FALSE, wt = wt)
        for (m in 1:n.dat) {
          if (n.dat > 1) {
            assign(paste("H.train", m, sep = ""), 
                   eval(parse(text = paste("fit.train$H[[", 
                                           m, "]]", sep = ""))))
          }
          else {
            assign("H.train", eval(parse(text = "fit.train$H")))
          }
        }
        XHt <- 0
        for (m in 1:n.dat) {
          if (n.dat > 1) {
            XHt <- XHt + sqrt(wt[m]) * d.test[[m]] %*% 
              t(eval(parse(text = paste("H.train", 
                                        m, sep = ""))))
          }
          else {
            XHt <- XHt + d.test[[m]] %*% t(eval(parse(text = "H.train")))
          }
        }
        HHt <- 0
        for (m in 1:n.dat) {
          if (n.dat > 1) {
            HHt <- HHt + sqrt(wt[m]) * eval(parse(text = paste("H.train", 
                                                               m, sep = ""))) %*% t(eval(parse(text = paste("H.train", 
                                                                                                            m, sep = ""))))
          }
          else {
            HHt <- HHt + eval(parse(text = "H.train")) %*% 
              t(eval(parse(text = "H.train")))
          }
        }
        W.predict <- XHt %*% MASS::ginv(HHt)
        W.predict <- W.predict + abs(min(W.predict))
        predicted.cluster.mem <- apply(W.predict, 1, 
                                       which.max)
        fit.test <- nmf.mnnals(dat = d.test, k = k, maxiter = maxiter, 
                               st.count = st.count, n.ini = 1, ini.nndsvd = FALSE, 
                               seed = FALSE, wt = wt) # Error in svd(X) : a dimension is zero
        computed.cluster.mem <- fit.test$clusters
        R.ind <- c(R.ind, mclust::adjustedRandIndex(predicted.cluster.mem, 
                                                    computed.cluster.mem))
        count <- count + 1
        if (progress & round(count/(n.runs * length(k.range) * 
                                    n.fold) * 100, 0) %in% seq(5, 100, by = 5)) {
          message(paste(round(count/(n.runs * length(k.range) * 
                                       n.fold) * 100, 0), "% complete", sep = ""))
          flush.console()
        }
      }
      CPI[k - 1, i] <- mean(R.ind)
    }
  }
  if (make.plot) {
    dev.new(width = 4, height = 5)
    plot(k.range, CPI[, 1], ylim = c(min(CPI), max(CPI)), 
         pch = 20, main = "", xlab = "k", ylab = "CPI")
    for (m in 2:n.runs) points(k.range, CPI[, m], pch = 20)
    lines(k.range, apply(CPI, 1, mean), col = "red", 
          lwd = 2)
    mtext("Optimum k", outer = TRUE, cex = 1, line = -2)
  }
  if (result) 
    return(CPI)
}

nmf.mnnals <- function (dat = dat, k = k, maxiter = 200, st.count = 20, n.ini = 30, 
                        ini.nndsvd = TRUE, seed = TRUE, wt = if (is.list(dat)) rep(1, 
                                                                                   length(dat)) else 1) 
{
  if (seed) 
    set.seed(12345)
  if (!is.list(dat) & length(wt) > 1) 
    stop("Weight is not applicable for single data")
  if (!is.list(dat)) 
    dat <- list(dat)
  n <- nrow(dat[[1]])
  n.dat <- length(dat)
  if (n.dat != length(wt)) 
    stop("Number of weights must match number of data")
  for (i in 1:n.dat) assign(paste("d", i, sep = ""), 
                            eval(parse(text = paste("dat[[", i, "]]", 
                                                    sep = ""))))
  for (i in 1:n.dat) {
    if (!all(eval(parse(text = paste("d", i, sep = ""))) >= 
             0)) 
      stop(paste("All values must be positive. There are -ve entries in dat", 
                 i, sep = ""))
  }
  min.f.WH <- NULL
  W.list <- NULL
  for (i in 1:length(dat)) assign(paste("H", i, ".list", 
                                        sep = ""), NULL)
  convergence.list <- NULL
  consensus.list <- NULL
  for (j in 1:n.ini) {
    abs.diff <- NA
    consensus <- matrix(0, n, n)
    if (j == 1 && ini.nndsvd) {
      for (i in 1:n.dat) assign(paste("tmp.H", i, 
                                      sep = ""), IntNMF:::.nndsvd.internal(dat[[i]], k, 
                                                                           flag = 0)$H)
      tmp.H <- NULL
      for (i in 1:n.dat) tmp.H <- c(tmp.H, list(eval(parse(text = paste("tmp.H", 
                                                                        i, sep = "")))))
      W <- W.fcnnls(x = tmp.H, y = dat, weight = wt)$coef
      W <- t(W)
      rm(tmp.H)
      rm(list = paste("tmp.H", 1:length(dat), sep = ""))
    }
    else {
      W <- matrix(runif(n * k, min = min(unlist(dat)), 
                        max = max(unlist(dat))), nrow = n)
    }
    connect.old <- matrix(0, nrow = n, ncol = n)
    count <- 1
    iter <- 1
    convergence <- NULL
    while ((iter < maxiter) && (count < st.count)) {
      W <- sweep(W, 2, pmax(sqrt(colSums(W^2)), .Machine$double.eps), 
                 "/")
      for (i in 1:n.dat) assign(paste("H", i, sep = ""), 
                                H.fcnnls(W, dat[[i]])$coef) # Error in svd(X) : a dimension is zero
      tmp.H <- NULL
      for (i in 1:n.dat) tmp.H <- c(tmp.H, list(eval(parse(text = paste("H", 
                                                                        i, sep = "")))))
      W <- W.fcnnls(x = tmp.H, y = dat, weight = wt)$coef
      W <- t(W)
      rm(tmp.H)
      if (iter == 1) {
        for (i in 1:n.dat) assign(paste("d", i, 
                                        ".old", sep = ""), W %*% eval(parse(text = paste("H", 
                                                                                         i, sep = ""))))
        f.WH <- 0
        for (i in 1:n.dat) f.WH <- f.WH + wt[i] * sum((eval(parse(text = paste("d", 
                                                                               i, sep = ""))) - eval(parse(text = paste("d", 
                                                                                                                        i, ".old", sep = ""))))^2)
      }
      else {
        for (i in 1:n.dat) assign(paste("d", i, 
                                        ".new", sep = ""), W %*% eval(parse(text = paste("H", 
                                                                                         i, sep = ""))))
        abs.diff <- 0
        for (i in 1:n.dat) abs.diff <- abs.diff + sum(abs(eval(parse(text = paste("d", 
                                                                                  i, ".new", sep = ""))) - eval(parse(text = paste("d", 
                                                                                                                                   i, ".old", sep = "")))))/sum(eval(parse(text = paste("d", 
                                                                                                                                                                                        i, ".old", sep = ""))))
        for (i in 1:n.dat) assign(paste("d", i, 
                                        ".old", sep = ""), eval(parse(text = paste("d", 
                                                                                   i, ".new", sep = ""))))
        f.WH <- 0
        for (i in 1:n.dat) f.WH <- f.WH + wt[i] * sum((eval(parse(text = paste("d", 
                                                                               i, sep = ""))) - eval(parse(text = paste("d", 
                                                                                                                        i, ".new", sep = ""))))^2)
      }
      f.WH <- sqrt(f.WH)
      clust.mem <- apply(W, 1, which.max)
      tmp1 <- matrix(rep(clust.mem, n), nrow = n, byrow = T)
      tmp2 <- matrix(rep(clust.mem, n), nrow = n, byrow = F)
      connect.new <- ifelse(tmp1 == tmp2, 1, 0)
      if (all(connect.new == connect.old)) 
        count <- count + 1
      else count <- 0
      convergence <- rbind(convergence, c(iter, count, 
                                          ifelse(all(connect.new == connect.old), 1, 0), 
                                          abs.diff, f.WH))
      consensus <- consensus + connect.new
      connect.old <- connect.new
      iter <- iter + 1
      tol <- abs.diff
      rm(tmp1, tmp2, clust.mem, connect.new)
    }
    colnames(convergence) <- c("iter", "count", 
                               "stability", "abs.diff", "f.WH")
    consensus <- consensus/(iter - 1)
    min.f.WH <- c(min.f.WH, f.WH)
    W.list <- c(W.list, list(W))
    for (i in 1:n.dat) assign(paste("H", i, ".list", 
                                    sep = ""), c(eval(parse(text = paste("H", 
                                                                         i, ".list", sep = ""))), list(eval(parse(text = paste("H", 
                                                                                                                               i, sep = ""))))))
    convergence.list <- c(convergence.list, list(convergence))
    consensus.list <- c(consensus.list, list(consensus))
    H.list <- NULL
    for (i in 1:n.dat) H.list <- c(H.list, list(eval(parse(text = paste("H", 
                                                                        i, ".list", sep = "")))))
    names(H.list) <- paste("H", 1:n.dat, ".list", 
                           sep = "")
    if (j > 1) {
      if (min(min.f.WH[-j]) < min.f.WH[j]) {
        W.list[-which.min(min.f.WH)] <- "Not an Optimum Solution"
        for (i in 1:n.dat) H.list[[i]][-which.min(min.f.WH)] <- "Not an Optimum Solution"
        convergence.list[-which.min(min.f.WH)] <- "Not an Optimum Solution"
        consensus.list[-which.min(min.f.WH)] <- "Not an Optimum Solution"
      }
      else {
        W.list[-j] <- "Not an Optimum Solution"
        for (i in 1:n.dat) H.list[[i]][-j] <- "Not an Optimum Solution"
        convergence.list[-j] <- "Not an Optimum Solution"
        consensus.list[-j] <- "Not an Optimum Solution"
      }
    }
  }
  W <- W.list[[which.min(min.f.WH)]]
  H <- NULL
  if (n.dat > 1) {
    for (i in 1:n.dat) {
      assign(paste("H", i, sep = ""), H.list[[i]][[which.min(min.f.WH)]])
      H <- c(H, list(eval(parse(text = paste("H", 
                                             i, sep = "")))))
    }
    names(H) <- paste("H", 1:n.dat, sep = "")
  }
  else {
    assign("H", H.list[[i]][[which.min(min.f.WH)]])
    names(H) <- "H"
  }
  consensus <- consensus.list[[which.min(min.f.WH)]]
  dimnames(consensus) <- list(rownames(W), rownames(W))
  convergence <- convergence.list[[which.min(min.f.WH)]]
  clusters <- apply(W, 1, which.max)
  output <- list(consensus = consensus, W = W, H = H, convergence = convergence, 
                 min.f.WH = min.f.WH, clusters = clusters)
  return(output)
}

H.fcnnls <- function (x = x, y = y, verbose = FALSE, pseudo = FALSE, eps = 0) 
{
  if (any(dim(y) == 0L)) {
    stop("Empty target matrix 'y' [", paste(dim(y), 
                                            collapse = " x "), "]")
  }
  if (any(dim(x) == 0L)) {
    stop("Empty regression variable matrix 'x' [", 
         paste(dim(x), collapse = " x "), "]")
  }
  C <- x
  A <- y
  nObs = nrow(C)
  lVar = ncol(C)
  if (nrow(A) != nObs) 
    stop("C and A have imcompatible sizes")
  pRHS = ncol(A)
  W = matrix(0, lVar, pRHS)
  iter = 0
  maxiter = 3 * lVar
  CtC = crossprod(C)
  CtA = crossprod(C, A)
  K = adj.cssls(CtC, CtA)
  Pset = K > 0
  K[!Pset] = 0
  D = K
  Fset = which(colSums(Pset) != lVar)
  oitr = 0
  while (length(Fset) > 0) {
    oitr = oitr + 1
    if (verbose && oitr > 5) 
      cat(sprintf("%d ", oitr))
    K[, Fset] = adj.cssls(CtC, CtA[, Fset, drop = FALSE], 
                          Pset[, Fset, drop = FALSE]) # Error in svd(X) : a dimension is zero
    Hset = Fset[colSums(K[, Fset, drop = FALSE] < eps) > 
                  0]
    if (length(Hset) > 0) {
      nHset = length(Hset)
      alpha = matrix(0, lVar, nHset)
      while (nHset > 0 && (iter < maxiter)) {
        iter = iter + 1
        alpha[, 1:nHset] = Inf
        ij = which(Pset[, Hset, drop = FALSE] & (K[, 
                                                   Hset, drop = FALSE] < eps), arr.ind = TRUE)
        i = ij[, 1]
        j = ij[, 2]
        if (length(i) == 0) 
          break
        hIdx = (j - 1) * lVar + i
        negIdx = (Hset[j] - 1) * lVar + i
        alpha[hIdx] = D[negIdx]/(D[negIdx] - K[negIdx])
        alpha.inf <- alpha[, 1:nHset, drop = FALSE]
        minIdx = max.col(-t(alpha.inf))
        alphaMin = alpha.inf[minIdx + (0:(nHset - 1) * 
                                         lVar)]
        alpha[, 1:nHset] = matrix(alphaMin, lVar, nHset, 
                                  byrow = TRUE)
        D[, Hset] = D[, Hset, drop = FALSE] - alpha[, 
                                                    1:nHset, drop = FALSE] * (D[, Hset, drop = FALSE] - 
                                                                                K[, Hset, drop = FALSE])
        idx2zero = (Hset - 1) * lVar + minIdx
        D[idx2zero] = 0
        Pset[idx2zero] = FALSE
        K[, Hset] = adj.cssls(CtC, CtA[, Hset, drop = FALSE], 
                              Pset[, Hset, drop = FALSE])
        Hset = which(colSums(K < eps) > 0)
        nHset = length(Hset)
      }
    }
    W[, Fset] = CtA[, Fset, drop = FALSE] - CtC %*% K[, Fset, drop = FALSE]
    Jset = which(colSums((ifelse(!(Pset[, Fset, drop = FALSE]), 
                                 1, 0) * W[, Fset, drop = FALSE]) > eps) == 0)
    Fset = setdiff(Fset, Fset[Jset])
    if (length(Fset) > 0) {
      mxidx = max.col(t(ifelse(!Pset[, Fset, drop = FALSE], 
                               1, 0) * W[, Fset, drop = FALSE]))
      Pset[(Fset - 1) * lVar + mxidx] = TRUE
      D[, Fset] = K[, Fset, drop = FALSE]
    }
  }
  list(coef = K, Pset = Pset)
}

adj.cssls <- function (CtC, CtA, Pset = NULL) 
{
  K = matrix(0, nrow(CtA), ncol(CtA))
  if (is.null(Pset) || length(Pset) == 0 || all(Pset)) {
    K = MASS::ginv(CtC) %*% CtA
  }
  else {
    lVar = nrow(Pset)
    pRHS = ncol(Pset)
    codedPset = as.numeric(2^(seq(lVar - 1, 0, -1)) %*% Pset)
    sortedPset = sort(codedPset)
    sortedEset = order(codedPset)
    breaks = diff(sortedPset)
    breakIdx = c(0, which(breaks > 0), pRHS)
    for (k in seq(1, length(breakIdx) - 1)) {
      cols2solve = sortedEset[seq(breakIdx[k] + 1, breakIdx[k + 1])]
      vars = Pset[, sortedEset[breakIdx[k] + 1]] # All FALSE
      if (any(vars)) {
        K[vars, cols2solve] = MASS::ginv(CtC[vars, vars]) %*% CtA[vars, cols2solve]
      }
    }
  }
  K
}

W.fcnnls <- function (x = x, y = y, weight = wt, verbose = FALSE, pseudo = FALSE, eps = 0) 
{
  wt = NULL
  for (i in 1:length(x)) {
    if (any(dim(y[[i]]) == 0L)) {
      stop("Empty target matrix 'y' [", paste(dim(y), 
                                              collapse = " x "), "]")
    }
    if (any(dim(x[[i]]) == 0L)) {
      stop("Empty regression variable matrix 'x' [", 
           paste(dim(x), collapse = " x "), "]")
    }
  }
  C <- t(x[[1]])
  A <- t(y[[1]])
  nObs = nrow(C)
  lVar = ncol(C)
  if (nrow(A) != nObs) 
    stop("C and A have imcompatible sizes")
  pRHS = ncol(A)
  W = matrix(0, lVar, pRHS)
  iter = 0
  maxiter = 3 * lVar
  CtC <- 0
  for (i in 1:length(x)) CtC <- CtC + sqrt(weight[i]) * t(x[[i]] %*% 
                                                            t(x[[i]]))
  CtA <- 0
  for (i in 1:length(x)) CtA <- CtA + sqrt(weight[i]) * t(y[[i]] %*% 
                                                            t(x[[i]]))
  K = adj.cssls(CtC, CtA)
  Pset = K > 0
  K[!Pset] = 0
  D = K
  Fset = which(colSums(Pset) != lVar)
  oitr = 0
  while (length(Fset) > 0) {
    oitr = oitr + 1
    if (verbose && oitr > 5) 
      cat(sprintf("%d ", oitr))
    K[, Fset] = adj.cssls(CtC, CtA[, Fset, drop = FALSE], 
                          Pset[, Fset, drop = FALSE])
    Hset = Fset[colSums(K[, Fset, drop = FALSE] < eps) > 
                  0]
    if (length(Hset) > 0) {
      nHset = length(Hset)
      alpha = matrix(0, lVar, nHset)
      while (nHset > 0 && (iter < maxiter)) {
        iter = iter + 1
        alpha[, 1:nHset] = Inf
        ij = which(Pset[, Hset, drop = FALSE] & (K[, 
                                                   Hset, drop = FALSE] < eps), arr.ind = TRUE)
        i = ij[, 1]
        j = ij[, 2]
        if (length(i) == 0) 
          break
        hIdx = (j - 1) * lVar + i
        negIdx = (Hset[j] - 1) * lVar + i
        alpha[hIdx] = D[negIdx]/(D[negIdx] - K[negIdx])
        alpha.inf <- alpha[, 1:nHset, drop = FALSE]
        minIdx = max.col(-t(alpha.inf))
        alphaMin = alpha.inf[minIdx + (0:(nHset - 1) * 
                                         lVar)]
        alpha[, 1:nHset] = matrix(alphaMin, lVar, nHset, 
                                  byrow = TRUE)
        D[, Hset] = D[, Hset, drop = FALSE] - alpha[, 
                                                    1:nHset, drop = FALSE] * (D[, Hset, drop = FALSE] - 
                                                                                K[, Hset, drop = FALSE])
        idx2zero = (Hset - 1) * lVar + minIdx
        D[idx2zero] = 0
        Pset[idx2zero] = FALSE
        K[, Hset] = adj.cssls(CtC, CtA[, Hset, drop = FALSE], 
                              Pset[, Hset, drop = FALSE])
        Hset = which(colSums(K < eps) > 0)
        nHset = length(Hset)
      }
    }
    W[, Fset] = CtA[, Fset, drop = FALSE] - CtC %*% K[, Fset, 
                                                      drop = FALSE]
    Jset = which(colSums((ifelse(!(Pset[, Fset, drop = FALSE]), 
                                 1, 0) * W[, Fset, drop = FALSE]) > eps) == 0)
    Fset = setdiff(Fset, Fset[Jset])
    if (length(Fset) > 0) {
      mxidx = max.col(t(ifelse(!Pset[, Fset, drop = FALSE], 
                               1, 0) * W[, Fset, drop = FALSE]))
      Pset[(Fset - 1) * lVar + mxidx] = TRUE
      D[, Fset] = K[, Fset, drop = FALSE]
    }
  }
  list(coef = K, Pset = Pset)
}