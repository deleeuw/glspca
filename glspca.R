library(RSpectra)

ulsPCA <- function(x, p) {
  s <- svd(x, nu = p, nv = p)
  a <- s$u
  b <- s$v %*% diag(s$d[1:p])
  return(list(a = a, b = b, ab = tcrossprod(a, b)))
}

ulsAdd <- function(x) {
  m <- mean(x)
  r <- apply(x, 1, mean) - m
  s <- apply(x, 2, mean) - m
  return(list(
    m = m,
    r = r,
    s = s,
    rs = outer(r, s, "+") + m
  ))
}

ulsBoth <- function(x, p) {
  h1 <- ulsAdd(x)
  h2 <- ulsPCA(x - h1$rs, p)
  return(list(
    m = h1$m,
    r = h1$r,
    s = h1$s,
    a = h2$a,
    b = h2$b,
    y = h1$rs + h2$ab
  ))
}

glsLoss <- function(x, y, u, v) {
  d <- x - y
  return(sum(v * crossprod(d, (u %*% d))))
}

glsAdd <- function(x,
                   u,
                   v,
                   type = 3,
                   p = 2,
                   itmax = 10000,
                   eps = 1e-6,
                   verbose = FALSE) {
  yold <- switch(type, 
                 ulsAdd(x)$rs,
                 ulsPCA(x, p)$ab,
                 ulsBoth(x, p)$y
  )
  sold <- glsLoss(x, yold, u, v)
  lbdm <- eigs_sym(u, 1)$values * eigs_sym(v, 1)$values
  itel <- 1
  repeat {
    d <- x - yold
    g <- yold + u %*% d %*% v / lbdm
    ynew <- switch(type, 
          ulsAdd(g)$rs,
          ulsPCA(g, p)$ab,
          ulsBoth(g, p)$y
    )
    snew <- glsLoss(x, ynew, u, v)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    sold <- snew
    yold <- ynew
  }
  return(list(y = ynew, loss = snew, itel = itel))
}