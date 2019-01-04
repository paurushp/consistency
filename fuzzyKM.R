##### Fuzzy Clustering K-Means
fuzzy_k_means=function (X, k, m, RS, startU){
    n = nrow(X)
    p = ncol(X)
    conv = 1e-09
    maxit = 1e+6
    if(is.null(rownames(X))) 
      rn = paste("Obj", 1:n, sep = " ")
    else 
      rn = rownames(X)
    if(is.null(colnames(X))) 
      cn = paste("Var", 1:p, sep = " ")
    else 
      cn = colnames(X)
    X = as.matrix(X)
    Xraw = X
    rownames(Xraw) = rownames(X)
    colnames(Xraw) = colnames(X)
    X = scale(X, center = TRUE, scale = TRUE)[, ]
    value = vector(length(RS), mode = "numeric")
    cput = vector(length(RS), mode = "numeric")
    it = vector(length(RS), mode = "numeric")
    check = 0
    func.opt = 10^10 * sum(X^2)
    for(rs in 1:RS){
        if((rs == 1) & (check != 1)) 
            U = startU
        else{
            set.seed(rs)
            U = matrix(runif(n * k, 0, 1), nrow = n, ncol = k)
            U = U/apply(U, 1, sum)
        }
        D = matrix(0, nrow = n, ncol = k)
        H = matrix(0, nrow = k, ncol = p)
        U.old = U + 1
        iter = 0
        cputime = system.time({
            while((sum(abs(U.old - U)) > conv) && (iter < maxit)) {
                iter = iter + 1
                U.old = U
                for(c in 1:k) H[c, ] = (t(U[, c]^m) %*% X)/sum(U[,c]^m)
                for(i in 1:n) {
                  for(c in 1:k) {
                    D[i, c] = sum((X[i, ] - H[c, ])^2)
                  }
                }
                for(i in 1:n) {
                  if(min(D[i, ]) == 0) {
                    U[i, ] = rep(0, k)
                    U[i, which.min(D[i, ])] = 1
                  }
                  else{
                    for(c in 1:k) {
                      U[i, c] = ((1/D[i, c])^(1/(m - 1)))/sum(((1/D[i,])^(1/(m - 1))))
                    }
                  }
                }
            }
        })
        func = sum((U^m) * D)
        cput[rs] = cputime[1]
        value[rs] = func
        it[rs] = iter
        if(func < func.opt) {
            U.opt = U
            H.opt = H
            func.opt = func
        }
    }
    rownames(H.opt) = paste("Clus", 1:k, sep = " ")
    colnames(H.opt) = cn
    rownames(U.opt) = rn
    colnames(U.opt) = rownames(H.opt)
    names(value) = paste("Start", 1:RS, sep = " ")
    names(cput) = names(value)
    names(it) = names(value)
    names(k) = c("Number of clusters")
    names(m) = c("Parameter of fuzziness")
    if(stand != 1) 
        stand = 0
    names(stand) = c("Standardization (1=Yes, 0=No)")
    clus = cl.memb(U.opt)
    out = list()
    out$U = U.opt
    out$H = H.opt
    out$F = NULL
    out$clus = clus
    out$medoid = NULL
    out$value = value
    out$cput = cput
    out$iter = it
    out$k = k
    out$m = m
    out$ent = NULL
    out$b = NULL
    out$vp = NULL
    out$delta = NULL
    out$stand = stand
    out$Xca = X
    out$X = Xraw
    out$call = match.call()
    class(out) = c("fclust")
    return(out)
}
