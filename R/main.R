
click.EM <- function(X, y = NULL, K, eps = 10^-10, r = 100, iter = 5, min.beta = 10^-3, min.gamma = 10^-3, scale.const = 1.0){

	if (K < 1) stop("Wrong number of mixture components K...\n")
	if (eps <= 0) stop("Wrong value of eps...\n")
	if (r < 1) stop("Wrong number of random restarts r...\n")
	if (iter < 1) stop("Wrong number of iterations iter...\n")
	if (min.beta < 0) stop("Wrong lower bound min.beta...\n")
	if (min.gamma < 0) stop("Wrong lower bound min.gamma...\n")
	if (scale.const <= 0) stop("Wrong value of scale.const...\n")

	xdims <- dim(X)
	p <- xdims[1]
	n <- xdims[3]

	x1 <- as.vector(X)
	id <- rep(0, n)

	alpha <- rep(0, K)
	beta <- rep(0, K*p)
	gamma <- rep(-1, p*p*K)
	z <- rep(0, n*K)
	l <- c(0, 0)

	if (!is.null(y)){

		y <- y - 1

		Q <- .C("runClickClust_", p1 = as.integer(p), K1 = as.integer(K), n1 = as.integer(n), x1 = as.integer(x1), y = as.integer(y), alpha = as.double(alpha), beta1 = as.double(beta), Pi1 = as.double(gamma), gamma1 = as.double(z), id = as.integer(id), e1 = as.double(eps), r1 = as.integer(r), shortem1 = as.integer(iter), l = as.double(l), lowbeta1 = as.double(min.beta), lowPi1 = as.double(min.gamma), scaleconst1 = as.double(scale.const), PACKAGE = "ClickClust")

	} else {

		Q <- .C("runClickClust", p1 = as.integer(p), K1 = as.integer(K), n1 = as.integer(n), x1 = as.integer(x1), alpha = as.double(alpha), Pi1 = as.double(gamma), gamma1 = as.double(z), id = as.integer(id), e1 = as.double(eps), r1 = as.integer(r), shortem1 = as.integer(iter), l = as.double(l), lowPi1 = as.double(min.gamma), scaleconst1 = as.double(scale.const), PACKAGE = "ClickClust")

	}

	a <- array(Q$Pi1, c(K, p, p))
	b <- array(NA, c(p, p, K))
	for (i in 1:K){
		b[,,i] <- t(a[i,,])
	}
	b[b == -1] <- NA

	if (!is.null(y)){

		return(list(z = t(matrix(Q$gamma1, nrow = K)), id = Q$id + 1, alpha = Q$alpha, beta = matrix(Q$beta, ncol = p, byrow = TRUE), gamma = b, logl = Q$l[1], BIC = Q$l[2]))

	} else {

		return(list(z = t(matrix(Q$gamma1, nrow = K)), id = Q$id + 1, alpha = Q$alpha, gamma = b, logl = Q$l[1], BIC = Q$l[2]))

	}

}
