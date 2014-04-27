
calc.props <- function(y){

	sums <- apply(y, 1, sum)
	sums <- ifelse(sums == 0, 1, sums)
	p <- sweep(y, 1, FUN = "/", sums)

	return(p)

}


calc.color <- function(colors, levels, x){

	z <- which(x <= levels)[1]-1

	return(colors[z])

}


click.plot <- function(X, y = NULL, file = NULL, id, margs = c(0.1, 0.1, 0.1, 0.1), font.cex = 2, font.col = "black", cell.cex = 1, cell.lwd = 1.3, cell.col = "red", sep.lwd = 1.3, sep.col = "magenta", obs.lwd = 0.5){

	K <- max(id)
	p <- dim(X)[1]
	n <- dim(X)[3]

	P1 <- X

	for (i in 1:n) P1[,,i] <- calc.props(X[,,i])

	Nk <- NULL
	P <- P1
	last <- 0
	if (!is.null(y)) y.new <- NULL

	# sort data according to group assignments

	for (k in 1:K){

		ind <- which(id == k)
		nk <- length(ind)
		Nk <- c(Nk, nk)
		P[,,(last+1):(last+nk)] <- P1[,,ind]
		last <- last + nk

		if (!is.null(y)){
			y.new <- c(y.new, y[ind])
		}

	}

	Nk.cum <- cumsum(Nk)


	col.levels <- 10
	levels = seq(0.1, 1.0, length.out = col.levels)
	color.back <- rgb(1, 0.980, 0.980)
	colors <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEC661", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026")
	
	col.vec <- rep(NA, p-1)
	levels <- c(0, levels)

	grid <- seq(0, 1, length.out = p + 1)
	grid.step <- grid[2] - grid[1]

	par(mar = margs)
	if (is.null(y)){
		plot( c(-grid.step/2, 1), c(0, 1 + grid.step/2), type = "n", xlab = "", ylab = "", axes = FALSE)
	} else {
		plot( c(-grid.step/2, 1 + grid.step), c(0, 1 + grid.step/2), type = "n", xlab = "", ylab = "", axes = FALSE)
	}
	box()

	# state numbers

	y1 <- 1 + grid.step / 2.5
	for (j in 1:p){
		x1 <- (grid[j] + grid[j+1]) / 2
		text(x1, y1, j, cex = font.cex, col = font.col)
	}

	x1 <- -grid.step / 2.5
	for (j in 1:p){
		y1 <- (grid[j] + grid[j+1]) / 2
		text(x1, y1, p+1-j, cex = font.cex, col = font.col)
	}

	# margin between cells
	eps <- grid.step / 20 / cell.cex


	# observation lines

	step <- (grid.step - 2 * eps) / n

	for (j in 1:p){

		for (i in 1:p){

			curr <- P[i,j,]

			for (h in 1:n){

				lines(c(grid[j]+eps*1.2, grid[j+1]-eps*1.2), c(grid[p-i+1]+eps+h*step, grid[p-i+1]+eps+h*step), col = calc.color(colors, levels, curr[h]), lwd = obs.lwd)

			}

			if (K != 1){

				for (k in 1:(K-1)){

					lines(c(grid[j]+eps*1.2, grid[j+1]-eps*1.2), c(grid[p-i+1]+eps+(Nk.cum[k]+0.5)*step, grid[p-i+1]+eps+(Nk.cum[k]+0.5)*step), col = sep.col, lwd = sep.lwd)

				}

			}

		}

	}


	# cell frames

	for (j in 1:p){
		x1 <- grid[j]
		y1 <- grid[j]
		for (i in 1:p){
			polygon(c(grid[j]+eps, grid[j+1]-eps, grid[j+1]-eps, grid[j]+eps),
				c(grid[p-i+1]+eps, grid[p-i+1]+eps, grid[p-i+2]-eps, grid[p-i+2]-eps),
				border = cell.col, lwd = cell.lwd)
		}

	}


	if (!is.null(y)){

		# additional column of cells to represent betas

		for (j in 1:p){

			for (i in 1:n){

				lines(c(max(grid) + 5*eps*1.2, max(grid) + grid.step / 2 + 5*eps*1.2), c(grid[p-j+1]+eps+i*step, grid[p-j+1]+eps+i*step), col = calc.color(colors, levels, y.new[i] == j), lwd = obs.lwd)

			}

			if (K != 1){

				for (k in 1:(K-1)){

					lines(c(max(grid) + 5*eps*1.2, max(grid) + grid.step / 2 + 5*eps*1.2), c(grid[p-j+1]+eps+(Nk.cum[k]+0.5)*step, grid[p-j+1]+eps+(Nk.cum[k]+0.5)*step), col = sep.col, lwd = sep.lwd)

				}	

			}

		}


		# cell frames

		for (j in 1:p){
			polygon(c(max(grid) + grid.step / 2 + 5*eps*1.2, max(grid) + 5*eps*1.2, max(grid) + 5*eps*1.2, max(grid) + grid.step / 2 + 5*eps*1.2), c(grid[p-j+1]+eps, grid[p-j+1]+eps, grid[p-j+2]-eps, grid[p-j+2]-eps), border = cell.col, lwd = cell.lwd)
		}

	}

	if (!is.null(file)) dev.copy2pdf(file = file)

}


