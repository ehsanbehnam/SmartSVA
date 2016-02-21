require(sva)
require(isva)

# Computing p-value from an F-test.
f.pval <- function (dat, orth11, orth01, y.norm, rss00, df00)  {
	
	n <- dim(dat)[2]
	
	df11 <- dim(orth11)[2]
	df01 <- dim(orth01)[2]
	
	prj11 <- dat %*% orth11
	prj01 <- dat %*% orth01
	
	rss11 <- y.norm - rowSums(prj11 * prj11)
	rss01 <- y.norm - rowSums(prj01 * prj01)
	
	fstats <- ((rss01 - rss11)/(df11 - df01))/(rss11/(n - df11))
	p1 <- 1 - pf(fstats, df1 = (df11 - df01), df2 = (n - df11))
	
	fstats <- ((rss00 - rss01)/(df01 - df00))/(rss01/(n - df01))
	p2 <- 1 - pf(fstats, df1 = (df01 - df00), df2 = (n - df01))
	
	return(list(p1=p1, p2=p2))
}

# Most arguments are the same as SVA
# Additional two parameters:
## alpha: Determines the initial point for optimization.
#         It affects the convergence speed. 
## epsilon (convergence threshold): the spearman correlation 
#          between posterior probabilities of consecutive runs.
#          epsilon=0.005 also gives very reasonable results;
#          however, we suggest epsilon=1e-3 as a conservative
#          threshold.
#
## B (iteration number): typically changes from 20 to 100.
smartsva <-  function(dat, mod, mod0 = NULL, n.sv, B = 100, 
		alpha=0.25, epsilon=1e-3, VERBOSE = F) {  
	if (is.null(mod0)) {
		mod0 <- mod[, 1]
	}
	
	qr.obj <- qr(mod)
	orth1 <- qr.Q(qr.obj)
	uu <- eigen(crossprod(dat - tcrossprod(dat %*% orth1, orth1)), 
			symmetric=T)$vectors[, 1:n.sv, drop=F]
	
	# Precompute the quantites
	y.norm <- rowSums(dat * dat)
	mod00 <- cbind(mod0)
	orth00 <- qr.Q(qr(mod00))
	prj00 <- dat %*% orth00
	rss00 <- y.norm - rowSums(prj00 * prj00)
	df00 <- dim(orth00)[2]
	
	if (VERBOSE)
	  cat(paste("Iteration (out of", B, "):"))
	
	i = 0
	rho = 0
	
	while (i < B && rho < 1 - epsilon) {
		i <- i + 1
		mod11 <- cbind(mod, uu)
		mod01 <- cbind(mod0, uu)
		
		orth11 <- qr.Q(qr(mod11))
		orth01<- qr.Q(qr(mod01))
		
		ptmp <- f.pval(dat, orth11, orth01, y.norm, rss00, df00)
		
		if (i == 1) {
			pprob.b <- (1 - sva:::edge.lfdr(ptmp[['p1']])^alpha)
		} else {
			pprob.b <- (1 - sva:::edge.lfdr(ptmp[['p1']]))
		}
		
		pprob.gam <- (1 - sva:::edge.lfdr(ptmp[['p2']]))
		pprob <- pprob.gam * (1 - pprob.b)
		
		uu <- eigen(crossprod(dat * pprob - rowMeans(dat * pprob)), 
				symmetric=T)$vectors[, 1:n.sv, drop=F]
		# Update spearman Rho.
		if (i > 1) {
			rho <- cor(x=pprob, y=p.prev, use="pairwise.complete.obs",
					method="spearman")
			p.prev <- pprob
		}else{
			p.prev <- pprob
		}
		if (VERBOSE)
		  cat(paste(i, " ", rho, "\n"))
	}
	
	sv <- uu[, 1:n.sv]
	retval <- list(sv = sv, rho = rho, iter = i, n.sv = n.sv)
	return(retval)
}

###################################################################
# A simple example

# Methylation M values (CpG by Sample)
Y <- matrix(rnorm(200*10000), 10000, 200)                 
df <- data.frame(pred=gl(2, 100))

# Determine the number of SVs
Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
# Add one to be conservative
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
mod <- model.matrix( ~ pred, df)
sv.obj <- smartsva(Y, mod, mod0=NULL, n.sv=n.sv)
