\name{synth}
\docType{data}
\alias{synth}
\title{Illustrative dataset: sequences of five states}

\description{ The data represents the synthetic dataset used as an
illustrative example in the Journal of Statistical Software paper
discussing the use of the package.\cr
There are 5 states denoted as \code{A}, \code{B}, \code{C}, \code{D}, and \code{E}. Categorical sequences have lengths varying from 10 to 50. 
}

\usage{data(synth)}
\format{$data contains a vector of 250 strings representing categorical sequences; $id is the original classification vector.}
\source{Melnykov, V. (2015)}
\references{
Melnykov, V. (2014) Model-based biclustering of clickstream data, accepted by Computational Statistics & Data Analysis.

Melnykov, V. (2015) Model-based biclustering of clickstream data, under review at Journal of Statistical Software.
}
\seealso{click.read}

\examples{

data(synth)
head(synth$data)

# FUNCTION THAT REPLACES CHARACTER STATES WITH NUMERIC VALUES
repl.levs <- function(x, ch.lev){
	for (j in 1:length(ch.lev)) x <- gsub(ch.levs[j], j, x)
	return(x)
}

# DETECT ALL STATES IN THE DATASET
d <- paste(synth$data, collapse = " ")
d <- strsplit(d, " ")[[1]]
ch.levs <- levels(as.factor(d))

# CONVERT DATA TO THE FORM USED BY click.read()
S <- strsplit(synth$data, " ")
S <- sapply(S, repl.levs, ch.levs)
S <- sapply(S, as.numeric)
head(S)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dataset}% __ONLY ONE__ keyword per line
