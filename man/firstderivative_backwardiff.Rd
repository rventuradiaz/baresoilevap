\name{firstderivative_backwardiff}
\alias{firstderivative_backwardiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Partial first derivative in backward difference approximation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the partial first derivative of any dependent vector to any independent vector. Returns a Vector
}
\usage{
firstderivative_backwardiff(j = 4, Y = matrix(rep(0, 4), nrow = 4, ncol = 1), X = matrix(rep(0, 4), nrow = 4, ncol = 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Klaus A Hoffmann & Steve T Chiang. Computational Fluid Dynamics Vol 1. Wichita, Kan: Engineering Eductation System, 2000.
Junsei Kondo and Nobuko Saigusa. Modelling the Evaporation from Base Soil with a Formula for Vaporization in the Soil Pores. Meteorological Society of Japan.
}
\author{
%%  ~~who you are~~
Rafael Ventura
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (j = 4, Y = matrix(rep(0, 4), nrow = 4, ncol = 1), X = matrix(rep(0,
    4), nrow = 4, ncol = 1))
{
    nrow <- dim(X)[[1]]
    jminus1 <- if (j > 1) {
        j - 1
    }
    else {
        1
    }
    paddingj <- matrix(rep(0, nrow - j), nrow = nrow - j, ncol = 1)
    paddingjminus1 <- matrix(rep(0, nrow - jminus1), nrow = nrow -
        jminus1, ncol = 1)
    xj <- rbind(paddingj, X[1:j, , drop = FALSE])
    xjminus1 <- rbind(paddingjminus1, X[1:xminus1, , drop = FALSE])
    yj <- rbind(paddingj, Y[1:j, , drop = FALSE])
    yjminus1 <- rbind(paddingjminus1, X[1:xminus1, , drop = FALSE])
    dxj <- xj - xjminus1
    dydx <- mapply(function(x, y) x/y, yj - yjminus1, dxj)
    return(dydx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
