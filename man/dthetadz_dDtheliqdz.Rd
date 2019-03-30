\name{dthetadz_dDtheliqdz}
\alias{dthetadz_dDtheliqdz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculation of producto of derivative of theta and D theta,liq
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function is defined to calculate mixed term derivative products using lagging formula and Newton formula for linearizing nonlinear term
}
\usage{
dthetadz_dDtheliqdz(i, j, X = matrix(NA, nrow = 4, ncol = 4), tim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{tim}{
%%     ~~Describe \code{tim} here~~
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
}
\author{
%%  ~~who you are~~
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
function (i, j, X = matrix(NA, nrow = 4, ncol = 4), tim)
{
    if (j < 3) {
        result <- 1
    }
    else {
        result <- 1
    }
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
