\name{fitMultiGEV}
\alias{fitMultiGEV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitMultiGEV(annMax, durLevels, start, lower, gevDl = durLevels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annMax}{
%%     ~~Describe \code{annMax} here~~
}
  \item{durLevels}{
%%     ~~Describe \code{durLevels} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{lower}{
%%     ~~Describe \code{lower} here~~
}
  \item{gevDl}{
%%     ~~Describe \code{gevDl} here~~
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
function (annMax, durLevels, start, lower, gevDl = durLevels) 
{
    stopifnot(length(start) == length(lower))
    nDl <- min(ncol(annMax), length(durLevels))
    optFun <- function(par) {
        attr(par, "gevDl") <- gevDl
        res <- -sum(sapply(1:nDl, function(cn) dMultiGEV(annMax[, 
            cn], durLevels[cn], par, log = TRUE)))
        if (is.infinite(res)) 
            return(1e+05)
        else res
    }
    optVal <- optim(start, optFun, lower = lower, method = "L-BFGS-B")
    res <- optVal$par
    attr(res, "optim") <- optVal
    attr(res, "gevDl") <- gevDl
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
