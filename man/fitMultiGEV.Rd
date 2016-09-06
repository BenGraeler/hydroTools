\name{fitMultiGEV}
\alias{fitMultiGEV}

\title{A fitting routine for the multiGEV}

\description{
A fitting routine for the multiGEV
}

\usage{
fitMultiGEV(annMax, durLevels, start,
            lower = rep(0, length(start)), gevDl = durLevels)
}

\arguments{
  \item{annMax}{A numeric matrix; holding the annual maxima (rows) per duration level (columns).}
  \item{durLevels}{A nuemric vector; the duration levels used during the estimation.}
  \item{start}{A numeric vector; starting values for the parameter optimisation (\code{\link{optim}}).}
  \item{lower}{A numeric vector; lower bounds used durng the optimisation (\code{\link{optim}}).}
  \item{gevDl}{A numeric vector; the duration levels used in the (\code{\link{multiGEV}}).}
}
\value{A list with slots for each GEV parameter (\code{loc}, \code{scale}, \code{shape}). Two attributes carry the number of \code{gevDl} and the output of \code{\link{optim}}. 
}
\author{
Benedikt Graeler}

\seealso{\code{\link{multiGEV}}}

\examples{
# set up initial (1) and increment (2) parameters of GEV for location sclae and shape
multiGEVpar <- c(50,5,5, # location
                 10,.5,.5, # scale
                 0,.5,.5) # shape

# set the base duration levels 
attr(multiGEVpar, "gevDl") <- c(1,6,24)

# what is the implicit correlation in the rMultiGEV?

smplRnd <- rMultiGEV(1000, c(1,3,6,12,24), multiGEVpar)
fitMultiGEV(smplRnd, c(1,3,6,12,24), multiGEVpar, gevDl = c(1,6,24))
}

\keyword{ distribution }
