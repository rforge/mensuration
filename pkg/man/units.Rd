\name{Units}
\alias{Units}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Unit conversions via the system GNU \dQuote{units} command
}

\description{
This routine will do unit conversions based on the GNU \dQuote{Units}
command (see References section). Any conversions that are legal there will work here.
}

\usage{
Units(x = c(3, 5), from = c("m"), to = "ft")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector of quantities in the \code{from} units to convert
	from}
  
  \item{from}{ a vector of units that apply to \code{x}; if scalar, it will be
	repeated to the length of \code{x} 
  }
  
  \item{to}{a vector of units that we want \code{x} converted to; it will be recycled as above
}
}


\details{ See examples below, just make sure things are compatible with
  the GNU \dQuote{Units} command. If any unit conversions are illegal,
  it will throw an error. It does not do much checking beyond this, so
  it will take whatever is given to it. This works flawlessly on Linux,
  and may be installed with \code{dnf} (Red Hat/Fedora and derivatives)
  or the system package manager. It requires downloading and installing
  GNU \dQuote{Units} on other operating systems (see below).  }

\value{
  A data frame with the first column the actual command that can be used in
  \code{system}, the second column contains the original quantities \code{x} to be
  converted, and the third column holds the conversions asked for.
}

\note{The routine uses a capital \dQuote{U} because it was conflicting
  with the \code{units} command in \R's base package.}

\references{
  Please see \url{http://www.gnu.org/software/units/units.html}
  information about \code{units} and its use (and installation on other systems).
}

\author{
JHG
}


\seealso{
\code{\link{system}}
}

\examples{
Units(1,'m','ft')                              #simple conversion factor
Units(c(25,50),'m^2/hectare','ft^2/acre')      #basal area quantities
Units(c(250,500),'ft^3/acre','m^3/hectare')    #volume quantities
Units(c(250,1,1),c('ft^3/acre','acre','hectare'),c('m^3/hectare','hectare','acre')) #mixture

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
