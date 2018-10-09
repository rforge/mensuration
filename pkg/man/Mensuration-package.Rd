\name{Mensuration-package}
\alias{Mensuration-package}
\alias{Mensuration}
\docType{package}
\title{
Miscellaneous Mensuration Routines
}
\description{
This package begain as a place to put a few handy mensuration routines
that I tend to use over and over. It now has been extended to include a
few other routines that incorporate the \pkg{shiny} package, such as one
to look at \dQuote{q} distributions.
}

\details{
\tabular{ll}{
Package: \tab Mensuration\cr
Type: \tab Package\cr
Version: \tab 1.1-3\cr
Date: \tab 2018-10-4\cr
License: \tab gpl\cr
LazyLoad: \tab yes\cr
}

See the individual package routines for more details. 
} %details

\section{List of functions\ldots}{

  \tabular{ll}{
	\code{\link{groupDBH}} \tab Places diameters in classes\cr
	\code{\link{dbhClassLimits}} \tab Sets class limits for \code{groupDBH}\cr
	\code{\link{standTable}} \tab Create a stand table using the above
	two routines\cr
	\code{\link{Units}} \tab A wrapper for the \code{units} shell
	command\cr
	\code{\link{findq}} \tab Determine \dQuote{q} from stand basal area
	and number of trees in the maximum dbh class\cr
	\code{\link{qDistnShiny}} \tab A shiny app for playing with \dQuote{q}\cr
	\code{\link{crDistnShiny}} \tab A shiny app for the Chapman-Richards
	distribution\cr
  }

} %section  


\author{JHG }

\references{
  Please see \url{http://shiny.rstudio.com/articles/} for more
  information about Shiny and its use.
}

%\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%}
