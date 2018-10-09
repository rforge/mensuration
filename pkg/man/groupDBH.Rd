\name{groupDBH}
\alias{groupDBH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Group DBH values into Classes.
}
\description{ This is used in conjunction with
  \code{\link{dbhClassLimits}} to simply group the raw diameters into
  class midpoints.
}
\usage{
groupDBH(dbh, classLimits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbh}{ a vector of diameter values to be classed. }
  \item{classLimits}{ A data frame with lower, mid and upper class
	limits from \code{\link{dbhClassLimits}}. }
}
\details{
This just computes the diameter classes based on the midpoints of the
class for all dbh values passed. If any diameters lie outside the chosen
limits, they are marked as \code{NA} and a warning is issued.
}
\value{ A list with\ldots
  \item{dbh.grp}{ A vector with midpoint diameter classes. }
  \item{idx}{ A vector of the integer index values of each diameter in
	the class limits data frame. }
}
\author{ JHG. }


\seealso{ \code{\link{dbhClassLimits}}, or for a complete stand table
  along the lines of the example below, \code{\link{standTable}}. }

\examples{
lims = dbhClassLimits()
dbh = rexp(40, .25)      #some diameters
gd = groupDBH(dbh, lims)
rbind(gd$dbh.grp, dbh)
table(gd$dbh.grp)
dbhcl = factor(groupDBH(dbh, lims)$dbh.grp, levels=lims[,2])
as.data.frame(table(dbhcl))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
