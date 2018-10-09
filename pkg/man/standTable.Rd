\name{standTable}
\alias{standTable}
\title{ Stand Table Calculation
}
\description{ Calculates a simple stand table from some tree diameters. }

\usage{
standTable(dbh,
           wantBA = TRUE,
           isEnglish = TRUE,
           runQuiet = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbh}{ A vector of tree diameters. }
  \item{wantBA}{ \code{TRUE}: calculate basal area too; \code{FALSE}:
	just frequency. }
  \item{isEnglish}{ \code{TRUE}: English units; \code{FALSE}: metric
	units. }
  \item{runQuiet}{ \code{TRUE}: No feedback; \code{FALSE}: a little
	report. }
  \item{\dots}{Arguments passed on to \code{\link{dbhClassLimits}}. }
}
\details{ This routine calculates an unweighted stand table based on the
  diameters passed. The diameters are classed using
  \code{\link{dbhClassLimits}} to define the limits, and
  \code{\link{groupDBH}} to do the grouping. Basal area is calculated
  based on the raw (ungrouped/non-midpoint) diameters passed.
}
\value{
  A data frame with the stand table composed of the dbh class limits,
  tree frequency, and basal area (if desired).
}
\author{JHG. }
\note{ Sometime we might want to add some functionality for point
  sampling, etc.
}


\seealso{\code{\link{groupDBH}}, \code{\link{dbhClassLimits}}
}
\examples{

dbh = rexp(100, 0.25)
st = standTable(dbh, runQuiet=FALSE) #should produce a warning

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
