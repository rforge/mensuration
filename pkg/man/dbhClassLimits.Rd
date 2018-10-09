\name{dbhClassLimits}
\alias{dbhClassLimits}
\title{ Make lower, mid and upper DBH class limits }
\description{
A very simple little routine to generate dbh class limits and the
midpoint of the range for each class.
}
\usage{
dbhClassLimits(dbh.low = 1,
               dbh.up = 20,
               dbh.width = 1,
               forceIntegerWidth = TRUE,
               adjustLowest = NULL,
               classNames = c('dbh.low','dbh.mid','dbh.hi'),
               runQuiet = TRUE,
               ... )
}
\arguments{
  \item{dbh.low}{ The lower \emph{midpoint} dbh.}
  \item{dbh.up}{ The upper \emph{midpoint} dbh.}
  \item{dbh.width}{ The DBH class width (should be integer) for all
	classes, except possibly the first.}
  \item{forceIntegerWidth}{TRUE: dbh class widths must be integer; FALSE: can be real.}
  \item{classNames}{ a vector of length 3 for the labels of the data frame returned; always
                     in the form lower, mid and upper.}
  \item{adjustLowest}{  NULL, NA, or <=0 means no adjustment; other positive means to
                        subtract this amount from the lower bound of the smallest class.}
  \item{runQuiet}{ TRUE: no feedback; FALSE: some feedback.} 
  \item{...}{ Arguments to be passed on or gobbled up.}
}
\details{
There's really not much too this. I just put it into a package because I
end up re-inventing it so much over the year and find it handy to now
have to call upon.
}
\value{ A data frame with the lower, mid, and upper limts for all
  diameter classes as request.
}
\author{
J.H.G.
}


\seealso{ \code{\link{groupDBH}} }

\examples{
dbhClassLimits()
\dontrun{
dbhClassLimits(1,10,3.2) # warning and error
}
dbhClassLimits(3.2/2,10,3.2, forceIntegerWidth=FALSE) # no warning or error
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
