\name{findq}
\alias{findq}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Find a \dQuote{q} value from stand parameters }

\description{
This routine will simply return a \dQuote{q} value for a stand from the
desired basal area and number of trees in the largest class.
}

\usage{
findq(dbh.low = 0.5,
      dbh.up = 19.5,
      dc.width = 1,
      B = 80,
      N.max = 2,
      units = c("English", "metric"),
      runQuiet = FALSE,
      ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dbh.low}{ The lower \emph{midpoint} dbh.}

  \item{dbh.up}{ The upper \emph{midpoint} dbh.}

  \item{dc.width}{The DBH class width (should be integer) for all classes. }

  \item{B}{ Basal area per acre/hectare. }
  \item{N.max}{ Number of trees in the largest dbh class. }

  \item{units}{ Measurement units (character): \dQuote{English} or
              \dQuote{metric}. }

  \item{runQuiet}{  TRUE: no feedback; FALSE: some feedback. }

  \item{\dots}{ Arguments to be passed on to \code{dbhClassLimits} or gobbled up.}
}

\details{
This is very simple, and solves the well-known relation between the
diminution coefficient and the stand parameters. Just pass it the
diameter ranges, the basal area and number of trees in the largest class
to find the \dQuote{q} value.

The results can be easily corroborated with \code{\link{qDistnShiny}};
e.g., one can generate a stand table from the \dQuote{q} value (to two
decimal places) with \code{qDistnShiny}.

Note that if \code{\link{uniroot}} fails (this is what is used to find
the \dQuote{q} value), then it is probably because your input parameters
do not make sense.

}

\value{
  A list with\ldots
  \item{qValue}{ The \dQuote{q} value.}
  \item{qTable}{ The stand table for this stand. }
}
%\references{
%% ~put references to the literature/web site here ~
%}

\author{
JHG.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dbhClassLimits}} for how the diameter classes are handled.
}

\examples{
findq(0.5, 13.5, 1, 80, 8.0744732)$qValue       #q=1.1999988
findq(2, 18, 4, 100, 23.7)                      #q=1.199835
qv = findq(1, 23, 2, 18, 57.026315, 'metric')   #q=1.1999988
qv$qTable
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
