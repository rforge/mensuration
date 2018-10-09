\name{qDistnShiny}
\alias{qDistnShiny}
\alias{qDistnUI}
\alias{qDistnServer}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simple Shiny Application for \dQuote{q} Distributions
}

\description{
This is a simple Shiny application that will display in a webpage. It
allows one to manipulate simple inputs to change the parameters of the
diameter distribution (e.g. dbh class width) and the \dQuote{q} value
illustrating the results in graphs and tabular output.
}

\usage{
qDistnShiny(units = c("English", "metric"),
            returnVar = NA,
            useDT = TRUE,
            shiftZero = FALSE,
            ...)
qDistnUI(conv, useDT, ...)
qDistnServer(returnVar, conv, useDT, shiftZero, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{units}{ \dQuote{English}  or \dQuote{metric} units.}

  \item{returnVar}{ This is the name of a variable that will hold a list
     of relevant stand parameters, including the data frame shown in the
     webpage; it will be created in the user's workspace (global
     environment). It must be a character string. Use it like a normal
     assignment from the return of a function when run from the
     workspace---more details below.  }

  \item{useDT}{ \code{TRUE}: the display uses the DataTable format (default---very
             nice); \code{FALSE}: uses a default html table format.}
	 
  \item{conv}{ A list of conversion factors, no need to be concerned
	about this since you should not call \code{qDistnServer} or
	\code{qDistnUI} directly (see below).
  }
	 
  \item{shiftZero}{ \code{TRUE}: shift the lower DBH limit down to zero
	 diameter; \code{FALSE}: no shift, just include small trees to the
	 normal lower bound of the first dbh class.}
	
  \item{\dots}{
     Just gobbled at present.
}
}

\details{
Running the \code{qDistnShiny} function will display the application in
your web browser. Its use is self-explanatory, there are several sliders
for input that change both the graphs and the tabular output essentially
instantaneously. The best thing to do is play with it.

Both \code{qDistnServer} and \code{qDistnUI} are meant to be called only
from \code{qDistnShiny}. So please do not run these yourself (I know you
will since I said this).

Now, because of the way a Shiny application is structured, there is no
way to return a variable in the usual manner through \R's \code{return}
statement when calling a function: Shiny is in control of both the user
interface (UI) and the server. Therefore, one has to kludge to get the
final result back. In this case, the argument \var{returnVar} can be
used to capture the tabular results on exit (actually all throughout the
application's existence). It always gets saved to the user's workspace,
so treat it like any other object with regard to the fact that it will
delete any existing object of the same name. (If requested, an
environment argument could be added to change where it is saved.)

Note that when you request \code{shiftZero=TRUE}, the slider for maximum
diameter will normally be larger than the upper limit to the largest
class shown. This is because the slider is set up for integers that may
not align with the shift. This should be obvious and does not affect the
results in any way.

Please use the \dQuote{Exit Application} button to stop the
application gracefully. Otherwise at the command line type control-c or
escape, whatever works on your system.
}

\value{
  Returns an object of class \dQuote{\code{shiny.appobj}}. If you want
  to capture this, simply assign it to a variable. Otherwise, the return
  value gets \dQuote{printed}, which runs the application.

  However, if you choose to save the results of the run using the
  \code{returnVar} argument, then a \code{list} is placed in the session
  \code{.GlobalEnv} with the following components\ldots

  \item{stParams}{ Final prameter values from the session. }
  \item{stSum}{Stand totals: numbers, basal area, quadratic MSD.}
  \item{df}{A data frame with the tabular results.}
}

\references{
  Please see \url{http://shiny.rstudio.com/articles/} for more
  information about Shiny and its use.
}

\author{
JHG.
}


\seealso{
The \pkg{\link{shiny}} package. Also see \code{\link{crDistnShiny}}.
}
\examples{
#
# cut and paste these examples if desired...
#
\dontrun{
qDistnShiny('metric', 'q.df')  #run immediately
qsh = qDistnShiny()            #assign to qsh
str(qsh)                       #look at it
qsh                            #run it
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
