\name{crDistnShiny}
\alias{crDistnShiny}
\alias{crDistnUI}
\alias{crDistnServer}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Simple Shiny application for the Chapman-Richards Distribution
}

\description{ This is a simple Shiny application that will display in a
webpage. It allows one to explore the range of shapes of the
numbers density using sliders to change the parameters of the diameter
distribution (i.e., the vital rate parameters) illustrating the results
in graphs and tabular output.  }

\usage{
crDistnShiny(R = 7,
             M = 0.051,
             recruitRange = c(1, 100),
             mortRange = c(0.001, 0.5),
             m = 0.4,
             gamma = 0.053,
             eta = 0.4,
             mRange = c(-2, 0.78),
             gammaRange = c(0.001, 0.1),
             etaRange = c(0, 1),
             d0 = 1e-08,
             units = c('English', 'metric'),
             returnVar = NA,
             useDT = TRUE,
             maxBA = 200,
             includeTiny = FALSE,
             ...)
crDistnUI(params.run,
          mortRange,
          recruitRange,
          mRange,
          gammaRange,
          etaRange,
          useDT,
          ...
         )
crDistnServer(growth,         #function
              recruitment,    #function
              mortality,      #function
              ssd.cr,         #function
              params.run,
              d0,
              returnVar,
              conv,
              useDT,
              maxBA,
              includeTiny,
              ...
             )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{R}{ The initial recruitment rate into the smallest class in terms of
    numbers. }

  \item{M}{ The initial per capita mortality rate. }

  \item{recruitRange}{ The slider range for \code{R}. }

  \item{mortRange}{  The slider range for \code{M}. }

  \item{m}{ The initial \eqn{m} parameter value for the Chapman-Richards model. }

  \item{gamma}{ The initial \eqn{\gamma}{gamma} parameter value for the
      Chapman-Richards model. }

  \item{eta}{ The initial \eqn{\eta}{eta} parameter value for the
      Chapman-Richards model. }

  \item{mRange}{ The slider range for \code{m}. }

  \item{gammaRange}{  The slider range for \code{gamma}. }

  \item{etaRange}{  The slider range for \code{eta}. }

  \item{d0}{ Lower diameter limit for calculations, this must be
        slightly larger than zero! }

  \item{units}{ \dQuote{English}  or \dQuote{metric} units.}

  \item{returnVar}{  This is the name of a variable that will hold a list
     of relevant stand parameters, including the data frame shown in the
     webpage; it will be created in the user's workspace (global
     environment). It must be a character string. Use it like a normal
     assignment from the return of a function when run from the
     workspace---more details below. }

  \item{useDT}{ \code{TRUE}: use the Shiny DataTable for nicer display
     of the data; \code{FALSE}: use a normal html table for the table display. }

  \item{maxBA}{ This is the maximum basal area stocking allowed, before
	the application shows a warning message. }

  \item{includeTiny}{ \code{TRUE}: include all trees down to \code{d0}
       diameter; \code{FALSE}: only include small trees to the lower
       bound of the first dbh class. When \code{TRUE}, the routine will
       remind you that the smallest class is wider than all other
       classes (except perhaps the largest, depending on the maximum DBH).}
	 
  \item{conv}{ A list of conversion factors, no need to be concerned
	about this since you should not call \code{crDistnServer} 
	directly (see below).
  }

  \item{growth}{ A growth function of the form \code{growth =
    function(d, params, ...)}, where \var{d} is diameter and \var{params}
    are as in \var{params.run} below. This argument should not be a
    concern as \code{crDistnServer} should only be accessed from
    \code{crDistnShiny}. 
  }

  \item{mortality}{ A mortality function of the form \code{mortality =
    function(d, params, ...)}, where \var{d} is diameter and \var{params}
    are as in \var{params.run} below.  This argument should not be a
    concern as \code{crDistnServer} should only be accessed from
    \code{crDistnShiny}. 
  }
  
  \item{recruitment}{ A recruitment function of the form
    \code{recruitment = function(d, params, ...)}, where \var{d} is
    diameter and \var{params} are as in \var{params.run} below.  This
    argument should not be a concern as \code{crDistnServer} should only
    be accessed from \code{crDistnShiny}.
  }

  \item{params.run}{ A vector of current CR growth and vital rate
	parameters.
  }

  \item{ssd.cr}{ The Chapman-Richards distribution function.  This
    argument should not be a concern as \code{crDistnServer} should only
	be accessed from \code{crDistnShiny}.
  }
  
	   
  \item{\dots}{ Passed on. }

}


\details{

Running the \code{crDistnShiny} function will display the application in
your web browser. Its use is self-explanatory, there are several sliders
for input that change both the graphs and the tabular output essentially
instantaneously. The best thing to do is play with it.

Both \code{crDistnServer} and \code{crDistnUI} are meant to be called only
from \code{crDistnShiny}. So please do not run these yourself (I know you
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

Please use the \dQuote{Exit Application} button to stop the
application gracefully. Otherwise at the command line type control-c or
escape, whatever works on your system.

All of the arguments mentioned above can be changed when starting the
\code{crDistnShiny} function. \strong{Please be careful}, there is no
checking in the code on the ultimate limits that are assigned to the
ranges for the parameters, and it is quite easy (even with the defaults)
to generate errors. The interplay of the parameters is such the defaults
are a reasonable compromise, without restricting the density too
much.

The default parameter values are for \dQuote{English} units, they will
generate a very small stand if you choose \dQuote{metric}. The best way
around this is to increase the largest diameter using the \var{eta}
parameter, decrease mortality, and increase recruitment. All these taken
together will get you in the ballpark for a more realistic metric
stand. \code{crDistnShiny} does not distiguish diameter units, so the
default parameters will give the same range in cm as in inches with no
conversion, this is why the default settings when used in metric will
generate a stand with very low basal area. The \dQuote{metric} setting
in the \var{units} argument is basically used only to convert the basal
area-size distribution to the correct units. This is why one must play
with the parameters a bit to get a realistic metric stand.

Note that the width of the largest diameter class will vary with each
change of parameters as it depends on the three Chapman-Richards growth
parameters. You may notice that the superimposed numbers or basal area
curves end lower or higher than the bars in this last class, this is due
to the width of the class being larger or smaller, respectively, than
the nominally chosen class width. This behavior is fine and is to be
expected; when \code{includeTiny=TRUE} the smallest class will extend to
\code{d0} and the curve will display similar behavior as the largest
class. The program now lets you know about the sizes of the lower and
upper class to remind you when they are not of nominal width.

Note that the numbers and basal area density bars, as well as the values
given in the tabular output are exact and come from integrate the
respective numbers density and basal area-size density. Midpoint
approximations are also calculated and returned if
\code{returnVar=TRUE}. If requested, the beta density curve is also
exact, but the tabular values and bar heights are from the midpoint
approximation (just to reduce calculations for quicker response time).
  
}


\value{
  Returns an object of class \dQuote{\code{shiny.appobj}}. If you want
  to capture this, simply assign it to a variable. Otherwise, the return
  value gets \dQuote{printed}, which runs the application.

  However, if you choose to save the results of the run using the
  \code{returnVar} argument, then a \code{list} is placed in the session
  \code{.GlobalEnv} with the following components\ldots

  \item{dcl }{ Diameter class parameters.}
  \item{dcWidths}{Dbh class widths (the last can be different than the
	others).}
  \item{totals}{Stand totals: numbers, basal area.}
  \item{params.run}{The final parameter values for the session.}
  \item{df}{A data frame with the tabular results (exact).}
  \item{df.approx}{A data frame with the midpoint approximations.}
}

\references{
  
\describe{
%
\item{}{Please see \url{http://shiny.rstudio.com/articles/} for more
  information about Shiny and its use.}
%
\item{}{Gove, J. H. 2015. A demographic study of the exponential
  distribution as applied to uneven-aged forests.
  \emph{Submitted}.}
%
\item{}{Gove, J. H. and T. B. Lynch. 2015. Notes on the Chapman-Richards size-structured
  distibution. \emph{Unpublished \R vignette} (Available from the author).}
} %describe
}

\author{
J. H. Gove
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The \pkg{\link{shiny}} package. Also see \code{\link{qDistnShiny}}.
}
\examples{
#
# cut and paste these examples if desired...
#
\dontrun{
crDistnShiny(recruitRange = c(1, 20))  #run immediately
crsh = crDistnShiny()                  #assign to crsh
str(crsh)                              #look at it
crsh                                   #run it
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
