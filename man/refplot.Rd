%    Copyright (c) 2017 Russell V. Lenth                                %

\name{refplot}
\alias{refplot}

\title{Dot plot of effects with a reference distribution}
\description{The \code{refplot} function constructs a dot plot of effects along with a rweference distribution (either normal or simulated) to help in visually identifying active effects in a screening experiment.}

\usage{
refplot(effects, ref = TRUE, half = TRUE, method = "Zahn", 
    col = half, guides = FALSE, ID = FALSE, pch = 16, xlab, xlim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effects}{Numeric vector of effects or contrasts to be explored.}
  \item{ref}{Logical or character value. If \code{TRUE}, a reference line or curve is added to the plot determined by \code{method}. If \code{FALSE}, no reference is added to the plot.

A character value matching \code{"normal"} or \code{"simulated"} is also permitted.
\code{"normal"} is equivalent to \code{ref = TRUE}. With \code{"simulated"}, 
, a kernel density estimate is displayed; it is obtained by scaling the simulated reference distribution for \code{"method"} by its observed PSE. }
  \item{half}{Logical value. If \code{TRUE}, a dot plot of the absolute effects is constructed. If \code{FALSE}, the original signed effects are plotted.}
  \item{method}{Character value. When \code{ref} is not false, the method to use in determining the PSE for scaling the reference curve (and also the method used for simulating the referemnce distribution when \code{ref = "simulated"}). This must be the name of a provided pseudo-standard-error method (see \code{\link{PSE}}), or a compatible user-supplied one.}
  \item{col}{Scalar or vector of colors; or a logical value. If logical, a value of \code{TRUE} colors the positive effects blue, the negative effects red, and any zeros as black. A logical value of \code{FALSE} colors them all black.}
  \item{guides}{Logical value. If \code{TRUE}, dotted lines are added that illustrate guide lines that could be used to draw the normal curve by hand.}
  \item{ID}{Logical value. If logical and \code{TRUE}, then after the plot is contructed, the \code{\link{dot.id}} invoked so that the user may click on points to be labeled on the plot. If a numeric value is given, it is used as a threshold by which all effects greater than \code{ID[1]} in absolute value are labeled.}
  \item{pch, xlab, xlim, \dots}{Additional graphical parameters (see \code{\link{par}}) passed to \code{\link{dot.plot}}.}
}

\details{
If the returned environment is saved, then \code{\link{dot.id}} or \code{\link{dot.mod}} may be used later as for \code{\link{dot.plot}} results.}
\value{
An environment that can be modified using \code{\link{dot.id}} or \code{\link{dot.mod}}.
}

%%\references{ }

\author{
Russell V. Lenth
}

%%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Other ways of assessing active effects include a half-normal plot (\code{\link{hnplot}}), a Pareto plot of effects (see \code{\link{parplot}}), and a tabular style of presenting effects and \eqn{P} values (see \code{eff.test}). For more information on PSEs and methods, see \code{\link{PSE}} and \code{\link{ref.dist}}.

}

\examples{
require(unrepx)

refplot(pdEff, ID = ME(pdEff))

\dontrun{
    # Batman lives!
    refplot(pdEff, ref = "sim", method = "Lenth", half = FALSE)
}
}

\keyword{design}
\keyword{hplot}
