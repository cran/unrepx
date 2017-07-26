%    Copyright (c) 2017 Russell V. Lenth                                %

\name{dot.plot}
\alias{dot.plot}
\alias{dot.id}
\alias{dot.mod}

\title{Resizeable dot plots}
\description{Flexible stacked-dot plots of relatively small samples, and provisions for identifying points. The dot plot is responsive to resizing of the plot window: a vertical resizing does not affect the vertical spacing of the dots, and a horizontal resizing may lead to re-binning to fit well within the range.}
\usage{
dot.plot(x, pch = 16, cex.dot = 1, spacing = 1, xlab, xlim = range(x), ...)

dot.id(env, height.id = 2, cex.id = 1, col.id = "black")
dot.mod(env, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric vector of values to be plotted.}
  \item{pch}{Plotting character for the dots. Default is a filled circle.}
  \item{cex.dot}{Size of dots relative to \code{par("cex")}.}
  \item{spacing}{Factor for adjusting the vertical spacing of stacked dots.}
  \item{xlab, xlim}{The usual graphical parameters (see \code{\link{par}}), but used by \code{dot.plot} or \code{dot.id}.}
  \item{env}{Environment returned by \code{dot.plot}.}
  \item{height.id}{Height (in character heights) above the horizontal axis for displaying identification labels.}
  \item{cex.id, col.id}{\code{cex} and \code{col} settings for labeling identified points.}
  \item{\dots}{In \code{dot.plot}, additional graphical parameters (see \code{\link{par}}) used in constructing the plot. In \code{dot.mod}, parameters to modify, from among \code{cex.dot}, \code{cex.id}, \code{col.id}, and \code{height.id}.}
}
\details{
\code{dot.id} and \code{dot.mod} work only with interactive graphics devices.

The \code{dot.id} works similarly to \code{\link{identify}}, but all dots in a particular stack are identified at once. The user should click along the bottom row of dots. No labels are displayed until the user exits identify mode (and then it may be necessary to refresh the plot by resizing it slightly). Also, the points that are labeled may change if the plot is resized horizontally, because only the x values are actually identified.

The user may call \code{dot.id(..., modify = FALSE)} more than once, if it is desired to change which values are identified.

After a call to \code{dot.id} or (especially) to \code{dot.mod}, the graph may need to be refreshed. There appears to be no standard way to do this, so the user may need to resize its window slightly.
}
\value{
\code{dot.plot} returns an environment that is used to hold information on points that are identified. \code{dot.id} modifies this environment; it (invisibly) returns the currently identified x values. (These are values from the original call to \code{dot.plot}, not their rounded values used in the plot.)}

%%\references{}

\author{
Russell V. Lenth
}

%%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{}

\examples{
require("unrepx")

educ = swiss$Education
names(educ) = abbreviate(row.names(swiss))
dp <- dot.plot(educ, xlab = "Percent Post-Primary Education (Switzerland, 1888)")

\dontrun{
    dot.id(dp)
    
    dot.mod(dp, height.id = 3, cex.id = .5)
}
}

\keyword{hplot}
\keyword{iplot}
