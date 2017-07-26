%    Copyright (c) 2017 Russell V. Lenth                                %

\name{hnplot}
\alias{hnplot}

\title{Half-normal plots (Daniel plots) of effects}
\description{The \code{hnplot} function constructs a (half-) normal plot of effects (see Daniel 1959) that is traditionally used to identify active effects in a screening experiment. Reference lines and various other options and extensions are supported.}

\usage{
hnplot(effects, ref = TRUE, half = TRUE, horiz = TRUE, method = "Zahn", 
       a = 0.375, col = half, pch = 16, ID = FALSE, alpha, ...)
}

\arguments{
  \item{effects}{Numeric vector of effects or contrasts to be explored.}
  \item{ref}{Logical value. If \code{TRUE}, a reference line is added to the plot determined by \code{method}. If \code{FALSE}, no reference is shown.
}
  \item{half}{Logical value. If \code{TRUE}, a plot based on the absolute effects is constructed. If \code{FALSE} (not recommended: see Details), the original signed effects are used.}
  \item{horiz}{Logical value. If \code{TRUE}, the (absolute) effects are plotted on the horizontal scale and the (half-) normal scores or labels are plotted on the vertical scale. If \code{FALSE}, these axes are reversed.}
  \item{method}{Character value. When \code{ref} is \code{TRUE}, the method to use in determining the reference line, curve, and/or critical values. This must be the name of a provided pseudo-standard-error method (see \code{\link{PSE}}), or a compatible user-supplied one.}
  \item{a}{The adjustment used in scaling and centering ranks in the interval (0, 1). The \eqn{i}th ordered (half-) normal score is computed as the \eqn{(i - a) / (n + 1 - 2*a)} quantile of the reference distribution.}
  \item{col}{Scalar or vector of colors; or a logical value. If logical, a value of \code{TRUE} colors the positive effects blue, the negative effects red, and any zeros as black. A logical value of \code{FALSE} colors them all black.}
  \item{pch}{Plotting character(s) to use.}
  \item{alpha}{Numeric value. If specified, a null reference distribution for \code{method} is used (see \code{\link{ref.dist}}) to determine a margin of error (labeled \sQuote{ME}) and simultaneous margin of error (labeled \sQuote{SME}) corresponding to a significance level of \code{alpha}, and reference lines are added to the plot at those positions as an aid to assessing the statistical significance of the effects. This is based on a suggestion in Mee (2015). These reference lines are omitted when \code{alpha} is left unspecified.}
  \item{ID}{Logical or numeric value. If logical and \code{TRUE}, then after the plot is contructed, the plot is put in \code{\link{identify}} mode, where the user may click on points to be labeled on the plot. If a numeric value is supplied, it is used as a threshold by which all effects greater than \code{ID[1]} in absolute value are labeled.}
  \item{\dots}{Additional graphical parameters (see \code{\link{par}}) used in constructing the plot.}
}

\details{
Use of \code{half = FALSE} is not recommended because it can be misleading to the user. Inactive effects are those that are close to zero, and a regular normal plot displays deviations from normality rather than deviations from zero.
}

%\value{}

\references{
Daniel, C (1959)
  Use of Half-Normal Plots in Interpreting Factorial Two-Level Experiments.
  \emph{Technometrics}, 1(4), 311-341
  
Mee, R (2015) 
  Discussion: Better, not Fewer, Plots. 
  \emph{Journal of Quality Technology}, 47(2), 107-109
}

\author{
Russell V. Lenth
}

%%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Other ways of assessing active effects include a dot plot with a reference curve (\code{\link{refplot}}), a pareto plot of effects (see \code{\link{parplot}}), and a tabular style of presenting effects and \eqn{P} values (see \code{eff.test}). For more information on methods, see \code{\link{PSE}} and \code{\link{ref.dist}}.
}

\examples{
require("unrepx")

hnplot(pdEff, ID = ME(pdEff))
}

\keyword{design}
\keyword{hplot}
