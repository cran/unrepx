%    Copyright (c) 2017 Russell V. Lenth

\name{ref.dist}
\alias{ref.dist}
\alias{eff.test}

\title{Analysis of effects from screening experiments}
\description{These functions facilitate making inferences based on effect estimates in an unreplicated experiment, withn an underlying effect-sparsity model.}
\usage{
ref.dist(method, n.effects, nsets, save = TRUE)

eff.test(effects, method = "Zahn", pareto = TRUE, refdist, save = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{The method to use in determining the reference line, curve, and/or critical values. This must be the name of a provided pseudo-standard-error method (see \code{\link{PSE}}), or a compatible user-supplied one.}
  \item{n.effects}{Integer number of effects estimated.}
  \item{nsets}{The number of complete-null samples of size \code{n.effects} to be simulated. If omitted, \code{nsets} is determined so that the total number of simulated effects is about 40,000.}
  \item{save}{Logical value. If \code{TRUE}, the simulated reference distribution is saved in the workspace under \code{.Last.ref.dist}. Other routines in this package try to avoid re-simulating a reference distribution if \code{.Last.refdist} exists and matches the current \code{method} and \code{n.effects}.}
  \item{effects}{Vector of observed effects to be tested against the reference distribution.}
  \item{pareto}{Logical value. If \code{TRUE}, the effects are presented in decreasing order of their absolute size.}
  \item{refdist}{A result of a previous call to \code{ref.dist}, in case the user wishes to manually supply a previously simulated reference distribution. Note however that \code{eff.test} will automatically reuse \code{.Last.ref.dist} if it is available and matches.}
}
\details{
\code{ref.dist} simulates samples of effects from the standard normal dstribution. For each sample, the pseudo standaerd error (PSE) of the effects (typically some kind of outlier-resistant estimate of the SD) is obtained via a call to \code{\link{PSE}} with specified method. The absolute \eqn{t} values are obtained as ratios of the simulated effects and the PSE, as well as and the maxima of these absolute \eqn{t} values. Quantiles and tail areas of these simulated distributions then form a reference for obtaining critical values and P values in testing an observed sample of effects.

\code{eff.test} performs a traditional-style analysis for an observed sample of effects. It outputs the effects, PSE, \eqn{t} ratios; and uses tail areas of the associated reference distribution to compute individual and simultaneous \eqn{P} values. The simultaneous \eqn{P} values implement a multiplicity correction for \emph{any} type-I errors occurring among the tests. 
}

\value{
  \code{ref.dist} returns an object of class \code{"eff_refdist"} -- structurally, a \code{list} with elements \code{abst} (the absolute values of the simulated \eqn{t} statistics), \code{max.abst} (the sample maxima of \code{abst}), and \code{sig} (a signature of the form \code{method_n.effects}). There is a \code{print} method for this class that displays a summary.

  \code{eff.test} returns a \code{data.frame} containing the estimates, \eqn{t ratios}, and estimated P values as tail areas of \code{abst} and \code{max.abst} from the reference distribution.}
  
%%\references{}

\author{
Russell V. Lenth
}

%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{}

\examples{
require("unrepx")

zahn15 <- ref.dist("Zahn", 15)
eff.test(pdEff, refdist = zahn15)
}

\keyword{design}
\keyword{htest}
