%    Copyright (c) 2017 Russell V. Lenth                                %

\name{parplot}
\alias{parplot}

\title{Pareto plot of effects}
\description{Constructs a bar plot of ordered effects, along with cutoff values for the margin of error (ME) and simultaneous margin of error (SME). Such a plot is suggested in Lenth (1989), but other methods may be used for obtaining the ME and SME.}

\usage{
parplot(effects, pareto = TRUE, absolute = TRUE, horiz = FALSE, col = absolute,
    critvals, method = "Zahn", alpha = .05, refdist, sim.opts, 
    ylab = "Estimated effects", top = n.effects, cex.annot = 0.75, ...)
}

\arguments{
  \item{effects}{Numeric vector of effects or contrasts to be explored.}
  \item{pareto}{Logical value. If \code{TRUE}, the effects are plotted in decreasing order of their absolute values.}
  \item{absolute}{Logical value. If \code{TRUE}, the absolute effects are plotted. If \code{FALSE}, the original signed effects are used, so that there are potentially positive- and negative-going bars in the plot.}
  \item{horiz}{Logical value. If \code{TRUE}, the bars are horizontal, and if \code{FALSE}, they are vertical.}
  \item{col}{A logical value, or valid color code(s) or names(s). If logical, a value of \code{TRUE} shades the positive effects in light blue, and the negative effects in pink. A logical value of \code{FALSE} colors them all light gray.}
  \item{critvals}{Numeric value(s). If these are provided, the first two elements of \code{critvals} are used as the ME and SME respectively (on the absolute scale of the effects). When \code{critvals} is specified, \code{method}, \code{alpha}, and \code{refdist} are ignored.}
  \item{method}{Character value designating the method to use in determining the margins of error displayed in the plot when \code{critvals} is not given. This must be the name of a provided pseudo-standard-error method (see \code{\link{PSE}}), or a compatible user-supplied one.}
  \item{alpha}{Numeric value. A null reference distribution for \code{method} is used (see \code{\link{ref.dist}}) to determine a margin of error (labeled \sQuote{ME}) and simultaneous margin of error (labeled \sQuote{SME}) corresponding to a significance level of \code{alpha}, and reference lines are added to the plot at those positions as an aid to assessing the statistical significance of the effects.}
  \item{refdist}{A result of \code{\link{ref.dist}}. If given, it is used to obtain critical values, rather than running a new simulation of the null distribution. The user should be careful that \code{refdist} indeed matches \code{method} and the number of effects.}
  \item{sim.opts}{A \code{list} containing arguments \code{nsets} and/or \code{save} to pass to \code{\link{ref.dist}} in case a new reference distribution needs to be simulated. See also details below.}
  \item{ylab}{Character axis label (overrides the default).}
  \item{top}{Numeric value giving the number of effects to display (this may help make all the important effect names visible). When \code{top} is less than the number of effects (\code{n.effects}), this forces \code{pareto = TRUE} and only the largest \code{top} effects are displayed. When this happens, an annotation is added to the plot to help clarify that not all effects are displayed.}
  \item{cex.annot}{Character magnification for annotations}
  \item{\dots}{Additional graphical parameters (see \code{\link{par}}) used in constructing the plot.}
}

\details{
The cutoff values displayed in the plot are labeled \dQuote{ME}, the margin of error, and \dQuote{SME}, the simultaneous margin of error. If not specified using \code{crtivals}, they are obtained using the 1-\code{alpha} quantiles of the reference distribution of absolute pseudo-\eqn{t} ratios. ME is based on the distribution of \eqn{|t|}. SME is based on the distribution of the maximum \eqn{|t|} for a whole set of null effects.

In determining cutoff values, \code{parplot} tries to avoid re-simulating the reference distribution. Specifically, if the global variable \code{.Last.ref.dist} exists, and its contents match the given \code{method} and number of effects, it is used as the reference distribution. Similarly, if \code{refdist} is supplied, it is used (without checking). If a suitable reference distribution is not found, then it is simulated via \code{\link{ref.dist}}, with any arguments from \code{sim.opts} added.

If \code{critvals} is supplied, the specified values are used as the ME and SME: no reference distribution is needed, and hence \code{method}, \code{alpha}, and \code{refdist} are ignored. 

The plot is scaled so that the ME cutoff always shows. The SME cutoff will only be visible if an observed effect is near or exceeds that boundary. The numeric values of the ME and SME are also shown in an annotation in the plot.
}

\value{
Invisibly, the vector of the ME and SME values.
}
\references{
Lenth, R (1989)
Quick and Easy Analysis of Unrelicated Factorials.
\emph{Technometrics} 31(4), 469-473

}

\author{
Russell V. Lenth
}

%%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For more details on PSEs and reference distributions, see \code{\link{PSE}} and \code{\link{ref.dist}}. Note that \code{parplot} produces in essence a graphical version of the information from \code{\link{eff.test}}, but the latter provides more resolution in terms of \eqn{P} values.

Other graphical ways of assessing active effects include a dot plot with a reference curve (\code{\link{refplot}}) and a half-normal plot (see \code{\link{hnplot}}).
}

\examples{
require("unrepx")

parplot(pdEff, top = 10)
}

\keyword{design}
\keyword{hplot}
\keyword{htest}
