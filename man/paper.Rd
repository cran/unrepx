%    Copyright (c) 2017 Russell V. Lenth                                %

\name{daniel.paper}
\alias{daniel.paper}

\title{Half-normal graph paper}
\description{The \code{daniel.paper} function draws a graph-paper grid suitable for manually constructing a half-normal plot of effects. The paper is customized to the number of effects to be plotted, thus making the task as simple as possible.}

\usage{
daniel.paper(n.effects = 15, linear = 0:40, a = 0.375, 
    rank.axis = c("y", "x"), type = c("half.normal", "normal"), 
    theme = "blue", lin.lab, rank.lab)
}

\arguments{
  \item{n.effects}{Numeric vector of effects or contrasts to be explored.}
  \item{linear}{sequence of numbers to use for the linear scale.}
 \item{a}{The adjustment used in scaling and centering ranks in the interval (0, 1). The \eqn{i}th ordered (half-) normal score is computed as the \eqn{(i - a) / (n + 1 - 2*a)} quantile of the reference distribution. The value should always be less than 1, and is recommended to be in [0,.5]. The scale becomes somewhat more nonlinear as \eqn{a} increases. }
  \item{rank.axis}{Character value of \code{"y"} or \code{"x"} selecting whether the ranks of effects should be plotted as the vertical or the horizontal coordinate.}
  \item{type}{Character value matching \code{"half.normal"} or \code{"normal"} specifying whether the scaling is for a half-normal or a normal plot. The latter is not recommended.}
  \item{theme}{Character giving the name of the theme to use for the grid lines. The default of \code{"blue"} plots the major grid lines in blue and the minor ones in cyan. Other built-in themes are \code{"bw"} (black majors, and dotted black minors), \code{"gray"} (black majors and gray minors), and \code{"bold"} (all lines black, but major lines are heavier). See details below for information on how to create one's own theme.}
  \item{lin.lab, rank.lab}{Labels for the linear and rank scales.}
}

\details{
Use of \code{type = "normal"} is discouraged, as is discussed in the Details section in \code{\link{hnplot}}.

A custom theme, say \code{"foo"}, may be created by creating a \code{list} named \code{foo_thm}, which should contain two elements named \code{major} and \code{minor}, and (optionally) \code{axis}. \code{major} must be a named \code{list} with graphical parameters \code{col}, \code{lty}, and \code{lwd} as described in \code{\link{par}}. \code{minor} only requires \code{col} and \code{lty}, as minor lines are always drawn using \code{lwd = 1}. \code{axis} must contain \code{col} and \code{lwd}, which are used for the tick marks and labels on the rank scale. If \code{axis} is not included, \code{major} is used.
}

\value{Nothing is returned.}

\references{
Daniel, C (1959)
  Use of Half-Normal Plots in Interpreting Factorial Two-Level Experiments.
  \emph{Technometrics}, 1(4), 311-341
}

\author{
Russell V. Lenth
}

\seealso{
\code{\link{hnplot}}
}

\examples{
require("unrepx")

old.par <- par(mar = c(5, 5, 1, 1) + .1)
daniel.paper(15)

neon_thm <- list(
    major = list(col = "green", lty = 1, lwd = 2),
    minor = list(col = "orange", lty = 1),
    axis = list(col = "magenta", lwd = 2)
)
daniel.paper(31, rank.axis = "x", theme = "neon")

par(old.par)
}

\keyword{design}
\keyword{hplot}
