% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_functions.R
\name{LLO}
\alias{LLO}
\title{Linear Log Odds (LLO) Recalibration Function}
\usage{
LLO(x, delta, gamma)
}
\arguments{
\item{x}{a numeric vector of predicted probabilities of an event. Must only
contain values in [0,1].}

\item{delta}{numeric, must be > 0, parameter \eqn{\delta} in LLO
recalibration function.}

\item{gamma}{numeric, parameter \eqn{\gamma} in LLO recalibration function.}
}
\value{
Vector of LLO-adjusted probabilities via specified \eqn{\delta} and
\eqn{\gamma}.
}
\description{
LLO-adjust predicted probabilities based on specified \eqn{\delta} and
\eqn{\gamma}.
}
\details{
The Linear Log Odds (LLO) recalibration function can be written as
\deqn{c(x_i;\delta, \gamma) = \frac{\delta x_i^\gamma}{\delta x_i^\gamma +
(1-x_i)^\gamma}} where \eqn{x_i} is a predicted probability,
\eqn{\delta > 0} and \eqn{\gamma \in \mathbb{R}}.  Then \eqn{c(x_i;\delta,
\gamma)} is the corresponding LLO-adjusted probability that has been shifted
by \eqn{\delta} and scaled by \eqn{\gamma} on the log odds scale.  When
\eqn{\delta = \gamma = 1}, there is no shifting or scaling imposed on \code{x}.
}
\examples{

# Vector of probability predictions from 0 to 1
x1 <- seq(0, 1, by=0.1)
x1

# LLO-adjusted predictions via delta = 2, gamma = 3
x1_llo23 <- LLO(x1, 2, 3)
x1_llo23  

# LLO-adjusted predictions via delta = 1, gamma = 1
x1_llo11 <- LLO(x1, 1, 1)
x1_llo11  # no change

# Create vector of 100 probability predictions
x2 <- runif(100)

# LLO-adjust via delta = 2, gamma = 3
x2_llo23 <- LLO(x2, 2, 3)

plot(x2, x2_llo23)
}
\references{
Turner, B., Steyvers, M., Merkle, E., Budescu, D., and Wallsten,
T. (2014) Forecast aggregation via recalibration, \emph{Machine Learning}
95, 261–289.

Gonzalez, R., and Wu, G. (1999), On the shape of probability weighting
function, \emph{Cognitive Psychology} 38, 129–66.
}
