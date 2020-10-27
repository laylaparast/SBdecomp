\name{pred.smooth}
\alias{pred.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Nonparametric conditional mean
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates expected outcome given a single predictor (conditional mean) using a nonparametric kernel estimator
}
\usage{
pred.smooth(zz, zi.one, bw = NULL, y1, weight = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zz}{
%%     ~~Describe \code{zz} here~~
predictor; conditional mean function is applied to this predictor to produce predicted outcome
}
  \item{zi.one}{
%%     ~~Describe \code{zi.one} here~~
predictor; used to estimate conditional mean function (often same as zz)
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
bandwidth
}
  \item{y1}{
%%     ~~Describe \code{y1} here~~
outcome
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
weight
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
expected outcome given predictor
}
\references{
%% ~put references to the literature/web site here ~
Parast, L and Griffin, BA (2020). "Quantifying the Bias due to Observed Individual Confounders in Causal Treatment Effect Estimates." Statistics in Medicine, In press.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\keyword{smooth}
\keyword{robust}
\keyword{nonparametric}
\keyword{internal}
